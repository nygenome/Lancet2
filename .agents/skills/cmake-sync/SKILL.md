---
name: cmake-sync
description: Use when adding a new source file, a new dependency, a new build target, a new build option, or when changing how an existing target links or compiles. Trigger on "I added a new .cpp file; do I need to update CMakeLists?", "add this new dep", "the test isn't being built", "what does LANCET_PROFILE_MODE do?", "the new layer needs its own library target", "this header isn't being installed". Walks the dependency declaration → target wiring → option propagation chain so all four pieces stay coordinated. Aware of the sanitizer build trees (LANCET_SANITIZE_BUILD), the profiling tree (LANCET_PROFILE_MODE), the static-link option (LANCET_BUILD_STATIC), and the test/benchmark inclusion options. Does NOT cover modifying the cmake/ subdirectory's superbuild fragments — those are protected and changes there go through normal review.
allowed-tools: Read, Glob, Grep, Edit, Write, Bash
---

# CMake sync on Lancet2

The Lancet2 build is structured around a top-level `CMakeLists.txt` (~450 lines) plus the per-layer `tests/CMakeLists.txt` and `benchmarks/CMakeLists.txt`. Dependencies are declared in `cmake/dependencies.cmake` (FetchContent-based, vendored under `cmake-build-*/_deps/`). Build options live near the top of the top-level file and propagate through the rest via standard CMake plumbing.

Most edits the agent will reasonably make are in the top-level `CMakeLists.txt` and the test/benchmark subdirectory CMakeLists. The `cmake/` modules — `dependencies.cmake`, `optimization_flags.cmake`, `platform_checks.cmake`, etc. — are write-blocked by the protected-paths hook. Changes there require explicit user override and go through normal review (they're high-blast-radius: a wrong toolchain fragment changes behavior across the codebase).

## Step 1 — Identify the change category

Six categories cover essentially all build-system changes:

1. **Adding a source file to an existing target.** A new `.cpp` in `src/lancet/<layer>/` needs a corresponding entry in the layer's source list. A new `_test.cpp` in `tests/<layer>/` needs an entry in `tests/CMakeLists.txt`.
2. **Adding a new layer or major module.** A new `src/lancet/<new-layer>/` requires its own subdirectory CMakeLists, library target declaration, and a `target_link_libraries` line on the dependent layers (per the layer-direction rule in `validate_layer_direction.py`).
3. **Adding a new dependency.** A new `FetchContent_Declare` in `cmake/dependencies.cmake`, plus `target_link_libraries(<consumer> ... <new-dep>)` wherever it's consumed. Touches `cmake/` so it's a protected-paths edit.
4. **Adding a new build option.** New `option(LANCET_X ...)` near the existing options, plus the conditional logic that reads the option, plus pixi tasks (in `pixi.toml`) that set it.
5. **Changing a target's link line.** Adding/removing a library from `target_link_libraries`, changing PRIVATE/PUBLIC visibility, adding `target_compile_definitions`. Common case for sanitizer support (`LANCET_SANITIZE_BUILD` excludes mimalloc-static).
6. **Changing a target's compile flags.** Adding/removing flags via `target_compile_options` or `CMAKE_CXX_FLAGS`. Rare; usually goes through `cmake/optimization_flags.cmake` which is protected.

Identify which category the user is in before editing.

## Step 2 — Find the coordinated locations

For each category, the touched files are:

**Adding a source file (category 1):**
- The `add_executable` or `add_library` block in `CMakeLists.txt` (or the layer subdirectory CMakeLists). Find the existing entries via `grep -n "<existing_file>.cpp" CMakeLists.txt`.
- For test files: `tests/CMakeLists.txt`'s `add_executable(TestLancet2 ...)` block lists every test source explicitly. Add the new file in the layer-grouped section (the file is organized by `# Layer 1: base`, `# Layer 2: hts`, etc.).
- For benchmark files: `benchmarks/CMakeLists.txt` similarly.

**Adding a new layer (category 2):**
- New subdirectory under `src/lancet/<new-layer>/` with its own `CMakeLists.txt` defining `add_library(lancet_<new-layer> ...)`.
- Top-level `CMakeLists.txt` adds an `add_subdirectory(src/lancet/<new-layer>)` and updates the layer chain in `validate_layer_direction.py` (in the bundle's `.Codex/hooks/`).
- The new library is added to `target_link_libraries` of any layer that depends on it.
- The new layer's path-scoped rule file: `.Codex/rules/<new-layer>.md` (per the rules README).

**Adding a new dependency (category 3, requires protected-path override):**
- `cmake/dependencies.cmake` adds `FetchContent_Declare(<dep> ...)` and `FetchContent_MakeAvailable(<dep>)`. Mirror the existing entries: GIT_REPOSITORY, GIT_TAG, SYSTEM keyword.
- `target_link_libraries(<consumer-target> PRIVATE <dep-target>)` wherever the dep is used.
- If the dep needs special flags or includes, those go in the consuming target's `target_*_directories` / `target_compile_options`.
- Add the dep to the `# mimalloc — high-performance allocator (replaces system malloc)` style header comment block at the top of `cmake/dependencies.cmake` so future readers know the role.

**Adding a build option (category 4):**
- `option(LANCET_X "..." "OFF")` near the existing options. The default is almost always `OFF` (opt-in).
- `set_property(CACHE LANCET_X PROPERTY STRINGS "OFF" "ON")` so `ccmake` shows it as a toggle.
- Conditional logic somewhere — often `if (LANCET_X)` in the same CMakeLists, or a `target_compile_definitions(... LANCET_X=1)` to surface it to source.
- Document the option in the comment block above where it's defined; this is where future readers look.
- If the option needs a pixi task to set it, add the task to `pixi.toml`.

**Changing link line (category 5):**
- Just the `target_link_libraries` block in `CMakeLists.txt` (or subdirectory). Simple, localized.

**Changing compile flags (category 6):**
- Usually `cmake/optimization_flags.cmake` (protected). Discuss with user before editing.

## Step 3 — Verify the wiring is correct

After every CMakeLists edit, run a clean reconfigure of the affected build tree:

```bash
pixi run configure-release
```

(or `configure-debug`, `configure-profile`, or one of the sanitizer configures depending on what's affected). A clean reconfigure surfaces:

- Missing source files (`Cannot find source file ...`).
- Missing libraries (`Target ... links to target ... but the target ...`).
- Option-mismatch errors (`<option> requires <other-option>`).
- Cyclic target dependencies.

If the reconfigure passes, build:

```bash
pixi run build-release
```

And test:

```bash
pixi run test    # depends on build-release; runs cmake-build-release/tests/TestLancet2
```

A green build + green tests is the gate for any CMake change. If the build passes but tests fail, the change probably affected linkage in a way that's correct at compile time but wrong at runtime — typical when adding a new dep that has a static initializer or a new target that's missing a runtime library.

## Step 4 — Update the layer rule if a new layer was added

If category 2 (new layer), the path-scoped rule file in `.Codex/rules/<new-layer>.md` needs to be created. Use one of the existing rules (`base.md`, `hts.md`, etc.) as a template. The frontmatter must have:

```yaml
---
description: <one-line description loaded into context at session start>
paths:
  - "src/lancet/<new-layer>/**"
---
```

The body documents the design rules, threading invariants, and bug classes the new layer prevents. Without this, agents editing the new layer get no scoped guidance — only AGENTS.md applies.

The rule file goes in the same commit as the CMakeLists changes. The bundle's `validate_layer_direction.py` hook also needs updating to know about the new layer in the layer chain — that's a `.Codex/hooks/` edit, which goes through the bundle's `/audit-bundle` review.

## Step 5 — Cross-validate with the architecture guide

If the change affects layer structure (category 2) or introduces a new major dependency (category 3), `docs/guides/architecture.md` describes the layer chain at a high level. Update it via the `doc-sync` skill in the same commit.

For categories 1, 4, 5, 6 the architecture guide usually doesn't need changes.

## Step 6 — Commit with the right scope

> **Note:** `lancet_assembly` below is a hypothetical example layer used to illustrate the sync pattern; the real Lancet2 layers are `base`, `hts`, `cbdg`, `caller`, `core`, `cli`.


The Lancet2 chglog filter accepts only `feat`, `fix`, `perf`, and `chore` as commit types; `build:`, `refactor:`, `docs:`, etc. parse but are silently dropped from `CHANGELOG.md`. Pure build-system changes (a new layer, a new build option, a CMakeLists refactor) take `chore:` because the user does not see them in pipeline output. A new layer that ships a user-visible feature alongside the build wiring takes `feat:`. A build change that fixes a broken build takes `fix:`.

```
chore: add new lancet_assembly layer

New library target lancet_assembly under src/lancet/assembly/.
Sits between cbdg and caller in the layer chain. Houses the
SPOA-driven multiple-sequence-alignment code that was inline
in caller previously.

CMakeLists.txt: new add_subdirectory and target_link_libraries
   updates on caller and dependents.
.Codex/rules/assembly.md: new path-scoped rule.
validate_layer_direction.py: layer chain updated to include
   assembly between cbdg and caller.
docs/guides/architecture.md: layer diagram updated.
```

For routine source-file additions (category 1), `feat:` or `fix:` with the source addition mentioned in the body is sufficient (the chglog header pattern does not support scopes, so write `fix: ...` rather than `fix(<layer>): ...`).

## When NOT to use this skill

Do not use this skill for:

- Editing `cmake/optimization_flags.cmake` or `cmake/platform_checks.cmake` directly (write-blocked by protected paths). User must override the protection and the change goes through normal review.
- Editing `pixi.toml` directly (write-blocked). New pixi tasks for new options go through normal review.
- Editing `cmake-build-*/CMakeCache.txt` or any `cmake-build-*/_deps/` content (write-blocked; build-system-managed).
- "Just rebuild the project" — that's a build invocation (`pixi run build-debug`), not a build-system change.

## When the change is a layer rename or a major restructuring

If the user proposes renaming a layer (`cbdg` → `graph`, etc.) or merging two layers, the change is invasive: every `target_link_libraries`, every `validate_layer_direction.py` entry, every layer-rule file, every layer-specific test directory name, and every `src/lancet/<layer>/` path needs to update. This is a larger surface than this skill covers; discuss the structure with the user, plan as a multi-commit refactor, and consider whether the rename earns its breakage cost.
