# AGENTS.md — Lancet2 project memory

This file is the canonical project-memory document for any agentic coding tool used on Lancet2 (Claude Code, GitHub Copilot, Codex, Cursor, etc.). It carries the conventions, build system, architecture, and workflow rules that apply to all changes. Rules that apply only to a specific layer or workflow live in the corresponding subagent, skill, or path-scoped rule under `.claude/rules/`, not here. See `.claude/cost-model.md` for why this file is intentionally short.

A thin wrapper at `CLAUDE.md` imports this file (`@AGENTS.md`) so Claude Code reads it through the conventional CLAUDE.md path; the canonical content lives here for cross-tool portability.

## Project at a glance

Lancet2 is a modern C++20 somatic variant caller (SNVs and InDels) by the New York Genome Center. It performs joint multi-sample localized colored de Bruijn graph assembly. The CLI entry point is `Lancet2 pipeline`. The `main` branch is the source of truth; the documentation site at `nygenome.github.io/Lancet2` typically lags `main`.

## Build and tooling

The build is driven by [pixi](https://pixi.sh) (`pixi task list` for the full catalog). The toolchain (gcc, clang, cmake, ninja, clang-tidy, IWYU, go) is pinned in `pixi.toml` and installed under `.pixi/`. Always prefer pixi tasks over hand-rolled CMake invocations so you pick up the locked toolchain rather than whatever is on `PATH`.

```
pixi run configure-release            # CMake Release configure → cmake-build-release/
pixi run build-release                # Release build (auto-runs configure-release first)
pixi run configure-debug              # Debug configure with tests + compile_commands.json
pixi run build-debug                  # Debug build (auto-runs configure-debug first)
pixi run test                         # Run Catch2 tests (auto-builds Debug first)
pixi run bench                        # Run benchmarks (configure with -DLANCET_BENCHMARKS=ON)
pixi run lint-all                     # fmt-check + lint-check + iwyu-check (CI gates)
pixi run fmt-fix                      # apply clang-format
pixi run iwyu-fix                     # apply IWYU then re-format
pixi run docs-build                   # strict, zero-warning
```

Sanitizer trees (`cmake-build-asan/`, `cmake-build-tsan/`) are configured manually by the `sanitizer-triage` skill because they require disabling the static mimalloc link. The profiling build (`pixi run configure-profile` / `pixi run build-profile`) is documented in the `profile-and-optimize` skill; reach for the skill rather than configuring profiling builds inline.

<important_if context="invoking pixi commands">
**Do NOT run `pixi run lint-fix`.** Clang-tidy auto-fixes have historically broken compilation. Always resolve `lint-check` warnings manually.
</important_if>

## Architecture: six-layer dependency chain

The codebase is organized as six static libraries with a strict dependency direction documented at the top of `CMakeLists.txt`:

```
lancet_base   → types, hashing, repeat/complexity scoring, logging, crash handler
lancet_hts    → HTSlib wrappers: BAM/CRAM/FASTA I/O, CIGAR, SAM flags, BGZF, URI utils
lancet_cbdg   → Colored de Bruijn graph: kmer/node/edge/path, max-flow,
                cycle finder, complexity, probe diagnostics, DOT serializer
lancet_caller → MSA (SPOA), variant bubble extraction, local+combined scorers,
                genotype likelihoods, posterior base quality, FORMAT assembly
lancet_core   → Window builder, active-region detector, read collector,
                variant builder/annotator/store, async worker, pipeline executor
lancet_cli    → CLI11 parsing, VCF header construction, pipeline_runner
Lancet2       → main.cpp, links lancet_cli + mimalloc
```

A file in layer N may depend on layers 1 through N. It must NOT depend on layers above it. The at-write hook in `.claude/hooks/validate_layer_direction.py` enforces this on every edit.

`CMakeLists.txt` documents each layer's intra-module data flow at the top of its `add_library(...)` block. Sources are listed in foundational → orchestrator order, not alphabetical. **When you add or rename files in a layer, update the data-flow comment at the top of that layer's `add_library(...)` block** — it goes stale silently.

Layer-specific design rules — the algorithmic invariants, threading contracts, numerical-stability requirements, and bug classes each layer prevents — live in path-scoped rule files under `.claude/rules/`. These auto-load only when Claude is editing files in the matching layer:

- `.claude/rules/base.md` — Welford numerical stability, AVX2/NEON SIMD primitives, Mann-Whitney effect-size derivation, polar-coords feature engineering, longdust GC-bias correction
- `.claude/rules/hts.md` — HTSlib RAII wrapping, the `bam1_t` non-owning-proxy lifetime contract, per-thread `Extractor`/`Iterator` invariant, `BgzfOstream` lifecycle
- `.claude/rules/cbdg.md` — k-mer canonicalization, BCALM 2 sign continuity, three-track Node read-support model, MaxFlow walk-tree arena, three-color cycle detection
- `.claude/rules/caller.md` — SPOA convex scoring with int16 SIMD lane discipline, Dirichlet-Multinomial PLs in VCF-standard ordering, allele-assignment scoring formula
- `.claude/rules/core.md` — 256-shard `VariantStore`, the chunked sorted VCF flush with `NUM_BUFFER_WINDOWS=100` lag, lock-free worker queues, jthread cooperative cancellation
- `.claude/rules/cli.md` — CLI11 flag surface as public API, `BgzfOstream` open→header→flush→Execute→Close ordering, `BuildVcfHeader` as the canonical FORMAT/INFO declaration site

When Claude is reasoning about correctness inside a layer, the relevant rule file loads automatically and provides design context the source files alone don't make explicit.

## Naming and code conventions

`docs_dev/style/` (entry point: `docs_dev/style/README.md`) is authoritative; `.clang-tidy` is the source of truth for the auto-enforced rules. The conventions below are the ones reviewers flag and are NOT all caught by clang-format or clang-tidy.

**Names.** `PascalCase` for types, methods, free functions. `mPascalCase` for member variables (the `m` prefix is required; trailing underscores are not used). `snake_case` for parameters, locals, and namespaces. `UPPER_CASE` for constants, enumerators, and macros. The `i8`/`i16`/`i32`/`i64`/`u8`/`u16`/`u32`/`u64`/`f32`/`f64`/`usize` aliases (defined in `src/lancet/base/types.h`) are project-wide and exempt from the snake_case rule; new code uses them rather than `int`/`unsigned`/`size_t`/`double`.

**East const, everywhere.** `int const x`, `auto const& ref` — never `const int x` or `const auto& ref`. Clang-format does not fix this; reviewers will.

**Quote vs angle-bracket includes.** `<>` is reserved for the C/C++ standard library only. Project headers, Abseil, spdlog, htslib, fmtlib, SPOA, minimap2, moodycamel, Catch2 — all use `""`. Clang-format will not fix this. C-only third-party headers (htslib, minimap2) are wrapped in `extern "C" {}` with a short comment noting they are POSIX/C headers.

**Member layout: descending alignment, with size annotations.** Place 8-byte members first, then 4-byte, then 2-byte, then 1-byte. Annotate members with `// 8B`, `// 4B`, `// 2B`, `// 1B` comments so padding is visible at a glance. Constructor initializer lists must match declaration order. Audited gold-standard examples: `cbdg::Read`, `ReadEvidence`, `RawVariant`, `VariantCall` (the `assembly-and-calling-expert` subagent's body has the full list with paths).

**Logging and formatting.** Use spdlog macros (`SPDLOG_INFO`, `SPDLOG_DEBUG`, etc.) and fmtlib (`fmt::format`, `fmt::print`). `std::format` and `std::print` are forbidden; the at-write hook hard-blocks both.

**Assertions and namespaces.** Use `LANCET_ASSERT` (defined in `src/lancet/base/assert.h`); bare `assert()` is forbidden. `using namespace std` is forbidden in headers. Both are hook-enforced.

**Algorithms over manual loops.** Prefer `<algorithm>`, `<ranges>`, and Abseil utility headers when intent is clearer. Exception: when the algorithm form would hurt readability, complicate the code, or degrade performance, keep the loop. Use judgment — the goal is clarity, not dogma.

<important_if intent="adding a NOLINT suppression of any kind">
**NOLINT discipline.** Suppressions must be scoped: `// NOLINTNEXTLINE(check-name)` for one line, paired `NOLINTBEGIN(check-name)` / `NOLINTEND(check-name)` for blocks. Bare `// NOLINT` and `// NOLINT(check-name)` (the inline same-line forms) are forbidden — they pollute the statement line and make scannability hard. The `validate_naming.py` hook hard-blocks bare-NOLINT additions inside `src/lancet/`. The first answer to "should I suppress this?" is always "no — fix the underlying issue." Reach for a suppression only when the check is genuinely wrong about this case, the diagnostic cannot be fixed without harming clarity or correctness, and you can articulate why. When refactoring, do not carry suppressions forward; remove all and re-add only after a documented simplification attempt.

**Scoped vs bare NOLINT — the structural reason.** The bare inline form (`int x; // NOLINT(check)`) is fragile under clang-format: when clang-format wraps long lines, it can move the comment away from the statement it suppresses, silently invalidating the suppression and re-introducing the warning. The scoped forms (`NOLINTNEXTLINE` and `NOLINTBEGIN`/`NOLINTEND`) anchor to surrounding lines, not the wrapped statement, so they survive any future formatting pass. This is not stylistic preference — it is the only form that's robust to the project's auto-format policy.

**Rationale placement: always on lines ABOVE the NOLINT directive, regardless of length.** Every suppression carries a WHY-rationale, and that rationale is written on one or more `//` comment lines **immediately above** the `NOLINTNEXTLINE` or `NOLINTBEGIN` line. Do not put the rationale inline on the NOLINT line itself (no trailing `-- <reason>` clause). This rule applies to short rationales as well as long ones — uniform placement keeps the directive line scannable and keeps the WHY adjacent to the directive without crowding it. The canonical form is:

```cpp
// <WHY this suppression is justified>
// NOLINTNEXTLINE(check-name)
<line that needs the suppression>
```

For block suppressions:

```cpp
// <WHY this block needs to be silenced>
// NOLINTBEGIN(check-name)
<block>
// NOLINTEND(check-name)
```

`NOLINTEND` carries no rationale (it's a closing token; the WHY belongs on the matching `NOLINTBEGIN`). The `validate_naming.py` hook enforces this layout for new NOLINT additions inside `src/lancet/` — it hard-blocks any new `NOLINTNEXTLINE`/`NOLINTBEGIN` that lacks a `//` comment on the line(s) immediately preceding it.

**Disclosure rule for suppressions.** Whenever Claude adds a NOLINT suppression of any form (NOLINTNEXTLINE or NOLINTBEGIN/NOLINTEND), it must surface the addition explicitly in the summary of work — what was suppressed, where, why, and what was tried first. The disclosure is not optional and not for the user's review only; it is the mechanism that prevents the suppression from being a quiet workaround for a real violation. The `fresh-reviewer` subagent flags any undisclosed new suppression at PR time as a finding.
</important_if>

**Function complexity ceilings** (from `.clang-tidy`): cognitive complexity ≤ 35, statement count ≤ 150, line count ≤ 200, parameter count ≤ 6, nesting depth ≤ 4. The Rule of 5 is enforced.

**Integer literal separators.** Five or more digits use `'`: `1'048'576` (decimal, groups of 3), `0xFF'00'FF` (hex, groups of 2), `0b1111'0000` (binary, groups of 4).

## Code comments

Comments serve a developer with the source open trying to understand *what* the code does, *why* it was designed this way, and *what would break* if it were done differently. Explain the why, not the what. Include math with representative values. State invariants and boundary conditions with a `CRITICAL:` prefix when non-obvious. Use ASCII diagrams (wrapped in `// clang-format off` / `// clang-format on`) for spatial code. No filler words, no anthropomorphization, no unexplained jargon.

The full code-comment convention with worked examples and the filler-word catalog lives in `docs_dev/style/code_comments.md`. The `semantic-audit` skill encodes the methodology that catches comment violations across a subdirectory; use it after a refactor that renamed a metric or field, and on a quarterly cadence.

## Downstream sync requirements

The following pairs go stale silently. Grep the codebase when changing any of them.

**VCF FORMAT/INFO/FILTER field names.** The field is defined in `src/lancet/caller/variant_call.{h,cpp}` and the value is computed in `src/lancet/caller/variant_support.cpp`. The VCF header `##FORMAT` / `##INFO` / `##FILTER` line is registered in `src/lancet/cli/vcf_header_builder.cpp` (the canonical site for header definitions). The user-facing description is in `docs/guides/vcf_output.md` and, for FORMAT fields, also in `docs/guides/alignment_annotations.md`. A change to the field name, the formula, or the value range must update all four locations. **Always invoke the `vcf-validator` subagent before merging FORMAT/INFO/FILTER changes.**

**CLI flag names and defaults.** Defined in `src/lancet/cli/cli_interface.cpp`. Documented in `docs/reference.md`. Cross-linked from `docs/guides/*.md` for any flag with behavioral implications.

**Architecture and data-flow comments.** Each layer's `add_library(...)` block in `CMakeLists.txt` carries a top-of-block comment. When files move between layers, are added, or are renamed, update it.

**Update-When table for documentation.**

| Page | Path | Update When |
|:-----|:-----|:------------|
| Pipeline Architecture | `docs/guides/architecture.md` | Any core component, pipeline step, or global tuning parameter changes |
| Feature Guides | `docs/guides/*.md` | Any mathematical metric, algorithmic heuristic, or statistical model changes |
| VCF Output Reference | `docs/guides/vcf_output.md` | Any INFO or FORMAT field added, modified, or removed |
| Alignment Annotations | `docs/guides/alignment_annotations.md` | Any FORMAT field whose computation involves CIGAR or alignment metrics changes |
| CLI Reference | `docs/reference.md` | Any CLI flag added, modified, or removed |
| MkDocs Navigation | `mkdocs.yml` | Any documentation page added or removed |

## Agent memory

Each subagent has a memory file under `.claude/agent-memory/<agent-name>.md` that persists across sessions. Five of the six are project-scoped (git-tracked, shared with the team); `perf-analyst.md` is user-scoped (gitignored) because performance observations are hardware-specific. See `.claude/agent-memory/README.md` for the format spec, knowledge kinds, and the decision rationale.

Project-scoped memory updates bundle into the next feat/fix/perf/chore commit — Claude does not create standalone "update agent memory" commits. The `pre_commit_summary.sh` hook surfaces the count of bundled memory files in the commit summary so the user knows when memory is implicitly traveling with their work. The `/audit-bundle` Pairing 12 reviews each project-scoped file for accuracy quarterly.

Because memory updates always ride inside typed commits (`feat:`, `fix:`, `perf:`, `chore:`), they do not produce standalone CHANGELOG entries — git-chglog already filters on commit type, so an agent-memory-only commit (which doesn't exist by design) wouldn't appear, and a typed commit that incidentally touches memory appears under its primary type. No `.chglog/config.yml` change is needed.

## Workflow rules

Use plan mode (Shift-Tab twice) for any change touching multiple files or the algorithmic core (`cbdg/`, `caller/`).

Use one of the four conventional-commit prefixes the project's `.chglog/config.yml` filters on: `feat:` (new features), `fix:` (bug fixes), `perf:` (performance), `chore:` (everything else). Other prefixes parse but are silently dropped by chglog. The header pattern does NOT support scopes, so write `fix: ...` rather than `fix(caller): ...`. Body shape (trivial vs. substantive, two-section format with file-list bullets, 72-char wrap) is enforced by the `validate_commit_message` hook against `.claude/commit-style.json`. Use the `/commit` slash command rather than writing the message by hand.

Use a worktree (`git worktree add`) for any change above thirty minutes or spanning multiple layers.

Use the `fresh-reviewer` subagent before merging any change to `main`. The fresh-context review catches what the writer rationalized away.

## Where to reach for what

A few specific routes worth keeping in mind, because the natural default is wrong:

- **Tests.** When adding a unit or integration test, use the `add-cpp-test` skill. It encodes the `tests/<layer>/<file>_test.cpp` layout, the project's Catch2 v3.14.0 conventions, and a complete cheat sheet of Catch2 idioms (assertions, sections, generators, fixtures, matchers, logging, skipping) plus a flexible-invocation reference for the test binary (tag/name filtering, `--section`, `--generator-index`, reporters, sharding, seed reproduction). Reach for the skill any time the test workflow needs more than `pixi run test`.
- **Test data.** Before picking a fixture for a sanitizer run, profile, or `/e2e`, consult the `test-data-reference` skill — it has a decision tree mapping workflow to fixture (germline NA12878 / chr1; somatic HCC1395 / chr4; truth comparisons; T2T).
- **Profiling and optimization.** When the user mentions slow runs, hotspots, or wants to optimize, use the `profile-and-optimize` skill. It composes with the `perf-analyst` subagent for analysis.
- **Variant-discovery code (cbdg/, caller/).** Delegate review and deep correctness reasoning to the `assembly-and-calling-expert` subagent. It carries the bi-directed de Bruijn graph mental model and the layer-specific correctness concerns. For VCF schema specifics, that agent delegates to `vcf-validator`.
- **Substantive features.** Use `/spec` first, before writing code. It interviews you to surface scope, acceptance criteria, and constraints, then writes both a spec and a kickoff prompt to `notes/<FEATURE>/` for a fresh execution session. Skipping the interview on substantive work is a known failure mode.
- **Probe tracking / missed-variant forensics.** Use the `probe-tracking` skill for the operational mechanics (running `truth_concordance.py`, the Lancet2 probe invocation, `analyze_probe_results.py`) and the `probe-interpreter` subagent for interpreting an existing report. The three slash commands `/probe-concordance`, `/probe-run`, `/probe-analyze` are independent and composable; run them in sequence for end-to-end, run any one alone if its inputs already exist. The skill+agent split is intentional: the skill is the playbook, the agent is the consultant.
- **Substantive design decisions or postmortems.** Use `/arch-decision-record` for an Architecture Decision Record (a deliberate, single-decision document future contributors need to understand the *why*). Use `/investigation` for a postmortem, debug archaeology, or performance investigation (an immutable snapshot of a moment of debugging). Both are interview-driven and produce drafts under `docs_dev/architecture/` or `docs_dev/investigations/` respectively, following the writing guides in those directories' READMEs. The failure mode without these commands is writing the rationale inline in the chat instead of producing the persistent document.

`/check` mid-session validates without waiting for the Stop hook. `/e2e` confirms pipeline-level changes did not break variant calling on either profile.

## Scratch space

When you need to experiment, write throwaway scripts, or jot analysis notes, use `notes/scratch/`. The protected-paths hook explicitly does not block writes there, and the path is gitignored. Do NOT scatter experimental files in `src/`, `tests/`, or the project root.

## What you must not do

Do not WRITE to `pixi.toml`, `pixi.lock`, the `cmake/` superbuild modules, `Dockerfile`, `.github/`, `data/`, `CHANGELOG.md`, `LICENSE`, or any path under `cmake-build-*/`, `_deps/`, `.pixi/`, `.worktrees/`, or `.chglog/`. The at-write hook blocks these paths. The full canonical list lives in `.claude/protected_paths.txt`. Note the asymmetry: these paths are write-blocked, not read-blocked. Source under `cmake-build-*/_deps/` (vendored htslib, abseil, spdlog, SPOA, minimap2, etc.) remains freely readable and is often the authoritative answer for deep "why does this dependency behave this way?" questions — read it directly rather than guessing.

Do not rename or remove user-facing CLI flags without explicit user approval. The flag surface is documented and downstream-keyed.

<important_if intent="proposing or accepting an optimization">
Do not propose optimizations whose correctness implications you cannot articulate. Lancet2's priorities are correctness first, performance second.
</important_if>

Do not modify VCF FORMAT/INFO/FILTER definitions without invoking `vcf-validator`. Schema drift breaks downstream pipelines silently.
