# AGENTS.md — Lancet2 project memory

This file is the canonical project-memory document for any agentic coding tool used on Lancet2 (Claude Code, GitHub Copilot, Codex, Cursor, etc.). It carries the conventions, build system, architecture, and workflow rules that apply to all changes. Rules that apply only to a specific layer or workflow live in the corresponding subagent, skill, or path-scoped rule under `.claude/rules/`, not here. See `.claude/cost-model.md` for why this file is intentionally short.

A thin wrapper at `CLAUDE.md` imports this file (`@AGENTS.md`) so Claude Code reads it through the conventional CLAUDE.md path; the canonical content lives here for cross-tool portability.

## Project at a glance

Lancet2 is a modern C++20 variant caller (SNVs and InDels) by the New York Genome Center. It performs joint multi-sample localized colored de Bruijn graph assembly. The CLI entry point is `Lancet2 pipeline`. The `main` branch is the source of truth; the documentation site at `nygenome.github.io/Lancet2` typically lags `main`.

## Build and tooling

The build is driven by [pixi](https://pixi.sh). Always prefer pixi tasks over hand-rolled CMake invocations so you pick up the locked toolchain (gcc, clang, cmake, ninja, clang-tidy, IWYU, go) rather than whatever is on `PATH`. Run `pixi task list` for the full catalog.

Sanitizer trees (`cmake-build-asan/`, `cmake-build-tsan/`, `cmake-build-msan/`, `cmake-build-ubsan/`) are configured by the `sanitizer-build-analysis` skill and the matching pixi tasks; they require disabling the static mimalloc link via the `LANCET_SANITIZE_BUILD` source guard. The profiling build (`pixi run configure-profile` / `pixi run build-profile`) is documented in the `profile-and-optimize` skill; reach for the skill rather than configuring profiling builds inline.

The `test` task uses the Release tree (`cmake-build-release/tests/TestLancet2`) because `lint-check` and `iwyu-check` already build that tree at commit time, and `tests/CMakeLists.txt` defines `LANCET_DEBUG_MODE` on the test target unconditionally so `LANCET_ASSERT` keeps firing under Release. The `build-debug` tree exists for explicit debug iteration but is not on the default flow.

`/fix-and-validate` is the canonical validation entry point. It runs `pixi run iwyu-fix` (mutating; auto-formats and fixes includes), then read-only `lint-check` and `test` against the Release tree. On full success it writes a key=value record to `.claude/cache/validation-state.txt` (gitignored) with four fields: `stash_hash` (a `git stash create` hash of the validated working tree), `head_sha`, `timestamp`, and `status` (pass/fail). The `pre_commit_gate.sh` PreToolUse hook reads `stash_hash` at `git commit` time and silently passes if it matches the about-to-be-committed working tree (~50ms). The `statusline.sh` Claude-Code statusline reads `status` and surfaces it as `[validation-check: pass|fail|stale]` on every prompt. If the marker is missing or stale, the gate blocks with "run /fix-and-validate first". Run `/fix-and-validate` whenever you've made non-trivial source changes; commits without a fresh validation will be blocked.

For clang-tidy violations, NOLINT discipline, and the project-specific reason `clang-tidy --fix` is forbidden, see the `clang-tidy-discipline` skill — it carries the procedural how-to. The canonical rule statements live in `docs_dev/style/cpp_style.md`.

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

A file in layer N may depend on layers 1 through N. It must NOT depend on layers above it. The at-write hook in `.claude/hooks/validate_layer_direction.py` enforces this on every edit; a companion hook `.claude/hooks/validate_cpp_identifiers.py` enforces naming conventions and NOLINT-suppression discipline (see the clang-tidy-discipline skill for details).

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

`docs_dev/style/cpp_style.md` is the source of truth for naming, member layout, include style, NOLINT discipline, complexity ceilings, logging/formatting (spdlog + fmtlib; `std::format` and `std::print` are forbidden), and assertions/namespaces (`LANCET_ASSERT` only; bare `assert()` and `using namespace std` in headers are forbidden). `.clang-tidy` is the source of truth for the auto-enforced rules.

**Function complexity ceilings** (from `.clang-tidy`): cognitive complexity ≤ 35, statement count ≤ 150, line count ≤ 200, parameter count ≤ 6, nesting depth ≤ 4. The Rule of 5 is enforced.

**Integer literal separators.** Five or more digits use `'`: `1'048'576` (decimal, groups of 3), `0xFF'00'FF` (hex, groups of 2), `0b1111'0000` (binary, groups of 4).

## Code comments

Comments serve a developer with the source open trying to understand *what* the code does, *why* it was designed this way, and *what would break* if it were done differently. Explain the why, not the what. Include math with representative values. State invariants and boundary conditions with a `CRITICAL:` prefix when non-obvious. Use ASCII diagrams (wrapped in `// clang-format off` / `// clang-format on`) for spatial code. No filler words, no anthropomorphization, no unexplained jargon.

The full code-comment convention with worked examples and the filler-word catalog lives in `docs_dev/style/code_comments.md`. The `doc-sync` skill encodes the methodology that catches comment violations across a subdirectory; use it after a refactor that renamed a metric or field, and on a quarterly cadence. Use Mode A for targeted post-change sync; Mode B for periodic semantic audits.

## Downstream sync requirements

The following pairs go stale silently. Grep the codebase when changing any of them.

**VCF FORMAT/INFO/FILTER field names.** The field is defined in `src/lancet/caller/variant_call.{h,cpp}` and the value is computed in `src/lancet/caller/variant_support.cpp`. The VCF header `##FORMAT` / `##INFO` / `##FILTER` line is registered in `src/lancet/cli/vcf_header_builder.cpp` (the canonical site for header definitions). The user-facing description is in `docs/guides/vcf_output.md` and, for FORMAT fields, also in `docs/guides/alignment_annotations.md`. A change to the field name, the formula, or the value range must update all four locations. **Always invoke the `vcf-validator` subagent before merging FORMAT/INFO/FILTER changes.** For any rename, cardinality change, removal, or silent-semantic-change, use the `external-interface-changes` skill.

**CLI flag names and defaults.** Defined in `src/lancet/cli/cli_interface.cpp`. Documented in `docs/reference.md`. Cross-linked from `docs/guides/*.md` for any flag with behavioral implications. CLI flag changes also go through `external-interface-changes`.

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

Use a worktree (`git worktree add`) for any change above thirty minutes or spanning multiple layers. The statusline shows `wt:<name>` when CWD is a linked worktree.

Use the `fresh-reviewer` subagent before merging any change to `main`. The fresh-context review catches what the writer rationalized away. When receiving its findings, verify before implementing, push back with technical reasoning when warranted, and implement in severity order — see `fresh-reviewer.md` § "Receiving review findings" for the full discipline.

When the user says "grill me," "push back hard," or otherwise asks for adversarial review, the responding session should hold its assumptions to a higher bar: question the framing, demand specific evidence for claims, and refuse to validate vague reasoning. The grill-me posture is opt-in by the user and applies to the immediate exchange, not the rest of the session.

When the user asks a binary or near-binary question and a single-syllable answer is the truth, give that answer. The "caveman acknowledgment" — "Yes." / "No." / "Don't know" / "Maybe" — is appropriate when the user is checking a fact, not when they want analysis. Wrapping a single-syllable answer in three sentences of qualification dilutes the signal.

For substantive features, use `/spec` first, before writing code. It interviews you to surface scope, acceptance criteria, and constraints, then writes a spec at `notes/<topic>/<YYYY-MM-DD>-spec.md` for a fresh execution session via `/execute-spec`. For exploratory work where the approach is unclear, use `/brainstorm` first to generate options, then `/spec` against the chosen option. Skipping the interview on substantive work is a known failure mode.

The `notes/<topic>/` directory naming convention applies to all interview-driven artifacts: `<YYYY-MM-DD>-spec.md` for `/spec`, `<YYYY-MM-DD>-options.md` for `/brainstorm`, `<YYYY-MM-DD>-investigation.md` for `/investigate`. The date prefix preserves history when the same topic is revisited later.

## Where to reach for what

A few specific routes worth keeping in mind, because the natural default is wrong:

- **Tests.** When adding a unit or integration test, use the `add-cpp-test` skill. It encodes the `tests/<layer>/<file>_test.cpp` layout, the project's Catch2 v3.14.0 conventions, and a complete cheat sheet of Catch2 idioms (assertions, sections, generators, fixtures, matchers, logging, skipping) plus a flexible-invocation reference for the test binary (tag/name filtering, `--section`, `--generator-index`, reporters, sharding, seed reproduction). Reach for the skill any time the test workflow needs more than `pixi run test`.
- **Test data.** Before picking a fixture for a sanitizer run, profile, or `/e2e-pipeline-test`, consult the `test-data-locations` skill — it has a decision tree mapping workflow to fixture (germline NA12878 / chr1; somatic HCC1395 / chr4; truth comparisons; T2T).
- **Profiling and optimization.** When the user mentions slow runs, hotspots, or wants to optimize, use the `profile-and-optimize` skill. It composes with the `perf-analyst` subagent for analysis.
- **Variant-discovery code (cbdg/, caller/).** Delegate review and deep correctness reasoning to the `assembly-and-calling-expert` subagent. It carries the bi-directed de Bruijn graph mental model and the layer-specific correctness concerns. For VCF schema specifics, that agent delegates to `vcf-validator`.
- **Substantive features.** Use `/spec` first, before writing code. Skipping the interview on substantive work is a known failure mode.
- **Open problem space, no chosen approach.** Use `/brainstorm` to generate options before committing to one.
- **Probe tracking / missed-variant forensics.** Use the `probe-tracking` skill for the operational mechanics (running `truth_concordance.py`, the Lancet2 probe invocation, `analyze_probe_results.py`) and the `probe-interpreter` subagent for interpreting an existing report. The three slash commands `/probe-concordance`, `/probe-run`, `/probe-analyze` are independent and composable; run them in sequence for end-to-end, run any one alone if its inputs already exist. The skill+agent split is intentional: the skill is the playbook, the agent is the consultant.
- **Substantive design decisions or postmortems.** Use `/arch-decision-record` for an Architecture Decision Record (a deliberate, single-decision document future contributors need to understand the *why*). Use `/investigate` for a postmortem, debug archaeology, or performance investigation (an immutable snapshot of a moment of debugging). Both are interview-driven and produce drafts under `docs_dev/architecture/` or `docs_dev/investigations/` respectively, following the writing guides in those directories' READMEs. The failure mode without these commands is writing the rationale inline in the chat instead of producing the persistent document.

`/fix-and-validate` mid-session validates without waiting for the Stop hook. `/e2e-pipeline-test` confirms pipeline-level changes did not break variant calling on either profile. `/wrap-branch` finalizes the branch (merge / PR / keep / discard) when work is done.

## Scratch space

When you need to experiment, write throwaway scripts, or jot analysis notes, use `notes/scratch/`. The protected-paths hook explicitly does not block writes there, and the path is gitignored. Do NOT scatter experimental files in `src/`, `tests/`, or the project root.

## What you must not do

Do not WRITE to `pixi.toml`, `pixi.lock`, the `cmake/` superbuild modules, `Dockerfile`, `.github/`, `data/`, `CHANGELOG.md`, `LICENSE`, or any path under `cmake-build-*/`, `_deps/`, `.pixi/`, `.worktrees/`, or `.chglog/`. The at-write hook blocks these paths. The full canonical list lives in `.claude/protected_paths.txt`. Note the asymmetry: these paths are write-blocked, not read-blocked. Source under `cmake-build-*/_deps/` (vendored htslib, abseil, spdlog, SPOA, minimap2, etc.) remains freely readable and is often the authoritative answer for deep "why does this dependency behave this way?" questions — read it directly rather than guessing.

Do not rename or remove user-facing CLI flags without explicit user approval. The flag surface is documented and downstream-keyed. Use `external-interface-changes` for any operation on the CLI flag surface.

<important_if intent="proposing or accepting an optimization">
Do not propose optimizations whose correctness implications you cannot articulate. Lancet2's priorities are correctness first, performance second.
</important_if>

Do not modify VCF FORMAT/INFO/FILTER definitions without invoking `vcf-validator`. Schema drift breaks downstream pipelines silently. Use `external-interface-changes` for any operation on the VCF schema.
