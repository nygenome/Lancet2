---
name: python-script-sync
description: Use when changing a Python script under scripts/, when changing C++ code that a Python script consumes, or when figuring out which Python script to reach for. Trigger on "I just changed analyze_profile.py", "the C++ probe-results format changed; does analyze_probe_results.py need an update?", "what does run_clang_tidy.py do?", "which script generates the walk palette?", "the export_docs.py is producing wrong output". Maintains the contract between the C++ pipeline and the Python tooling that consumes its output, runs lint-style checks on script edits, and surfaces drift between scripts and the project conventions they implement (clang-tidy invocation, format checks, IWYU). Scoped to scripts/*.py and scripts/*.sh; does NOT cover Python under tests/ (those are C++ test fixtures, not pipeline tooling) or Python under docs/ (those are mkdocs plugins).
allowed-tools: Read, Glob, Grep, Edit, Write, Bash
---

# Python script sync on Lancet2

The `scripts/` directory holds 15 files (12 Python, 3 shell, 1 Jinja template) that fall into three roles:

**Pipeline-output analysis** — these consume artifacts produced by the C++ Lancet2 binary or its profiling output. Drift here is silent: the C++ side changes a TSV column, the Python side keeps reading by index, and the analysis is wrong without anyone noticing.

- `analyze_debug_run.py` — reads pipeline debug logs and surfaces structured findings.
- `analyze_probe_results.py` — the 27-level cascade attribution engine for probe tracking. Reads `probe_results.tsv`, `missed_variants.txt`, optionally `concordance_details.txt` and `lancet2_debug.log`.
- `analyze_profile.py` — gperftools profile analysis wrapper. Reads `*.bin` profile files via pprof.
- `gen_walk_palette.py` — generates color palettes for graph visualization output.
- `truth_concordance.py` — compares Lancet2 VCF against truth VCFs (GIAB, Manta) and produces `missed_variants.txt`.
- `profile_report.html.j2` — Jinja2 template consumed by `analyze_profile.py` for HTML report rendering.

**Build/release tooling** — these wrap or extend the project's build system.

- `build_conda_local.sh` — builds the local conda package.
- `build_push_image.sh` — builds and pushes the Docker image.
- `bump_version.py` — increments the version in `CMakeLists.txt`, `pixi.toml`, and elsewhere.
- `export_docs.py` — concatenates all docs into a single Markdown file (consumed by the `docs-export` pixi task).
- `prep_stm_viz.sh` — preps state-machine visualization output.
- `update_changelog.sh` — regenerates `CHANGELOG.md` via git-chglog.

**Lint and IWYU enforcement** — these implement the project's check/fix tasks.

- `run_clang_format.py` — clang-format wrapper used by `pixi run fmt-check` / `fmt-fix`.
- `run_clang_tidy.py` — clang-tidy wrapper used by `pixi run lint-check`. Auto-fix mode (`--fix`) is intentionally NOT supported in this project; clang-tidy auto-fix has historically broken compilation, so violations are resolved by editing source manually. See AGENTS.md "invoking clang-tidy" callout.
- `run_iwyu.py` — IWYU wrapper used by `pixi run iwyu-check` / `iwyu-fix`.

This skill walks the sync between these scripts and the rest of the codebase.

## Step 1 — Identify what the script consumes or produces

Before editing a script, identify its inputs and outputs from the source itself. Read the docstring at the top, then walk the argument parser and the I/O calls. Don't guess from the filename — `analyze_debug_run.py` and `analyze_profile.py` sound similar but consume entirely different artifacts.

For pipeline-output analysis scripts specifically, find the C++ code that produces the consumed artifact:

- `analyze_probe_results.py` consumes `probe_results.tsv` produced by `src/lancet/caller/probe_diagnostics.cpp` (or whichever file owns the probe output emission).
- `analyze_profile.py` consumes `*.bin` files produced by gperftools, with environment baked in `src/lancet/cli/pipeline_runner.cpp`.
- `truth_concordance.py` consumes Lancet2 VCFs and produces `missed_variants.txt` consumed in turn by the C++ `--probe-variants` flag (a round-trip).

The C++ side is the source of truth for the data format. If the script is wrong, fix the script. If the C++ format genuinely needs to change, the change goes through `external-interface-changes` (for VCF) or normal review (for TSV/log formats), and the script update is part of the same commit.

## Step 2 — When the C++ side changes, find the affected scripts

If the user just changed C++ code that emits a TSV column, a log line, or a VCF field, scan `scripts/` for consumers:

```bash
grep -l "<column_name>" scripts/*.py
```

For schema-affected fields, the search includes all the places the field name might appear: as a TSV header, as a literal string in a parser, as a Python dict key. Be thorough — partial coverage is worse than none because it gives the appearance of having checked.

For a VCF field change, `truth_concordance.py` and `analyze_probe_results.py` are the primary consumers. For a log line change, `analyze_debug_run.py` is usually the only consumer.

## Step 3 — When the script is the source of truth, find the C++ consumers

The lint/IWYU/format scripts are the project's enforcement mechanism. If `run_clang_tidy.py` changes its CLI surface, the pixi task `lint-check` that invokes it needs to update too. Same for `run_clang_format.py` (`fmt-check`/`fmt-fix`) and `run_iwyu.py` (`iwyu-check`/`iwyu-fix`).

If `export_docs.py` changes its output format, the `docs-export` pixi task and any downstream tooling that consumes the exported Markdown need to update.

For the build/release scripts (`bump_version.py`, `update_changelog.sh`, `build_conda_local.sh`, `build_push_image.sh`), the consumers are the release workflow itself — the user invoking them, plus possibly CI. Changes there are higher-impact and should be discussed before editing.

## Step 4 — Apply the change with surrounding awareness

Edit the script. Run a smoke test if there's an obvious one available:

- For `run_clang_format.py` / `run_clang_tidy.py` / `run_iwyu.py`: invoke against a small input (`pixi run fmt-check` on a single file) to confirm the wrapper still works.
- For `analyze_*.py`: have a known-good input artifact and confirm the script produces a sensible report against it. If you don't have one, either run the pipeline against `LANCET_TEST_*_REGION_SMALL` to generate one, or skip the smoke test and rely on review.
- For `bump_version.py` / `update_changelog.sh`: the dry-run path is `git diff --stat` after running; revert with `git restore` if anything looks wrong. Do not push until the user has confirmed.

## Step 5 — Update documentation if the script's user contract changed

If the change altered:

- The script's CLI surface (new flag, renamed flag, changed default)
- The script's output format (new column, renamed column, different file extension)
- The script's input expectations (new required argument, changed file format)

Then the documentation that describes the script needs updating too. The `doc-sync` skill handles the procedural part. Specifically:

- Profile-related script changes → `notes/scratch/profiling-runs/` and the `profile-and-optimize` skill.
- Probe-related script changes → `.Codex/skills/probe-tracking/SKILL.md` and the `probe-tracking` skill.
- Lint/format script changes → likely no user-facing docs (they're invoked through pixi tasks); update the pixi task descriptions in `pixi.toml` if behavior changed.

If the change was internal-only (refactor, performance, comment edits), no doc update needed.

## Step 6 — Commit with the right type

The Lancet2 chglog filter accepts only `feat`, `fix`, `perf`, and `chore` as commit types; `refactor:`, `docs:`, `build:`, etc. parse but are silently dropped from `CHANGELOG.md`. The chglog header pattern also does not support scopes — `fix(scripts): ...` parses as if the type were the literal string `fix(scripts)`, which is not in the filter list and gets dropped. Always write the type with no scope and mention the script name in the subject text or body instead:

- Behavior-changing → `feat: ...` or `fix: ...` naming the affected script in the subject (e.g., `fix: handle missing concordance_details in analyze_probe_results.py`).
- Refactor, comment, or docstring-only changes → `chore:` (these are internal-facing and `chore` is the catch-all the filter recognizes).

Mention the consumer side in the body if there's a coordinated update:

```
fix: handle missing concordance_details.txt in analyze_probe_results.py

Previously raised KeyError when --concordance-details was omitted.
Now skips the §1 scorecard enrichment cleanly and logs a warning.

C++ side unchanged; this is purely a Python consumer fix.
```

## When NOT to use this skill

Do not use this skill for:

- Python under `tests/` — those are not pipeline tooling; they're test fixtures or test harnesses for C++ code.
- Python in mkdocs config (`docs/overrides/`, mkdocs plugins) — those are documentation infrastructure, handled by `doc-sync` if relevant.
- One-off ad-hoc analysis scripts in `notes/scratch/` — those are throwaway; no contract to maintain.

## When a new script should be added

If the user proposes adding a new file under `scripts/`, walk through:

1. **Does an existing script already do this, or could be extended?** A new script is a maintenance commitment; first check `analyze_debug_run.py`, `analyze_profile.py`, `analyze_probe_results.py` — these are the analysis-tool home.
2. **Does it need a pixi task?** If the script is meant to be invoked routinely, define a pixi task for it (in `pixi.toml`). One-shot scripts without a pixi task should live in `notes/scratch/` instead.
3. **Does it consume a C++-produced artifact?** Document the contract in the script's module docstring, citing the C++ file that owns the format.

A new script with a clear long-term role and a documented C++-side contract is fine. A new script that's "I just need this for one thing" goes in `notes/scratch/`.
