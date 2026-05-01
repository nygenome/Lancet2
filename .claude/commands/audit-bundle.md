# /audit-bundle — Comprehensive drift audit against project source

Walk every claim the bundle makes about the Lancet2 source code, compare against the actual source, and report drift. This is the source-of-truth refresh ritual that complements `/sync-cost-model` (which refreshes against Anthropic's external docs); this command refreshes against the project's own source files.

## When to run this

After a substantive refactor (renamed structs, moved files, restructured layers). After a release that changes the FORMAT/INFO/FILTER schema or adds CLI flags. On a quarterly cadence as a maintenance sweep, paired with `/sync-cost-model`. When something in the bundle has felt stale and you want a systematic check rather than an ad-hoc one.

For a focused check after a VCF schema change specifically, use `/audit-vcf-schema` instead — it's faster and the report is more targeted.

## Procedure

The audit walks twelve source-of-truth pairings. Run them all, collect findings, then produce a single structured report at the end. Do not fix anything inline; the report is the deliverable, and the user decides what to fix afterward.

For each pairing, the procedure is the same: read the source, read the bundle's claims, diff them, record findings with `file:line` references on both sides.

### Pairing 1 — Pixi tasks

**Source of truth:** `pixi.toml` `[tasks]` block.

Read `pixi.toml` and extract every task name. The names that should exist (from the project's history plus the bundle's integration additions) are `configure-release`, `build-release`, `configure-debug`, `build-debug`, `configure-profile`, `build-profile`, `test`, `bench`, `lint-all`, `lint-check`, `fmt-check`, `fmt-fix`, `iwyu-check`, `iwyu-fix`, `docs-build`, `docs-serve`, plus the 12 sanitizer tasks the bundle adds (`configure-asan`, `build-asan`, `test-asan` and the parallel triples for `msan`, `tsan`, `ubsan`). The `lint-fix` task was deliberately removed (clang-tidy auto-fix is forbidden in this project; see AGENTS.md "invoking clang-tidy" callout) — its absence is correct, its presence is a finding. If any other task has been renamed, added, or removed, that is a finding.

**Bundle locations to check:** `AGENTS.md` (build-and-tooling section), `.claude/commands/check.md`, `.claude/commands/e2e.md`, `.claude/commands/README.md`, `.claude/hooks/pre_commit_gate.sh`, `.claude/hooks/stop_reminder.sh`, `.claude/skills/profile-and-optimize/SKILL.md`, `.claude/skills/sanitizer-triage/SKILL.md`, `.claude/agents/perf-analyst.md`. Grep for `pixi run` in each.

Report any task name in the bundle that does not appear in `pixi.toml`.

### Pairing 2 — Six-layer dependency chain

**Source of truth:** `CMakeLists.txt` `target_link_libraries` declarations.

Parse the layer chain from `CMakeLists.txt`. The expected chain is `lancet_base` → `lancet_hts` → `lancet_cbdg` → `lancet_caller` → `lancet_core` → `lancet_cli`, with `Lancet2` linking `lancet_cli` and `mimalloc`. If a new layer has been added, an existing layer renamed, or the dependency chain reordered, that is a finding.

**Bundle locations to check:** `AGENTS.md` (architecture section, the box diagram), `.claude/hooks/validate_layer_direction.py` (the layer-order list inside the script), `.claude/agents/assembly-and-calling-expert.md` (Dependencies section), `.claude/agents/README.md` (any references to the chain).

Report any layer name or dependency that does not match `CMakeLists.txt`.

### Pairing 3 — CLI flag surface

**Source of truth:** `src/lancet/cli/cli_interface.cpp`.

Grep `cli_interface.cpp` for `AddOpt` calls. Each call defines a flag with name, target, and required-ness. The bundle makes claims about specific flags being present and their required-ness; verify each claim.

Specific things to verify:
- `--tumor` and `--normal` both exist; `--normal` is `required=true` (third arg to `AddOpt`).
- `--reference`, `--region`, `--num-threads`, `--out-vcfgz` all exist.
- The `--probe-variants` and `--probe-results` pair, if mentioned in any docs, is mutually-required.
- Any flag mentioned in `/e2e` (in the help-flag verification loop) actually exists.

**Bundle locations to check:** `AGENTS.md` (build-and-tooling section), `.claude/commands/check.md`, `.claude/commands/e2e.md`, `.claude/commands/README.md`, `.claude/hooks/pre_commit_gate.sh`, `.claude/hooks/stop_reminder.sh`, `.claude/skills/profile-and-optimize/SKILL.md`, `.claude/skills/sanitizer-triage/SKILL.md`, `.claude/agents/perf-analyst.md`. Grep for `pixi run` in each.

Report any flag the bundle relies on that is no longer present, plus any flag whose required-ness has changed in a way that affects the bundle's invocations.

### Pairing 4 — Test-data fixtures

**Source of truth:** `data/download_test_data.sh` (which lists what gets downloaded) plus `gs://lancet2-test-datasets/test_harness_data/` (the upstream bucket). For practical purposes, only the script needs to be read; the bucket can be sampled if the script is unclear.

Read `download_test_data.sh` and extract the actual file names and the destination directory. The current convention is files land in `data/` (not `data/test_harness_data/`).

**Bundle locations to check:** `.claude/settings.local.json.example` (every path), `.claude/skills/test-data-reference/SKILL.md` (decision tree, env-var table, "what's in the dataset" section), `.claude/commands/e2e.md` (download instruction line), `AGENTS.md` (any path references).

Report any path or filename that has drifted from what `download_test_data.sh` actually produces.

### Pairing 5 — VCF FORMAT/INFO/FILTER schema

**Source of truth:** `src/lancet/cli/vcf_header_builder.cpp`.

Grep `vcf_header_builder.cpp` for `##FORMAT=`, `##INFO=`, `##FILTER=` lines. Each line defines a field with `ID`, `Number`, `Type`, and `Description`.

For brevity in the comprehensive audit, only flag schema changes that affect bundle claims (i.e., specific field names mentioned in agents, skills, or AGENTS.md). Do not enumerate every field. For a full schema walk, the user should run `/audit-vcf-schema` instead.

**Bundle locations to check:** `.claude/agents/vcf-validator.md` (the "Schema invariants you must enforce" section, which is the canonical home for the project's invariants list and which other agents and skills delegate to), `AGENTS.md` downstream-sync section. `.claude/agents/assembly-and-calling-expert.md` should mention VCF schema only as a delegation pointer to `vcf-validator`; if that agent has accumulated its own invariants list, that is a finding (duplication risk).

Report any specific field name the bundle mentions that is no longer in the header (or whose Number/Type has changed in a way the bundle's narrative depends on).

### Pairing 6 — Audited gold-standard structs

**Source of truth:** the actual struct definitions in source.

The bundle calls out specific structs as "gold-standard" memory layout examples: `cbdg::Read` (`src/lancet/cbdg/read.h`), `ReadEvidence` (`src/lancet/caller/variant_support.h`), `ReadAlleleAssignment` and `Mm2AlnResult` (`src/lancet/caller/genotyper.h`), `RawVariant` (`src/lancet/caller/raw_variant.h`), `VariantCall` (`src/lancet/caller/variant_call.h`).

For each struct, verify: the file path is correct (the file still exists at that path); the struct still exists in the file (grep for `struct <Name>` or `class <Name>`); the struct still has the descending-alignment + size-annotation pattern (a quick read of the struct body should show `// 8B`, `// 4B`, etc. comments).

**Bundle locations to check:** `.claude/agents/assembly-and-calling-expert.md` (gold-standard structs table), `AGENTS.md` (member-layout convention paragraph).

Report any struct that has moved, been renamed, or no longer follows the audited pattern.

### Pairing 7 — Code-style ceilings

**Source of truth:** `.clang-tidy` configuration values.

Read `.clang-tidy` and extract the configured complexity ceilings: `cognitive-complexity`, `statement-count`, `line-count`, `parameter-count`, `nesting-depth`. The bundle's AGENTS.md states specific values (35, 150, 200, 6, 4 are the current values).

**Bundle locations to check:** `AGENTS.md` (function-complexity-ceilings paragraph), `.claude/hooks/validate_naming.py` (if any of the rules it checks have moved into clang-tidy or vice versa).

Report any ceiling value that has changed in `.clang-tidy` but not in `AGENTS.md`.

### Pairing 8 — Protected paths source-of-truth sync

**Source of truth:** `.claude/protected_paths.txt`.

This file is the canonical list of paths the agent must not edit without explicit user override. It is consumed by two surfaces that must stay in sync with it:

1. **`.claude/hooks/block_protected_paths.py`** — reads the file directly at hook execution time. The hook should NOT carry a hardcoded list of its own; if it does, that is a finding (regression to the duplicated-source posture this consolidation eliminated).

2. **`.claude/settings.json` `permissions.deny`** — declarative deny rules that Claude Code evaluates alongside the hook. Per Claude Code's permission model, `Edit(...)` rules apply only to file-writing tools (Edit/Write/MultiEdit), not to Read/Glob/Grep. Each pattern in `protected_paths.txt` (except `.git/`, which Claude Code special-cases) should have a corresponding `Edit(<glob>)` rule in `permissions.deny`. Glob conventions: directory patterns end in `/**` (e.g., `cmake-build-*/**`), single-file patterns are bare (e.g., `LICENSE`, `pixi.lock`). The `.git/` pattern is handled by Claude Code itself and needs no deny rule.

The audit also verifies the read/write asymmetry is preserved: `protected_paths.txt` should NOT be referenced by any `Read(...)` deny rule in `settings.json`. The intent is that vendored dependency source under `cmake-build-*/_deps/` remains freely readable so the agent can answer deep questions by reading dependency source directly. A `Read(...)` deny rule covering any of these paths would defeat that and is a finding.

For each pattern in `protected_paths.txt`, verify a matching `Edit(...)` deny rule exists in `settings.json` (allowing for trailing `/` and `**` glob differences). For each `Edit(...)` deny rule in `settings.json`, verify the corresponding pattern exists in `protected_paths.txt`. Report any mismatch.

### Pairing 9 — Bundle-internal counts

**Source of truth:** the bundle's actual file counts.

The bundle's READMEs and cost-model document make claims about how many subagents, skills, slash commands, and hooks the bundle has. These claims drift silently every time something is added or removed. Compute the current counts:

```bash
ls .claude/agents/ | grep -v README | grep '\.md$' | wc -l       # subagent count
find .claude/skills -maxdepth 1 -mindepth 1 -type d | wc -l       # skill count (folders)
ls .claude/commands/*.md | grep -v README | wc -l                # slash command count
ls .claude/hooks/*.{py,sh} 2>/dev/null | wc -l                   # hook count
wc -l CLAUDE.md AGENTS.md 2>/dev/null                            # CLAUDE.md / AGENTS.md size
```

**Bundle locations to check** for stale counts:

- `README.md` Layout block (the comment after each `.claude/<dir>/` line) — names a count
- `README.md` Design philosophy section — should mention the actual command count
- `.claude/cost-model.md` "How the bundle's choices map to this" — names all four counts
- `.claude/agents/README.md`, `.claude/skills/README.md`, `.claude/commands/README.md`, `.claude/hooks/README.md` — each may name a per-directory count

For each location, compare the named count to the computed count. Report any mismatch.

### Pairing 10 — Path-scoped rules layer alignment

**Source of truth:** `src/lancet/<layer>/` directory listing in the project + `validate_layer_direction.py`'s layer-order list.

The bundle ships path-scoped rule files at `.claude/rules/<layer>.md` whose YAML frontmatter `paths:` glob targets `src/lancet/<layer>/**`. These rules carry layer-specific design logic that auto-loads when Claude edits files in the matching layer. The expected layer set is `base, hts, cbdg, caller, core, cli` — same as the dependency chain in pairing 2.

For each rule file in `.claude/rules/` (excluding `README.md`):

1. The `paths:` glob in frontmatter must match an existing `src/lancet/<layer>/` directory in the project. Stale rules pointing at a removed/renamed layer are a finding.
2. The layer name must appear in `validate_layer_direction.py`'s layer list. A rule that names a layer the hook doesn't enforce is a finding.

For each layer in `validate_layer_direction.py`'s layer list:

3. There should be a corresponding `.claude/rules/<layer>.md`. A new layer added without a rule file is a finding (the layer's design context is missing).

4. The `.claude/rules/README.md` table should list the rule file. Drift between the table and the actual file set is a finding.

**Bundle locations to check:** `.claude/rules/*.md` (frontmatter `paths:` field), `.claude/rules/README.md` (the "What lives here" table), `.claude/hooks/validate_layer_direction.py` (the layer list), `AGENTS.md` (the architecture-section bullet list pointing into `.claude/rules/`).

### Pairing 11 — Rule source-claim verification

**Source of truth:** the named source-code constants and class names cited inside `.claude/rules/*.md` rule body content.

The path-scoped rules under `.claude/rules/` are intentionally self-contained — they don't point at source files for rationale; they restate the rationale inline. This eliminates one drift mode (broken pointers to source) but introduces another: the inline restatement can drift away from the source over time. A rule that says "`MSA_MATCH_SCORE` is 0 to fix the SPOA tie-breaking bug" stays correct only as long as `msa_builder.h` actually defines `MSA_MATCH_SCORE = 0`.

This pairing grep-checks each named constant and class in the rules against the source. It fires only on names that would matter — specific numeric thresholds, named struct fields, magic constants — not on every word in the rules. The audit script grep-searches for each cited name in the source tree and verifies presence and (where relevant) value.

Specific checks the audit performs, each grounded in what the rules actually claim:

1. In `caller.md`, the SPOA scoring constants block: `MSA_MATCH_SCORE`, `MSA_MISMATCH_SCORE`, `MSA_OPEN1_SCORE`, `MSA_EXTEND1_SCORE`, `MSA_OPEN2_SCORE`, `MSA_EXTEND2_SCORE`. Each must appear in `src/lancet/caller/msa_builder.h` with the value the rule cites. A change to any value in the source without a corresponding rule update is a finding; a rule citing a constant that no longer exists in the source is a finding.

2. In `caller.md`, the Dirichlet-Multinomial constants: `DM_BACKGROUND_ERROR`, `DM_OVERDISPERSION`, `DM_ALPHA_FLOOR`. Each must appear in `src/lancet/caller/genotype_likelihood.cpp` with the cited values.

3. In `core.md`, the variant-store and pipeline-executor constants the rule names: `NUM_SHARDS` (must appear in `src/lancet/core/variant_store.h`), `NUM_BUFFER_WINDOWS` (must appear in `src/lancet/core/pipeline_executor.cpp` or `pipeline_executor.h`), `BATCH_SIZE` (must appear in `src/lancet/core/window_builder.h`), `BATCH_THRESHOLD` (must appear in `src/lancet/core/pipeline_executor.cpp`), the window-builder defaults (`DEFAULT_WINDOW_LENGTH`, `DEFAULT_PCT_OVERLAP`, `DEFAULT_REGION_PADDING` — all in `src/lancet/core/window_builder.h`), and the read-collector default (`DEFAULT_MAX_WINDOW_COVERAGE` in `src/lancet/core/read_collector.h`, where it actually lives — `core.md` cites it as `ReadCollector::DEFAULT_MAX_WINDOW_COVERAGE`).

4. In `cli.md`, the VCF-header constants the rule names: `FORMAT_STR_HEADER` and `CASE_CTRL_INFO_HDR_LINES` (both must appear in `src/lancet/cli/vcf_header_builder.cpp`), the placeholder names (`CONTIG_HDR_LINES`, `CONDITIONAL_INFO_LINES`, `ANNOTATION_INFO_LINES`, `RUN_TIMESTAMP`, `FULL_VERSION_TAG`, `FULL_COMMAND_USED`, `REFERENCE_PATH`), the CLI option groups (`GRP_DATASETS`, `GRP_REQUIRED`, `GRP_REGIONS`, `GRP_PARAMETERS`, `GRP_FLAGS`, `GRP_OPTIONAL`) which appear in `src/lancet/cli/cli_interface.cpp`, and `CONTIGS_BUFFER_SIZE` if cited in the rule.

5. In `cli.md`, the FORMAT/INFO/FILTER field IDs the rule cites by name (e.g., AD, DP, RMQ, NPBQ, SB, RPCD, CMLOD, etc.). Each must appear as a `##FORMAT=<ID=...>` or `##INFO=<ID=...>` line in `src/lancet/cli/vcf_header_builder.cpp`.

6. In `hts.md`, the HTSlib RAII wrapper class names the rule names: `Extractor` (the file/handle wrapper), `Iterator` (the position-based traversal handle), `Alignment` (the non-owning per-record proxy over `bam1_t`), `BgzfOstream` (the BGZF write stream), and `Reference` (the FASTA wrapper). The internal deleter types `HtsFileDeleter` and `SamHdrDeleter` live in the `extractor.h::detail` namespace. Each must exist with the contracts the rule describes in the corresponding `src/lancet/hts/<name>.h` header.

7. In `base.md`, the macros and primitives the rule names: `LANCET_ASSERT` (must appear in `src/lancet/base/assert.h`), `LANCET_DEBUG_MODE`, `LANCET_VERBOSE_MODE`, the `compute_stats.h::OnlineStats` class with Welford recurrence members.

If any cited name is absent from source, or any cited value differs from source, that's a finding. The user updates either the rule (if the source change was intentional and the rule is now stale) or the source (if the rule was right and the source drifted). Names in the rules that are NOT in the lists above are not part of this pairing's standing checks; the audit grows as new pairings are added.

**Bundle locations to check:** every `.claude/rules/<layer>.md` body that names a source-code constant or class.

### Pairing 12 — Agent memory accuracy

**Source of truth:** the project's current source code, against which each project-scoped agent-memory file is checked.

Project-scoped memory files at `.claude/agent-memory/<agent>.md` (excluding `perf-analyst.md` which is user-scoped) accumulate active knowledge, decision logs, and rejected-decisions lists over time. Without periodic review they bloat and drift: active knowledge entries refer to renamed classes, decision-log entries cite line numbers that have moved, rejected-decisions list ideas that the project has since adopted.

For each project-scoped memory file:

1. **Stale active knowledge.** Entries should reference current source. An entry that says "the genotyper's `EvaluateAllele` handles edge case X" should still match an `EvaluateAllele` member in the genotyper (or be flagged for update if the function was renamed `ScoreAllele`). Walk each active-knowledge entry; spot-check at least three entries per file against current source.

2. **Decision-log entries older than 12 months.** These are candidates for either promotion (if the decision is now bedrock and belongs in AGENTS.md or a layer rule), retirement (if the decision is no longer relevant), or compaction (if the entry has been superseded by later entries). Old decisions that aren't being applied should be removed; old decisions that ARE being applied should be promoted somewhere durable.

3. **REJECTED decisions that have since been adopted.** If a memory file says "we tried X and rejected it because Y" and the current source contains X, either the rejection was reversed (mark as such and move to active knowledge) or X is being snuck in without revisiting the rationale (a real finding).

4. **Bloat.** Memory files over 200 lines should be compacted. Long files mean the agent reads through more material every invocation — the per-invocation token cost rises linearly. Compaction means: merge similar entries, drop superseded ones, summarize verbose entries.

5. **Format consistency.** Each file should have the three top-level sections (Active knowledge / Decision log / REJECTED decisions). Files missing a section or with custom sections are findings (drift away from the documented format makes the agent's own access to memory less reliable).

The user-scoped `perf-analyst.md` is exempt from this audit — it's per-developer, not shared. The agent itself is encouraged to audit its own file periodically.

**Bundle locations to check:** every `.claude/agent-memory/<agent>.md` except `perf-analyst.md`.

## Report format

After walking all twelve pairings, produce a single report with this shape:

```
## /audit-bundle drift report

### Summary
<N total findings across <K of 12> pairings>

### Pairing 1 — Pixi tasks
Status: <CLEAN | DRIFT FOUND>
Findings:
  - bundle/<file>:<line> claims `<old>`; pixi.toml has `<new>`
  - ...

### Pairing 2 — Six-layer chain
Status: ...

[... etc ...]

### Suggested fixes
<for each finding, a proposed edit; one bullet per file>

### Manual review needed
<any finding the script could not auto-resolve, e.g., a layer chain restructure>
```

After producing the report, ask the user how to proceed: apply all auto-fixable findings, apply selectively, or do nothing (report-only mode). Do not apply any fixes without explicit approval.

## When NOT to use this command

Do not use this command for ad-hoc "is X still right?" questions — the comprehensive walk is overkill. Just read the relevant source file and ask Claude.

Do not use this command if `pixi.toml`, `CMakeLists.txt`, or any of the source-of-truth files are mid-refactor (uncommitted changes, broken syntax). The audit relies on the source being parseable; a half-finished refactor will produce noisy false-positive findings.

Do not use this command as a substitute for actually testing changes. The audit catches documentation drift; it does not catch behavioral regression.

## Maintenance

This command's pairing list grows when the bundle adds new claims about source. When you add a new agent or skill that cites specific files, struct names, or values from `pixi.toml`/`CMakeLists.txt`/`.clang-tidy`, add the new bundle locations to the relevant pairing here. When the bundle stops claiming something (e.g., a deleted agent's claims are gone), remove from the pairing.

Pairings 5 and 6 (VCF schema, audited structs) are the most likely to grow, since they're the surfaces where new claims accumulate fastest.
