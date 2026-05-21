---
name: external-interface-changes
description: Use when proposing to add, remove, rename, change cardinality of, or silently change semantics of a VCF FORMAT, INFO, or FILTER field, OR when changing a user-facing CLI flag (rename, default, semantics). Trigger on "add a FORMAT field", "rename SB to STRBIAS", "change AD from R to A", "deprecate this field", "the field still exists but means something different now", "rename the --tumor flag", "change the default for --num-threads". Walks the user through a 5-operation matrix (add / rename / cardinality-change / remove / silent-semantic-change), produces a coordinated edit set across vcf_header_builder.cpp + caller code + tests + docs, contributes a structured paragraph to the commit body, and invokes the vcf-validator subagent before finalizing. Does NOT cover internal struct fields or non-VCF data formats — those go through normal review.
allowed-tools: Read, Glob, Grep, Edit, Write, Bash
---

# External interface changes on Lancet2

Lancet2's VCF schema and CLI flag set are external interfaces that downstream pipelines depend on. Most VCF FORMAT/INFO fields are read by automated tooling that doesn't tolerate silent renames, cardinality changes, or semantic shifts. CLI flags are the same — wrapper scripts break when `--tumor` becomes `--case-sample` or `--num-threads`'s default changes.

This skill formalizes the procedure so external-interface changes are explicit, coordinated across all the places that need to know, and reviewed by `vcf-validator` before merge.

## Self-check before proceeding

Before reading further, verify the assumptions this skill rests on:

1. **The user named a specific FORMAT/INFO/FILTER field or CLI flag.** If they only said "the schema needs updating" without a specific field, ask which one. The procedure is per-field, not per-spec-version.
2. **The change is to an external interface.** Internal struct fields, internal class APIs, and non-VCF data formats (probe TSVs, debug logs) are NOT external interfaces and go through normal review, not this skill.
3. **The file paths cited in this skill exist in the current source.** This skill names specific files (`src/lancet/cli/vcf_header_builder.cpp`, `src/lancet/caller/variant_call.cpp`). Confirm they exist via `Glob` before relying on the file paths. If a path is stale, surface that as a finding rather than blindly editing.
4. **The user has not already started editing.** If they have, read the partial diff first and integrate it into the proposed edit set rather than ignoring their work.

## The five operations

Every external-interface change falls into one of these:

1. **Add** — new FORMAT/INFO/FILTER field or new CLI flag. Lowest-risk; downstream consumers either start using the new field or keep ignoring it. Still needs documentation.
2. **Rename** — an existing field gets a new name with the same semantics. High-risk: every downstream parser breaks. Requires deprecation period and dual-write if widely consumed.
3. **Cardinality change** — a FORMAT field's `Number=` attribute changes (e.g., `Number=R` to `Number=A`, `Number=1` to `Number=A`). Subtle: the new file parses successfully but downstream code that indexed by allele count gets wrong values. Highest-risk class because the breakage is silent.
4. **Remove** — a FORMAT/INFO/FILTER definition or CLI flag is deleted. Downstream code that referenced it now fails or gets default values. Requires deprecation period.
5. **Silent semantic change** — the field's name and cardinality stay the same but the underlying meaning changes (e.g., `SB` was Phred-scaled strand bias, now is log odds ratio). The most dangerous because no parser breaks; downstream interpretations are silently wrong. **ALWAYS rename the field instead** of changing semantics — even a small rename like `SB` → `SB2`.

## Step 1 — Interview before editing

Do NOT begin editing on first mention. Ask the user:

1. **Which operation is this?** Show the five-operation list and ask them to pick. If they say "we're just changing what SB means," walk through why that's operation 5 and propose a rename instead.
2. **Which field, and what's the change?** Get the exact field name. The current `##FORMAT`/`##INFO`/`##FILTER` line lives in `src/lancet/cli/vcf_header_builder.cpp` (FORMAT defs in `FORMAT_STR_HEADER`; case/control INFO defs in `CASE_CTRL_INFO_HDR_LINES`). Read the file to see exact line spans rather than assuming.
3. **What downstream consumers will be affected?** Names of pipelines, scripts, notebooks. The user may not know all of them; document what they do know in the commit body.
4. **For renames and removes: what is the deprecation strategy?** Either dual-write for N releases, or a hard cutover with a CHANGELOG warning. Hard cutover is rare and should be justified.
5. **For cardinality changes and silent-semantic changes: are you sure?** These are the highest-risk classes. Push back if the user hasn't fully reasoned about downstream impact.

If any answer is unclear, stop and ask. Don't proceed until the operation type and rationale are explicit.

## Step 2 — Find every coordinated location

A FORMAT/INFO/FILTER change touches at least:

- `src/lancet/cli/vcf_header_builder.cpp` — the canonical declaration site for FORMAT/INFO/FILTER lines.
- `src/lancet/caller/` — the code that populates the field (typically `variant_call.cpp` for FORMAT field emission, `variant_support.cpp` for value computation).
- The relevant test under `tests/caller/` — for the FORMAT-field-emission unit tests.
- `docs/guides/vcf_output.md` — the user-facing schema documentation.
- `.chglog/CHANGELOG.tpl.md` and the upcoming CHANGELOG entry — schema changes are a release note category.

A CLI flag change touches:

- `src/lancet/cli/cli_interface.cpp` and the relevant `pipeline_runner.cpp` — flag definition and handling.
- **No `tests/cli/` directory exists** — the CLI is exercised end-to-end via `.Codex/scripts/e2e_pipeline_test.sh` and the broader `pixi run test` task. There are no unit tests for CLI parsing directly. If the CLI flag change has subtle parsing implications (e.g., a new alias, a default change), surface a recommendation to the user that an end-to-end test in `e2e_pipeline_test.sh` may need updating.
- `docs/index.md`, `docs/guides/wgs_analysis.md`, `docs/guides/targeted_analysis.md` — user-facing CLI usage.

Use `Glob` and `Grep` to find every occurrence of the field/flag name across the codebase before editing. Build the full list before touching any file.

## Step 3 — Propose the coordinated edit set

Before applying any change, present the user with the full file list and the proposed edit at each. The exact format depends on the change; for the actual file inventory, walk what Step 2 produced. Do NOT use a stock template — the file list is empirical, not assumed.

Wait for explicit confirmation. If the user says "refine," accept their corrections and re-present.

## Step 4 — Apply the edits

Apply the edit set as a single coordinated change. Run `pixi run lint-check` and `pixi run test` after each major edit to surface compilation or test failures early. The protected-paths hook will not block any of these (none are on the protected list).

## Step 5 — Invoke the vcf-validator subagent

Before finalizing, hand the change to the `vcf-validator` subagent:

```
Use the vcf-validator subagent to review this <operation> on FORMAT/INFO/FILTER <field>.

Edits applied: <list from Step 3>

Verify the schema invariants in vcf-validator's frame and surface any
breakage (Number= cardinality, Description completeness, Type validity,
multi-sample consistency, missing-value semantics).
```

The subagent reads the actual VCF output (use a small `pixi run test-asan` invocation, or run the pipeline against `LANCET_TEST_*_REGION_SMALL`) and checks every invariant. Address findings before merge.

## Step 6 — Contribute to the commit body

External-interface changes get a structured paragraph in the commit body:

```
feat: <subject describing the change>

Operation: <add/rename/cardinality-change/remove/silent-semantic-change>
Deprecation: <dual-write strategy, removal version, or "hard cutover">
Rationale: <why this change>

Downstream impact: <which consumers affected, and how they should migrate>

vcf-validator: <clean / findings>
```

The structured form makes the CHANGELOG entry mechanical to generate (`scripts/update_changelog.sh` via git-chglog) and gives downstream maintainers the impact summary they need.

Lancet2's `.chglog/config.yml` does not support scopes; mention the layer in the subject text instead of `feat(caller):`.

## When NOT to use this skill

Do not use this skill for changes to internal struct fields, internal class APIs, or non-VCF data formats (probe TSV columns, debug log lines). Those go through normal code review. Do not use it for purely additive Description changes (typo fixes in an existing `##FORMAT` line) — those are spelling fixes, not schema changes.

## When the change is operation 5 (silent semantic change)

Strongly resist this. Walk through with the user why a rename is the correct response — even a small one like `SB` → `SBNEW`. Silent semantic changes are how downstream pipelines break in production six months later. If the user insists, document in the commit body that this was a deliberate operation-5 with full awareness of the downstream risk, and add a CHANGELOG warning at the top of the entry.
