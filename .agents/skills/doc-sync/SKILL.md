---
name: doc-sync
description: Use when keeping in-source comments, dev docs in docs_dev/, and user docs in docs/ in sync with each other and with code. Two modes — Mode A (targeted sync after a code change: find every documentation reference to the changed code, propose a coordinated edit set, gate on a consistency check) and Mode B (periodic semantic audit: read a subdirectory file-by-file applying a five-pass methodology — structural compliance, comment quality, language compliance, header-vs-implementation split, cross-file synchronization — and produce a tiered backlog of findings). Trigger on "I just changed X, did I miss the docs?", "stale comment audit", "rigorous docs check", "the architecture guide says Y but the code does Z", "we renamed FOO to BAR; what else needs updating?", "comment audit on src/lancet/cbdg/", "quarterly audit on the caller layer", "the docstring on Foo doesn't match the implementation". Does NOT cover docs_dev/investigations/ (immutable by design), agent-memory files, .Codex/ READMEs, AGENTS.md, or CHANGELOG.md.
allowed-tools: Read, Glob, Grep, Edit, Write, Bash
---

# Doc sync on Lancet2

Lancet2 has three documentation layers that must stay synchronized with source code: in-source comments, developer documentation under `docs_dev/`, and user-facing documentation under `docs/`. The synchronization rule is mechanical: every claim about the codebase made in any layer must remain consistent with the actual code as the codebase evolves. The authoritative wording of the rule is in `docs_dev/style/sync_and_verification.md`.

This skill operates in two modes. **Mode A** is targeted sync — the user changed code (or a doc) and wants the rest of the layers to follow. **Mode B** is a periodic semantic audit — a file-by-file walk of a subdirectory checking for accumulated drift, comment-quality issues, and cross-file inconsistencies.

The choice of mode is determined by the trigger. "I just renamed X" or "the docs are out of date for this subsystem" is Mode A. "Audit `src/lancet/cbdg/`" or "quarterly comment-quality sweep" is Mode B.

The exception that crosses both modes is `docs_dev/investigations/`. Investigations capture a moment in time and are intentionally NOT maintained as the codebase evolves. If understanding later changes, write a new investigation that supersedes the old one rather than editing the old one in place.

## Mode A — Targeted sync after a code change

### Step A1: Identify the direction(s) of sync

Three possible starting directions:

- **Code-to-docs.** The user changed code. Find every `docs_dev/` and `docs/` location describing the changed code. Most cases are this direction.
- **docs_dev-to-code-and-user-docs.** The user updated a developer doc; find code comments and user-facing docs that describe the same thing.
- **user-docs-to-code-and-dev-docs.** The user updated a user-facing doc; find code comments describing the implementation and dev docs describing design history.

### Step A2: Run the cross-layer audit

```bash
# Replace OLD_NAME with the term being renamed/changed
grep -rn "OLD_NAME" \
    src/**/*.h src/**/*.cpp \
    docs/**/*.md \
    docs_dev/**/*.md \
    tests/**/*.cpp
```

High-value paths to always include for VCF-schema-adjacent changes:

- `src/lancet/cli/vcf_header_builder.cpp` — VCF FORMAT/INFO/FILTER definitions (canonical declaration site).
- `src/lancet/caller/variant_call.{h,cpp}` — FORMAT field-name comments and per-record emission code.
- `src/lancet/caller/variant_support.cpp` — value computation for the metric.
- `docs/guides/vcf_output.md` — user-facing schema documentation.

For CLI flag changes: `src/lancet/cli/cli_interface.cpp`, `pipeline_runner.cpp`, `docs/index.md`, the relevant `docs/guides/*.md`.

For subsystem changes: implementing `src/lancet/<layer>/`, `docs_dev/subsystems/<topic>.md`, `docs_dev/workflows/<topic>.md`, `docs/guides/<topic>.md`.

For algorithm changes: implementing source files, possibly `docs_dev/architecture/` if cross-cutting consequences, and the user-facing `docs/guides/<topic>.md`.

### Step A3: Surface every divergence before editing

Before applying any change, list every location found and present them:

```
Code change: renamed <old-name> to <new-name>.

Drift found:
  Code comments:
    - src/lancet/<layer>/<file>:<line> — <what references the old name>
    - ...
  docs_dev:
    - docs_dev/<sub>/<file>:<line> — <reference>
    - ...
  docs:
    - docs/<sub>/<file>:<line> — <reference>
    - ...

Proposed updates: rename in all <N> locations.

Confirm to proceed (yes/no/refine)?
```

Wait for confirmation. If "refine", accept corrections and re-present. The "show before applying" discipline is the same as `external-interface-changes` — the user catches mistakes (or scope misjudgments) before edits propagate.

### Step A4: Apply the coordinated edit set

A single coordinated change set across all three layers wherever drift was found. Run the relevant verification after each major edit:

- For `docs/` changes: `pixi run docs-build` (mkdocs `--strict` mode; broken cross-links fail the build).
- For `docs_dev/` changes: re-run the Step A2 audit with the OLD name and confirm zero matches outside `docs_dev/investigations/`.
- For code comment edits: review is the only verification (a comment-only edit cannot break compilation).

### Step A5: Note the doc updates in the commit

Lancet2's `.chglog/config.yml` accepts only `feat`, `fix`, `perf`, and `chore` as commit types; `docs:`, `refactor:`, etc. parse but are silently dropped. The chglog header pattern also does not support scopes — `refactor(caller): ...` is dropped.

For documentation-only commits, use `chore:`. For mixed code+docs changes, the primary code commit type wins; mention the doc updates in the body:

```
chore: rename <old-name> to <new-name>

<rationale for the rename>

Doc updates:
  Code comments: <files>
  docs_dev:      <files>
  docs:          <files>
```

## Mode B — Periodic semantic audit

A file-by-file, line-by-line read of a target subdirectory's source, tests, and documentation to verify that comments and docs say true things, explain the right things, use the right words, are in the right place, and stay synchronized across code and docs.

Process namespaces in dependency order: `src/lancet/base/` → `src/lancet/hts/` → `src/lancet/cbdg/` → `src/lancet/caller/` → `src/lancet/core/` → `src/lancet/cli/`. Then tests/benchmarks (only Pass 2 applies — they don't need algorithmic commentary). Then `docs/guides/*.md` and `docs/reference.md` (only Passes 3 and 5 apply).

A skilled auditor catches all five categories of finding in one pass — do not split into five sequential reads of the same file.

### Pass 1 — Structural compliance

For each `.h` and `.cpp`:

- **Include ordering**: four-tier convention (main header → project → third-party → C++ stdlib → C stdlib), blank lines between tiers.
- **Include style**: angle brackets only for C/C++ stdlib; quoted includes for everything else (Abseil, spdlog, htslib, minimap2, spoa, moodycamel, Catch2). Clang-format does not catch this.
- **`extern "C"` wrapping**: C-only third-party headers (htslib, minimap2) are wrapped with a brief comment.
- **Member layout**: structs/classes declare members in descending alignment (8B → 4B → 2B → 1B) with size annotations. Constructor initializer lists match declaration order.
- **NOLINT discipline**: scoped `NOLINTNEXTLINE(check-name)` or `NOLINTBEGIN/END(check-name)` only — never bare, never inline, always with a rationale on the line(s) above. Defer to `clang-tidy-discipline` skill for the full discipline.

### Pass 2 — Comment quality

For each comment block, ask four questions:

**Does it explain why, not what?** The code shows what; the comment must explain why this approach was chosen and what would break if it were done differently. `++mAligned; // increment aligned count` restates the code. The replacement explains rationale: `// Record the weakest base quality in the variant region — the confidence in an indel observation is bounded by the least confident base (weakest-link).`

**Is every claim factually correct?** Read the code the comment describes. Common drift patterns: comment says "returns X" but code returns Y after a refactor; comment cites parameter Z that was renamed or removed; comment claims O(N) but the actual complexity is O(N log N) due to a sort.

**Is the math correct and annotated?** For any formula, verify each variable corresponds to a real program variable, check that representative values are correct (e.g. `Q=30 → 0.999`), and confirm the formula produces the stated result.

**Are invariants and boundary conditions stated?** For any function with non-obvious preconditions, the coordinate system (0-based vs 1-based, absolute vs relative), edge cases (empty input, zero coverage, single-element groups), and sentinel return values (`nullopt`, `0.0`, `SIZE_MAX`) should be explicit.

### Pass 3 — Language compliance

This pass is not a checklist exercise. The rules below catch common violations, but new filler words appear in future code that no table can anticipate. The core skill is **semantic reading of English**: reading each sentence and asking *does every word carry technical meaning*?

**The deletion test is the primary tool.** Remove the suspect word; re-read the sentence. If the technical meaning is unchanged, the word was filler. If the sentence loses a real distinction, the word earns its place.

Filler words previously caught: `natively` (always delete), `organically` (always delete; code is deterministic, not organic), `seamlessly` (always delete; never true for code readers debugging an issue), `fundamentally` (acceptable when describing structural/topological difference; delete when "very" or "importantly"), `inherently` (acceptable for a mathematical property; delete when "obviously"), `intrinsically` (same test as `inherently`).

Filler adverbs previously caught: `mathematically` (keep when describing a derivation step; delete when restating arithmetic), `explicitly` (keep when contrasting with implicit behavior; delete when emphasizing), `structurally` (keep when describing data-structure topology; delete when modifying a verb for emphasis), `aggressively` (keep when quantifying a tradeoff with a number; delete as vague intensifier), `dynamically` (keep when distinguishing runtime from compile-time; delete when redundant), `deliberately` (keep when contrasting with accidental; delete when redundant with a `CRITICAL:` note), `unconditionally` (keep for guard clauses; delete when redundant with "always"), `affirmatively` (never keep; replace with specifics), `maliciously` (never keep; anthropomorphizes software).

**Watch for adverb clusters.** Two or more adverbs modifying the same clause is a strong signal that at least one is filler.

**Watch for anthropomorphization.** Code does not "want", "try", "refuse", "fight", or "maliciously win". If a comment attributes intent to software, rephrase to describe the mechanism.

**Watch for hedging.** "Essentially", "basically", "more or less", "in a sense" — if the claim is true, state it directly. If it is an approximation, quantify the approximation.

**Unexplained jargon.** Apply the biologist test: would a biologist or clinician have to stop and look it up? If yes, either replace with plain language ("independent" not "orthogonal") or define on first occurrence with a brief parenthetical.

### Pass 4 — Header vs implementation split

Headers (`.h`) document the public API: what each method does, what types it accepts and returns, what invariants it maintains. They should NOT contain algorithmic walkthroughs. Implementation files (`.cpp`) document the algorithm: step-by-step logic, formula derivations, ASCII diagrams, coordinate-system handling. They should NOT duplicate the API contract. Massive comment blocks in headers explaining algorithms are a Tier 2 finding; the explanation belongs in the `.cpp`.

### Pass 5 — Cross-file synchronization

For each concept the current file documents:

- If the file mentions a VCF FORMAT/INFO/FILTER field, the field name matches `variant_call.h`, the registration in `vcf_header_builder.cpp`, and the user-facing description in `docs/guides/vcf_output.md`.
- If the file documents a metric (NPBQ, SB, RPCD, etc.), the formula and description match the implementation in `variant_support.cpp` and the corresponding doc in `docs/guides/`.
- If the file mentions another file, function, or class by name, that name still exists.

The patterns that go stale most reliably: architecture diagram comments referencing FORMAT field lists, metric descriptions in tangentially related functions, `ReadEvidence` member comments referencing FORMAT field names, operating-mode tables in `vcf_output.md` referencing QUAL computation methods.

### Findings classification (Mode B)

Every finding is classified into exactly one tier:

- **Tier 1 — Convention violation, must fix.** A codified rule from `docs_dev/style/` or `.clang-tidy` is broken: angle-bracket include on a third-party header, banned filler word with no technical meaning, members declared out of alignment order, inline `NOLINT`, etc.
- **Tier 2 — Comment quality issue, should fix.** The comment is technically not rule-breaking but fails to serve its purpose: comment restates the code instead of explaining why, non-obvious logic has no rationale, formula in comment does not match the code, missing boundary-condition documentation.
- **Tier 3 — Documentation sync gap.** A web documentation page or cross-file reference is stale or inconsistent: `docs/guides/architecture.md` describes a metric that was renamed in code, a VCF FORMAT field comment uses a different formula than `variant_support.cpp`.

### Report template (Mode B)

Structure the report per-tier, per-file. Each finding has:

```markdown
### [path/to/file.cpp:NNN] — <Rule Category>

**Current**: <verbatim offending comment text>

**Rule**: <quote the specific convention rule being violated>

**Fix**: <the exact replacement text>
```

End the report with an explicit list of files that passed the audit with no findings. The pass list proves the audit was exhaustive — every file in scope appears either in findings or in the pass list — and identifies gold-standard files that new contributors should study.

### Anti-patterns to avoid during audit (Mode B)

**Grep-only scanning.** Grep finds words, not meaning. `grep "natively"` catches the word but not "the DM model handles multi-allelic sites" (a sentence that restates the class docstring 20 lines above). Read the code.

**Checking comments in isolation.** A comment is correct only in the context of the code it describes. Reading a comment without reading the adjacent code cannot detect factual inaccuracy.

**Applying rules mechanically.** "fundamentally" is sometimes correct (structural topology) and sometimes filler (emphasis). Every instance requires judgment.

**Reporting style preferences as violations.** The audit checks rules and factual accuracy, not personal writing style. "Subtracting `local_raw_score` carves a hole" is fine — it is a metaphor that aids understanding. "Subtracting `local_raw_score` natively carves an exact algebraic hole" is a violation — three filler words.

### Audit cadence (Mode B)

The audit fits into the project at three frequencies:

- **Per-PR (writer self-audit, then reviewer).** The developer applies Passes 1 through 4 to their own changes before opening the PR.
- **Quarterly (full-scope sweep).** A full-namespace audit (one namespace at a time, in dependency order). Findings committed under `notes/audits/YYYY-QN/`.
- **Post-refactor (focused sync).** Any rename, restructure, or metric change triggers a focused Pass-5 sync audit on every file that referenced the changed concept.

## When NOT to use this skill

- `docs_dev/investigations/` — immutable by design.
- AGENTS.md, AGENTS.md, `.Codex/` READMEs — bundle-internal, updated via `/audit-bundle`.
- Agent-memory files — append-only logs.
- CHANGELOG.md — generated from commits via git-chglog.
- Vendored deps under `cmake-build-*/_deps/` — read-only.
- Pre-merge review of a specific diff — that's `fresh-reviewer`'s job. The two are complementary: `fresh-reviewer` asks "did this diff introduce a new violation?"; this skill (Mode B) asks "what's the existing comment-and-docs backlog across this subdirectory?"

## When the change is large enough to need a new doc

If the code change introduces a new subsystem, substantive feature, or fundamentally changes a workflow, the right move is often a new document rather than incremental edits:

- New subsystem → new `docs_dev/subsystems/<topic>.md` (and possibly a paired `docs/guides/<topic>.md`).
- New recurring workflow → new `docs_dev/workflows/<topic>.md`.
- Substantive design decision → use `/arch-decision-record`.
- Postmortem or debug archaeology → use `/investigate`.
