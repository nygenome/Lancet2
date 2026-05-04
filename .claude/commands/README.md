# Slash commands

Sixteen slash commands live here. Each one is a workflow the user wants to start deliberately rather than have inferred from context. Skills fire automatically when Claude judges them relevant; slash commands fire only when the user types them. The discriminator: does explicit user intent matter for this workflow?

## What's here

Operational commands:

- **`/fix-and-validate`** — fix-mode validation: runs `pixi run iwyu-fix` (which chains clang-format), `pixi run lint-check` (read-only; clang-tidy auto-fix is forbidden), then `pixi run test` against the Release tree. Mirrors the scope of the `pre_commit_gate.sh` PreToolUse hook but applies fixes where the gate runs read-only equivalents. Use to make the diff clean before staging; the gate is the safety net at `git commit` time.
- **`/e2e-pipeline-test`** — runs the full Lancet2 pipeline twice: germline (NA12878 / chr1) then somatic (HCC1395 tumor + HCC1395BL normal / chr4). Confirms exit-0 and a reasonable variant count for each. Truth-set comparison is a separate harness and is not done here.
- **`/commit`** — composes a commit message in the project's two-section style (subject + context paragraphs + optional file-list bullets), grounded in `docs_dev/style/cpp_style.md` § Git commit messages and `.chglog/config.yml`. The validate_commit_message hook independently double-checks at the moment of `git commit`.

Feature kickoff workflow (three composable steps; `/brainstorm` → `/spec` → `/execute-spec`):

- **`/brainstorm`** — option exploration when the user has a problem but no chosen approach. Surfaces 2-4 distinct approaches with explicit tradeoffs, picks one, hands off to `/spec`. Output: `notes/<topic>/<YYYY-MM-DD>-options.md`.
- **`/spec`** — interview-driven feature specification. Asks the user 3-7 questions via `AskUserQuestion`, drafts a spec at `notes/<FEATURE>/spec.md`, and writes a kickoff prompt at `notes/<FEATURE>/kickoff.md` for a fresh execution session. Use for substantive features (above ~50 lines or touching multiple files).
- **`/execute-spec`** — picks up a spec written by `/spec` and executes it as a sequence of TodoWrite-tracked tasks with per-task verification and explicit stop conditions. The spec drives; the executing session follows.

End-of-branch:

- **`/wrap-branch`** — end-of-branch cleanup. Pre-checks the tree, then asks the user to pick: merge to main locally, push + PR, keep the branch, or discard. Handles worktree cleanup based on the choice.

Probe tracking workflow (three composable steps; run any subset):

- **`/probe-concordance`** — Step 1: runs `scripts/truth_concordance.py` to compare a truth VCF against a Lancet2 output VCF. Produces `missed_variants.txt` (input for step 2) and `concordance_details.txt` (optional enrichment for step 3).
- **`/probe-run`** — Step 2: runs Lancet2 with `--probe-variants` and `--probe-results` to collect forensic data on missed variants. Includes a build-staleness check (clean tree + embedded SHA matches HEAD; rebuild via `pixi run build-release` if not).
- **`/probe-analyze`** — Step 3: runs `scripts/analyze_probe_results.py` to attribute each missed variant to a pipeline stage via the 27-level cascade, then hands off to the `probe-interpreter` subagent for focused recommendations.

Drift-and-refresh commands:

- **`/sync-cost-model`** — refreshes `cost-model.md` against current Anthropic Claude Code documentation, current community practice, and the bundle's current contents. Quarterly maintenance ritual.
- **`/audit-bundle`** — comprehensive drift check against project source. Walks twelve source-of-truth pairings (pixi tasks, layer chain, CLI flags, test data, VCF schema, audited structs, code-style ceilings, protected paths sync, bundle-internal counts, path-scoped rules layer alignment, rule source-claim verification, agent memory accuracy) and reports findings with file:line citations. Run quarterly or after a substantive refactor. Does not auto-apply fixes; produces a report and asks the user how to proceed.
- **`/audit-vcf-schema`** — focused VCF schema drift check. Walks every `##FORMAT`/`##INFO`/`##FILTER` line in `vcf_header_builder.cpp` and compares against the bundle's schema claims. Run after every release that touches VCF emission.
- **`/audit-probe-pipeline`** — focused drift check for the probe tracking infrastructure. Walks `docs_dev/subsystems/probe_tracking.md` and the bundle's probe-tracking artifacts (skill, agent, three slash commands) against the actual `src/lancet/cbdg/probe_*.cpp`, `src/lancet/core/probe_diagnostics.cpp`, and the two scripts.

Document scaffolding commands:

- **`/arch-decision-record`** — interview-driven scaffolding for an Architecture Decision Record. Walks the four required sections (Context, Decision, Alternatives Considered, Consequences) via `AskUserQuestion`, assigns the next sequential ADR number, and produces a draft at `docs_dev/architecture/<NNNN>-<slug>.md` following the writing guide in `docs_dev/architecture/README.md`. Use when making a substantive architectural decision worth preserving.
- **`/investigate`** — interview-driven scaffolding for an investigation document (postmortem, debug archaeology, or performance investigation). Asks one disambiguating question for the type, then walks the seven required sections, producing a draft at `docs_dev/investigations/<YYYY-MM-DD>-<slug>.md` following the writing guide in `docs_dev/investigations/README.md`. Use after debugging a non-trivial issue worth preserving as an immutable snapshot.

## Drift-and-refresh commands at a glance

The bundle has four "drift" commands that look different at a glance but serve the same maintenance purpose: keeping the bundle's claims accurate.

- **`/sync-cost-model`** drifts against *external* docs — Anthropic's Claude Code documentation. Run quarterly.
- **`/audit-bundle`** drifts against *project* source — the twelve source-of-truth pairings. Run quarterly.
- **`/audit-vcf-schema`** is a focused subset of `/audit-bundle`. Run more often, since the VCF schema is the surface that drifts most frequently.
- **`/audit-probe-pipeline`** drifts against the probe tracking surface specifically — wide enough (8 C++ files + 2 scripts + 1 doc) that bundling into `/audit-bundle` would dilute focus.

If you want one mnemonic: the bundle is correct against four surfaces (Anthropic docs, project source, VCF schema, probe-tracking pipeline), one command per surface.

## Why slash commands and not skills?

Slash commands are user-typed, so they fire only with explicit intent. Skills can fire from context, which is the right pattern for "consult on demand" workflows but the wrong pattern for "run this specific procedure right now."

The sixteen commands here all benefit from explicit intent:

- `/fix-and-validate` — running it is a deliberate "I want to verify this is clean now" decision; it takes minutes, not seconds.
- `/e2e-pipeline-test` — runs the pipeline against tens of GB of test data; takes 10-30+ minutes; explicit intent matters.
- `/commit` — the user wants the commit composed at a specific moment, not whenever Claude infers it might be appropriate.
- `/spec` — the interview phase is a deliberate setup for a substantive feature, not something Claude should infer from a passing remark.
- `/brainstorm` — option exploration before commitment. The interview clarifies the underlying problem, then 2-4 distinct approaches are surfaced with explicit tradeoffs before the user picks. Hands off to `/spec` once a direction is chosen.
- `/execute-spec` — TodoWrite-tracked execution of a `/spec` document, with per-task verification and explicit stop conditions. The spec drives; the executing session follows.
- `/wrap-branch` — end-of-branch workflow with four reasonable paths (merge, push + PR, keep, discard). Makes the choice explicit and the steps consistent rather than recalling the right git incantation under "I just finished work and want to ship it" conditions.
- `/sync-cost-model`, `/audit-bundle`, `/audit-vcf-schema`, `/audit-probe-pipeline` — maintenance rituals, deliberate by design. Auto-firing them on context would mean refreshing/auditing during normal work, which is the wrong time.
- `/probe-concordance`, `/probe-run`, `/probe-analyze` — each step is independently expensive (especially step 2's 10-30 minute Lancet2 run) and the user often wants to run only one or two of them. Forcing a single monolithic command would obscure what each step contributes.
- `/arch-decision-record`, `/investigate` — document creation is a deliberate act. The failure mode without explicit intent is Claude writing the architectural rationale or postmortem inline in the chat instead of producing the persistent document the user wants on disk.

Workflows that don't need explicit intent — "while you're working on X, also keep Y in mind" — are skills, not slash commands.

## Why not push everything into slash commands?

Slash commands force the user to remember the command name and type it. For workflows the user does often (like reviewing a diff, or following the project's TDD discipline), the friction of recall and typing is real. Skills with good triggers fire when the situation matches; slash commands only fire when invoked.

The right default is "skill if Claude can infer the trigger; slash command if the user wants explicit control."

## Slash command structure conventions

- Slash commands are markdown files at `.claude/commands/<name>.md`. The file name (minus extension) is the command name.
- The body of the file is the prompt that gets injected into the conversation when the user types `/<name>`. It can include shell preprocessing via the `!` syntax (rare; only when fresh data is needed at invocation time).
- Most commands follow a structure of: rationale (why this command exists), procedure (the steps Claude should walk through), and "When NOT to use" section.
- For commands that produce artifacts (like `/spec` writing files to `notes/<feature>/`), explicitly document the output paths. Claude needs to know where to write.

## Maintenance lifecycle

### Adding a new slash command

A slash command earns its place when:

1. The workflow is deliberate (user explicitly chooses to invoke).
2. The workflow is repeatable across many sessions (otherwise just type the prompt directly).
3. The procedure is multi-step or requires specific framing — otherwise the user might as well just describe the task in chat.

Procedure: write the markdown file, name it after the trigger (kebab-case, short). Write the body as a prompt directed at Claude — first person, instructional tone. Document the procedure as numbered steps when the order matters; as goals when it doesn't. Include a "When NOT to use" section. Add to this README's "What's here" section.

Cost reasoning: slash command bodies cost nothing baseline (loaded only when typed). The cost is essentially free; the constraint is whether the user will remember the command exists. See `../cost-model.md`.

### Deleting a slash command

Delete a slash command when:

- You haven't typed it in the last quarter.
- The workflow it encoded has migrated elsewhere (a hook, a skill, an AGENTS.md note).
- A different command or skill covers the same ground better.

The bias should be toward keeping slash commands you do use — they cost nothing baseline and adding them is real friction. But unused commands clutter the `/<tab>` menu and contribute to "I don't know what's available" cognitive load.

Procedure: remove the markdown file, update this README's "What's here" section, add a brief micro-changelog entry.

### Refactoring a slash command

Common patterns:

**Tightening the procedure.** When a command fires but Claude consistently does the wrong thing inside it, the prompt body needs work. Read what Claude did wrong, identify which step in the procedure was ambiguous, and rewrite that step. The most common fix is replacing prescriptive instructions with explicit goals.

**Adding shell preprocessing.** When a command needs fresh data at invocation time (like the current branch name, or a recent commit hash), embed `!command` in the body. The output replaces the placeholder before Claude sees the prompt.

**Splitting a command.** When a command tries to do two things (`/check-and-commit`), users end up wanting one or the other independently. Splitting into two commands is usually a win.

### Reviewing the slash command set

Quarterly. Walk all sixteen commands:

- Have you typed each one? When? Was the result useful?
- Does each command's body still match the codebase? Have any of the cited paths or commands moved?
- Are there workflows you did manually this quarter that should be commands? (One quarter of evidence is usually enough to commit.)

### Retiring a slash command

Retire when on the deletion candidate list two quarters in a row. Slash commands are cheap to keep, so the bias toward retention is reasonable; but commands you don't use erode the discoverability of the ones you do.

## Cost model

See `../cost-model.md` for the full breakdown. Briefly: slash command names are listed but bodies are not loaded into context until typed. Per-session cost is essentially zero; per-invocation cost is the body length added to the main session.
