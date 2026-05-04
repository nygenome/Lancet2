# /arch-decision-record — Draft a new Architecture Decision Record

Interview the user to surface the four required sections of an ADR (Context, Decision, Alternatives Considered, Consequences), assign the next sequential ADR number, and write the document into `docs_dev/architecture/`.

This command exists because the ADR format has a specific shape — four required sections, present-tense Decision, past-tense Context, honest Alternatives section — that is easy to half-do without an interview. The interview is the quality gate. The writing standards live in `docs_dev/architecture/README.md`; read that document if any of the prompts below seem ambiguous.

## When to use

When you've made (or are about to make) an architectural decision that will constrain how the codebase is structured. Examples: choosing one library over another, adopting or removing a layer in the dependency chain, committing to a specific algorithm or data structure where alternatives existed, deciding what *not* to do.

If you're just documenting how an existing system works (no decision involved), you don't want an ADR — you want either a subsystem deep-dive (`docs_dev/subsystems/`) or an architecture overview document (`docs_dev/architecture/<name>.md`, no number). See `docs_dev/architecture/README.md` for the disambiguation.

If you're capturing what was learned during a debugging session, you don't want an ADR — you want an investigation in `docs_dev/investigations/`.

## Procedure

### Phase 1 — Verify this is the right format

Ask the user to confirm in one sentence what decision they're capturing. Then check whether ADR is the right format:

- If the user describes a *decision between alternatives*, ADR is right.
- If the user describes *how a system works*, redirect to subsystem deep-dive or architecture overview.
- If the user describes *what was learned during debugging*, redirect to `/investigate`.
- If the user can't describe the alternatives that were considered, ADR is wrong — they need either an overview (the system is what it is) or no document at all.

Do not proceed to drafting without a clear yes that this is a decision-with-alternatives.

### Phase 2 — Determine the ADR number

```bash
# Find the highest existing ADR number; increment by one.
ls docs_dev/architecture/ 2>/dev/null \
    | grep -E '^[0-9]{4}-' \
    | sort -V \
    | tail -1
```

The number is zero-padded to four digits. If no ADRs exist yet, this is `0001`. The slug is assigned later from the user's Decision section.

### Phase 3 — Interview the user

Use `AskUserQuestion` to walk through each of the four required sections in order. The questions are:

1. **Context.** "What was true about the system, the constraints, the team, or the requirements that forced this decision? Two or three paragraphs in past tense. Be concrete: cite file paths, version numbers, measured numbers, or external constraints — not aspirations."

2. **Decision.** "What did you decide? Write a single present-tense sentence stating the chosen position, then a paragraph elaborating what adopting it looks like in concrete terms."

3. **Alternatives considered.** "What else was on the table? For each alternative, name it, say why it was considered, and say why it lost. Be specific — 'performance' is not a reason; '~15% wall-clock cost on the standard somatic fixture' is. List 2–4 alternatives."

4. **Consequences.** "What does this decision constrain? What follow-on work or non-obvious trade-offs does it create? An ADR with no listed consequences is a red flag — every architectural decision closes some doors. Name them."

If the user gives one-line answers to any of these, push back. The Context section needs *what was true*, not *what we wanted*. The Decision needs the position stated up front, not buried under reasoning. The Alternatives section needs concrete reasons each lost, not a bullet list of names. The Consequences section needs honest trade-offs, not just upsides.

### Phase 4 — Pick the slug

After the Decision is drafted, derive a 2–5 word kebab-case slug from it. Examples: `pixi-over-conda`, `six-layer-dependency-rule`, `spoa-int16-simd-path`. The slug should match the verb tense of the Decision — present-tense statement of the chosen position.

Confirm with the user. The slug locks in the filename: `0007-pixi-over-conda.md`.

### Phase 5 — Write the document

Write to `docs_dev/architecture/<NNNN>-<slug>.md` using the template in `docs_dev/architecture/README.md` § "ADR template (copy-paste)." Required fields in the header:

```markdown
# ADR <NNNN>: <Title — 2-5 words, capitalized>

**Status:** Proposed
**Date:** YYYY-MM-DD
**Last reviewed:** YYYY-MM-DD
```

Both dates start as today.

### Phase 6 — Update the index

Append an entry to the `# Existing documents` section in `docs_dev/architecture/README.md`:

```markdown
- **ADR NNNN:** [<Title>](NNNN-<slug>.md) — <one-sentence summary>
```

### Phase 7 — Surface the result

End the response with:

> ADR <NNNN> drafted at `docs_dev/architecture/<NNNN>-<slug>.md`. Status: Proposed.
> When the decision is reviewed and accepted, change the Status field to "Accepted."
> If a future ADR supersedes this one, add a Note line at the top of this document
> linking forward.

Do not commit the file or open a PR — that's the user's choice.

## When NOT to use

- **For a decision under 50 lines of explanation.** If the decision is genuinely small, a comment in the source plus a `chore:` commit may be enough. Don't pad an ADR to look substantial.
- **For documenting how something works.** That's not a decision; it's an explanation. Use `docs_dev/architecture/<name>.md` (overview) or `docs_dev/subsystems/<name>.md` (deep-dive).
- **For postmortems or debug archaeology.** Those have a different lifecycle (immutable snapshots) and live in `docs_dev/investigations/`. Use `/investigate` instead.
- **For tentative or speculative decisions.** ADRs document decisions that are in force or proposed for force. A "we might want to do X someday" note belongs in a design doc, not an ADR.

## Maintenance

The ADR format is defined in `docs_dev/architecture/README.md` § "Architecture Decision Records (ADRs)". If that document changes, this command's interview prompts may need to update. The `/audit-bundle` command does not currently audit ADR contents.
