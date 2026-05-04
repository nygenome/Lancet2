# /investigate — Capture a postmortem, debug archaeology, or performance investigation

Interview the user to capture an investigation document — a snapshot of what was learned during a hard-to-find bug, a regression postmortem, or a performance puzzle. Write the result into `docs_dev/investigations/`.

This command exists because investigation documents have a strict immutable lifecycle and a fixed seven-section structure that benefits from being walked through in order. The writing standards live in `docs_dev/investigations/README.md`; read that document if any of the prompts below seem ambiguous.

## When to use

When you've just finished tracking down a hard bug or completing a performance investigation, AND the situation meets at least two of these criteria:

- The investigation took more than a working day and the path was non-obvious.
- The bug or issue was caused by a subtle interaction, a violated invariant, or a misunderstood library behaviour — something a future contributor could plausibly hit again.
- A future maintainer working in the same area would benefit from knowing what was already tried and ruled out.
- The fix is in production but the *reasoning* is not obvious from the diff or commit message.

If the bug was found in five minutes with a one-line fix, do not write an investigation. The commit message is enough.

If the issue was "user reported flakiness, can't reproduce, closing," there's nothing to capture — do not write an investigation.

## Procedure

### Phase 1 — Verify this is the right format

Ask the user to summarize the situation in one sentence and the resolution in one sentence. Then check whether an investigation document is the right output:

- If the user describes a debugging session that took serious time and ended with a finding worth preserving, investigation is right.
- If the user describes a decision they made about how to structure something, redirect to `/arch-decision-record`.
- If the user describes how a system *currently works* (no debugging involved), redirect to `docs_dev/subsystems/` or `docs_dev/architecture/`.
- If the user can summarize the entire investigation in two sentences with no false leads or surprises, an investigation document is probably overkill — confirm before proceeding.

### Phase 2 — Pick the type

Use `AskUserQuestion` to clarify which of the three types fits:

- **Postmortem** — incident or regression with user-visible impact. Primary value: why we didn't catch it earlier; what process or test changed.
- **Debug archaeology** — subtle bug that took serious effort to find. Primary value: the path through the symptoms, the false leads, the eventual root cause.
- **Performance investigation** — "why was X slow" that required real measurement work. Primary value: what was measured, what was tried, what worked, why.

The types blur. Pick the one that best describes the situation; the document body can mention overlaps.

### Phase 3 — Determine the filename

Use either a specific date or a quarter prefix:

```bash
# Specific date — for investigations centred on a known day
date +%Y-%m-%d
# 2026-04-29

# Quarter — for investigations spanning longer periods
# Use YYYY-QN format manually based on user input
```

The slug is descriptive lowercase kebab-case, specific enough that a future grep can find it. *"2026-04-15-msa-shifted-coordinate-bug.md"* — not *"2026-04-15-bug.md"*.

### Phase 4 — Interview through the seven sections

Use `AskUserQuestion` to walk each of the seven required sections in order. The questions:

1. **Summary.** "Write one paragraph (~5 sentences) covering: what was the problem, what was found, what was done about it. A reader of just this paragraph should know whether to keep reading."

2. **What we observed.** "What were the concrete symptoms? Logs, error messages, performance numbers, screenshots, test failures. State what was *seen*, not what was suspected at the time."

3. **Investigation path.** "How did you get from symptom to root cause? Be honest about the false leads — 'we initially suspected X because of Y, but ruled it out when Z' is more useful than a clean retroactive narrative."

4. **Root cause.** "What was actually wrong? Cite specific code, specific commits, specific assumptions that turned out to be false. Be concrete."

5. **Resolution.** "What changed? Link to the fix commit or PR. If the fix was non-trivial, briefly explain the reasoning — but do not duplicate the commit message."

6. **What we learned.** "Lessons that would benefit a future contributor. What invariant was violated? What test gap allowed this through? What documentation was missing or wrong? What was a surprising behaviour of a library or tool? Be specific enough that someone hitting similar symptoms can recognize the pattern."

7. **What we did not change.** "What *should* arguably have changed but didn't, and why? If a class of bugs is now known to exist but you chose not to fix them all, say so. This section prevents readers from later asking 'why didn't they fix this properly?'"

If any answer is short or vague, push back. The Summary should be five sentences, not one. The Investigation path should include the false leads. The Root cause should cite specific code. "What we did not change" should not be empty — every investigation has trade-offs that weren't pursued.

### Phase 5 — Write the document

Write to `docs_dev/investigations/<date>-<slug>.md` using the template in `docs_dev/investigations/README.md` § "Investigation template (copy-paste)." Required header:

```markdown
# <Type>: <Specific descriptive title>

**Type:** Postmortem | Debug Archaeology | Performance Investigation
**Captured on:** YYYY-MM-DD
**Status:** Final
**Author:** <user's git handle or name>
```

The author's name comes from `git config user.name` or `git config user.email` if it's available; otherwise ask the user.

### Phase 6 — Update the index

Append an entry to the `# Existing investigations` section in `docs_dev/investigations/README.md`. **Most-recent at top.** Format:

```markdown
- **<date>** — [<title>](<filename>) (<Type>) — <one-sentence summary>
```

### Phase 7 — Surface the result and the immutability rule

End with:

> Investigation captured at `docs_dev/investigations/<date>-<slug>.md`.
>
> Reminder: investigations are **immutable** once merged. Do not edit this document
> later to reflect new understanding. If understanding changes, write a new
> investigation that supersedes this one — see `docs_dev/investigations/README.md`
> § "Superseded investigations" for the procedure.

Do not commit or open a PR — that's the user's choice.

## When NOT to use

- **For a bug found in five minutes.** Use the commit message.
- **For a flaky test that resolved itself.** Use a comment in the test file or close the issue.
- **For decisions about how to structure code.** Use `/arch-decision-record` instead — that's an architecture decision, not an investigation.
- **For ongoing or unresolved issues.** Investigations capture *what was learned*, not *what we're still figuring out*. If the answer isn't known yet, the investigation isn't ready to write. Capture it after the resolution lands.

## Maintenance

The investigation format is defined in `docs_dev/investigations/README.md`. If that document changes, this command's interview prompts may need to update. The `/audit-bundle` command does not audit investigation contents (and shouldn't — they're immutable).

The "Existing investigations" index in the README must be appended to (most recent at top); the slash command handles this in Phase 6.
