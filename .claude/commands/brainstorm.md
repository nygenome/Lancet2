# /brainstorm — Generate options before committing to a direction

Open-ended exploration of approaches to a problem. The output is a structured options document under `notes/<topic>/<YYYY-MM-DD>-options.md` that the user can review and pick from. The next step after `/brainstorm` is usually `/spec` against the chosen option.

## Why this command exists

When the user has a problem but no chosen approach, the failure mode is to pick the first reasonable solution and execute it. That's fine when the problem space is small. For substantive features, the right move is to surface 2-4 distinct options with explicit tradeoffs, then pick from a position of comparison rather than from first-fit.

`/brainstorm` formalizes that exploration: an interview clarifies the problem, the options are generated and traded off, and the user picks before committing to a spec.

## Procedure

The command has four phases. Do not skip the interview, do not generate options before clarifying the problem, and do not write the document before the options are stable.

### Phase 1 — Clarify the problem

Use `AskUserQuestion` to interview the user. Plan two to four questions covering:

**The actual problem.** What outcome is the user trying to achieve? "I want to add X" is a solution, not a problem. The problem is "I want users to be able to Y." Drive to the underlying outcome.

**The constraints.** What can the solution NOT do? Existing schema can't change; layer direction must be respected; the benchmark cannot regress; downstream pipelines depend on field FOO. Constraints shape the option space as much as the goal.

**The space.** Is this a small problem with one obvious shape, or a larger one with multiple genuine approaches? If the user can already describe the answer in one sentence, the right command is `/spec` directly, not `/brainstorm`.

**The scope of exploration.** "Quick — give me 2 options" or "deep — walk through 4 distinct approaches"? The user's appetite for exploration shapes the document length.

Ask the questions all together in a single tool call when they are independent; in sequence only when an answer changes which downstream question to ask. This matches the convention `/spec` uses.

### Phase 2 — Generate options

Generate 2-4 distinct options. For each option:

- **Approach** — one or two sentences describing the shape of the solution.
- **How it works** — a short paragraph or bullet list of the moving parts.
- **Tradeoffs** — what this option is good at, and what it sacrifices.
- **Good fit when** — the conditions under which this is the right pick.
- **Bad fit when** — the conditions under which it's the wrong pick.

The options must be **distinct in shape**, not minor variations. "Use a hash map" vs "use a hash map with a different hash function" is one option. "Use a hash map" vs "use a sorted vector with binary search" vs "use a probabilistic structure (Bloom filter)" is three distinct options.

If you cannot generate 2 distinct options, tell the user — the problem may not have a meaningful option space, in which case `/spec` directly is correct.

### Phase 3 — Write the document

Write `notes/<topic>/<YYYY-MM-DD>-options.md` where `<topic>` is a short slug derived from the problem (kebab-case, ≤30 characters) and `<YYYY-MM-DD>` is today's date. If the directory `notes/<topic>/` doesn't exist, create it. The datetime convention matches `/spec` (`<YYYY-MM-DD>-spec.md`) and `/investigate` (`<YYYY-MM-DD>-investigation.md`).

Document template:

```markdown
# Brainstorm: <problem statement>

**Date:** <YYYY-MM-DD>
**Goal:** <restated outcome>
**Constraints:** <list>

## Option 1: <approach name>

<approach paragraph>

**How it works:**
- <bullet>
- <bullet>

**Tradeoffs:** <good at / sacrifices>

**Good fit when:** <conditions>
**Bad fit when:** <conditions>

## Option 2: <approach name>

[same shape]

## Option 3: <approach name>

[same shape]

## Recommendation

<your read of which option fits this user's stated constraints best, with one sentence of rationale. Don't restate the option; just say "Option N because <reason>". The user makes the final call.>
```

After writing the document, summarize the options to the user in chat (one-line per option). Do not paste the full document back — the file path is the deliverable.

### Phase 4 — Hand off

Use `AskUserQuestion` to ask which option to take forward:

```
question: Which option do you want to take forward?
options:
  - Option 1 — <short label>
  - Option 2 — <short label>
  - Option 3 — <short label>
  - Brainstorm more — none of these fit
  - Hold off — review the document first, decide later
```

If the user picks an option, hand off to `/spec` with the chosen option as input: tell them the brainstorm document at `notes/<topic>/<YYYY-MM-DD>-options.md` is ready, suggest running `/spec` with that document as input. Do NOT auto-invoke `/spec` — the handoff is explicit; the user types `/spec` themselves once they've picked.

If "brainstorm more": ask what specifically to explore further (a new constraint, a different angle on the same problem, etc.) and re-do Phases 2-3.

If "hold off": tell the user the document is ready at `notes/<topic>/<YYYY-MM-DD>-options.md` and stop.

## When NOT to use this command

Do not use `/brainstorm` for problems where the answer is obvious in one sentence. The interview overhead is wasted. Use `/spec` directly.

Do not use it for trivial decisions (typo fix, dependency bump, single-line behavior change). The exploration overhead is wasted on changes whose answer is in the first reasonable thought.

Do not use it as procrastination. If you've been brainstorming for two days and still have no preferred option, the bottleneck is decision-making, not options.
