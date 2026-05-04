# /spec — Interview-driven feature specification

Capture a precise spec for a substantive feature or change before any code is written. The output is two files in `notes/<topic>/` that a fresh Claude Code session can pick up to execute the work.

## Why this command exists

Substantive feature work (anything above ~50 lines or touching multiple files) goes badly when the spec is implicit. Claude infers gaps from context, reviewers infer gaps differently, and the result is rework or the wrong feature. The cost of writing a good spec up front is real but bounded; the cost of executing the wrong spec compounds.

This command makes the spec writing structured: an interview drives the spec, the spec lands in a known location, and a fresh session executes the spec without the conversational overhead that drove it.

## Procedure

The command has three phases. Do not skip the interview, and do not start writing the spec before the interview is complete.

### Phase 1 — Interview

Use the `AskUserQuestion` tool to interview the user. Plan three to seven questions covering the dimensions below; ask the most important ones first so the user can shape the spec before fatigue sets in. Ask the questions all together in a single tool call when they are independent; ask in sequence only when an answer changes which downstream question to ask.

The dimensions to cover, in priority order:

**Goal.** What is the user-facing or correctness outcome? "Add a CMLOD field" is too thin — "verify CMLOD computation against a reference implementation and emit it as a per-ALT FORMAT field" is closer. The goal should be testable; if you cannot describe a test that would pass when the goal is met and fail when it is not, the goal is too vague.

**Scope.** Which directories, files, or layers does this touch? Which does it deliberately not touch? An explicit "not touching" list shapes the scope as much as the "touching" list — it pre-commits to discipline that prevents scope creep.

**Acceptance criteria.** What concrete checks must pass before the work is considered done? Unit tests, end-to-end VCF concordance, a benchmark not regressing, lint/format clean, a particular edge case handled. List each criterion as a separate bullet so they can be verified individually.

**Constraints.** What can the implementation not do? Must not change existing FORMAT/INFO/FILTER semantics. Must not introduce upward layer dependencies. Must not regress benchmark X. Must not change the VCF column order. The constraint list is what makes a spec defensive.

**Open questions.** What does the user not yet know that needs research before implementation can begin? Items here are flagged as TODO in the spec; the executing session will need to resolve them before writing code.

**Test data.** Which test data fixtures will be used to validate the change? The `test-data-locations` skill catalogs the available datasets; if the answer is "I don't know," consult that skill before continuing.

### Phase 2 — Draft the spec

After the interview, write `notes/<topic>/<YYYY-MM-DD>-spec.md` where `<FEATURE>` is a short slug derived from the goal (kebab-case, ≤30 characters). If the directory does not exist, create it.

The spec template is:

```markdown
# <Feature title>

## Goal
<1–3 sentence summary of the outcome>

## Scope
**Touching:**
- <file or directory>: <what changes there>
- ...

**Not touching:**
- <file or directory>: <reason>
- ...

## Acceptance criteria
- [ ] <criterion 1>
- [ ] <criterion 2>
- ...

## Constraints
- <constraint 1>
- <constraint 2>
- ...

## Open questions
- <question 1> — <plan for resolving>
- ...

## Test data
<which fixtures, with paths or env-var names>

## Related context
<any prior commits, issues, docs, or code references that the executing
session should read first>
```

After writing the spec, summarize it back to the user in chat and ask for explicit approval. Do not proceed to Phase 3 until the user has said yes (or has corrected the spec — in which case loop back to draft, re-summarize, re-ask).

### Phase 3 — Write the kickoff prompt

After approval, write `notes/<topic>/kickoff.md`. This is the prompt the user will paste into a fresh Claude Code session to execute the work. The kickoff prompt:

```markdown
# Kickoff: <Feature title>

You are picking up implementation of a feature whose spec lives at
`notes/<topic>/<YYYY-MM-DD>-spec.md`. Read that file first, in full, before doing
anything else.

## Execution discipline

1. Read the spec carefully. If anything in the "Open questions" section is
   unresolved, resolve it before writing code — propose your resolution
   and confirm with the user.
2. Plan the work as a sequence of small commits, each one independently
   reviewable. Share the plan with the user before starting.
3. Implement one acceptance criterion at a time. Run the relevant lint /
   build / test gauntlet after each.
4. Use the layer-direction hook, naming hook, and commit-message hook —
   they will catch violations the executor would otherwise miss.
5. When all acceptance criteria are checked, run `/fix-and-validate` for the full
   gauntlet, then `/e2e-pipeline-test` if the change touches the pipeline runtime,
   then `/commit` to compose the merge commit.

## Skills and agents that may apply
<list relevant skills/agents from the bundle, with one-line reason for each>
```

The "skills and agents that may apply" section is filled in based on the spec content — for example, a feature touching the caller layer flags `assembly-and-calling-expert`, a feature with VCF schema implications flags `vcf-validator`, etc.

### Phase 4 — Offer to set up a feature worktree

After writing both files, ask the user whether they want a feature worktree set up. Use `AskUserQuestion`:

```
question: Set up a feature worktree under .worktrees/<FEATURE>/ so the
implementation session works in an isolated branch?
options:
  - Yes — create the worktree
  - No — work in the main checkout
```

The worktree path `.worktrees/<FEATURE>/` is the project's convention; the bundle's `protected_paths.txt` prevents the agent from writing into `.worktrees/` (so an implementation session can't accidentally pollute it from outside), and `.gitignore` excludes the directory.

If the user picks "Yes":

```bash
git worktree add .worktrees/<FEATURE> -b <FEATURE>
```

The branch name matches the worktree directory name and the FEATURE slug used for the spec. After creation, surface the next-step instruction:

```
Worktree created at .worktrees/<FEATURE>/ on branch <FEATURE>.

To execute the spec in a fresh session:
  cd .worktrees/<FEATURE>
  claude  # start a fresh Claude Code session in the worktree

Then paste the contents of notes/<topic>/kickoff.md as the first
message. Note that notes/<topic>/<YYYY-MM-DD>-spec.md is in the main checkout
but visible to the worktree session through the kickoff prompt's
explicit path reference.
```

If the user picks "No," skip the worktree setup. The implementation session runs in the main checkout. This is fine for short features; the worktree is most valuable when the change is multi-day or interleaves with other work.

After both files are written and the worktree decision is resolved, tell the user the spec and kickoff are ready. Do not start the implementation in the current session.

## When NOT to use this command

Do not use `/spec` for trivial changes (typos, formatting, dependency bumps, single-line bug fixes). The interview overhead is wasted on changes whose spec fits in a commit subject. Use it for changes that involve judgment about scope, multiple files, or non-obvious acceptance criteria.

Do not use it as procrastination. If you can describe the change in one sentence and you know what success looks like, just do the change.

## Maintenance

This command's interview list is the highest-leverage thing to tune. If you find that a particular dimension consistently produces unhelpful answers (or that a critical dimension is missing from your typical features), edit the dimensions list above. The template can also evolve — keep it minimal but add sections that consistently turn out to be needed.
