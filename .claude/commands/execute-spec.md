# /execute-spec — Execute a spec from a fresh session

Pick up implementation of a feature whose spec was captured by `/spec`. Reads the spec, plans the work as a sequence of tasks, tracks them via TodoWrite, executes them in order with per-task verification, and stops when the spec is complete.

## Why this command exists

A `/spec` document is a contract between the user and a future implementation session. The implementation is most reliable when the spec is read in full, the work is decomposed into tasks before any code is written, and progress is tracked explicitly so a long-running session doesn't lose its place.

`/execute-spec` formalizes that discipline. The spec drives; the executing session follows.

## Usage

```
/execute-spec <path-to-spec>
```

The path is **required**. Example:

```
/execute-spec notes/cmlod-format-field/2026-04-30-spec.md
```

If no path is provided, error out with usage instructions. Do NOT default to "find the most recent spec" — explicit is better than implicit, and a wrong default produces wrong work.

## Procedure

### Phase 1 — Read the spec

Read the entire spec file in full. Pay particular attention to:

- **Goal** — the testable outcome.
- **Scope (Touching / Not touching)** — the explicit boundaries.
- **Acceptance criteria** — the checks that must pass.
- **Constraints** — what the implementation must NOT do.
- **Open questions** — items flagged as TODO that need resolution before code can be written.
- **Test data** — which fixtures will validate the work (germline NA12878 / chr1, somatic HCC1395 / chr4, etc.).

If the **Open questions** section has unresolved items, **stop** before any code. Propose a resolution to the user for each open question and confirm before proceeding. Implementing against an unresolved spec produces work that needs to be redone.

### Phase 2 — Plan as tasks

Decompose the work into a sequence of tasks via `TodoWrite`. Each task should be:

- **Independently reviewable** — a reviewer can read the resulting commit and assess it in isolation.
- **Small enough to verify** — at most a few hundred lines of change before a verification pass.
- **Ordered by dependency** — earlier tasks unblock later ones; the test added in task 3 fails until task 4 implements the behavior.

Share the task list with the user before starting. Wait for explicit approval. If the user requests changes to the plan, accept them and re-share the updated list before any code.

### Phase 3 — Execute one task at a time

For each task:

1. **Mark the task as in-progress** via `TodoWrite`.
2. **Implement the task** — make the smallest change that satisfies the task's stated outcome.
3. **Verify per the spec.** The spec specifies what verification each task needs. Some tasks need `pixi run test`; some need `/fix-and-validate`; some need a specific Catch2 test (`./cmake-build-release/tests/TestLancet2 "[layer]"`); some need a `pixi run docs-build` for documentation pages; some need `pixi run e2e-pipeline-test` for pipeline-level changes. Follow the spec's instruction. If the spec doesn't specify, run `pixi run test` as a minimum.
4. **Mark the task as completed** via `TodoWrite`.
5. **Commit the task's work.** Use `/commit` to compose the commit message; do not write commit messages by hand inside an executing session.

If a verification fails, do not move to the next task. Diagnose the failure, fix it, re-run verification, and only then mark the task complete. The task list is not a milestone-tracker — incomplete verification is incomplete task.

### Phase 4 — Stop conditions

Stop the execution and ask the user when any of these occur:

- **A blocker.** A task cannot be completed because of a missing dependency, a contradictory constraint, or a code-base reality that contradicts the spec.
- **A gap in the plan.** Mid-execution, you realize the spec doesn't cover a case the implementation must handle. Surface the gap; ask the user whether to add a task to the plan or to revise the spec.
- **An unclear instruction.** The spec says X; X has two reasonable interpretations; the choice between them affects the result. Ask which interpretation is intended before guessing.
- **Repeated verification failure.** The same verification has failed three times in a row despite different fix attempts. The diagnosis is wrong; ask the user.

Do NOT continue past a stop condition. Surface the issue, share what you understand, and wait for direction: when in doubt, stop.

### Phase 5 — Final validation

After all tasks are completed and verified:

1. Run `/fix-and-validate` once more for the full gauntlet (`iwyu-fix` + `lint-check` + `test`). This catches any drift accumulated across the task sequence — a refactor in task 3 may not have been picked up by task 5's narrower verification.
2. If `/fix-and-validate` passes cleanly, **suggest** to the user that `/wrap-branch` is the next step — but do NOT auto-invoke it. The user must explicitly type `/wrap-branch` to trigger the finalization flow.

If `/fix-and-validate` fails, surface the failures and stop. Do not loop on fixes silently — the user should know the branch isn't merge-ready.

## What this command does NOT do

`/execute-spec` does not write the spec — that's `/spec`'s job. It does not finalize the branch — that's `/wrap-branch`'s job (suggested but not auto-invoked). It does not interpret the spec creatively — if the spec is ambiguous, stop and ask.

The discipline is: spec drives, executing session follows, stop conditions are real, finalization is explicit.

## When NOT to use this command

Do not use `/execute-spec` if there is no spec. Run `/spec` first.

Do not use it for trivial changes whose execution doesn't benefit from task tracking. A typo fix doesn't need `/execute-spec`; just fix it.

Do not use it for spec documents that haven't been reviewed by the user. The spec is a contract; executing an unreviewed contract produces work that gets rejected at PR time.
