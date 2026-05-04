# /wrap-branch — Finish a development branch

Drive the cleanup workflow at the end of feature work. Runs a pre-check, then asks the user how to finalize the branch — merge locally, open a PR, keep the branch for later, or discard. Handles worktree cleanup based on the choice.

## Why this command exists

The end-of-branch workflow has four reasonable paths and each has its own steps (merge command, push + PR, cleanup vs no-cleanup of worktrees). Doing it manually means recalling the right git incantation under the conditions of "I just finished work and want to ship it." `/wrap-branch` makes the choice explicit and the steps consistent.

## Procedure

### Phase 0 — Sanity checks

Refuse to run if the current branch is `main`. The workflow assumes a feature branch:

```bash
branch=$(git rev-parse --abbrev-ref HEAD)
if [ "$branch" = "main" ]; then
    echo "ERROR: /wrap-branch is for feature branches; you're on main."
    exit 1
fi
```

### Phase 1 — Pre-check

Run `pixi run test` as a quick sanity check that the branch is in a working state:

```bash
pixi run test
```

This is a narrower check than `/fix-and-validate` (which adds lint and IWYU); the goal here is "do tests pass" — a faster gate before asking the user about finalization.

If `pixi run test` fails, surface the failure and **STOP**. Do not proceed to Phase 2; the user needs to fix tests before deciding how to finalize. After they fix, suggest running `/fix-and-validate` for the full pre-merge gauntlet, then re-running `/wrap-branch`.

If `pixi run test` passes, recommend running `/fix-and-validate` next for the full validation (`iwyu-fix` + `lint-check` + `test`). The decision to actually run it is the user's — `/wrap-branch` will continue to Phase 2 either way, since the pre-check has cleared the basic gate.

### Phase 2 — Choose the finalization path

Use `AskUserQuestion` to present four options:

```
question: How do you want to finalize this branch?
options:
  - Merge to main locally
  - Push branch and create PR
  - Keep branch (handle later)
  - Discard this work
```

The user picks one. Each option triggers a different downstream flow.

### Phase 3a — Merge to main locally

```bash
# Capture the current branch name and commit count
branch=$(git branch --show-current)
ahead=$(git rev-list --count main..HEAD)

echo "About to merge $branch ($ahead commit(s) ahead of main) into main."

# Switch to main, pull, merge
git checkout main
git pull --ff-only
git merge --no-ff "$branch" -m "Merge branch '$branch' into main"

# Push to remote
git push origin main
```

After the merge, ask: "Delete the local branch `$branch`?" If yes, `git branch -d "$branch"`. If the branch had a worktree at `.worktrees/$branch/`, also clean up the worktree:

```bash
if [ -d ".worktrees/$branch" ]; then
    git worktree remove ".worktrees/$branch"
fi
```

### Phase 3b — Push branch and create PR

```bash
branch=$(git branch --show-current)

# Push to remote (set upstream if needed)
git push -u origin "$branch"
```

Compose a PR body from the diff. Use a mini-version of `/commit`'s logic: walk the changed files, group by purpose, produce a brief summary. Then:

```bash
# Title from the most recent commit subject
title=$(git log -1 --pretty=format:%s)

# Body from the composed summary, written to a temp file
gh pr create --title "$title" --body-file /tmp/pr_body.md
```

If `gh` is not available, surface the URL where the user can manually create the PR (`https://github.com/<owner>/<repo>/pull/new/$branch`) and copy the composed body to clipboard or print it for manual paste.

Do NOT clean up the worktree — the user is still iterating on the branch via the PR. Cleanup happens after the PR merges, separately.

### Phase 3c — Keep branch (handle later)

No action. Tell the user the branch is preserved as-is. Do NOT clean up the worktree.

### Phase 3d — Discard this work

This is the destructive path. Require typed confirmation using the branch name:

```
This will permanently delete:
  - The branch <branch-name>
  - The worktree at .worktrees/<branch-name>/ (if it exists)
  - All commits on this branch that aren't on main

Type the branch name (<branch-name>) to confirm:
```

Wait for the user to type the exact branch name. If they type anything else (including "discard", "yes", or a typo'd version of the branch name), abort the discard and tell them no action was taken. The branch-name confirmation is deliberately friction-heavy because the action is irreversible — auto-completing on "yes" is too easy to do by accident.

If they type the exact branch name:

```bash
branch=$(git branch --show-current)

# Switch off the branch first
git checkout main

# Delete the branch (force-delete because it's not merged)
git branch -D "$branch"

# Remove the worktree if it existed
if [ -d ".worktrees/$branch" ]; then
    git worktree remove --force ".worktrees/$branch"
fi
```

After the discard, confirm to the user what was deleted.

## Worktree handling summary

| Option | Worktree cleanup |
|---|---|
| 1. Merge locally | Yes (after merge succeeds) |
| 2. Push + PR | No (user is still iterating) |
| 3. Keep branch | No |
| 4. Discard | Yes (part of the discard) |

This matches the following pattern: cleanup on options that *finalize* the work (merge, discard); no cleanup on options that *defer* it (PR-in-progress, keep-for-later).

## When NOT to use this command

Do not use `/wrap-branch` if the branch you're on is `main`. The command checks for this and aborts.

Do not use it before the work is actually done. The pre-check catches "tests fail" but doesn't catch "the feature is half-implemented." Run `/fix-and-validate` first if you want a more comprehensive pre-flight; run the relevant `vcf-validator` or `assembly-and-calling-expert` review if the change is in their wheelhouse; only then `/wrap-branch`.
