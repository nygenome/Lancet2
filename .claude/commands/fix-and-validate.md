---
description: Run the Lancet2 fix-mode validation suite (iwyu-fix + lint-check + test). Mirrors the pre-commit gate's scope but applies auto-fixes; pre-commit gate runs read-only equivalents.
allowed-tools: Bash
---

# /fix-and-validate — fix-mode validation

Run the developer-facing "fix everything I can, surface the rest" pass. The scope mirrors the pre-commit gate (`pre_commit_gate.sh`), but `/fix-and-validate` invokes the **fix** variants where they exist, while the pre-commit gate is strictly read-only. Use `/fix-and-validate` whenever you want the tooling to make the diff clean before you stage and commit; the gate is the safety net that runs at `git commit` time.

The relationship between the two is intentional. `/fix-and-validate` mutates the tree (clang-format, IWYU rewrites). The pre-commit gate does not — what you stage at `git commit` time is exactly what lands. If `/fix-and-validate` mutates files, you re-stage the modified files and commit; the gate then runs the read-only equivalents and confirms the tree is clean.

Use this command in three situations. First, after completing a meaningful chunk of work mid-session, when you want the formatter and IWYU to catch up before you stage. Second, before invoking `fresh-reviewer` or pushing to a remote, as a sanity check that the change is in a green state. Third, when the pre-commit gate has just failed and you want the auto-fixable parts handled in one pass.

## What `/fix-and-validate` actually runs (in order)

1. **`pixi run iwyu-fix`** — runs IWYU's `--fix` over the Release build tree, then re-applies clang-format to the rewritten files. Because IWYU's fix mode chains a format pass, this single step covers both the include-hygiene and the formatting concerns. Mutates files when violations are found.
2. **`pixi run lint-check`** — clang-tidy in read-only mode, against the same Release tree. Reports violations; does not auto-fix. Clang-tidy auto-fix has historically broken compilation in this project (see AGENTS.md "invoking clang-tidy" callout); resolve violations by editing source by hand.
3. **`pixi run test`** — runs `cmake-build-release/tests/TestLancet2` (the Release-tree test binary, with `LANCET_DEBUG_MODE` defined on the test target so `LANCET_ASSERT` keeps firing). Reports pass/fail.

All three steps share the `cmake-build-release` build tree, so the cost is one Release configure+compile (warm: ~30s; cold: ~3–7 min) plus the per-step analysis time.

Output streams live to stdout/stderr — no per-step log capture, no `tail -N` truncation. On failure the tool's full diagnostic output is already on screen; we just print a one-line "step X failed" marker so the boundary between steps is visible. A failed `iwyu-fix` halts the sequence (the lint and test passes won't be meaningful against a tree IWYU couldn't normalize); a failed `lint-check` or `test` reports but does not halt the next step.

After all three steps pass, `/fix-and-validate` writes a marker at `.claude/cache/validation-state.txt` containing `git stash create` output for the current working-tree state. The `pre_commit_gate.sh` hook reads this marker at `git commit` time and silently passes if the marker matches the about-to-be-committed tree — that's how the gate stays at ~50ms instead of re-running the 3-7 min validation cycle on every commit. If any step fails, the marker is removed (or never written) so subsequent commits stay blocked until `/fix-and-validate` passes again. See `.claude/hooks/pre_commit_gate.sh` for the verification mechanism.

## What this command does

```bash
echo "─── /fix-and-validate fix-mode validation ───────────────────────────"

# Stream live: each tool's full output goes to the terminal as it runs.
# This avoids capturing to a log and tailing N lines — IWYU and clang-tidy
# can produce hundreds of lines of recommendations/diagnostics, and a
# tail-based view risks truncating the actual failure. The trade-off is
# more output during successful runs (which is fine — pixi --quiet keeps
# the underlying tools' verbosity tame).

# Track aggregate success — we only write the gate marker if ALL steps
# pass. A partial success (lint clean, tests fail, etc.) leaves the
# previous marker stale or removed so the pre-commit gate stays red.
all_ok=1

echo ""
echo "1/3 iwyu-fix ..."
if ! pixi run --quiet iwyu-fix; then
  echo "────────────────────────────────────────────────────────────"
  echo "❌ iwyu-fix failed (see full output above). Halting."
  rm -f .claude/cache/validation-state.txt
  exit 1
fi
if [ -n "$(git status --porcelain -- '*.cpp' '*.h' '*.hpp' '*.cc' '*.cxx' '*.hxx')" ]; then
  echo "⚠ iwyu-fix rewrote includes/format — review and stage:"
  git status --short -- '*.cpp' '*.h' '*.hpp' '*.cc' '*.cxx' '*.hxx'
else
  echo "✓ iwyu-fix ok (no changes)"
fi

echo ""
echo "2/3 lint-check ..."
if ! pixi run --quiet lint-check; then
  echo "❌ clang-tidy violations above (resolve manually — no auto-fix)."
  all_ok=0
else
  echo "✓ lint-check ok"
fi

echo ""
echo "3/3 test ..."
if ! pixi run --quiet test; then
  echo "❌ tests failed (see Catch2 output above)."
  all_ok=0
else
  echo "✓ tests ok"
fi

if [ -n "$(git status --porcelain)" ]; then
  echo ""
  echo "⚠ uncommitted changes:"
  git status --short
fi

# Write/clear the gate marker based on aggregate success.
# We hash the TREE of git stash create rather than the commit object
# itself: stash-create produces a commit with a current-time timestamp,
# so the commit hash differs between calls on the same working tree —
# but the resolved tree is deterministic. The pre_commit_gate hook
# computes the same tree hash and compares.
mkdir -p .claude/cache
if [ "$all_ok" = "1" ]; then
  stash_commit=$(git stash create 2>/dev/null || true)
  if [ -n "$stash_commit" ]; then
    state=$(git rev-parse "${stash_commit}^{tree}" 2>/dev/null || echo "")
  fi
  if [ -z "$state" ]; then
    state=$(git rev-parse "HEAD^{tree}" 2>/dev/null || echo "no-tree")
  fi
  echo "$state" > .claude/cache/validation-state.txt
  echo "✓ pre-commit gate marker updated (.claude/cache/validation-state.txt)"
else
  rm -f .claude/cache/validation-state.txt
  echo "✗ pre-commit gate marker cleared (validation failed; commit will be blocked)"
fi

echo "────────────────────────────────────────────────────────────"
```

## When you want narrower test runs

`/fix-and-validate` runs the full gauntlet against the full test set. It's the right tool for "is this change merge-ready?" but it's overkill for "did my one change break the test I just wrote?". For narrower invocations, drive the test binary directly:

```bash
# Just the test you're iterating on
./cmake-build-release/tests/TestLancet2 "test name"

# Just one layer
./cmake-build-release/tests/TestLancet2 "[caller]"

# Reproduce a random-order failure (seed printed in the failing /fix-and-validate output)
./cmake-build-release/tests/TestLancet2 --rng-seed <reported-seed>
```

The `add-cpp-test` skill has the full Catch2 invocation reference (tag filtering, sections, generators, reporters, sharding, seed reproduction). Use `/fix-and-validate` for the validation gauntlet; use the binary directly for everything else.
