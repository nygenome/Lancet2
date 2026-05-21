#!/usr/bin/env bash
# Lancet2 PreToolUse hook: pre_commit_gate
#
# Verifies that `/fix-and-validate` has passed against the EXACT working-tree state
# about to be committed. Does NOT re-run validation — that would put the
# full Release-build + lint + IWYU + test cycle (3-7 min warm) on the
# critical commit path.
#
# How: `/fix-and-validate` writes a 4-field key=value record to
# `.claude/cache/validation-state.txt` on completion:
#   stash_hash=<git stash create hash at validation time>
#   head_sha=<HEAD commit sha at validation time>
#   timestamp=<unix seconds>
#   status=pass|fail
# This hook reads `stash_hash` and compares it with the about-to-be-
# committed working tree. Match -> silent exit 0 (~50ms). Mismatch or
# missing marker -> exit 2 with remediation pointing at /fix-and-validate.
# (statusline.sh consumes the `status` field independently.)

set -e

PROJECT_DIR="${CLAUDE_PROJECT_DIR:-$PWD}"
MARKER_FILE="$PROJECT_DIR/.claude/cache/validation-state.txt"

# ── Read JSON payload from stdin and extract the command ──────────────────
payload=$(cat)
command=$(printf '%s' "$payload" | python3 -c '
import json, sys
try:
    data = json.load(sys.stdin)
    print(data.get("tool_input", {}).get("command", ""))
except Exception:
    pass
' 2>/dev/null)

# ── Fast path: only fire on actual git commit invocations ─────────────────
case "$command" in
    *"git commit"*) ;;
    *) exit 0 ;;
esac

# ── --amend without -m / --message reuses the prior subject; skip ─────────
if printf '%s' "$command" | grep -qE -- '--amend' \
   && ! printf '%s' "$command" | grep -qE -- '(-m\b|-am\b|--message\b|--message=)'; then
    exit 0
fi

# ── Fast-path skip: nothing relevant staged ───────────────────────────────
relevant=$(git diff --cached --name-only -- \
    '*.cpp' '*.h' '*.hpp' '*.cc' '*.cxx' '*.hxx' \
    'CMakeLists.txt' '**/CMakeLists.txt' '*.cmake' '**/*.cmake' \
    'scripts/*.py' 'scripts/*.sh' \
    '.clang-tidy' '.clang-format' \
    2>/dev/null)
if [ -z "$relevant" ]; then
    exit 0
fi

# ── Compute current tree hash and compare against marker ──────────────────
# `git stash create` produces a commit object capturing the full working-tree
# state (staged + unstaged). We resolve to its tree (`commit^{tree}`) to get
# a deterministic content-only hash; the bare commit object includes a
# timestamp and would differ between calls on the same tree state. On a
# clean tree, stash-create returns empty — fall back to HEAD's tree.
stash_commit=$(git stash create 2>/dev/null || true)
if [ -n "$stash_commit" ]; then
    current=$(git rev-parse "${stash_commit}^{tree}" 2>/dev/null || echo "")
fi
if [ -z "$current" ]; then
    current=$(git rev-parse "HEAD^{tree}" 2>/dev/null || echo "no-tree")
fi

last=""
if [ -f "$MARKER_FILE" ]; then
    last=$(cat "$MARKER_FILE" 2>/dev/null)
fi

if [ -n "$last" ] && [ "$current" = "$last" ]; then
    # /fix-and-validate has been run against this exact tree state and passed.
    exit 0
fi

# ── Block: marker missing or stale ────────────────────────────────────────
{
    echo ""
    echo "─── pre_commit_gate ──────────────────────────────────────"
    echo "❌ /fix-and-validate has not been run against the current tree state."
    echo ""
    if [ -z "$last" ]; then
        echo "   No marker found at $MARKER_FILE."
    else
        echo "   Marker is stale:"
        echo "     current tree:  $current"
        echo "     last /fix-and-validate:   $last"
    fi
    echo ""
    echo "Run /fix-and-validate to validate (auto-fixes formatting + includes,"
    echo "surfaces lint and test issues). After /fix-and-validate passes, re-stage"
    echo "any files mutated by iwyu-fix, then re-attempt the commit."
    echo "──────────────────────────────────────────────────────────"
} >&2

exit 2
