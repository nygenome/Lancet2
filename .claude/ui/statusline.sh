#!/usr/bin/env bash
# Lancet2 statusline — printed at the top of Claude Code's UI on every prompt.
#
# Format: <branch>[|wt:<name>] [layers: <a>+<b>] [validation-check: <state>]
#
# Field semantics:
#   <branch>            : current git branch name. Always present (in a repo).
#   |wt:<name>          : appended only when CWD is a git linked worktree
#                         (under .worktrees/). The name is the basename of CWD.
#   [layers: a+b]       : shown only when git diff (working+staged) touches
#                         more than one src/lancet/<layer>/ directory. Layer
#                         names are sorted alphabetically and +-separated.
#   [validation-check]  : reads .claude/cache/validation-state.txt's status field.
#                         "pass" / "fail" / "stale" / absent.
#
# Fallback: if not in a git repo, prints "(not a git repo)" and exits 0.
#
# Performance: this runs on every prompt; must complete in <100ms.
# Implementation: single `git diff` call combined with `git diff --cached`,
# both filtered through a single grep. No subshell explosions.

set -e

# ── Detect git repo or fallback ──────────────────────────────────────────────
if ! git rev-parse --git-dir >/dev/null 2>&1; then
    echo "(not a git repo)"
    exit 0
fi

PROJECT_DIR="${CLAUDE_PROJECT_DIR:-$(git rev-parse --show-toplevel 2>/dev/null || echo "")}"

# ── Field 1: branch [|wt:name] ───────────────────────────────────────────────
branch=$(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo "?")

# Detect worktree: if the git common-dir differs from the git-dir, we're in
# a linked worktree. The worktree directory name is the basename of the
# current working directory (Lancet2 convention: .worktrees/<feature>).
git_dir=$(git rev-parse --git-dir 2>/dev/null || echo "")
common_dir=$(git rev-parse --git-common-dir 2>/dev/null || echo "")

field_branch="$branch"
if [ -n "$git_dir" ] && [ "$git_dir" != "$common_dir" ]; then
    wt_name=$(basename "$(pwd)")
    field_branch="${branch}|wt:${wt_name}"
fi

# ── Field 2: [layers: a+b] (only if multi-layer diff) ────────────────────────
field_layers=""
changed_files=$(
    {
        git -C "$PROJECT_DIR" diff --name-only HEAD 2>/dev/null
        git -C "$PROJECT_DIR" diff --cached --name-only 2>/dev/null
    } | sort -u
)

layers=$(
    echo "$changed_files" \
        | grep -oE '^src/lancet/[a-z]+/' \
        | sort -u \
        | sed 's|src/lancet/||;s|/||' \
        | paste -sd+ -
)
layer_count=$(echo "$layers" | tr '+' '\n' | grep -c .)

if [ "$layer_count" -gt 1 ]; then
    field_layers=" [layers: ${layers}]"
fi

# ── Field 3: [validation-check: <state>] ─────────────────────────────────────
field_check=""
state_file="$PROJECT_DIR/.claude/cache/validation-state.txt"
if [ -f "$state_file" ]; then
    status=$(grep -E '^status=' "$state_file" 2>/dev/null | cut -d'=' -f2- | head -1)
    case "$status" in
        pass) field_check=" [validation-check: pass]" ;;
        fail) field_check=" [validation-check: fail]" ;;
        *)    field_check=" [validation-check: stale]" ;;
    esac
fi

# ── Compose ──────────────────────────────────────────────────────────────────
echo "${field_branch}${field_layers}${field_check}"
