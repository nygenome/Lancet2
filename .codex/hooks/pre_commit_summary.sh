#!/usr/bin/env bash
# Lancet2 PreToolUse hook: pre_commit_summary
#
# Prints a one-paragraph summary of the staged change before a `git commit`
# executes. Informational only — never blocks. The summary appears in stderr
# (which Claude sees), giving the agent and user one last chance to notice
# "wait, that does not look right" before the commit lands.
#
# Why a hook rather than a slash command: this needs to fire automatically
# on every commit attempt without the user having to remember to run a
# review step. The cost is small (one Python invocation) and the value is
# real (catches the "I committed the wrong staged set" mistake that is
# common in solo workflows).
#
# Why PreToolUse rather than PostToolUse: PreToolUse fires before the commit
# executes, so the agent has a chance to abort if the summary reveals an
# error.
#
# Why this is separate from validate_commit_message: that hook validates the
# message against commit-style.json. This hook summarizes the diff. Two
# different concerns; keeping them separate makes each easier to maintain.
#
# Architecture: this script is a presentation layer. The diff classification
# itself lives in .claude/hooks/lib/diff_analysis.py — a shared Python module
# also imported by agent-memory commit-bundling, external-interface-changes interview,
# and any other Lancet2 hook/skill that needs to reason about staged
# changes. Keeping the heuristics in one place avoids the diverged-
# duplicate-implementations bug class.

set -e

# ── Read JSON payload from stdin and forward to the analyzer ──────────────
# diff_analysis.py reads the same payload format and emits a single JSON
# object on stdout. We pipe stdin → analyzer → jq-equivalent extraction.
payload=$(cat)

# Fast path: only fire on actual git commit invocations.
command=$(printf '%s' "$payload" | python3 -c '
import json, sys
try:
    data = json.load(sys.stdin)
    print(data.get("tool_input", {}).get("command", ""))
except Exception:
    pass
' 2>/dev/null)

case "$command" in
    *"git commit"*) ;;
    *) exit 0 ;;
esac

# ── Run the analyzer ──────────────────────────────────────────────────────
# Fail open: if the analyzer fails for any reason, exit silently rather
# than blocking the commit on a missing measurement.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
analysis=$(printf '%s' "$payload" | python3 "${SCRIPT_DIR}/lib/diff_analysis.py" 2>/dev/null)
if [ -z "$analysis" ]; then
    exit 0
fi

# ── Extract fields from analysis JSON via Python (one-shot) ───────────────
# We could call jq, but jq is not guaranteed to be installed in every
# Lancet2 dev environment. Python is part of the pixi env and is already
# used by other hooks.
read -r is_git_repo file_count has_all_flag has_schema has_nolints \
        cross_layer needs_test_warn needs_docs_warn \
        staged_short unstaged_short layers kinds \
        files_head schema_lines nolint_count <<< "$(printf '%s' "$analysis" | python3 -c '
import json, sys
d = json.load(sys.stdin)

def s(v):  # safe scalar
    if isinstance(v, bool):
        return "1" if v else "0"
    return str(v)

def joined(xs):
    return ",".join(xs) if xs else "-"

def head(xs, n=8):
    if not xs:
        return "-"
    sliced = xs[:n]
    if len(xs) > n:
        sliced.append(f"... ({len(xs) - n} more)")
    return "|".join(sliced)

print(s(d.get("is_git_repo")), end=" ")
print(s(d.get("file_count")), end=" ")
print(s(d.get("has_all_flag")), end=" ")
print(s(d.get("has_schema_changes")), end=" ")
print(s(d.get("has_nolint_additions")), end=" ")
print(s(d.get("is_cross_layer")), end=" ")
print(s(d.get("needs_test_warning")), end=" ")
print(s(d.get("needs_docs_warning")), end=" ")
# Multi-word fields use a single-quote-wrapped form to survive read
print(repr(d.get("staged_shortstat", "")), end=" ")
print(repr(d.get("unstaged_shortstat", "")), end=" ")
print(repr(joined(d.get("layers_touched", []))), end=" ")
print(repr(joined(d.get("kinds", []))), end=" ")
print(repr(head(d.get("files", []))), end=" ")
print(repr(head(d.get("schema_lines_added", []), 4)), end=" ")
print(s(len(d.get("nolint_additions", []))))
' 2>/dev/null)"

# Bail if not in a git repo.
if [ "$is_git_repo" != "1" ]; then
    exit 0
fi

# Strip Python repr() quotes from the multi-word fields.
strip_quotes() {
    local v="$1"
    # Remove leading/trailing single quotes if present
    v="${v#\'}"
    v="${v%\'}"
    # Unescape \\n if any
    printf '%s' "$v"
}
staged_short=$(strip_quotes "$staged_short")
unstaged_short=$(strip_quotes "$unstaged_short")
layers=$(strip_quotes "$layers")
kinds=$(strip_quotes "$kinds")
files_head=$(strip_quotes "$files_head")
schema_lines=$(strip_quotes "$schema_lines")

# ── Emit the summary ──────────────────────────────────────────────────────
{
    echo "──────────────────────────────────────────────────────────────────────"
    echo "ℹ commit-summary: about to commit"
    echo ""

    if [ "$has_all_flag" = "1" ]; then
        echo "  Files (staged + -a unstaged): $file_count"
        [ -n "$staged_short" ] && echo "  Staged:  $staged_short"
        [ -n "$unstaged_short" ] && echo "  -a will add: $unstaged_short"
    else
        echo "  Files (staged): $file_count"
        [ -n "$staged_short" ] && echo " $staged_short"
    fi

    if [ "$layers" != "-" ]; then
        echo "  Lancet2 layers touched: ${layers//,/ }"
    fi

    if [ "$kinds" != "-" ]; then
        echo "  Kinds: $kinds"
    fi

    echo ""
    echo "  Paths:"
    if [ "$files_head" != "-" ]; then
        printf '%s\n' "$files_head" | tr '|' '\n' | sed 's/^/    /'
    fi

    # Heuristic warnings.

    if [ "$needs_docs_warn" = "1" ]; then
        echo ""
        echo "  ⚠ Multiple layers touched without a docs change. Consider"
        echo "    whether docs/guides/architecture.md or layer guides need an update."
    fi

    if [ "$needs_test_warn" = "1" ]; then
        echo ""
        echo "  ⚠ Source change without a test/benchmark change. If this is a fix"
        echo "    or feature (not a refactor or chore), a regression test is expected."
    fi

    if [ "$has_schema" = "1" ]; then
        echo ""
        echo "  ⚠ VCF header line(s) appear to be added/changed:"
        if [ "$schema_lines" != "-" ]; then
            printf '%s\n' "$schema_lines" | tr '|' '\n' | sed 's/^/      /'
        fi
        echo "    Run vcf-validator subagent before merging if this is a schema change."
    fi

    if [ "$has_nolints" = "1" ]; then
        echo ""
        echo "  ⚠ ${nolint_count} NOLINT suppression line(s) added in this commit."
        echo "    Per AGENTS.md NOLINT discipline, each new suppression must be"
        echo "    disclosed in the summary of work with a rationale. fresh-reviewer"
        echo "    will flag undisclosed suppressions at PR time."
    fi

    # Agent memory bundling note. If memory files were touched in this
    # commit, surface that the user is implicitly bundling memory updates
    # with this work commit (the design decision in AGENTS.md / cost-model).
    memory_count=$(printf '%s' "$analysis" | python3 -c '
import json, sys
d = json.load(sys.stdin)
print(len(d.get("agent_memory_files_touched", [])))
' 2>/dev/null)
    if [ -n "$memory_count" ] && [ "$memory_count" -gt 0 ]; then
        echo ""
        echo "  ℹ ${memory_count} agent-memory file(s) included in this commit."
        echo "    These are bundled with the primary change per the project"
        echo "    convention; chglog filters on commit type rather than file"
        echo "    path, so the memory edits inherit the host commit's type"
        echo "    (feat/fix/perf/chore) and do not produce a separate"
        echo "    CHANGELOG entry."
    fi

    echo "──────────────────────────────────────────────────────────────────────"
} >&2

# Always exit 0. Hook is informational; abort decisions are made by the
# agent or user based on the summary content.
exit 0
