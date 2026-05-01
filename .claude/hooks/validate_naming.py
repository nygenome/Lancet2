#!/usr/bin/env python3
"""
Lancet2 PreToolUse hook: validate_naming

Catches common Lancet2 naming-convention violations at write time, before
clang-tidy runs at Stop. The full source of truth is .clang-tidy; this hook
catches only the highest-confidence patterns to give fast feedback.

What it catches (only inside src/lancet/<layer>/ files):
  - Member variables not prefixed with 'm' (e.g., `int counter_;` should be `mCounter`)
  - using namespace std; in headers (forbidden)
  - std::format or std::print usage (project uses spdlog/fmtlib instead)
  - bare assert() instead of LANCET_ASSERT
  - bare `// NOLINT` (inline) suppressions — only scoped NOLINTNEXTLINE /
    NOLINTBEGIN+NOLINTEND forms are allowed
  - Trailing underscore convention (some C++ codebases use trailing _, Lancet2 uses leading m)
  - NEW (Phase 5+): scoped `NOLINTNEXTLINE`/`NOLINTBEGIN` additions without
    a `//` rationale comment on the line(s) immediately above the directive
    are hard-blocked. Every suppression must be accompanied by an explicit
    WHY written on the lines directly preceding the NOLINT directive — the
    inline `-- <reason>` form is no longer accepted, regardless of brevity.
    `NOLINTEND` is exempt (its rationale lives on the matching NOLINTBEGIN).

What it deliberately does NOT check:
  - Function/class case (clang-tidy catches reliably)
  - Local variable case (too many false positives at write time)
  - Constant case (depends on storage class context)

═════════════════════════════════════════════════════════════════════════
Hook calibration philosophy — read this BEFORE adding a new pattern
═════════════════════════════════════════════════════════════════════════

This hook (and every Lancet2 PreToolUse hook) is calibrated to a
near-zero-false-positive bar. The reason is empirical: if a hook produces
false positives at any meaningful rate, both the human and Claude learn
to ignore the hook output. Once an enforcement mechanism is being
routinely overridden, it has negative value — it adds noise to the
session's stderr without preventing anything.

The bar is operationalized as three rules:

  1. Hard blocks (exit 2) require near-zero false positives. A pattern
     earns hard-block status only if (a) the violation is explicit and
     unambiguous in the source text, (b) reviewers have flagged it
     repeatedly in the past, and (c) the cost of an erroneous block is
     trivially recoverable (the user adds an inline comment override
     and re-runs).

  2. Soft warnings (exit 0 with stderr) are appropriate when the
     pattern is high-signal but ambiguous in some real cases. Soft
     warnings MUST be acknowledged: when Claude sees a soft warning,
     the next assistant turn should briefly reference what was warned
     about (one phrase is enough — "noted, the m-prefix warning is on
     a struct field, leaving it as-is intentional"). The acknowledgment
     keeps the warning's signal alive instead of letting it become
     wallpaper.

  3. Calibration is validated quarterly. The /audit-bundle command
     reviews each hook's stderr output from the prior quarter (logged
     via the InstructionsLoaded observability hooks) and checks that
     hard-block invocations were genuine and soft warnings were
     either fixed or acknowledged. A hook with a >5% false-positive
     rate is downgraded to soft warning or removed.

When proposing a new check, classify it against this bar before
choosing the exit code. If you're not sure whether the pattern is FP-
near-zero, default to soft warning — promoting later is easy; demoting
a hard block after Claude has been trained to ignore it is much
harder.

The doctrine is also documented in .claude/hooks/README.md — keep the
two in sync. The hook calibration philosophy is project-wide, not
per-hook.

The hook runs in well under 100ms.
"""
import json
import re
import sys

# ── Hard blocks: usage that violates explicit project conventions ───────────

HARD_PATTERNS = [
    (
        re.compile(r"\bstd::format\s*\("),
        "std::format is forbidden — use spdlog or fmtlib (fmt::format) instead. "
        "See .clang-tidy: -modernize-use-std-format is disabled because the project standardizes on fmtlib.",
    ),
    (
        re.compile(r"\bstd::print\s*\("),
        "std::print is forbidden — use spdlog (SPDLOG_INFO, SPDLOG_DEBUG, etc.) instead.",
    ),
    (
        re.compile(r"^\s*using\s+namespace\s+std\s*;", re.MULTILINE),
        "`using namespace std;` is forbidden in headers (and discouraged in .cpp). Qualify uses explicitly.",
    ),
    (
        re.compile(r"(?<!LANCET_)\bassert\s*\("),
        "Bare assert() is forbidden — use LANCET_ASSERT (defined in src/lancet/base/assert.h). "
        ".clang-tidy is configured to recognize LANCET_ASSERT for side-effect detection.",
    ),
    # Bare NOLINT (inline form, with or without a check name) is forbidden;
    # it pollutes the statement line and hurts scannability. Scoped forms
    # only: NOLINTNEXTLINE(check) for one line, NOLINTBEGIN(check)/NOLINTEND(check)
    # for blocks. Regex matches NOLINT as a whole token (negative lookahead
    # excludes NOLINTNEXTLINE/BEGIN/END).
    (
        re.compile(r"//\s*NOLINT(?!NEXTLINE|BEGIN|END)\b"),
        "Bare `// NOLINT` (inline form) is forbidden. Use a scoped form with\n"
        "  the rationale on a comment line ABOVE the directive:\n"
        "    // <WHY this suppression is justified>\n"
        "    // NOLINTNEXTLINE(check-name)\n"
        "    <line that needs the suppression>\n"
        "  or, for a multi-line block:\n"
        "    // <WHY this block needs to be silenced>\n"
        "    // NOLINTBEGIN(check-name)\n"
        "    <block>\n"
        "    // NOLINTEND(check-name)\n"
        "  The bare inline form pollutes the statement line and makes\n"
        "  scannability hard. The choice of scoped form is independent of\n"
        "  whether the suppression is justified — use the scoped form\n"
        "  regardless. Also, Claude must disclose any newly added suppression\n"
        "  in its summary of work; see AGENTS.md `NOLINT discipline`.",
    ),
]
# The Phase-4+ inline `--` rationale check has been replaced by a line-pair
# walk in `check_nolint_rationale_above` below — the rationale must live on
# the comment line(s) immediately above the directive, regardless of length.

# ── Soft warnings: suspected violations of the mPascalCase convention ──────

# Looks for class/struct member variable declarations without the m prefix.
# Heuristic only — false positives are tolerable since this is non-blocking.
# Matches lines like: `int counter;` `std::string mName;` etc. inside what looks
# like a class body. We only flag if the identifier starts with lowercase AND
# is not preceded by 'm' followed by uppercase.
SUSPECT_MEMBER_RE = re.compile(
    r"^\s+(?:[A-Z]\w*::)?[\w<>:,*&\s]+?\s+([a-z]\w*)\s*[;{=]",
    re.MULTILINE,
)

# Files where the soft warning makes sense (class-heavy headers, .cpp files).
# The "src/lancet/" pattern matches both relative paths (like
# "src/lancet/caller/genotyper.cpp") and absolute paths (like
# "/repo/src/lancet/caller/genotyper.cpp"); Claude Code typically passes
# absolute paths but the matcher does not require it.
_LANCET_SOURCE_PREFIX = re.compile(r"(?:^|/)src/lancet/")


def is_lancet_source(file_path: str) -> bool:
    return bool(_LANCET_SOURCE_PREFIX.search(file_path)) and (
        file_path.endswith(".cpp")
        or file_path.endswith(".cc")
        or file_path.endswith(".h")
        or file_path.endswith(".hpp")
    )


def is_lancet_header(file_path: str) -> bool:
    return bool(_LANCET_SOURCE_PREFIX.search(file_path)) and (
        file_path.endswith(".h") or file_path.endswith(".hpp")
    )


# Detector for `// NOLINTNEXTLINE(...)` or `// NOLINTBEGIN(...)` lines.
# Anchored to start-of-line (with optional indent) so it doesn't false-fire on
# the literal token appearing inside a string or longer comment.
_NOLINT_DIRECTIVE_RE = re.compile(
    r"^\s*//\s*NOLINT(?:NEXTLINE|BEGIN)\b",
)

# Detector for any `//` comment line. The rationale-above rule accepts any
# non-empty text after the `//` — we deliberately don't try to validate the
# rationale content; reviewers and the disclosure rule cover quality.
_COMMENT_LINE_RE = re.compile(r"^\s*//\s*\S")


def check_nolint_rationale_above(new_text: str) -> list[str]:
    """Return a list of human-readable error messages for NOLINTNEXTLINE /
    NOLINTBEGIN directives that lack a `//` comment on the line(s) immediately
    above them.

    Walks new_text line-by-line. For each directive line, looks at the
    previous line; if it is a non-empty `//` comment, the directive is
    rationaled. Otherwise it is flagged.

    Note: the inline `// NOLINTNEXTLINE(...) -- <reason>` form is NOT accepted
    by this check (the inline `--` clause is no longer the rationale token).
    """
    lines = new_text.splitlines()
    errors: list[str] = []
    for idx, line in enumerate(lines):
        if not _NOLINT_DIRECTIVE_RE.match(line):
            continue
        # Look at the previous line (line above the directive).
        if idx == 0:
            errors.append(
                f"line 1: `{line.strip()}` has no rationale above it (file starts with the directive)."
            )
            continue
        prev = lines[idx - 1]
        if _COMMENT_LINE_RE.match(prev):
            continue
        errors.append(
            f"line {idx + 1}: `{line.strip()}` lacks a `//` rationale comment on the line immediately above."
        )
    return errors


def collect_new_text(payload: dict) -> str:
    tool_name = payload.get("tool_name", "")
    tool_input = payload.get("tool_input", {})
    if tool_name == "Write":
        return tool_input.get("content", "") or ""
    if tool_name == "Edit":
        return tool_input.get("new_string", "") or ""
    if tool_name == "MultiEdit":
        edits = tool_input.get("edits", []) or []
        return "\n".join(e.get("new_string", "") for e in edits)
    return ""


def main() -> int:
    try:
        payload = json.load(sys.stdin)
    except (json.JSONDecodeError, ValueError):
        return 0

    file_path = payload.get("tool_input", {}).get("file_path", "") or ""
    if not is_lancet_source(file_path):
        return 0

    new_text = collect_new_text(payload)
    if not new_text:
        return 0

    # Hard blocks
    hard_violations: list[str] = []
    for pattern, message in HARD_PATTERNS:
        # `using namespace std` is only blocked in headers
        if "using namespace std" in message and not is_lancet_header(file_path):
            continue
        if pattern.search(new_text):
            hard_violations.append(f"  - {message}")

    nolint_errors = check_nolint_rationale_above(new_text)
    if nolint_errors:
        rationale_msg = (
            "NOLINT suppression added without a rationale on the line(s) ABOVE the directive.\n"
            "  Required form (Phase 5+):\n"
            "      // <WHY this suppression is justified>\n"
            "      // NOLINTNEXTLINE(check-name)\n"
            "      <line that needs the suppression>\n"
            "  The inline `-- <reason>` form is no longer accepted, regardless of brevity.\n"
            "  See AGENTS.md `NOLINT discipline` for the convention.\n"
            "  Offending directives:\n"
            + "\n".join(f"      {err}" for err in nolint_errors)
        )
        hard_violations.append(f"  - {rationale_msg}")

    if hard_violations:
        print(f"BLOCKED: naming/convention violation in {file_path}", file=sys.stderr)
        print("\n".join(hard_violations), file=sys.stderr)
        return 2

    # Soft warnings — print but allow
    if is_lancet_header(file_path):
        suspects: list[str] = []
        for m in SUSPECT_MEMBER_RE.finditer(new_text):
            identifier = m.group(1)
            # Skip if it looks like a function/method (followed by `(`)
            full_match = m.group(0)
            if "(" in full_match:
                continue
            # Skip common false positives: typedef/using, return, etc.
            if any(kw in full_match for kw in ("return ", "typedef ", "using ", "if ", "for ", "while ", "switch ")):
                continue
            # Skip identifiers that already follow mPascalCase
            if identifier.startswith("m") and len(identifier) > 1 and identifier[1].isupper():
                continue
            suspects.append(identifier)

        if suspects:
            unique = sorted(set(suspects))[:5]  # cap output
            print(
                f"⚠ Possible member-variable naming issues in {file_path}: {', '.join(unique)}\n"
                f"  Lancet2 convention is mPascalCase (e.g., mKmer, mEdges, mCurrK).\n"
                f"  This is a soft warning; clang-tidy will catch real violations at Stop.",
                file=sys.stderr,
            )
            # exit 0 — non-blocking

    return 0


if __name__ == "__main__":
    sys.exit(main())
