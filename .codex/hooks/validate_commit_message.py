#!/usr/bin/env python3
"""
Lancet2 PreToolUse hook: validate_commit_message

Intercepts `git commit -m "..."` (and `-am`, `--message`, `--message=`) bash
commands and validates the message against the rules in .claude/commit-style.json.

The rules are grounded in the project's actual .chglog/config.yml. The chglog
header pattern is `^(\\w*)\\:\\s(.*)$` with pattern_maps [Type, Subject] — type
and subject only, no scope group, exactly one space after the colon. The filter
list is exactly [feat, fix, perf, chore]; any other type silently fails to
appear in the generated CHANGELOG, so this hook rejects it pre-commit.

Hard blocks (exit 2):
  - subject does not match the chglog header pattern (no colon, wrong spacing,
    or scope syntax which chglog does not support)
  - type is not in allowed_types (would not appear in CHANGELOG)
  - type is not lowercase (chglog filter is case-sensitive)
  - subject too long, too short, ends with period, or starts with uppercase
  - substantive change (>30 staged lines) without a body

Skips entirely (exit 0 silently):
  - any non-commit bash command
  - `git commit` without -m (opens $EDITOR, out of scope; use a real
    prepare-commit-msg git hook for that surface)
  - `git commit --amend` reusing the prior message
  - merge and revert commits matching configured prefixes
  - fixup!/squash! markers (intentional; rebase will clean them up)

This hook is BLOCK-OR-PASS: it never prints advisory ("soft") warnings.
Earlier versions emitted soft-warns for non-imperative mood, body-line
length, body-shape, and trivial-pattern-with-large-diff. Those produced
noise on commits that passed anyway and were removed.
"""
import json
import os
import re
import shlex
import sys
from pathlib import Path


def load_config(project_dir: str) -> dict | None:
    """Load commit-style config from project root, then user home."""
    candidates = [
        Path(project_dir) / ".claude" / "commit-style.json",
        Path.home() / ".claude" / "commit-style.json",
    ]
    for path in candidates:
        if path.is_file():
            try:
                with path.open() as f:
                    return json.load(f)
            except (json.JSONDecodeError, OSError):
                # Don't block on a broken config; let the user see the error elsewhere.
                return None
    return None


# Match a HEREDOC command substitution of the canonical Claude form
# `$(cat <<'TAG'\n<body>\nTAG\n)`. The TAG may be single-quoted, double-quoted,
# or unquoted, and `<<-` (tab-stripping) is tolerated. The body capture is
# lazy so the closing `\nTAG\s*)` matches the *first* delimiter that follows
# the body, not the last. This is the most common pattern Claude uses to pass
# multi-line commit messages with embedded blank lines, and it is what bash
# would otherwise expand before invoking git — but this hook runs before
# bash, so we have to expand it ourselves to see the real message text.
_HEREDOC_RE = re.compile(
    r"^\$\(\s*cat\s+<<-?\s*['\"]?(?P<tag>\w+)['\"]?\s*\n(?P<body>.*?)\n(?P=tag)\s*\)$",
    re.DOTALL,
)


def _resolve_heredoc(value: str) -> str:
    """Return the inner body of a `$(cat <<'TAG' ... TAG)` HEREDOC, or `value`
    unchanged if no such pattern is found."""
    m = _HEREDOC_RE.match(value.strip())
    return m.group("body") if m else value


def _split_at_blank_lines(value: str) -> list[str]:
    """Split a multi-line message into paragraphs at blank-line boundaries.

    `git commit -m "subj\n\nbody1\n\nbody2"` is semantically equivalent to
    `git commit -m subj -m body1 -m body2` — git uses the first paragraph
    as the subject and remaining paragraphs as body sections. Returning a
    paragraph list lets the rest of this hook treat both invocation styles
    uniformly. A value with no blank lines returns a single-element list.
    """
    paragraphs: list[str] = []
    current: list[str] = []
    for line in value.split("\n"):
        if line.strip() == "":
            if current:
                paragraphs.append("\n".join(current).rstrip())
                current = []
        else:
            current.append(line)
    if current:
        paragraphs.append("\n".join(current).rstrip())
    return paragraphs or [value]


def extract_commit_messages(command: str) -> list[str] | None:
    """
    Parse a bash command string and return the list of -m / --message values
    if it is a `git commit` invocation. Returns None if not a git commit at all,
    or [] if it is `git commit` without -m (opens $EDITOR — out of scope).

    Multiple -m flags concatenate as separate paragraphs (git's documented
    behavior). The first -m is the subject; subsequent -m flags become the body.
    Compound shell commands like `cd /repo && git commit -m "..."` are handled
    by scanning for `git ... commit` and stopping at shell separators.
    """
    try:
        tokens = shlex.split(command, posix=True)
    except ValueError:
        # Unclosed quote or similar — let the actual git invocation surface it.
        return None
    if not tokens:
        return None

    messages: list[str] = []
    found_commit = False
    SHELL_SEP = ("&&", "||", ";", "|", "&")

    i = 0
    while i < len(tokens):
        if tokens[i] == "git":
            j = i + 1
            # Walk forward looking for `commit`, allowing `git -C dir commit` etc.
            while j < len(tokens) and tokens[j] not in SHELL_SEP:
                if tokens[j] == "commit":
                    found_commit = True
                    k = j + 1
                    while k < len(tokens) and tokens[k] not in SHELL_SEP:
                        ct = tokens[k]
                        if ct in ("-m", "--message"):
                            if k + 1 < len(tokens):
                                messages.append(tokens[k + 1])
                                k += 2
                                continue
                        elif ct.startswith("--message="):
                            messages.append(ct[len("--message="):])
                        elif ct in ("-am", "-ma"):
                            # Combined -a + -m
                            if k + 1 < len(tokens):
                                messages.append(tokens[k + 1])
                                k += 2
                                continue
                        k += 1
                    break
                j += 1
        i += 1

    if not found_commit:
        return None
    return messages


def is_auto_generated(subject: str, allowed_prefixes: list[str]) -> bool:
    return any(subject.startswith(p) for p in allowed_prefixes)


def is_fixup_or_squash(subject: str) -> bool:
    return subject.startswith("fixup!") or subject.startswith("squash!")


def validate_message(messages: list[str], cfg: dict, command: str = "") -> list[str]:
    """
    Validate the parsed commit messages against the config.
    Returns the list of hard-error strings; an empty list means the commit
    passes. The hook is block-or-pass; there are no soft-warning paths.
    """
    hard: list[str] = []

    if not messages:
        return hard

    subject = messages[0].strip()
    body_paragraphs = [m.strip() for m in messages[1:] if m.strip()]

    # Skip auto-generated merge/revert commits entirely.
    auto_prefixes = cfg.get("auto_generated_subject_prefixes_allowed", [])
    if is_auto_generated(subject, auto_prefixes):
        return hard

    # fixup!/squash! commits are intentional rebase fodder; pass through.
    if is_fixup_or_squash(subject):
        return hard

    # ── Header parsing — match chglog's pattern exactly ──────────────────
    pattern_str = cfg.get("header_pattern", r"^(?P<type>\w+): (?P<subject>.+)$")
    try:
        header_re = re.compile(pattern_str)
    except re.error as e:
        # Misconfigured pattern; fail open with a stderr note.
        print(f"⚠ commit-style: invalid header_pattern in config: {e}", file=sys.stderr)
        return hard

    m = header_re.match(subject)
    if not m:
        hard.append(
            f"Subject does not match the project header pattern from .chglog/config.yml.\n"
            f"  Pattern: {pattern_str}\n"
            f"  Got:     {subject!r}\n"
            f"  Required form: '<type>: <subject>' with exactly one space after the colon.\n"
            f"  Note: scopes like 'fix(caller): ...' are NOT supported by .chglog/config.yml.\n"
            f"  Examples of correct subjects:\n"
            + "\n".join(f"    {ex}" for ex in cfg.get("examples_good", [])[:3])
        )
        return hard

    # Pull captures by group name when available, falling back to position.
    # The default config uses named groups; a user-overridden pattern from
    # .chglog/config.yml ported verbatim might use numbered groups instead.
    try:
        ctype = m.group("type")
        subj_text = m.group("subject")
    except IndexError:
        groups = m.groups()
        ctype = groups[0] if len(groups) >= 1 else ""
        subj_text = groups[1] if len(groups) >= 2 else ""

    # ── Type validation ──────────────────────────────────────────────────
    allowed_types = cfg.get("allowed_types", [])
    if ctype != ctype.lower():
        hard.append(
            f"Type must be lowercase. Got '{ctype}'. The chglog filter is "
            f"case-sensitive, so '{ctype}' would silently fail to be grouped "
            f"into the changelog."
        )
    elif allowed_types and ctype not in allowed_types:
        sections = cfg.get("type_to_changelog_section", {})
        section_hints = ", ".join(
            f"{t}={sections.get(t, '?')!r}" for t in allowed_types
        )
        hard.append(
            f"Type '{ctype}' is not in .chglog/config.yml's filter list.\n"
            f"  Allowed types and their changelog sections: {section_hints}.\n"
            f"  A '{ctype}:' commit would parse but never appear in the generated "
            f"CHANGELOG.md. If this is a refactor, use 'chore:'. If it is a doc "
            f"or test or build change that genuinely should not appear in the "
            f"changelog, use 'chore:' or extend allowed_types in commit-style.json."
        )

    # ── Subject text validation ──────────────────────────────────────────
    max_len = cfg.get("max_subject_length", 72)
    if len(subject) > max_len:
        hard.append(
            f"Subject line is {len(subject)} characters; limit is {max_len}. "
            f"Shorten the subject and move detail to the body."
        )

    min_len = cfg.get("min_subject_length", 10)
    if len(subj_text) < min_len:
        hard.append(
            f"Subject text is {len(subj_text)} characters after the type; "
            f"minimum is {min_len}. Be specific about what changed."
        )

    if cfg.get("subject_must_not_end_with_period", True) and subj_text.endswith("."):
        hard.append("Subject must not end with a period.")

    if cfg.get("subject_must_be_lowercase_first_letter", False) and subj_text:
        first = subj_text[0]
        if first.isupper():
            hard.append(
                f"Subject text must start with a lowercase letter. Got '{first}' "
                f"({subj_text[:30]!r}...)."
            )

    # ── Body-required gate ───────────────────────────────────────────────
    # Per docs_dev/style/cpp_style.md § Git commit messages, trivial commits
    # (typos, exec-bit, formatting, dep bumps) ship with no body; substantial
    # commits use a two-section body. The only body-related rule the hook
    # enforces as a HARD block is: a non-trivial subject with a staged diff
    # exceeding `small_diff_line_threshold` lines must carry a body.
    #
    # Diff size unknowable (no git, not in a repo, etc.) → fail open.
    is_trivial = _matches_trivial_pattern(subject, cfg)
    threshold = cfg.get("small_diff_line_threshold", 30)

    if not body_paragraphs and not is_trivial:
        diff_size = _get_diff_size(command)
        if diff_size is not None and diff_size > threshold:
            hard.append(
                f"Substantive change ({diff_size} lines) without a body.\n"
                f"  Subject does not match a trivial pattern (typos, exec-bit, format,\n"
                f"  dep bumps), and the staged diff exceeds the {threshold}-line threshold\n"
                f"  beyond which the subject alone is unlikely to convey what changed.\n"
                f"  Add a body with a context paragraph and file-list bullets:\n"
                f"    git commit -m \"{subject}\" \\\n"
                f"      -m \"Context paragraph explaining why...\" \\\n"
                f"      -m \"- src/lancet/<layer>/file.cpp: what changed there\"\n"
                f"  If the change really is trivial, use a subject pattern that\n"
                f"  matches (e.g. 'chore: format ...', 'chore: bump ...', 'fix: typo ...')."
            )

    return hard


def _matches_trivial_pattern(subject: str, cfg: dict) -> bool:
    """True if subject matches any configured trivial-commit pattern."""
    patterns = cfg.get("trivial_commit_subject_patterns", [])
    for pat in patterns:
        try:
            if re.match(pat, subject, re.IGNORECASE):
                return True
        except re.error:
            continue  # skip malformed pattern silently
    return False


def _get_diff_size(command: str) -> int | None:
    """
    Return the staged diff size as insertions + deletions, suitable for the
    body-required threshold check.

    For `git commit -m`, this is `git diff --cached --shortstat` — the changes
    already in the index. For `git commit -am` (or -ma, --all), the `-a` flag
    has not yet executed at hook time; we add `git diff --shortstat` (working
    tree vs index) so the count reflects what `-a` would stage. Untracked
    files are not included in either case (they are not picked up by `-a`).

    Returns None if the diff size cannot be determined (no git binary, not
    in a repo, command timed out). Callers should treat None as "unknown,
    fail open" — a missing measurement should not block a commit.
    """
    import subprocess

    is_all_flag = bool(re.search(r"(?:^|\s)(?:-a|-am|-ma|--all)(?:\s|$)", command))

    def _shortstat_to_lines(stdout: str) -> int:
        """Parse output like ' 3 files changed, 42 insertions(+), 5 deletions(-)'."""
        if not stdout.strip():
            return 0
        ins = 0
        dels = 0
        m = re.search(r"(\d+)\s+insertions?\(\+\)", stdout)
        if m:
            ins = int(m.group(1))
        m = re.search(r"(\d+)\s+deletions?\(-\)", stdout)
        if m:
            dels = int(m.group(1))
        return ins + dels

    try:
        cached = subprocess.run(
            ["git", "diff", "--cached", "--shortstat"],
            capture_output=True, text=True, timeout=2, check=False,
        )
        if cached.returncode != 0:
            return None
        total = _shortstat_to_lines(cached.stdout)

        if is_all_flag:
            unstaged = subprocess.run(
                ["git", "diff", "--shortstat"],
                capture_output=True, text=True, timeout=2, check=False,
            )
            if unstaged.returncode == 0:
                total += _shortstat_to_lines(unstaged.stdout)
            # If unstaged probe fails, return what we have (cached only) rather than None.
        return total
    except (subprocess.SubprocessError, OSError, FileNotFoundError):
        return None


def main() -> int:
    try:
        payload = json.load(sys.stdin)
    except (json.JSONDecodeError, ValueError):
        # Hook protocol changed or input malformed; fail open.
        return 0

    command = payload.get("tool_input", {}).get("command", "") or ""
    if not command:
        return 0

    messages = extract_commit_messages(command)
    if messages is None:
        # Not a `git commit` command at all.
        return 0
    if not messages:
        # `git commit` with no -m — opens $EDITOR. Out of scope.
        return 0

    # Expand HEREDOC command substitutions and split multi-line single-flag
    # messages at blank lines so subject + body live as distinct entries
    # (matching the layout the rest of the hook expects from multi-`-m`
    # invocations). If a -m value contains a substitution we can't expand
    # (anything other than the canonical `$(cat <<TAG ... TAG)` HEREDOC),
    # we cannot see what git will actually commit; skip validation with a
    # soft warning rather than block on a literal we know is incorrect.
    #
    # Important: the substitution check applies to the RAW value, not the
    # expanded content. After successful HEREDOC expansion, the inner text
    # is the literal commit message — and that text legitimately contains
    # backticks (formatting identifier names) and may contain `$(` (citing
    # shell snippets in a commit body). Those are NOT shell substitutions
    # at the time git receives the -m value; bash already finished expanding
    # before invoking git.
    expanded: list[str] = []
    unresolvable = False
    for raw in messages:
        resolved = _resolve_heredoc(raw)
        if resolved == raw and ("$(" in raw or "`" in raw):
            unresolvable = True
            break
        expanded.extend(_split_at_blank_lines(resolved))
    if unresolvable:
        print(
            "⚠ commit-style: detected unresolvable shell substitution in a -m "
            "value; skipping validation. The actual git invocation will use the "
            "substituted text — please ensure it complies with Lancet2 commit "
            "style (see .claude/commit-style.json).",
            file=sys.stderr,
        )
        return 0
    messages = expanded

    project_dir = os.environ.get("CLAUDE_PROJECT_DIR", os.getcwd())
    cfg = load_config(project_dir)
    if cfg is None:
        # No config; nothing to enforce. Fail open.
        return 0

    hard = validate_message(messages, cfg, command=command)

    if hard:
        print(
            "BLOCKED: commit message does not satisfy Lancet2 commit-style rules.",
            file=sys.stderr,
        )
        for e in hard:
            print(f"  - {e}", file=sys.stderr)
        print("", file=sys.stderr)
        print(
            "Rules live in .claude/commit-style.json (grounded in .chglog/config.yml).",
            file=sys.stderr,
        )
        print(
            "Run /commit to have Claude compose a conformant message from the staged diff.",
            file=sys.stderr,
        )
        return 2

    return 0


if __name__ == "__main__":
    sys.exit(main())
