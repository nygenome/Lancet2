#!/usr/bin/env python3
"""
Lancet2 PreToolUse hook: block_protected_paths

Blocks edits to paths that should never be modified by the agent without
explicit human review. This is the first line of defense against
agentic-edit accidents in build, dependency, CI, and binary-asset
directories.

The list of protected paths is read from .claude/protected_paths.txt at
the project root. That file is the single source of truth; both this
hook AND settings.json's permissions.deny block reference it. The
/audit-bundle slash command verifies the two consumers stay in sync.

Protocol:
- Reads tool input as JSON on stdin.
- Resolves the proposed `file_path` relative to `project_root()` and substring-
  matches each protected pattern against the relative form. Edits inside a
  `.worktrees/<branch>/` checkout — the project's recommended workflow per
  `AGENTS.md` — therefore do not false-positive on the `.worktrees/` pattern,
  while edits from the main checkout into a sibling worktree still match
  (`.worktrees/<branch>/...` is the relative form there). Paths outside
  `project_root()` fall through to the absolute form so the hook stays
  conservative on unexpected inputs.
- Exits 2 (deny) if a pattern matches.
- Exits 0 (allow) otherwise, including on protocol errors (fail-open).
"""
import json
import os
import sys
from pathlib import Path


def project_root() -> Path:
    """Resolve the Lancet2 project root.

    Prefer CLAUDE_PROJECT_DIR (set by Claude Code at hook execution
    time). Fall back to walking up from this script's location, then to
    the current working directory.
    """
    env_dir = os.environ.get("CLAUDE_PROJECT_DIR", "").strip()
    if env_dir:
        return Path(env_dir)
    script_path = Path(__file__).resolve()
    # .claude/hooks/block_protected_paths.py → project root is parents[2]
    return script_path.parents[2]


def load_protected_patterns() -> list[str]:
    """Read patterns from .claude/protected_paths.txt.

    Returns an empty list if the file is missing — the hook then
    fails-open (allows everything), which is the correct posture for a
    missing source-of-truth file. The /audit-bundle command will catch
    the missing file separately.
    """
    candidates = [
        project_root() / ".claude" / "protected_paths.txt",
        Path.cwd() / ".claude" / "protected_paths.txt",
    ]
    for candidate in candidates:
        if candidate.is_file():
            return _parse_patterns_file(candidate)
    return []


def _parse_patterns_file(path: Path) -> list[str]:
    patterns: list[str] = []
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        patterns.append(stripped)
    return patterns


def main() -> int:
    try:
        payload = json.load(sys.stdin)
    except (json.JSONDecodeError, ValueError):
        # If the hook protocol changes, do not silently block edits.
        return 0

    file_path = payload.get("tool_input", {}).get("file_path", "") or ""
    if not file_path:
        return 0

    protected = load_protected_patterns()
    if not protected:
        # Source-of-truth file missing or empty — fail open. The audit
        # command catches the missing file as a separate concern.
        return 0

    # Match against the path RELATIVE to project_root so that operating *inside*
    # a `.worktrees/<branch>/` checkout doesn’t false-positive on the
    # `.worktrees/` pattern. That pattern is meant to block reaching across
    # worktrees from the MAIN checkout, not to block edits to a worktree’s
    # own files. When the current project root IS a worktree, files under it
    # have no `.worktrees/` segment in their relative form — that is the
    # principled match. Files outside project_root fall through to the
    # absolute path so the hook stays conservative on unexpected inputs.
    root = project_root().resolve()
    try:
        rel_path = str(Path(file_path).resolve().relative_to(root))
    except ValueError:
        rel_path = file_path

    for pattern in protected:
        if pattern in rel_path:
            print(
                f"BLOCKED: {file_path} matches protected pattern "
                f"'{pattern}'.\n"
                f"  This path is excluded from agent edits because "
                f"mistakes there\n"
                f"  are subtle and propagate broadly. The pattern is "
                f"declared in\n"
                f"  .claude/protected_paths.txt — edit that file (and "
                f"the matching\n"
                f"  permissions.deny rule in settings.json) if you "
                f"want to remove\n"
                f"  the protection.",
                file=sys.stderr,
            )
            return 2

    return 0


if __name__ == "__main__":
    sys.exit(main())
