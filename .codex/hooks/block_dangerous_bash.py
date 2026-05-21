#!/usr/bin/env python3
"""
Lancet2 PreToolUse hook: block_dangerous_bash

Validates bash commands before they execute. Permission-rule wildcards
in settings.json do not see compound commands like `cd && rm -rf`, so
this hook is the real defense.

Reads tool input as JSON on stdin. Exits 2 (deny) if the command would
cause irrecoverable damage. Exits 0 (allow) otherwise — for deletion
commands, the hook also prints a safety verdict to stderr so the
user-approval prompt that follows surfaces the categorization.

Three layers of policy:

  1. BLANKET-BANNED PATTERNS — `git push`, `git reset --hard origin`,
     `dd if=`, `mkfs.<fs>`, the fork-bomb literal, `chmod -R 777 /`,
     `chown -R`, and `xargs rm`. These have no safe variant; the
     maintainer handles them by hand outside Claude Code. `xargs rm`
     is in this list because the path arguments come from stdin (a
     previous command's output) and cannot be statically validated by
     this hook.

  2. PATH-VALIDATED DELETIONS — `rm`, `rmdir`, `shred`, `find … -delete`,
     `find … -exec rm`. These are NOT blanket-banned. Instead, the hook
     extracts every path argument and validates it. A path passes
     validation iff:
        - it is absolute (starts with `/`),
        - it contains no shell metacharacters that would expand to a
          different value at execution time (`$`, backtick, `*`, `?`,
          `[`, `]`, `{`, `}`),
        - it resolves cleanly (no `..` escape, no symlink to outside
          the allowlist),
        - it is not the bare allowlist root (e.g. refuse `rm /tmp`
          itself, allow `rm /tmp/foo`),
        - the resolved path lands under one of:
            * `/tmp/`, `/var/tmp/`, `/scratch/` (POSIX/system temp), OR
            * inside `CLAUDE_PROJECT_DIR` AND `git check-ignore` accepts
              the path (i.e., the path is gitignored).
     If any path fails validation, the whole command is blocked. If
     every path passes, the hook exits 0 with a stderr line explaining
     the verdict (e.g. "safe: under /tmp/", "safe: gitignored under
     project root").

  3. CLAUDE CODE PERMISSION SYSTEM (downstream of this hook) — the
     bundle's `settings.json` does NOT list `Bash(rm:*)` /
     `Bash(rmdir:*)` / `Bash(shred:*)` in `permissions.deny` and does
     NOT list them in `permissions.allow` either. With no rule
     matching, Claude Code prompts the user per-invocation. The hook
     is the safety net (any deletion that reaches the prompt has
     already cleared path validation); the user is the final
     authority (each delete is explicitly approved by hand).

Together: the hook says "this path is safe to delete," and the user
says "yes, delete it now." Either alone is insufficient; both must
agree. The "approve every deletion by hand" discipline is preserved
without forcing Claude to leak its own scratch dirs.

The hook runs in well under 100ms (one shlex parse, at most a couple
of `git check-ignore` invocations with a 2-second timeout each).
"""
import json
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Optional


# ── Layer 1: blanket-banned patterns ─────────────────────────────────────
#
# These have no safe variant, so they're matched as regex over the raw
# command string and rejected unconditionally.

_BLANKET_DEFS = [
    # ALL `git push` invocations are blocked. Pushing to a remote is a
    # human review gate: the maintainer reads the local commits and
    # pushes by hand. The agent should never issue a push, force or
    # otherwise; a non-force push of a wrong commit is just as
    # disruptive on a shared branch as a force push, just slower to
    # notice.
    (r"\bgit\s+push\b", "git push",
     "remote push must be reviewed and issued by hand"),
    # `git reset --hard origin/...` wipes local commits.
    (r"\bgit\s+reset\s+--hard\s+origin\b", "git reset --hard origin",
     "wipes local commits"),
    # Raw disk writes.
    (r"\bdd\s+[^\n]*\bif=", "dd if=…",
     "raw disk write"),
    # Filesystem creation: `mkfs.<fs>`.
    (r"\bmkfs\.[a-zA-Z0-9]+\b", "mkfs.<fs>",
     "filesystem creation"),
    # Classic fork bomb literal.
    (r":\(\)\{\s*:\|:&\s*\};:", ":(){ :|:& };:",
     "fork bomb"),
    # Permission chaos at root (`chmod -R 777 /`).
    (r"\bchmod\s+-R\s+777\s+/(?![A-Za-z0-9._-])", "chmod -R 777 /",
     "permission chaos at root"),
    # Recursive ownership change (often used in destructive recovery
    # scripts).
    (r"\bchown\s+-R\b", "chown -R",
     "destructive recovery patterns"),
    # `xargs rm` — paths come from stdin (a previous command's output),
    # so this hook cannot statically validate them. Block entirely; the
    # maintainer can split the work into an explicit enumeration if
    # needed.
    (r"\bxargs\b[^\n]*\brm\b", "xargs rm",
     "paths come from stdin and cannot be statically validated"),
]

BLANKET = [(re.compile(p), label, why) for p, label, why in _BLANKET_DEFS]


# ── Layer 2: path-validated deletion commands ───────────────────────────

ALLOWED_TEMP_ROOTS = ("/tmp", "/var/tmp", "/scratch")
SHELL_SEPS = {"&&", "||", ";", "|", "&"}
DESTRUCTIVE_BASENAMES = {"rm", "rmdir", "shred"}
SHELL_METACHARS = set("$`*?[]{}")


def basename_is(token: str, name: str) -> bool:
    """True if token's basename equals name. Handles `/usr/bin/rm`, `rm`."""
    return os.path.basename(token) == name


def basename_in(token: str, names: set) -> Optional[str]:
    """Return the matching name if token's basename is in names, else None."""
    base = os.path.basename(token)
    return base if base in names else None


def validate_delete_path(
    path_str: str, project_root: Path
) -> tuple:
    """Return (is_safe, reason). reason is human-readable for both branches."""
    # Require absolute path. Relative paths depend on cwd at execution
    # time, which the hook cannot reliably predict; the safest stance is
    # to require Claude to spell paths out fully.
    if not path_str.startswith("/"):
        return False, f"path '{path_str}' is not absolute"

    # Reject shell metacharacters that would expand differently at exec
    # time. The hook sees the literal command string before bash expands
    # `$VAR`, globs `*`/`?`/`[...]`, brace `{a,b}`, or command
    # substitution `` ` ``. If Claude wants to delete a glob match, it
    # must enumerate the literal paths first.
    bad = sorted(set(path_str) & SHELL_METACHARS)
    if bad:
        return False, (
            f"path contains shell metacharacters {bad}; "
            f"enumerate the literal paths instead"
        )

    # Resolve `..` and follow symlinks. If a symlinked path under /tmp
    # actually points at /etc, this catches it.
    try:
        resolved = Path(path_str).resolve(strict=False)
    except (OSError, RuntimeError) as exc:
        return False, f"path resolution failed: {exc}"

    resolved_str = str(resolved)

    # Refuse to delete an allowlist root itself (e.g. `rm -rf /tmp`).
    for root in ALLOWED_TEMP_ROOTS:
        if resolved_str == root:
            return False, f"refusing to delete the temp root '{root}' itself"

    # Allowlist 1: under a temp root.
    for root in ALLOWED_TEMP_ROOTS:
        if resolved_str.startswith(root + "/"):
            return True, f"safe: under {root}/"

    # Allowlist 2: inside project root AND gitignored.
    try:
        resolved.relative_to(project_root)
        in_project = True
    except ValueError:
        in_project = False

    if in_project:
        try:
            result = subprocess.run(
                ["git", "check-ignore", "-q", str(resolved)],
                cwd=project_root,
                capture_output=True,
                timeout=2,
                check=False,
            )
        except (subprocess.SubprocessError, OSError, FileNotFoundError) as exc:
            return False, f"could not consult git check-ignore: {exc}"

        if result.returncode == 0:
            return True, "safe: gitignored under project root"
        if result.returncode == 1:
            return False, (
                "under project root but NOT gitignored — "
                "refuse to delete tracked content"
            )
        return False, (
            f"git check-ignore returned exit {result.returncode}; "
            f"refuse to delete on uncertain git state"
        )

    return False, (
        "path is not under any approved deletion root "
        "(/tmp/, /var/tmp/, /scratch/, or gitignored project paths)"
    )


def normalize_separators(command: str) -> str:
    """
    Insert spaces around shell separators (`&&`, `||`, `;`, `|`, `&`)
    OUTSIDE of single- or double-quoted regions, so that the subsequent
    `shlex.split` produces them as standalone tokens. Without this,
    `rm /tmp/foo;rm /tmp/bar` would parse as a single rm with three
    args, not as two separate rm invocations — which would let chained
    deletions slip past the per-invocation validation.

    Quoted regions are passed through unchanged so that `echo
    "hello; world"` stays a single argument.
    """
    out = []
    i = 0
    n = len(command)
    quote = None
    while i < n:
        c = command[i]
        if quote:
            out.append(c)
            if c == "\\" and i + 1 < n:
                out.append(command[i + 1])
                i += 2
                continue
            if c == quote:
                quote = None
            i += 1
            continue
        if c in ("'", '"'):
            quote = c
            out.append(c)
            i += 1
            continue
        # Two-character separators first.
        if i + 1 < n and command[i:i + 2] in ("&&", "||"):
            out.append(" ")
            out.append(command[i:i + 2])
            out.append(" ")
            i += 2
            continue
        if c in (";", "|", "&"):
            out.append(" ")
            out.append(c)
            out.append(" ")
            i += 1
            continue
        out.append(c)
        i += 1
    return "".join(out)


def extract_destructive_invocations(command: str) -> list:
    """
    Parse the command and return a list of (label, [path_args]) tuples
    for every destructive invocation: rm, rmdir, shred, find … -delete,
    find … -exec rm.

    Raises ValueError on shell parse error.
    """
    try:
        tokens = shlex.split(normalize_separators(command), posix=True)
    except ValueError as exc:
        raise ValueError(f"shell parse error: {exc}") from exc

    invocations = []
    n = len(tokens)
    i = 0
    at_command_pos = True

    while i < n:
        tok = tokens[i]
        if tok in SHELL_SEPS:
            at_command_pos = True
            i += 1
            continue

        if at_command_pos:
            # rm / rmdir / shred
            cmd_name = basename_in(tok, DESTRUCTIVE_BASENAMES)
            if cmd_name is not None:
                paths = []
                j = i + 1
                end_of_flags = False
                while j < n and tokens[j] not in SHELL_SEPS:
                    arg = tokens[j]
                    if not end_of_flags:
                        if arg == "--":
                            end_of_flags = True
                            j += 1
                            continue
                        if arg.startswith("-"):
                            j += 1
                            continue
                    paths.append(arg)
                    j += 1
                invocations.append((cmd_name, paths))
                i = j
                continue

            # find … -delete or find … -exec rm
            if basename_is(tok, "find"):
                find_args = []
                j = i + 1
                while j < n and tokens[j] not in SHELL_SEPS:
                    find_args.append(tokens[j])
                    j += 1

                is_destructive = False
                label = "find"
                if "-delete" in find_args:
                    is_destructive = True
                    label = "find … -delete"
                for k in range(len(find_args) - 1):
                    if find_args[k] == "-exec" and basename_is(
                        find_args[k + 1], "rm"
                    ):
                        is_destructive = True
                        label = "find … -exec rm"
                        break

                if is_destructive:
                    # find search roots: every leading arg that does not
                    # start with `-` (find's predicates begin with `-`).
                    roots = []
                    for arg in find_args:
                        if arg.startswith("-"):
                            break
                        roots.append(arg)
                    if not roots:
                        # find's default search root is `.` — relative,
                        # so it will fail validation. Surface it
                        # explicitly so the error message is clear.
                        roots = ["."]
                    invocations.append((label, roots))
                i = j
                continue

        at_command_pos = False
        i += 1

    return invocations


def main() -> int:
    try:
        payload = json.load(sys.stdin)
    except (json.JSONDecodeError, ValueError):
        return 0

    command = payload.get("tool_input", {}).get("command", "") or ""
    if not command:
        return 0

    # Layer 1: blanket-banned patterns.
    for pattern, label, why in BLANKET:
        if pattern.search(command):
            print(
                f"BLOCKED: command contains dangerous pattern '{label}'.\n"
                f"  Command: {command}\n"
                f"  Reason : {why}.\n"
                f"  This pattern has no safe variant; run it directly in your\n"
                f"  shell outside Claude Code if you genuinely need it.",
                file=sys.stderr,
            )
            return 2

    # Layer 2: extract any destructive invocations and validate paths.
    try:
        invocations = extract_destructive_invocations(command)
    except ValueError as exc:
        print(
            f"BLOCKED: cannot parse command for path validation.\n"
            f"  Command: {command}\n"
            f"  Reason : {exc}.",
            file=sys.stderr,
        )
        return 2

    if not invocations:
        # No deletion command in this bash invocation.
        return 0

    project_root = Path(
        os.environ.get("CLAUDE_PROJECT_DIR", os.getcwd())
    ).resolve()

    safe_messages = []
    for cmd_label, paths in invocations:
        if not paths:
            print(
                f"BLOCKED: {cmd_label} has no path arguments.\n"
                f"  Command: {command}\n"
                f"  Reason : refuse to run a deletion command without explicit\n"
                f"           path arguments.",
                file=sys.stderr,
            )
            return 2

        for path_str in paths:
            ok, reason = validate_delete_path(path_str, project_root)
            if not ok:
                print(
                    f"BLOCKED: {cmd_label} path is not in an approved deletion root.\n"
                    f"  Path  : {path_str}\n"
                    f"  Reason: {reason}.\n"
                    f"  Approved roots:\n"
                    f"    - /tmp/, /var/tmp/, /scratch/ (POSIX/system temp)\n"
                    f"    - any path inside the project root that is gitignored\n"
                    f"  Path arguments must be absolute (start with '/') and must\n"
                    f"  not contain shell metacharacters ($, `, *, ?, [, ], {{, }}).\n"
                    f"  If you want to delete a glob match, enumerate the matched\n"
                    f"  paths first.\n"
                    f"  Command: {command}",
                    file=sys.stderr,
                )
                return 2
            safe_messages.append(f"  - {cmd_label} {path_str} → {reason}")

    # All paths cleared validation. Surface the verdict so the user-approval
    # prompt that follows in Claude Code's permission system shows the
    # categorization the hook applied. The hook still exits 0 (allow);
    # Claude Code's permission system is what asks the user.
    print(
        "block_dangerous_bash: deletion paths validated; awaiting user approval.\n"
        + "\n".join(safe_messages),
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
