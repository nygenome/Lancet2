#!/usr/bin/env python3
"""Run clang-tidy on Lancet2 source files (never touches third-party deps).

Usage:
    ./scripts/run_clang_tidy.py                              # check-only
    ./scripts/run_clang_tidy.py --build-dir build            # custom build dir

Auto-fix mode is intentionally NOT supported. Clang-tidy's auto-fix has
historically broken compilation and produced unreadable code in this
project; every clang-tidy diagnostic must be resolved by reading the
warning, understanding the root cause, and editing the source by hand.
See AGENTS.md for the project-wide rule.

Prerequisites:
    - Build directory with compile_commands.json
    - pixi environment with clang-tidy and run-clang-tidy
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path


def ensure_pixi() -> None:
    """Install pixi if it is not already on PATH."""
    if shutil.which("pixi") is not None:
        return
    print("pixi not found, installing...")
    subprocess.run(
        ["bash", "-c", "curl -fsSL https://pixi.sh/install.sh | bash"],
        check=True,
    )


REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_BUILD_DIR = "cmake-build-release"

# Regex filter applied to compile_commands.json entries.
# Only files whose absolute path matches this regex are processed.
# This prevents clang-tidy from touching _deps/, third-party, or build artifacts.
FILE_FILTER = r"(src/lancet|tests|benchmarks)/.*\.(cpp|h)$"

# Only report warnings for Lancet headers (suppresses noise from system/dep headers).
HEADER_FILTER = r"src/lancet/.*"


# ──────────────────────────────────────────────────────────────────────────────
# Output parser
# ──────────────────────────────────────────────────────────────────────────────
_RE_DIAG = re.compile(
    r"^(?P<file>.+):(?P<line>\d+):(?P<col>\d+): "
    r"(?P<level>warning|error): "
    r"(?P<msg>.+?) "
    r"\[(?P<check>[^\]]+)\]$"
)

# Strip ANSI escape sequences (colors, bold, reset) from clang-tidy output.
# clang-tidy emits these when UseColor: true is set in .clang-tidy, regardless
# of the -use-color flag passed to run-clang-tidy.
_RE_ANSI = re.compile(r"\x1b\[[0-9;]*m")


@dataclass
class Diagnostic:
    """A single clang-tidy diagnostic."""

    file: str
    line: int
    col: int
    level: str
    msg: str
    check: str


def parse_clang_tidy_output(
    output: str, repo_root: Path
) -> list[Diagnostic]:
    """Parse clang-tidy output into structured diagnostics."""
    diagnostics: list[Diagnostic] = []
    repo_prefix = str(repo_root) + "/"

    for line in output.splitlines():
        clean_line = _RE_ANSI.sub("", line)
        match = _RE_DIAG.match(clean_line)
        if not match:
            continue

        filepath = match.group("file")
        # Only include diagnostics from our source tree
        if not filepath.startswith(repo_prefix):
            continue

        diagnostics.append(
            Diagnostic(
                file=filepath[len(repo_prefix) :],
                line=int(match.group("line")),
                col=int(match.group("col")),
                level=match.group("level"),
                msg=match.group("msg"),
                check=match.group("check"),
            )
        )

    return diagnostics


# ──────────────────────────────────────────────────────────────────────────────
# Summary printing
# ──────────────────────────────────────────────────────────────────────────────
_YELLOW = "\033[33m"
_RED = "\033[31m"
_GREEN = "\033[32m"
_CYAN = "\033[36m"
_DIM = "\033[2m"
_BOLD = "\033[1m"
_RESET = "\033[0m"


def _color(text: str, code: str) -> str:
    if not sys.stdout.isatty():
        return text
    return f"{code}{text}{_RESET}"


def print_summary(diagnostics: list[Diagnostic]) -> None:
    """Print a concise, color-coded summary of clang-tidy findings."""
    errors = [d for d in diagnostics if d.level == "error"]
    warnings = [d for d in diagnostics if d.level == "warning"]

    # Group by check name for the category breakdown
    by_check: dict[str, int] = defaultdict(int)
    for d in diagnostics:
        by_check[d.check] += 1

    # Group by file
    by_file: dict[str, list[Diagnostic]] = defaultdict(list)
    for d in diagnostics:
        by_file[d.file].append(d)

    print()
    print(_color("─" * 80, _CYAN))
    print(_color("  clang-tidy Summary", _BOLD))
    print(_color("─" * 80, _CYAN))
    print(f"  Total diagnostics: {len(diagnostics)}")
    print(f"  Errors:            {_color(str(len(errors)), _RED)}")
    print(f"  Warnings:          {_color(str(len(warnings)), _YELLOW)}")
    print(f"  Files affected:    {len(by_file)}")
    print(f"  Check categories:  {len(by_check)}")
    print(_color("─" * 80, _CYAN))

    if not diagnostics:
        print(_color("\n  ✓ No clang-tidy warnings or errors.\n", _GREEN))
        return

    # Print check breakdown sorted by frequency
    print(_color("\n  By check:", _BOLD))
    for check, count in sorted(by_check.items(), key=lambda x: -x[1]):
        print(f"    {count:3d}  {check}")

    # Print per-file diagnostics
    print(_color("\n  By file:", _BOLD))
    for filepath in sorted(by_file):
        file_diags = by_file[filepath]
        print(f"\n  {_color(filepath, _BOLD)}")
        for d in sorted(file_diags, key=lambda x: x.line):
            level_color = _RED if d.level == "error" else _YELLOW
            loc = _color(f":{d.line}:{d.col}", _DIM)
            check = _color(f"[{d.check}]", _DIM)
            print(f"    {_color(d.level, level_color)} {loc} {d.msg} {check}")

    print()
    print(
        _color(
            f"  ✗ {len(diagnostics)} diagnostic(s) in {len(by_file)} file(s). "
            f"Fix the issues above to pass this check.",
            _YELLOW,
        )
    )
    print()


# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────
def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run clang-tidy on Lancet2 source files (check-only)."
    )
    parser.add_argument(
        "--build-dir",
        default=DEFAULT_BUILD_DIR,
        help=f"Path to CMake build directory (default: {DEFAULT_BUILD_DIR})",
    )
    args = parser.parse_args()

    build_dir = REPO_ROOT / args.build_dir
    compile_db = build_dir / "compile_commands.json"

    if not compile_db.exists():
        print(
            f"ERROR: {compile_db.relative_to(REPO_ROOT)} not found.\n"
            f"Run: pixi run build",
            file=sys.stderr,
        )
        return 1

    ensure_pixi()
    cmd: list[str] = [
        "pixi", "run", "run-clang-tidy",
        "-warnings-as-errors", "*",
        "-hide-progress",
        "-p", str(build_dir),
        "-quiet", "-use-color",
        f"-header-filter={HEADER_FILTER}",
    ]

    print("==> Running clang-tidy check on Lancet sources (no fixes)...")

    print(f"    Build dir: {build_dir.relative_to(REPO_ROOT)}")
    print(f"    File filter: {FILE_FILTER}")
    print()

    cmd.append(FILE_FILTER)

    # Stream output live while capturing it for the summary.
    # run-clang-tidy forks child clang-tidy processes — their diagnostics go to
    # stderr, so we must merge stderr into stdout to capture everything.
    output_lines: list[str] = []
    with subprocess.Popen(
        cmd, cwd=REPO_ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
    ) as proc:
        assert proc.stdout is not None
        for line in proc.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            output_lines.append(line)
        proc.wait()

    output = "".join(output_lines)
    diagnostics = parse_clang_tidy_output(output, REPO_ROOT)
    print_summary(diagnostics)

    if diagnostics:
        return 1

    # Fallback: if the parser missed diagnostics but run-clang-tidy itself
    # returned non-zero, propagate that exit code.
    if proc.returncode != 0:
        return proc.returncode

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
