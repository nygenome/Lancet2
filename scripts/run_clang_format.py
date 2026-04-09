#!/usr/bin/env python3
"""Run clang-format on Lancet2 source files.

Usage:
    ./scripts/run_clang_format.py                    # dry-run on all project sources
    ./scripts/run_clang_format.py --fix              # apply formatting in-place
    ./scripts/run_clang_format.py src/lancet/cbdg    # check a specific directory
    ./scripts/run_clang_format.py --fix src/lancet/caller/genotyper.cpp  # fix one file
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_SOURCE_DIRS = ["src/lancet", "tests", "benchmarks"]
SOURCE_EXTENSIONS = {".cpp", ".h"}


def collect_sources(paths: list[Path]) -> list[Path]:
    """Recursively collect all C++ source files under the given paths."""
    sources: list[Path] = []
    for path in paths:
        if not path.exists():
            continue
        if path.is_file() and path.suffix in SOURCE_EXTENSIONS:
            sources.append(path)
        elif path.is_dir():
            sources.extend(
                f for f in sorted(path.rglob("*")) if f.suffix in SOURCE_EXTENSIONS
            )
    return sources


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run clang-format on Lancet2 source files."
    )
    parser.add_argument(
        "--fix",
        action="store_true",
        help="Apply formatting changes in-place (default: dry-run check)",
    )
    parser.add_argument(
        "sources",
        nargs="*",
        default=None,
        help=(
            "Files or directories to format. "
            f"Default: {', '.join(DEFAULT_SOURCE_DIRS)}"
        ),
    )
    args = parser.parse_args()

    # Resolve source paths relative to repo root
    if args.sources:
        paths = [REPO_ROOT / s for s in args.sources]
    else:
        paths = [REPO_ROOT / d for d in DEFAULT_SOURCE_DIRS]

    files = collect_sources(paths)
    if not files:
        print("No source files found.")
        return 0

    print(f"==> {'Formatting' if args.fix else 'Checking'} {len(files)} file(s)...")

    if args.fix:
        cmd = ["pixi", "run", "clang-format", "-i", *[str(f) for f in files]]
        result = subprocess.run(cmd, cwd=REPO_ROOT)
        if result.returncode == 0:
            print(f"==> Done. All {len(files)} files formatted.")
        return result.returncode

    # Dry-run: check each file individually to report which ones need formatting
    failures: list[Path] = []
    for f in files:
        result = subprocess.run(
            ["pixi", "run", "clang-format", "--dry-run", "--Werror", str(f)],
            cwd=REPO_ROOT,
            capture_output=True,
        )
        if result.returncode != 0:
            failures.append(f)
            # Print the diff output from clang-format
            if result.stderr:
                sys.stderr.buffer.write(result.stderr)

    if failures:
        print(f"\n==> {len(failures)} file(s) need formatting:")
        for f in failures:
            print(f"    {f.relative_to(REPO_ROOT)}")
        print("\nRun with --fix to apply changes.")
        return 1

    print(f"==> All {len(files)} files are correctly formatted.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
