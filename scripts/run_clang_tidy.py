#!/usr/bin/env python3
"""Run clang-tidy on Lancet2 source files (never touches third-party deps).

Usage:
    ./scripts/run_clang_tidy.py                              # check-only
    ./scripts/run_clang_tidy.py --fix                        # apply auto-fixes
    ./scripts/run_clang_tidy.py --build-dir build            # custom build dir
    ./scripts/run_clang_tidy.py --fix --build-dir build      # fix with custom dir

Prerequisites:
    - Build directory with compile_commands.json
    - pixi environment with clang-tidy and run-clang-tidy
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
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


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run clang-tidy on Lancet2 source files."
    )
    parser.add_argument(
        "--fix",
        action="store_true",
        help="Apply clang-tidy auto-fixes (default: check-only)",
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
            f"Run: pixi run cmake -B {args.build_dir} -DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
            file=sys.stderr,
        )
        return 1

    ensure_pixi()
    cmd: list[str] = [
        "pixi", "run", "run-clang-tidy",
        "-p", str(build_dir),
        "-quiet", "-use-color",
        f"-header-filter={HEADER_FILTER}",
    ]

    if args.fix:
        cmd.append("-fix")
        print("==> Running clang-tidy with auto-fix on Lancet sources...")
        print("    WARNING: Review changes carefully — some auto-fixes may break compilation.")
        print("    Known problematic checks (disabled in .clang-tidy):")
        print("      - modernize-use-trailing-return-type (-> auto on lambdas)")
        print("      - modernize-use-ranges (std::sort -> ranges::sort breaks w/o <=>)")
        print("      - misc-include-cleaner (reorders includes)")
    else:
        print("==> Running clang-tidy check on Lancet sources (no fixes)...")

    print(f"    Build dir: {build_dir.relative_to(REPO_ROOT)}")
    print(f"    File filter: {FILE_FILTER}")
    print()

    cmd.append(FILE_FILTER)

    result = subprocess.run(cmd, cwd=REPO_ROOT)
    print("\n==> Done.")
    return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())
