#!/usr/bin/env python3
"""
diff_analysis — shared diff classification for Lancet2 hooks and skills.

Single source of truth for the heuristics that classify a staged git
diff: which Lancet2 layers it touches, what file kinds it spans, what
schema-affecting lines were added, and what NOLINT suppressions were
added. Originally inlined in pre_commit_summary.sh; extracted here so
external-interface-changes, agent-memory bundling, and future hooks can share
the same logic without divergent re-implementations.

Two entry points:

1. **Importable Python API** — call `analyze(repo_root)` to get a
   `DiffAnalysis` dataclass with all the classification fields. Use
   from agent-memory commit-bundling, external-interface-changes interview
   skills, etc.

2. **CLI fallback** — `python3 diff_analysis.py [--json]` emits the
   same fields as JSON, suitable for shell hooks that prefer to keep
   their bash structure (pre_commit_summary.sh) and just call out for
   the analysis.

The module deliberately fails open: if git is unavailable, not a repo,
or any subprocess errors, it returns an empty analysis rather than
raising. Callers can check `analysis.git_available` if they need to
distinguish "no diff" from "couldn't check".
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Iterable


# Fixed layer set, mirrors CMakeLists.txt and validate_layer_direction.py.
LANCET_LAYERS: tuple[str, ...] = ("base", "hts", "cbdg", "caller", "core", "cli")

# Source files (compiled C++).
_SRC_RE = re.compile(r"^src/.*\.(h|cpp|h\.inc)$")

# Schema-affecting VCF header line markers. Substring-anchored to
# `##FORMAT=`, `##INFO=`, `##FILTER=` so we catch additions inside the
# raw-string templates in vcf_header_builder.cpp as well as anywhere
# else they appear.
_SCHEMA_LINE_RE = re.compile(r"##(FORMAT|INFO|FILTER)=")

# NOLINT additions in the staged diff. Matches NOLINT, NOLINTNEXTLINE,
# NOLINTBEGIN, NOLINTEND. The leading `+` distinguishes additions from
# diff context lines.
_NOLINT_ADDITION_RE = re.compile(r"^\+[^+].*\bNOLINT(NEXTLINE|BEGIN|END)?\b")


@dataclass
class DiffAnalysis:
    """Result of analyzing the current git diff state.

    All fields default to empty/false when git is unavailable or no
    changes are staged, so consumers can call dataclass methods
    unconditionally.
    """

    git_available: bool = False
    is_git_repo: bool = False
    has_all_flag: bool = False  # whether `git commit -a/-am/--all` is in play
    file_count: int = 0
    staged_shortstat: str = ""
    unstaged_shortstat: str = ""
    files: list[str] = field(default_factory=list)
    layers_touched: list[str] = field(default_factory=list)
    has_source: bool = False
    has_tests: bool = False
    has_benchmarks: bool = False
    has_docs: bool = False
    has_build: bool = False
    has_config: bool = False
    has_other: bool = False
    schema_lines_added: list[str] = field(default_factory=list)
    nolint_additions: list[str] = field(default_factory=list)
    agent_memory_files_touched: list[str] = field(default_factory=list)

    @property
    def kinds(self) -> list[str]:
        labels = []
        if self.has_source:
            labels.append("source")
        if self.has_tests:
            labels.append("tests")
        if self.has_benchmarks:
            labels.append("benchmarks")
        if self.has_docs:
            labels.append("docs")
        if self.has_build:
            labels.append("build")
        if self.has_config:
            labels.append("config")
        if self.has_other:
            labels.append("other")
        return labels

    @property
    def has_schema_changes(self) -> bool:
        return bool(self.schema_lines_added)

    @property
    def has_nolint_additions(self) -> bool:
        return bool(self.nolint_additions)

    @property
    def is_cross_layer(self) -> bool:
        return len(self.layers_touched) > 1

    @property
    def needs_test_warning(self) -> bool:
        """Source change without a test/benchmark change."""
        return self.has_source and not (self.has_tests or self.has_benchmarks)

    @property
    def needs_docs_warning(self) -> bool:
        """Cross-layer change without a docs change."""
        return self.is_cross_layer and not self.has_docs

    def to_json(self) -> str:
        d = asdict(self)
        d["kinds"] = self.kinds
        d["has_schema_changes"] = self.has_schema_changes
        d["has_nolint_additions"] = self.has_nolint_additions
        d["is_cross_layer"] = self.is_cross_layer
        d["needs_test_warning"] = self.needs_test_warning
        d["needs_docs_warning"] = self.needs_docs_warning
        return json.dumps(d, indent=2)


def _run(args: list[str], cwd: Path | None = None) -> tuple[int, str]:
    """Run a subprocess, capturing stdout. Returns (rc, stdout). Empty
    on failure rather than raising, to match the fail-open contract."""
    try:
        result = subprocess.run(
            args,
            cwd=cwd,
            capture_output=True,
            text=True,
            check=False,
        )
        return result.returncode, result.stdout
    except (FileNotFoundError, OSError):
        return -1, ""


def _detect_layers(files: Iterable[str]) -> list[str]:
    found = []
    for layer in LANCET_LAYERS:
        prefix = f"src/lancet/{layer}/"
        if any(f.startswith(prefix) for f in files):
            found.append(layer)
    return found


def _classify_file(path: str) -> str | None:
    """Return one of: source, tests, benchmarks, docs, build, config,
    other. Returns None for empty input."""
    if not path:
        return None
    if _SRC_RE.match(path):
        return "source"
    if path.startswith("tests/"):
        return "tests"
    if path.startswith("benchmarks/"):
        return "benchmarks"
    if path.endswith(".md") or path.startswith("docs/") or path.startswith("mkdocs"):
        return "docs"
    if path in {"CMakeLists.txt", "pixi.toml", "pixi.lock"}:
        return "build"
    if path.startswith("cmake/") or path.startswith("conda/"):
        return "build"
    if path in {".clang-tidy", ".clang-format"}:
        return "config"
    if path.startswith(".chglog/") or path.startswith(".github/") or path.startswith(".claude/"):
        return "config"
    return "other"


def _detect_all_flag(commit_command: str) -> bool:
    """Detect whether `git commit -a` / `-am` / `--all` is in the
    bash command. Mirrors the substring check from
    validate_commit_message.py."""
    if not commit_command:
        return False
    # Check for --all or -a / -am. Avoid false-positive on -am inside
    # other words by checking for a space or end boundary.
    return bool(
        re.search(r"\b(--all|--all=)", commit_command)
        or re.search(r"(^|\s)-a(m)?(\s|$)", commit_command)
    )


def analyze(
    repo_root: Path | None = None,
    commit_command: str = "",
) -> DiffAnalysis:
    """Analyze the current git diff state.

    Args:
        repo_root: directory to run git commands in. Defaults to
            CLAUDE_PROJECT_DIR if set, otherwise CWD.
        commit_command: full bash command being executed (from the
            PreToolUse hook payload). Used to detect `-a` flag.

    Returns:
        DiffAnalysis dataclass with every field populated. Fail-open
        on errors: returns analysis with `git_available=False`.
    """
    cwd = repo_root or Path(os.environ.get("CLAUDE_PROJECT_DIR", "")) or Path.cwd()
    if isinstance(cwd, str):
        cwd = Path(cwd)
    if not cwd.is_dir():
        cwd = Path.cwd()

    rc, _ = _run(["git", "--version"])
    if rc != 0:
        return DiffAnalysis(git_available=False)

    rc, _ = _run(["git", "rev-parse", "--is-inside-work-tree"], cwd=cwd)
    if rc != 0:
        return DiffAnalysis(git_available=True, is_git_repo=False)

    has_all_flag = _detect_all_flag(commit_command)

    # Staged shortstat.
    _, staged_shortstat_raw = _run(["git", "diff", "--cached", "--shortstat"], cwd=cwd)
    staged_shortstat = staged_shortstat_raw.strip()

    # Unstaged shortstat (only meaningful when -a is in play).
    unstaged_shortstat = ""
    if has_all_flag:
        _, unstaged_raw = _run(["git", "diff", "--shortstat"], cwd=cwd)
        unstaged_shortstat = unstaged_raw.strip()

    # File lists.
    _, staged_files_raw = _run(["git", "diff", "--cached", "--name-only"], cwd=cwd)
    files = [f for f in staged_files_raw.splitlines() if f]
    if has_all_flag:
        _, unstaged_files_raw = _run(["git", "diff", "--name-only"], cwd=cwd)
        files.extend(f for f in unstaged_files_raw.splitlines() if f)
    # Deduplicate while preserving order.
    seen = set()
    unique_files = []
    for f in files:
        if f not in seen:
            seen.add(f)
            unique_files.append(f)
    files = unique_files

    layers_touched = _detect_layers(files)

    kinds = {k: False for k in ("source", "tests", "benchmarks", "docs", "build", "config", "other")}
    for f in files:
        kind = _classify_file(f)
        if kind:
            kinds[kind] = True

    # Schema-affecting lines added in the staged diff (and unstaged if
    # -a is in play). Use unified=0 to focus on actual changes.
    diff_args = ["git", "diff", "--cached", "-U0"]
    if has_all_flag:
        # Use HEAD as reference for combined diff.
        diff_args = ["git", "diff", "HEAD", "-U0"]
    _, diff_raw = _run(diff_args, cwd=cwd)

    schema_lines_added = []
    nolint_additions = []
    for line in diff_raw.splitlines():
        if line.startswith("+") and not line.startswith("+++"):
            if _SCHEMA_LINE_RE.search(line):
                schema_lines_added.append(line.rstrip())
            if _NOLINT_ADDITION_RE.match(line):
                nolint_additions.append(line.rstrip())

    # Agent memory files touched (under .claude/agent-memory/) — used by
    # the agent-memory bundling logic to decide whether the commit
    # already absorbs memory updates or needs a separate concern note.
    agent_memory_files_touched = [
        f for f in files if f.startswith(".claude/agent-memory/")
    ]

    return DiffAnalysis(
        git_available=True,
        is_git_repo=True,
        has_all_flag=has_all_flag,
        file_count=len(files),
        staged_shortstat=staged_shortstat,
        unstaged_shortstat=unstaged_shortstat,
        files=files,
        layers_touched=layers_touched,
        has_source=kinds["source"],
        has_tests=kinds["tests"],
        has_benchmarks=kinds["benchmarks"],
        has_docs=kinds["docs"],
        has_build=kinds["build"],
        has_config=kinds["config"],
        has_other=kinds["other"],
        schema_lines_added=schema_lines_added,
        nolint_additions=nolint_additions,
        agent_memory_files_touched=agent_memory_files_touched,
    )


def main(argv: list[str]) -> int:
    """CLI entrypoint. Reads the bash command from the JSON payload
    on stdin (matching the hook protocol) and emits a JSON analysis
    on stdout. Fail-open with `{}` on protocol errors."""
    try:
        payload = json.load(sys.stdin) if not sys.stdin.isatty() else {}
    except (json.JSONDecodeError, ValueError):
        payload = {}
    commit_command = payload.get("tool_input", {}).get("command", "")
    analysis = analyze(commit_command=commit_command)
    print(analysis.to_json())
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
