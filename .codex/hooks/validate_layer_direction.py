#!/usr/bin/env python3
"""
Lancet2 PreToolUse hook: validate_layer_direction

Enforces the six-layer dependency direction documented at the top of
CMakeLists.txt. A file in layer N may include from layers 1..N but never
from layers N+1..6. Violations are blocked at write time.

The layer order (lowest to highest) is:
    1. base
    2. hts
    3. cbdg
    4. caller
    5. core
    6. cli

Files outside src/lancet/ (tests, benchmarks, main.cpp) are exempt.
This hook only inspects #include directives in the proposed write content.
"""
import json
import re
import sys
from pathlib import Path

LAYER_ORDER = ["base", "hts", "cbdg", "caller", "core", "cli"]
LAYER_INDEX = {name: i for i, name in enumerate(LAYER_ORDER)}

# Match: #include "lancet/<layer>/..." or #include <lancet/<layer>/...>
INCLUDE_RE = re.compile(
    r'^\s*#\s*include\s*[<"]lancet/(' + "|".join(LAYER_ORDER) + r')/',
    re.MULTILINE,
)

# Match the layer of the file being edited: src/lancet/<layer>/...
FILE_LAYER_RE = re.compile(r"src/lancet/(" + "|".join(LAYER_ORDER) + r")/")


def file_layer(file_path: str) -> str | None:
    """Return the layer name a file belongs to, or None if outside the layered tree."""
    m = FILE_LAYER_RE.search(file_path)
    return m.group(1) if m else None


def collect_new_text(payload: dict) -> str:
    """Return the new content that will be written, concatenated for inspection."""
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
    if not file_path:
        return 0

    own_layer = file_layer(file_path)
    if own_layer is None:
        # File is outside src/lancet/<layer>/ — exempt.
        return 0

    own_idx = LAYER_INDEX[own_layer]
    new_text = collect_new_text(payload)
    if not new_text:
        return 0

    violations: list[tuple[str, str]] = []
    for match in INCLUDE_RE.finditer(new_text):
        included_layer = match.group(1)
        included_idx = LAYER_INDEX[included_layer]
        if included_idx > own_idx:
            full_line = match.group(0).strip()
            violations.append((included_layer, full_line))

    if violations:
        msg = [
            f"BLOCKED: layer-direction violation in {file_path}",
            f"  This file is in layer '{own_layer}' (index {own_idx}).",
            f"  It may include from layers: {', '.join(LAYER_ORDER[:own_idx + 1])}.",
            f"  It MUST NOT include from layers: {', '.join(LAYER_ORDER[own_idx + 1:])}.",
            "",
            "  Violations found in proposed change:",
        ]
        for layer, line in violations:
            msg.append(f"    - includes '{layer}': {line}")
        msg.extend([
            "",
            "  See the data-flow comments at the top of CMakeLists.txt for the",
            "  full architecture rationale. If you genuinely need an upward",
            "  dependency, you are probably looking at a design problem — consider",
            "  extracting an interface to the lower layer instead.",
        ])
        print("\n".join(msg), file=sys.stderr)
        return 2

    return 0


if __name__ == "__main__":
    sys.exit(main())
