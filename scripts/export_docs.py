#!/usr/bin/env python3
"""Combine all Lancet2 documentation into a single Markdown file.

Dynamically reads the navigation structure from mkdocs.yml so the output
stays in sync even as pages are added, removed, or reordered.

Features:
  - Images are embedded as base64 data URIs (self-contained document)
  - Internal cross-page links are rewritten to in-document anchors
  - Git version provenance and UTC timestamp in the header

Usage:
    python3 scripts/export_docs.py                          # writes to stdout
    python3 scripts/export_docs.py > Lancet2_Docs.md        # redirect to file
    python3 scripts/export_docs.py -o Lancet2_Docs.md       # write to file directly
"""

from __future__ import annotations

import argparse
import base64
import mimetypes
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

try:
    import yaml
except ImportError:
    sys.exit("PyYAML is required: pip install pyyaml")

ROOT = Path(__file__).resolve().parent.parent
DOCS_DIR = ROOT / "docs"
MKDOCS_YML = ROOT / "mkdocs.yml"


def get_version_info() -> dict[str, str]:
    """Extract version metadata from VERSION.txt and git.

    Format: VERSION_TAG_BRANCH_COMMIT[-dirty]
    Example: 2.9.0_main_a3b4c5d6e7
    """
    def git(*args: str) -> str:
        try:
            return subprocess.check_output(
                ["git", *args], cwd=ROOT, stderr=subprocess.DEVNULL
            ).decode().strip()
        except (subprocess.CalledProcessError, FileNotFoundError):
            return "unknown"

    version_file = ROOT / "VERSION.txt"
    version_tag = version_file.read_text().strip() if version_file.exists() else "unknown"

    branch = git("rev-parse", "--abbrev-ref", "HEAD")
    commit = git("rev-parse", "--short=10", "--verify", "HEAD")

    try:
        subprocess.check_call(
            ["git", "diff", "--quiet"], cwd=ROOT, stderr=subprocess.DEVNULL
        )
        dirty = ""
    except (subprocess.CalledProcessError, FileNotFoundError):
        dirty = "-dirty"

    build_tag = f"{version_tag}_{branch}_{commit}{dirty}"
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    return {"build_tag": build_tag, "timestamp": timestamp, "commit": commit}


# ---------------------------------------------------------------------------
# Nav entry: (title, filepath | None, depth)
#   filepath=None means a section header with no content file.
# ---------------------------------------------------------------------------
NavEntry = tuple[str, str | None, int]


def parse_nav(nav_list: list, depth: int = 1) -> list[NavEntry]:
    """Recursively parse the mkdocs.yml nav list into flat (title, path, depth) entries."""
    entries: list[NavEntry] = []
    for item in nav_list:
        if isinstance(item, str):
            entries.append((Path(item).stem.replace("_", " ").title(), item, depth))
        elif isinstance(item, dict):
            for title, value in item.items():
                if isinstance(value, str):
                    entries.append((title, value, depth))
                elif isinstance(value, list):
                    entries.append((title, None, depth))
                    entries.extend(parse_nav(value, depth + 1))
    return entries


def strip_frontmatter(content: str) -> str:
    """Remove YAML front matter (---...---) from markdown content."""
    return re.sub(r"^---\s*\n.*?\n---\s*\n", "", content, count=1, flags=re.DOTALL)


def bump_headings(content: str, level: int) -> str:
    """Bump all markdown headings by `level` levels so page content nests under the TOC."""
    def replacer(m: re.Match) -> str:
        return "#" * (len(m.group(1)) + level) + m.group(2)
    return re.sub(r"^(#{1,6})(.*)", replacer, content, flags=re.MULTILINE)


def make_anchor(title: str) -> str:
    """Convert a title to a GitHub-style markdown anchor."""
    anchor = title.lower()
    anchor = re.sub(r"[^a-z0-9 -]", "", anchor)
    anchor = anchor.strip().replace(" ", "-")
    anchor = re.sub(r"-+", "-", anchor)
    return anchor


# ---------------------------------------------------------------------------
# Image handling: embed as base64 or flatten to co-located filenames
# ---------------------------------------------------------------------------
def process_images(content: str, source_dir: Path, *, embed: bool,
                   referenced: set[str] | None = None) -> str:
    """Process image references in markdown content.

    When embed=True:  replace with inline base64 data URIs (self-contained).
    When embed=False: rewrite paths to flat filenames (assume co-located).

    Handles:
      ![alt](../assets/foo.png)
      ![alt](../assets/foo.png#only-light)
      ![alt](../assets/foo.png){ align=left, loading=lazy }
    """
    def replacer(m: re.Match) -> str:
        alt = m.group(1)
        raw_path = m.group(2)

        # Strip MkDocs Material URL fragments (#only-light, #only-dark)
        clean_path = raw_path.split("#")[0]

        # Skip external URLs
        if clean_path.startswith(("http://", "https://", "data:")):
            return m.group(0)

        img_path = (source_dir / clean_path).resolve()
        if not img_path.exists():
            print(f"  ⚠ Image not found: {img_path}", file=sys.stderr)
            return m.group(0)

        if embed:
            mime, _ = mimetypes.guess_type(str(img_path))
            if mime is None:
                mime = "application/octet-stream"
            data = base64.b64encode(img_path.read_bytes()).decode("ascii")
            return f"![{alt}](data:{mime};base64,{data})"

        # Flatten: use just the filename, track for manifest
        filename = img_path.name
        if referenced is not None:
            referenced.add(str(img_path))
        return f"![{alt}]({filename})"

    # Strip MkDocs Material { attrs } — they don't work in vanilla markdown
    return re.sub(r"!\[(.*?)\]\((.*?)\)(\{.*?\})?", replacer, content)


# ---------------------------------------------------------------------------
# Internal link rewriting: convert cross-page .md links to #anchor refs
# ---------------------------------------------------------------------------
def build_file_to_title_map(entries: list[NavEntry]) -> dict[str, str]:
    """Build a mapping from doc filepath -> nav title for anchor resolution."""
    mapping: dict[str, str] = {}
    for title, filepath, _ in entries:
        if filepath is not None:
            # Store multiple lookup keys for the same file:
            # "guides/architecture.md", "architecture.md", "../reference.md"
            mapping[filepath] = title
            mapping[Path(filepath).name] = title
            # For guides linking to root pages via ../
            if "/" not in filepath:
                mapping[f"../{filepath}"] = title
    return mapping


def rewrite_internal_links(content: str, file_map: dict[str, str]) -> str:
    """Rewrite internal .md links to in-document #anchor links."""
    def replacer(m: re.Match) -> str:
        link_text = m.group(1)
        raw_target = m.group(2)

        # Skip external, anchor-only, and non-md links
        if raw_target.startswith(("http://", "https://", "mailto:", "#", "data:")):
            return m.group(0)
        if ".md" not in raw_target:
            return m.group(0)

        # Split path and fragment: "architecture.md#section" -> ("architecture.md", "section")
        if "#" in raw_target:
            md_path, fragment = raw_target.split("#", 1)
        else:
            md_path, fragment = raw_target, ""

        # Also try stripping "guides/" prefix for root-level lookups
        title = file_map.get(md_path) or file_map.get(Path(md_path).name)

        if title is None:
            # Can't resolve — leave the link as-is
            return m.group(0)

        if fragment:
            return f"[{link_text}](#{fragment})"
        else:
            return f"[{link_text}](#{make_anchor(title)})"

    return re.sub(r"\[(.*?)\]\((.*?)\)", replacer, content)


def build_document(entries: list[NavEntry], *, embed_images: bool = True) -> tuple[str, set[str]]:
    """Assemble the combined markdown document from nav entries.

    Returns (document_text, set_of_referenced_image_absolute_paths).
    """
    parts: list[str] = []
    version = get_version_info()
    file_map = build_file_to_title_map(entries)
    referenced_images: set[str] = set()

    # Title page
    parts.append("# Lancet2 — Complete Documentation\n")
    parts.append(f"> **Version:** `{version['build_tag']}`\n")
    parts.append(f"> **Generated:** {version['timestamp']}\n")
    parts.append("> **Source:** <https://github.com/nygenome/Lancet2>\n")

    # Table of contents
    parts.append("\n---\n\n## Table of Contents\n")
    for title, filepath, depth in entries:
        indent = "  " * (depth - 1)
        if filepath is None:
            parts.append(f"{indent}- **{title}**")
        else:
            parts.append(f"{indent}- [{title}](#{make_anchor(title)})")
    parts.append("\n---\n")

    # Content
    for title, filepath, depth in entries:
        if filepath is None:
            parts.append(f"\n---\n\n# {title}\n")
            continue

        src = DOCS_DIR / filepath
        if not src.exists():
            parts.append(f"\n## {title}\n\n> *File not found: {filepath}*\n")
            continue

        raw = src.read_text(encoding="utf-8")
        content = strip_frontmatter(raw)

        # Process images (embed as base64 or flatten to co-located filenames)
        content = process_images(content, src.parent, embed=embed_images,
                                 referenced=referenced_images)

        # Rewrite internal .md links to #anchors
        content = rewrite_internal_links(content, file_map)

        # Strip the first H1 from the file — we use our own section header
        content = re.sub(r"^#\s+.*\n", "", content, count=1)
        content = content.strip()

        if depth == 1:
            parts.append(f"\n---\n\n# {title}\n")
        else:
            parts.append(f"\n---\n\n## {title}\n")

        bumped = bump_headings(content, depth)
        parts.append(bumped)
        parts.append("")

    return "\n".join(parts), referenced_images


def main() -> None:
    parser = argparse.ArgumentParser(description="Export Lancet2 docs to a single Markdown file")
    parser.add_argument("-o", "--output", default=None,
                        help="Output file path (default: stdout)")
    parser.add_argument("--no-embed-images", action="store_true",
                        help="Don't embed images as base64. Instead, rewrite paths to flat "
                             "filenames and print the required image manifest to stderr.")
    args = parser.parse_args()

    if not MKDOCS_YML.exists():
        sys.exit(f"mkdocs.yml not found at {MKDOCS_YML}")

    with open(MKDOCS_YML, encoding="utf-8") as fh:
        config = yaml.safe_load(fh)

    nav_list = config.get("nav")
    if not nav_list:
        sys.exit("No 'nav' key found in mkdocs.yml")

    entries = parse_nav(nav_list)
    document, referenced = build_document(entries, embed_images=not args.no_embed_images)

    if args.output:
        output_path = Path(args.output)
        output_path.write_text(document, encoding="utf-8")
        size_kb = output_path.stat().st_size / 1024
        print(f"✓ Written {output_path.name} ({size_kb:.0f} KB, {len(entries)} sections)", file=sys.stderr)
    else:
        sys.stdout.write(document)

    if args.no_embed_images and referenced:
        print(f"\n📎 {len(referenced)} image(s) expected in the same directory as the markdown file:", file=sys.stderr)
        for img in sorted(referenced):
            print(f"   {img}", file=sys.stderr)


if __name__ == "__main__":
    main()
