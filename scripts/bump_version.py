#!/usr/bin/env python3
import subprocess
import sys
import re
import os
import shutil

def run_cmd(cmd, check=True, capture=False):
    print(f"Running: {' '.join(cmd)}")
    res = subprocess.run(cmd, check=check, capture_output=capture, text=True)
    return res.stdout.strip() if capture else None



def get_latest_tag():
    try:
        # Get tags sorted by version descending
        out = run_cmd(["git", "tag", "-l", "--sort=-v:refname"], capture=True)
        if not out:
            return "v0.0.0"
            
        tags = [t for t in out.splitlines() if t]
        for tag in tags:
            # Find the highest semantic version tag
            if re.match(r'^v?\d+\.\d+\.\d+$', tag):
                return tag
        return "v0.0.0"
    except subprocess.CalledProcessError:
        return "v0.0.0"

def bump_version(tag, bump_kind):
    prefix = ""
    if tag.startswith('v'):
        prefix = "v"
        tag = tag[1:]
    
    parts = tag.split('.')
    major, minor, patch = int(parts[0]), int(parts[1]), int(parts[2])
    
    if bump_kind == "major":
        major += 1
        minor = 0
        patch = 0
    elif bump_kind == "minor":
        minor += 1
        patch = 0
    elif bump_kind == "patch":
        patch += 1
        
    return f"{prefix}{major}.{minor}.{patch}"

def main():
    bump_kind = "patch"
    if len(sys.argv) > 1:
        bump_kind = sys.argv[1].lower()
        
    if bump_kind not in ("major", "minor", "patch"):
        print("Usage: bump_version.py [major|minor|patch]")
        sys.exit(1)
        
    print(f"BUMP_KIND: {bump_kind}")
    
    # 1. Push upstream changes first
    run_cmd(["git", "push"])
    
    # 2. Compute incremented tag
    latest_tag = get_latest_tag()
    new_tag = bump_version(latest_tag, bump_kind)
    print(f"Bumping version: {latest_tag} -> {new_tag}")
    
    # 3. Create annotated local tag (bypasses nano interactive prompt)
    run_cmd(["git", "tag", "-m", f"Release {new_tag}", new_tag])
    
    # 4. Push commits and tags to remote
    run_cmd(["git", "push"])
    run_cmd(["git", "push", "--tags"])
    
    print(f"Successfully bumped version to {new_tag}.")

if __name__ == "__main__":
    main()
