#!/usr/bin/env bash
# Build a local Lancet2 conda package using rattler-build.
#
# This script mirrors the CI workflow (deploy_prefix.yml) by computing
# version metadata from git tags and injecting them as environment
# variables that the recipe.yaml reads via env.get().
#
# Usage:
#   ./scripts/build_conda_local.sh                  # output to ./output/
#   ./scripts/build_conda_local.sh /path/to/outdir  # output to custom dir

set -euo pipefail

# Ensure pixi is available, install if missing
if ! command -v pixi &>/dev/null; then
  echo "pixi not found, installing..."
  curl -fsSL https://pixi.sh/install.sh | bash
  export PATH="${HOME}/.pixi/bin:${PATH}"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
RECIPE="${REPO_ROOT}/conda/recipe.yaml"
OUTPUT_DIR="${1:-${REPO_ROOT}/output}"

# Read base version from VERSION.txt file (single source of truth)
BASE_VERSION=$(cat "${REPO_ROOT}/VERSION.txt" | tr -d '[:space:]')

# Augment with git metadata when available
BRANCH=$(git -C "${REPO_ROOT}" rev-parse --abbrev-ref HEAD 2>/dev/null || echo "")
HASH=$(git -C "${REPO_ROOT}" rev-parse --short=10 HEAD 2>/dev/null || echo "")
TAG=$(git -C "${REPO_ROOT}" describe --abbrev=0 --tags 2>/dev/null || echo "")
DISTANCE=$(git -C "${REPO_ROOT}" rev-list "${TAG}..HEAD" --count 2>/dev/null || echo "0")

# Stable release: version as-is (e.g. "2.10.0")
# Dev build: version_branch_hash (e.g. "2.10.0_main_441d081927")
if [ -z "${BRANCH}" ] || [ "${DISTANCE}" = "0" ]; then
  export LANCET_VERSION="${BASE_VERSION}"
else
  export LANCET_VERSION="${BASE_VERSION}_${BRANCH}_${HASH}"
fi
export LANCET_BUILD_STRING="${BRANCH:-local}_${HASH:-0000000000}_0"

echo "LANCET_VERSION:      ${LANCET_VERSION}"
echo "LANCET_BUILD_STRING: ${LANCET_BUILD_STRING}"
echo "RECIPE:              ${RECIPE}"
echo "OUTPUT_DIR:          ${OUTPUT_DIR}"
echo ""

pixi exec rattler-build build \
  --recipe "${RECIPE}" \
  --output-dir "${OUTPUT_DIR}"
