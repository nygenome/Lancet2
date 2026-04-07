#!/usr/bin/env bash
# shellcheck disable=SC2155

set -e

function ensure_deps() {
  local PROJECT_DIR=$(git rev-parse --show-toplevel)
  echo "Building git-chglog from source..."
  (cd "${PROJECT_DIR}/cmake-build-release/git-chglog" && go build -o "${GOBIN:-$HOME/go/bin}/git-chglog" ./cmd/git-chglog)
}

function main() {
  local PROJECT_DIR=$(git rev-parse --show-toplevel)

  echo "Updating changelog..."
  git-chglog --config "${PROJECT_DIR}/.chglog/config.yml" -o "${PROJECT_DIR}"/CHANGELOG.md
  
  git add "${PROJECT_DIR}"/CHANGELOG.md 
  git commit -vsm "chore: Update release notes in changelog"
  git push
  
  echo "Changelog updated and pushed successfully."
}

export PATH="$PATH:${GOBIN:-$HOME/go/bin}"
ensure_deps && main "$@"
