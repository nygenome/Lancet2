#!/usr/bin/env bash
# shellcheck disable=SC2155

set -e

function ensure_deps() {
  go install github.com/flaticols/semtag@latest
  go install github.com/git-chglog/git-chglog/cmd/git-chglog@latest
}

function main() {
  local PROJECT_DIR=$(git rev-parse --show-toplevel)
  local BUMP_KIND="${1:-patch}"

  echo "BUMP_KIND: ${BUMP_KIND}"
  git push && semtag --brave --no-tty "${BUMP_KIND}" && git-chglog --config "${PROJECT_DIR}/.chglog/config.yml" -o "${PROJECT_DIR}"/CHANGELOG.md && \
  git add "${PROJECT_DIR}"/CHANGELOG.md && git commit -vsm "chore: Update release notes in changelog" && git push --tags
}

export PATH=$PATH:$GOBIN
ensure_deps && main "$1"
