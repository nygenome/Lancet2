#!/usr/bin/env bash
# shellcheck disable=SC2155

set -e

function ensure_deps() {
  go install github.com/threecommaio/sembump@latest
  go install github.com/git-chglog/git-chglog/cmd/git-chglog@latest
}

function main() {
  local PROJECT_DIR=$(git rev-parse --show-toplevel)
  local BUMP_KIND="${1:-patch}"

  local MOST_RECENT_TAG=$(git describe --abbrev=0 --tags | cut -d'.' -f1-3)
  MOST_RECENT_TAG="${MOST_RECENT_TAG:-0.0.0}"
  echo "MOST_RECENT_TAG VERSION: ${MOST_RECENT_TAG}"

  echo "BUMP_KIND: ${BUMP_KIND}"
  echo "sembump ${MOST_RECENT_TAG} ${BUMP_KIND}"
  local NEW_VERSION=$(sembump "${MOST_RECENT_TAG}"  "${BUMP_KIND}")

  if [[ -z ${NEW_VERSION} ]]; then
      exit 1
  fi

  echo "Bumping version from ${MOST_RECENT_TAG} to ${NEW_VERSION}"
  git-chglog --next-tag "${NEW_VERSION}" -o "${PROJECT_DIR}"/CHANGELOG.md
  git add "${PROJECT_DIR}"/CHANGELOG.md
  git commit -vsm "chore: Bump version to ${NEW_VERSION}"
  git tag -m "${NEW_VERSION}" -sa "${NEW_VERSION}"
}

ensure_deps && main "$1"
