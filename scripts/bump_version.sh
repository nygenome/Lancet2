#!/usr/bin/env bash

set -e

function ensure_deps() {
  go get -u github.com/jfinken/sembump
  go get -u github.com/git-chglog/git-chglog/cmd/git-chglog
}

function main() {
  # shellcheck disable=SC2155
  local readonly PROJECT_DIR=$(git rev-parse --show-toplevel)
  local readonly BUMP_KIND="${1:-patch}"

  # shellcheck disable=SC2155
  local MOST_RECENT_TAG=$(git describe --abbrev=0 --tags)
  MOST_RECENT_TAG="${MOST_RECENT_TAG:-0.0.0}"
  echo "MOST_RECENT_TAG: ${MOST_RECENT_TAG}"

  echo "BUMP_KIND: ${BUMP_KIND}"
  echo "sembump --kind ${BUMP_KIND} ${MOST_RECENT_TAG}"
  # shellcheck disable=SC2155
  # shellcheck disable=SC2034
  local readonly NEW_VERSION=$(sembump --kind "${BUMP_KIND}" "${MOST_RECENT_TAG}")

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
