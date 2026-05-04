#!/usr/bin/env bash
# shellcheck disable=SC2155

set -e

function generate_tag() {
  local SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  local REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
  local VERSION_TAG=$(cat "${REPO_ROOT}/VERSION.txt" | tr -d '[:space:]')
  local GIT_BRANCH=$(git rev-parse --abbrev-ref HEAD)
  local GIT_COMMIT=$(git rev-parse --short=10 --verify HEAD)
  local GIT_STATE=$(git diff --quiet || echo '-dirty')
  local GIT_TAG=$(git describe --abbrev=0 --tags 2>/dev/null || echo "")
  local DISTANCE=$(git rev-list "${GIT_TAG}..HEAD" --count 2>/dev/null || echo "0")

  # Stable release: clean semver (e.g. "2.9.0") when HEAD is at a tag
  # and the working tree is clean. Mirrors deploy_prefix.yml logic.
  # Dev build: version_branch_hash (e.g. "2.9.0_main_441d081927")
  if [ "${DISTANCE}" = "0" ] && [ -z "${GIT_STATE}" ]; then
    echo "${VERSION_TAG}"
  else
    echo "${VERSION_TAG}_${GIT_BRANCH}_${GIT_COMMIT}${GIT_STATE}"
  fi
}

gcloud builds submit \
  --project nygc-comp-s-fd4e \
  --machine-type n1-highcpu-32 \
  --timeout 60m \
  --gcs-log-dir="gs://nygc-comp-s-fd4e_cloudbuild/logs" \
  --service-account "projects/nygc-comp-s-fd4e/serviceAccounts/lancet2-build@nygc-comp-s-fd4e.iam.gserviceaccount.com" \
  --tag "us-central1-docker.pkg.dev/nygc-app-c-148c/lancet-public/lancet:$(generate_tag)"
