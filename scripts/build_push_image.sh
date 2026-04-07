#!/usr/bin/env bash
# shellcheck disable=SC2155

set -e

function generate_tag() {
  local VERSION_TAG=$(git describe --abbrev=0 --tags | cut -d'.' -f1-3)
  local GIT_BRANCH=$(git rev-parse --abbrev-ref HEAD)
  local GIT_COMMIT=$(git rev-parse --short=10 --verify HEAD)
  local GIT_STATE=$(git diff --quiet || echo '-dirty')
  local BUILD_TAG="${VERSION_TAG}-${GIT_BRANCH}-${GIT_COMMIT}${GIT_STATE}"
  echo "${BUILD_TAG}"
}

gcloud builds submit \
  --project nygc-comp-s-fd4e \
  --machine-type n1-highcpu-32 \
  --timeout 60m \
  --gcs-log-dir="gs://nygc-comp-s-fd4e_cloudbuild/logs" \
  --service-account "projects/nygc-comp-s-fd4e/serviceAccounts/lancet2-build@nygc-comp-s-fd4e.iam.gserviceaccount.com" \
  --tag "us-central1-docker.pkg.dev/nygc-app-c-148c/lancet-public/lancet:$(generate_tag)"
