#!/bin/bash

set -euxo pipefail

readonly SCRIPT_DIR=$(dirname -- "$(readlink -f -- "${0}")")
readonly GCS_OAUTH_TOKEN="$(gcloud auth print-access-token)"

export GCS_OAUTH_TOKEN && \
gcloud storage cp "gs://lancet2-test-datasets/SEQC2_chr4_test_data/**" "${SCRIPT_DIR}"
