#!/usr/bin/env bash

set -euxo pipefail

readonly ROOT_DIR="${1}"
readonly C_COMPILER="${2}"
readonly PROCESSOR="${3:-x86_64}"

readonly ZLIBNG_BUILD_DIR=$(realpath "${ROOT_DIR}/../zlib-ng-build")
readonly LANCET_OPT_FLAGS=$(grep 'LANCET_OPT_FLAGS' "${ROOT_DIR}/../../CMakeCache.txt" | cut -d= -f2-)

echo "MINIMAP2 BUILD ROOT_DIR  : ${ROOT_DIR}"
echo "MINIMAP2 BUILD C_COMPILER  : ${C_COMPILER}"
echo "MINIMAP2 BUILD PROCESSOR : ${PROCESSOR}"

MAKE_ARGS=""
if [[ "${PROCESSOR}" == "arm64" || "${PROCESSOR}" == "aarch64" ]]; then
  MAKE_ARGS="arm_neon=1 aarch64=1"
fi

cd "${ROOT_DIR}" && env CC="${C_COMPILER}" \
  make libminimap2.a ${MAKE_ARGS} CFLAGS="${LANCET_OPT_FLAGS} -I${ZLIBNG_BUILD_DIR} -L${ZLIBNG_BUILD_DIR}"
