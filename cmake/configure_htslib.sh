#!/usr/bin/env bash

set -euxo pipefail

readonly ROOT_DIR="${1}"
readonly C_COMPILER="${2}"

readonly ZLIBNG_BUILD_DIR=$(realpath "${ROOT_DIR}/../zlib-ng-build")
readonly LIBDEFLATE_INC_DIR=$(realpath "${ROOT_DIR}/../libdeflate-src")
readonly LIBDEFLATE_LIB_DIR=$(realpath "${ROOT_DIR}/../libdeflate-build")

readonly LANCET_OPT_FLAGS=$(grep 'LANCET_OPT_FLAGS' "${ROOT_DIR}/../../CMakeCache.txt" | cut -d= -f2-)

echo "HTSLIB CONFIGURE ROOT_DIR : ${ROOT_DIR}"
echo "HTSLIB CONFIGURE C_COMPILER : ${C_COMPILER}"

env CC="${C_COMPILER}" CFLAGS="${LANCET_OPT_FLAGS}" \
  "${ROOT_DIR}"/configure --prefix="${ROOT_DIR}" --disable-libcurl --disable-plugins --with-libdeflate \
  CPPFLAGS="-I${ZLIBNG_BUILD_DIR} -I${LIBDEFLATE_INC_DIR}" LDFLAGS="-L${ZLIBNG_BUILD_DIR} -L${LIBDEFLATE_LIB_DIR}"
