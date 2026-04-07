#!/usr/bin/env bash

set -euxo pipefail

readonly ROOT_DIR="${1}"
readonly C_COMPILER="${2}"
readonly ENABLE_CLOUD_IO="${3:-OFF}"

readonly ZLIBNG_BUILD_DIR=$(realpath "${ROOT_DIR}/../zlib-ng-build")
readonly LIBDEFLATE_INC_DIR=$(realpath "${ROOT_DIR}/../libdeflate-src")
readonly LIBDEFLATE_LIB_DIR=$(realpath "${ROOT_DIR}/../libdeflate-build")

readonly LANCET_OPT_FLAGS=$(grep 'LANCET_OPT_FLAGS' "${ROOT_DIR}/../../CMakeCache.txt" | cut -d= -f2-)

echo "HTSLIB CONFIGURE ROOT_DIR : ${ROOT_DIR}"
echo "HTSLIB CONFIGURE C_COMPILER : ${C_COMPILER}"
echo "HTSLIB CONFIGURE ENABLE_CLOUD_IO : ${ENABLE_CLOUD_IO}"

CLOUD_FLAGS="--disable-libcurl"
if [ "${ENABLE_CLOUD_IO}" = "ON" ]; then
    CLOUD_FLAGS="--enable-libcurl --enable-s3 --enable-gcs"
fi

# Note: Use ${CPPFLAGS:-} and ${LDFLAGS:-} to safely inherit existing flags if they are set. 
# This is explicitly required for Conda/Rattler builds to correctly dynamically inject downstream 
# dependency headers/paths (e.g. libcurl, openssl within $PREFIX) without throwing `set -u` unbound locally.
env CC="${C_COMPILER}" CFLAGS="${LANCET_OPT_FLAGS:-}" \
  "${ROOT_DIR}"/configure --prefix="${ROOT_DIR}" ${CLOUD_FLAGS} --disable-plugins --with-libdeflate \
  CPPFLAGS="${CPPFLAGS:-} -I${ZLIBNG_BUILD_DIR} -I${LIBDEFLATE_INC_DIR}" \
  LDFLAGS="${LDFLAGS:-} -L${ZLIBNG_BUILD_DIR} -L${LIBDEFLATE_LIB_DIR}"
