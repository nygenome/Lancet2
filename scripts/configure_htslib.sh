#!/usr/bin/env bash

set -e

readonly HTSLIB_ROOT_DIR="${1}"
readonly C_COMPILER="${2}"
readonly LIBDEFLATE_ROOT_DIR=$(realpath "${HTSLIB_ROOT_DIR}"/../libdeflate)

echo "  HTSLIB CONFIGURE ROOT_DIR : ${HTSLIB_ROOT_DIR}"
echo "HTSLIB CONFIGURE C_COMPILER : ${C_COMPILER}"

env CC="${C_COMPILER}" CFLAGS="-O3 -DNDEBUG" "${HTSLIB_ROOT_DIR}"/configure \
    --prefix="${HTSLIB_ROOT_DIR}" --enable-libcurl --enable-s3 --enable-gcs --disable-plugins \
    CPPFLAGS="-I${LIBDEFLATE_ROOT_DIR}" LDFLAGS="$LDFLAGS -L${LIBDEFLATE_ROOT_DIR}" --with-libdeflate
