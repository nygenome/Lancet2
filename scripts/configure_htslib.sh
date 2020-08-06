#!/usr/bin/env bash

set -e

readonly ROOT_DIR="${1}"
readonly C_COMPILER="${2}"
readonly CFG_CFLAGS="-O3 -DNDEBUG"

echo "  HTSLIB CONFIGURE ROOT_DIR : ${ROOT_DIR}"
echo "HTSLIB CONFIGURE C_COMPILER : ${C_COMPILER}"

env CC="${C_COMPILER}" CFLAGS="${CFG_CFLAGS}" "${ROOT_DIR}"/configure \
    --prefix="${ROOT_DIR}" --enable-libcurl --enable-s3 --enable-gcs --disable-plugins
