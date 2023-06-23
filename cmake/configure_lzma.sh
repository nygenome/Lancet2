#!/usr/bin/env bash

set -euxo pipefail

readonly ROOT_DIR="${1}"
readonly C_COMPILER="${2}"

readonly LANCET_OPT_FLAGS=$(grep 'LANCET_OPT_FLAGS' "${ROOT_DIR}/../../CMakeCache.txt" | cut -d= -f2-)

echo "LZMA CONFIGURE ROOT_DIR : ${ROOT_DIR}"
echo "LZMA CONFIGURE C_COMPILER : ${C_COMPILER}"

env CC="${C_COMPILER}" CFLAGS="${LANCET_OPT_FLAGS}" "${ROOT_DIR}"/configure --prefix="${ROOT_DIR}" \
  --enable-static --disable-shared --disable-doc --disable-scripts --disable-xz \
  --disable-xzdec --disable-lzmadec --disable-lzmainfo --disable-lzma-links
