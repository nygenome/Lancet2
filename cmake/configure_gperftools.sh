#!/usr/bin/env bash

set -euxo pipefail

readonly ROOT_DIR="${1}"
readonly C_COMPILER="${2}"
readonly CXX_COMPILER="${3}"

readonly LANCET_OPT_FLAGS=$(grep 'LANCET_OPT_FLAGS' "${ROOT_DIR}/../../CMakeCache.txt" | cut -d= -f2-)
readonly CFLAGS="${LANCET_OPT_FLAGS}"

echo "GPERFTOOLS CONFIGURE ROOT_DIR : ${ROOT_DIR}"
echo "GPERFTOOLS CONFIGURE C_COMPILER : ${C_COMPILER}"
echo "GPERFTOOLS CONFIGURE CXX_COMPILER : ${CXX_COMPILER}"

env CC="${C_COMPILER}" CXX="${CXX_COMPILER}" CFLAGS="${CFLAGS}" "${ROOT_DIR}"/configure --prefix="${ROOT_DIR}" \
  --enable-static --disable-shared --disable-heap-profiler --disable-heap-checker --disable-debugalloc
