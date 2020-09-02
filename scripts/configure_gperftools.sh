#!/usr/bin/env bash

set -e

readonly ROOT_DIR="${1}"
readonly C_COMPILER="${2}"
readonly CXX_COMPILER="${3}"
readonly CFG_CFLAGS="-O3 -DNDEBUG"

echo "  GPERFTOOLS CONFIGURE ROOT_DIR : ${ROOT_DIR}"
echo "GPREFTOOLS CONFIGURE C_COMPILER : ${C_COMPILER}"

env CC="${C_COMPILER}" CXX="${CXX_COMPILER}" CFLAGS="${CFG_CFLAGS}" "${ROOT_DIR}"/configure --prefix="${ROOT_DIR}"
