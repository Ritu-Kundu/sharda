#!/usr/bin/env bash
set -euxo pipefail

cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
  -DCMAKE_PREFIX_PATH="${PREFIX}" \
  -DHTSLIB_DIR="${PREFIX}" \
  -DHIGHS_DIR="${PREFIX}" \
  -DSHARDA_BUILD_TESTS=OFF \
  -DSHARDA_FETCH_DEPS=OFF

cmake --build build -j"${CPU_COUNT}"
cmake --install build
