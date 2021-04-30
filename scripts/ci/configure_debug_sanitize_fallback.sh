#!/usr/bin/env bash
set -vex

BUILD_DIR=$1

# configure
# '--wrap-mode nofallback' prevents meson from downloading
# stuff from the internet or using subprojects.
meson \
  --default-library shared \
  --libdir lib \
  --unity "${ENABLED_UNITY_BUILD:-off}" \
  --prefix "${PREFIX_ARG:-/usr/local}" \
  --buildtype=debugoptimized \
  -Db_sanitize=address \
  -Dtests="${ENABLED_TESTS:-false}" \
  -Dgpu-cuda="${ENABLED_GPU_CUDA:-false}" \
  "${BUILD_DIR}" .
