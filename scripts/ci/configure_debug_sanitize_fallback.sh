#!/usr/bin/env bash
set -vex

# configure
# '--wrap-mode nofallback' prevents meson from downloading
# stuff from the internet or using subprojects.
meson \
  --buildtype=debugoptimized \
  --default-library shared \
  --libdir lib \
  --unity "${ENABLED_UNITY_BUILD:-off}" \
  --prefix "${PREFIX_ARG:-/usr/local}" \
  -Db_sanitize=address,undefined \
  -Dtests="${ENABLED_TESTS:-false}" \
  -Dgpu-cuda="${ENABLED_GPU_CUDA:-false}" \
  -Dsse41="${ENABLED_SSE41:-false}" \
  "${CURRENT_BUILD_DIR:-build}" .
