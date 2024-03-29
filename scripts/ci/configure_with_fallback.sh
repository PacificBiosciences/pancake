#!/usr/bin/env bash
set -vex

# configure
# '--wrap-mode nofallback' prevents meson from downloading
# stuff from the internet or using subprojects.
meson \
  --buildtype release \
  --default-library shared \
  --libdir lib \
  --unity "${ENABLED_UNITY_BUILD:-off}" \
  --prefix "${PREFIX_ARG:-/usr/local}" \
  -Dtests="${ENABLED_TESTS:-false}" \
  -Dgpu-cuda="${ENABLED_GPU_CUDA:-false}" \
  -Dgpu-cuda-gw="${ENABLED_GPU_CUDA_GW:-false}" \
  -Dsse41="${ENABLED_SSE41:-false}" \
  "${CURRENT_BUILD_DIR:-build}" .
#  --werror \
