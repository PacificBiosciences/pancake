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
  -Dtests="${ENABLED_TESTS:-false}" \
  -Dsse41="${ENABLED_SSE41:-false}" \
  "${BUILD_DIR}" .
