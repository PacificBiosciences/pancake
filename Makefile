CURRENT_BUILD_DIR?=build
CURRENT_DEBUG_BUILD_DIR?=build-debug
CURRENT_DEBUG_BUILD_DIR_SANITIZE?=build-debug-sanitize
ENABLED_TESTS?=true
export ENABLED_TESTS CURRENT_BUILD_DIR

.PHONY: all build conf conf-debug unit cram modules check-formatting build-debug build-debug2 conf-debug2 debug2



################
### Release. ###
################
build:
	mkdir -p ${CURRENT_BUILD_DIR} && ninja -C "${CURRENT_BUILD_DIR}" -v test

conf:
	rm -rf "${CURRENT_BUILD_DIR}"
	ENABLED_GPU_CUDA=false ENABLED_GPU_CUDA_GW=false ENABLED_SSE41=true ENABLED_TESTS=true bash -vex scripts/ci/configure_with_fallback.sh

conf-gpu:
	rm -rf "${CURRENT_BUILD_DIR}"
	ENABLED_GPU_CUDA=true ENABLED_GPU_CUDA_GW=false ENABLED_SSE41=true ENABLED_TESTS=true bash -vex scripts/ci/configure_with_fallback.sh

conf-gcc:
	rm -rf "${CURRENT_BUILD_DIR}"
	ENABLED_GPU_CUDA=false ENABLED_GPU_CUDA_GW=false ENABLED_SSE41=true ENABLED_TESTS=true CC=gcc-10 CXX=g++-10 PREFIX_ARG=$(CURDIR)/INSTALL bash -vex scripts/ci/configure_with_fallback.sh

all: conf build cram check-formatting
all-gpu: conf-gpu build cram check-formatting

################
### Debug.   ###
################
build-debug:
	mkdir -p ${CURRENT_DEBUG_BUILD_DIR} && ninja -C "${CURRENT_DEBUG_BUILD_DIR}" -v test

conf-debug:
	rm -rf "${CURRENT_DEBUG_BUILD_DIR}"
	bash -vex scripts/ci/configure_debug_fallback.sh ${CURRENT_DEBUG_BUILD_DIR}

debug: conf-debug build-debug



build-debug2:
	mkdir -p ${CURRENT_DEBUG_BUILD_DIR_SANITIZE} && ninja -C "${CURRENT_DEBUG_BUILD_DIR_SANITIZE}" -v test

conf-debug2:
	rm -rf "${CURRENT_DEBUG_BUILD_DIR_SANITIZE}"
	ENABLED_SSE41=true ENABLED_TESTS=true CURRENT_BUILD_DIR="${CURRENT_DEBUG_BUILD_DIR_SANITIZE}" bash -vex scripts/ci/configure_debug_sanitize_fallback.sh

debug2: conf-debug2 build-debug2



build-debug3:
	mkdir -p ${CURRENT_DEBUG_BUILD_DIR} && ninja -C "${CURRENT_DEBUG_BUILD_DIR}" -v test

conf-debug3:
	rm -rf "${CURRENT_DEBUG_BUILD_DIR}"
	meson --buildtype=debugoptimized --default-library shared --libdir lib -Dtests="true" -Dgpu-cuda="false" -Dsse41="true" -Dcpp_debugstl="true" "${CURRENT_DEBUG_BUILD_DIR}"

debug3: conf-debug3 build-debug3

##############
### Tests. ###
##############
unit:
	ninja -C "${CURRENT_BUILD_DIR}" -v test

cram: modules
	scripts/cram tests/cram/*.t

check-formatting:
	tools/check-formatting --all

##############
### Other. ###
##############
modules:
	git submodule update --init --recursive
