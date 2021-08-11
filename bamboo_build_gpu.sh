#!/usr/bin/env bash

type module >& /dev/null || . /mnt/software/Modules/current/init/bash

set -e

################
# DEPENDENCIES #
################

## Load modules
set +vx
source scripts/ci/modules.sh
set -vx

export CC="ccache gcc"
export CXX="ccache g++"
export CCACHE_BASEDIR="${PWD}"

if [[ -z ${bamboo_planRepository_branchName+x} ]]; then
  : #pass
elif [[ ! -d /pbi/flash/bamboo/ccachedir ]]; then
  echo "[WARNING] /pbi/flash/bamboo/ccachedir is missing"
elif [[ $bamboo_planRepository_branchName == develop ]]; then
  export CCACHE_DIR=/pbi/flash/bamboo/ccachedir/${bamboo_shortPlanKey}.${bamboo_shortJobKey}.develop
  export CCACHE_TEMPDIR=/scratch/bamboo.ccache_tempdir
elif [[ $bamboo_planRepository_branchName == master ]]; then
  export CCACHE_DIR=/pbi/flash/bamboo/ccachedir/${bamboo_shortPlanKey}.${bamboo_shortJobKey}.master
  export CCACHE_TEMPDIR=/scratch/bamboo.ccache_tempdir
elif [[ $USER == bamboo ]]; then
  _shortPlanKey=$(echo ${bamboo_shortPlanKey}|sed -e 's/[0-9]*$//')
  export CCACHE_DIR=/pbi/flash/bamboo/ccachedir/${bamboo_shortPlanKey}.${bamboo_shortJobKey}
  if [[ -d /pbi/flash/bamboo/ccachedir/${_shortPlanKey}.${bamboo_shortJobKey}.develop ]]; then
    ( cd /pbi/flash/bamboo/ccachedir/
      cp -a ${_shortPlanKey}.${bamboo_shortJobKey}.develop $CCACHE_DIR
    )
  fi
  export CCACHE_TEMPDIR=/scratch/bamboo.ccache_tempdir
fi

case "${bamboo_planRepository_branchName}" in
  develop|master)
    export PREFIX_ARG="/mnt/software/p/pancake/${bamboo_planRepository_branchName}"
    export BUILD_NUMBER="${bamboo_globalBuildNumber:-0}"
    ;;
  *)
    export BUILD_NUMBER="0"
    ;;
esac

####################################
# Get a custom Meson for CUDA.
git clone https://github.com/SoapGentoo/meson.git -b cuda-fixes
pushd meson
ln -sf meson.py meson
export PATH=$(pwd):$PATH
popd

module load python/3 ninja

export ENABLED_GPU_CUDA="true"
####################################

# call the main build+test scripts
export ENABLED_TESTS="true"
export ENABLED_INTERNAL_TESTS="${bamboo_ENABLED_INTERNAL_TESTS}"
# export LDFLAGS="-static-libstdc++ -static-libgcc"

source scripts/ci/configure_with_fallback.sh
source scripts/ci/build.sh
source scripts/ci/test.sh
 
echo "Not installing anything (branch: ${bamboo_planRepository_branchName}), exiting."

# if [[ -z ${PREFIX_ARG+x} ]]; then
#   echo "Not installing anything (branch: ${bamboo_planRepository_branchName}), exiting."
#   exit 0
# fi

# source scripts/ci/install.sh

# Also install racon, since it lacks its own PB repo.
# export PREFIX_ARG="/mnt/software/r/racon/${bamboo_planRepository_branchName}"
# set +vx
# source scripts/ci/racon.modules.sh
# set -vx
# rm -rf racon-v1.4.13/build-meson # TEMP fix for my mistake
# make -C scripts/ci all
