#!/usr/bin/env bash

module purge

module load gcc
module load ccache

module load meson
module load ninja

module load boost
module load cram/0.7-py3
module load python3
module load gtest

module load htslib
module load pbbam

module load samtools

module load cuda

case "${bamboo_planRepository_branchName}" in
  master)
    module load pbcopper/master
    ;;
  *)
    module load pbcopper/develop
    ;;
esac
