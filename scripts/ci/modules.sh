#!/usr/bin/env bash

module purge

module load gcc
module load ccache

module load meson
module load ninja

module load boost
module load cram
module load gtest

module load htslib
module load pbbam

module load samtools

module load cuda/11.2.0_460.27.04

case "${bamboo_planRepository_branchName}" in
  master)
    module load pbcopper/master
    ;;
  *)
    module load pbcopper/develop
    ;;
esac
