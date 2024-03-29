# Copyright (c) 2019, Pacific Biosciences of California, Inc.
# All rights reserved.
# See LICENSE.txt.
#
# Contributions from NVIDIA are Copyright (c) 2021, NVIDIA Corporation.
# All rights reserved.
# SPDX-License-Identifier: BSD-3-Clause-Clear

################
# version info #
################
pancake_ver_arr = meson.project_version().split('.')
pancake_major_version = pancake_ver_arr[0]
pancake_minor_version = pancake_ver_arr[1]
pancake_patch_version = pancake_ver_arr[2]

pancake_config = configuration_data()
pancake_config.set('PANCAKE_VERSION', meson.project_version())
pancake_config.set('PANCAKE_VERSION_MAJOR', pancake_major_version)
pancake_config.set('PANCAKE_VERSION_MINOR', pancake_minor_version)
pancake_config.set('PANCAKE_VERSION_PATCH', pancake_patch_version)

pancake_version_sources = [
  vcs_tag(
    input : 'LibraryGitHash.hpp.in',
    output : 'LibraryGitHash.hpp',
    command : ['git', 'describe', '--always', '--dirty=*'],
    replace_string : '@PANCAKE_GIT_SHA1@'),
  configure_file(
    input : 'LibraryVersion.hpp.in',
    output : 'LibraryVersion.hpp',
    configuration : pancake_config),

]

###########
# library #
###########

pancake_cpp_sources = files([
    'LibraryInfo.cpp',
    'main/dbfilter/DBFilterSettings.cpp',
    'main/dbfilter/DBFilterWorkflow.cpp',
    'main/mapclr/MapCLRWorkflow.cpp',
    'main/mapclr/MapCLRSettings.cpp',
    'main/overlaphifi/OverlapHifiSettings.cpp',
    'main/overlaphifi/OverlapHifiWorkflow.cpp',
    'main/seeddb/SeedDBSettings.cpp',
    'main/seeddb/SeedDBWorkflow.cpp',
    'main/seqdb/SeqDBSettings.cpp',
    'main/seqdb/SeqDBWorkflow.cpp',
    'main/seqdbdump/SeqDBDumpSettings.cpp',
    'main/seqdbdump/SeqDBDumpWorkflow.cpp',
    'main/seqdbinfo/SeqDBInfoSettings.cpp',
    'main/seqdbinfo/SeqDBInfoWorkflow.cpp',
    'main/seqfetch/SeqFetchSettings.cpp',
    'main/seqfetch/SeqFetchWorkflow.cpp',
    'pancake/AdapterFinder.cpp',
    'pancake/AlignmentTools.cpp',
    'pancake/AlignerBase.cpp',
    'pancake/AlignerBatchCPU.cpp',
    'pancake/AlignerKSW2.cpp',
    'pancake/AlignerEdlib.cpp',
    'pancake/AlignerSES1.cpp',
    'pancake/AlignerSES2.cpp',
    'pancake/AlignerFactory.cpp',
    'pancake/AlignmentSeeded.cpp',
    'pancake/CCSAlignmentResult.cpp',
    'pancake/CompressedSequence.cpp',
    'pancake/DPChain.cpp',
    'pancake/FastaSequenceId.cpp',
    'pancake/MapperBase.cpp',
    'pancake/MapperBatchBase.cpp',
    'pancake/MapperBatchCPU.cpp',
    'pancake/MapperBatchUtility.cpp',
    'pancake/MapperCLR.cpp',
    'pancake/MapperHiFi.cpp',
    'pancake/MapperUtility.cpp',
    'pancake/Minimizers.cpp',
    'pancake/Overlap.cpp',
    'pancake/OverlapWriterBase.cpp',
    'pancake/OverlapWriterFactory.cpp',
    'pancake/OverlapWriterIPAOvl.cpp',
    'pancake/OverlapWriterM4.cpp',
    'pancake/OverlapWriterPAF.cpp',
    'pancake/OverlapWriterSAM.cpp',
    'pancake/Secondary.cpp',
    'pancake/SeedHit.cpp',
    'pancake/SeedHitWriter.cpp',
    'pancake/SeedDBIndexCache.cpp',
    'pancake/SeedDBReader.cpp',
    'pancake/SeedDBReaderCached.cpp',
    'pancake/SeedDBReaderCachedBlock.cpp',
    'pancake/SeedDBReaderRawBlock.cpp',
    'pancake/SeedDBWriter.cpp',
    'pancake/SeedIndex.cpp',
    'pancake/SeqDBIndexCache.cpp',
    'pancake/SeqDBReader.cpp',
    'pancake/SeqDBReaderCached.cpp',
    'pancake/SeqDBReaderCachedBlock.cpp',
    'pancake/SeqDBWriter.cpp',
    'pancake/SequenceSeeds.cpp',
    'pancake/SequenceSeedsCached.cpp',
    'pancake/SesDistanceBanded.cpp',
    'pancake/Twobit.cpp',
    'pancake/util/FileIO.cpp',
    'pancake/util/RunLengthEncoding.cpp',
    'pancake/util/SeqLengthStats.cpp',
    'pancake/util/TicToc.cpp',
])

opt_cuda = (opt_cuda_edelweiss or opt_cuda_genomeworks)

if opt_cuda
    cuda = import('unstable-cuda')
    pancake_cuda_flags = cuda.nvcc_arch_flags(nvcc, 'Auto')
    arch_readable = cuda.nvcc_arch_readable(nvcc, 'Auto')
    driver_version = cuda.min_driver_version('>=11.3.0')
    message('NVCC version:   ' + nvcc.version())
    message('NVCC flags:     ' + ' '.join(pancake_cuda_flags))
    message('NVCC readable:  ' + ' '.join(arch_readable))
    message('Driver version: >=' + driver_version)

    if (get_option('b_ndebug') == 'true') or ((get_option('b_ndebug') == 'if-release') and (get_option('buildtype') == 'release'))
      pancake_cuda_flags += '-DNDEBUG'
    endif

    pancake_cpp_sources += files([
        'pancake/MapperBatchGPUUtility.cpp',
    ])
endif

if opt_cuda_edelweiss
    pancake_cpp_sources += files([
        'pancake/AlignerBatchGPUEdelweiss.cpp',
        'pancake/MapperBatchGPUEdelweiss.cpp',
    ])
endif

if opt_cuda_genomeworks
    pancake_cpp_sources += files([
        'pancake/AlignerBatchGPU.cpp',
        'pancake/MapperBatchGPU.cpp',
    ])
    pancake_cu_sources = files([
        'pancake/AlignerBatchGPUGenomeWorksInterface.cu',
    ])
    pancake_cpp_sources += pancake_cu_sources
endif

if opt_sse41
    pancake_cpp_sources += files([
    'pancake/DPChainSimd.cpp',
    ])
endif

pancake_cpp_sources += pancake_version_sources

###########
# Library #
###########
# install library if
# - either running as a proper project
# - or using shared libraries
pancake_lib_install = (not meson.is_subproject()) or (get_option('default_library') == 'shared')

pancake_exe_dir = meson.current_source_dir()

if opt_cuda
    pancake_lib = library(
      'pancake',
      pancake_cpp_sources,
      soversion : meson.project_version(),
      version : meson.project_version(),
      install : pancake_lib_install,
      dependencies : pancake_lib_deps,
      include_directories : [pancake_include_directories, pancake_third_party_include_directories],
      cpp_args : pancake_flags,
      cuda_args : pancake_cuda_flags,
      override_options : [
        'cuda_std=c++17',
        'warning_level=2'])
else
    pancake_lib = library(
      'pancake',
      pancake_cpp_sources,
      soversion : meson.project_version(),
      version : meson.project_version(),
      install : pancake_lib_install,
      dependencies : pancake_lib_deps,
      include_directories : [pancake_include_directories, pancake_third_party_include_directories],
      cpp_args : pancake_flags)
endif

##############
# executable #
##############

if not meson.is_subproject()
  pancake_main_exe = executable(
    'pancake', files([
      'main/main.cpp',
    ]),
    install : true,
    dependencies : pancake_lib_deps,
    include_directories : [pancake_include_directories, pancake_third_party_include_directories],
    link_with : pancake_lib,
    cpp_args : pancake_flags)
endif
