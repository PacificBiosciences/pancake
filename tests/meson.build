pancake_test_cpp_sources = files([
  'src/test_AlignerBatchCPU.cpp',
  'src/test_AlignmentSeeded.cpp',
  'src/test_AlignmentTools.cpp',
  'src/test_DPChain.cpp',
  'src/test_FastaSequenceCachedStore.cpp',
  'src/test_FileIO.cpp',
  'src/test_LIS.cpp',
  'src/test_MapperBatchCPU.cpp',
  'src/test_MapperCLR.cpp',
  'src/test_MapperHiFi.cpp',
  'src/test_MapperUtility.cpp',
  'src/test_Minimizers.cpp',
  'src/test_Overlap.cpp',
  'src/test_Pancake.cpp',
  'src/test_RangeTools.cpp',
  'src/test_RunLengthEncoding.cpp',
  'src/test_Secondary.cpp',
  'src/test_SeedIndex.cpp',
  'src/test_SeedDBReader.cpp',
  'src/test_SeedDBReaderCached.cpp',
  'src/test_SeedDBReaderCachedBlock.cpp',
  'src/test_SeedDBReaderRawBlock.cpp',
  'src/test_SeqDBReader.cpp',
  'src/test_SeqDBReaderCached.cpp',
  'src/test_SeqDBReaderCachedBlock.cpp',
  'src/test_SesAlignBanded.cpp',
  'src/test_SesDistanceBanded.cpp',
  'src/test_Ses2AlignBanded.cpp',
  'src/test_Ses2DistanceBanded.cpp',
  'src/test_Twobit.cpp',
  'src/test_Util.cpp',
  'src/TestHelperUtils.cpp',
])

if opt_cuda
    pancake_test_cpp_sources += files([
      'src/test_AlignerBatchGPU.cpp',
      'src/test_MapperBatchGPU.cpp',
    ])
endif
if opt_sse41
    pancake_test_cpp_sources += files([
        'src/test_DPChainSimd.cpp',
    ])
endif

pancake_test_cram_sources = files([
  'cram/test_bugfixes.t',
  'cram/test_hifi_ovl.t',
  'cram/test_hifi_ovl_mapping.t',
  'cram/test_hifi_ovl_multiple_hits_per_target.t',
  'cram/test_robustness.t',
  'cram/test_seqdb_dump.t',
  'cram/test_seqdb_info.t',
  'cram/test_seqdb_filter.t',
  'cram/test_seqdb_writer.t',
  'cram/test_seeddb_writer.t',
  'cram/test_seqfetch.t',
])

# find GoogleTest
pancake_gtest_dep = dependency('gtest_main', fallback : ['gtest', 'gtest_main_dep'])

pancake_pancakeTestData_h_config = configuration_data()
pancake_pancakeTestData_h_config.set('pancake_TestsDir', meson.current_source_dir())
pancake_pancakeTestData_h_config.set('CMAKE_CURRENT_BINARY_DIR', meson.current_build_dir())
pancake_pancakeTestData_h_config.set('GeneratedDir', meson.current_build_dir())
pancake_pancakeTestData_h_config.set('GeneratedTestDataDir', meson.current_build_dir())

pancake_pancakeTestData_h = configure_file(
  input : files('include/PancakeTestData.h.in'),
  output : 'PancakeTestData.h',
  configuration : pancake_pancakeTestData_h_config)
pancake_test_cpp_sources += pancake_pancakeTestData_h

pancake_test = executable(
  'pancake_test', [
    pancake_test_cpp_sources],
  dependencies : [pancake_gtest_dep] + pancake_lib_deps,
  include_directories : [pancake_include_directories, include_directories('include')],
  link_with : [pancake_lib],
  cpp_args : pancake_warning_flags,
  install : false)

#########
# tests #
#########

test(
  'pancake formatting check',
  pancake_clang_formatter,
  args : [
    '--all'],
  workdir : meson.source_root())

test(
  'pancake gtest unittests',
  pancake_test,
  args : [
    '--gtest_output=xml:' + join_paths(meson.build_root(), 'pancake-gtest-unittests.xml')],
  env : [
    'ARGS=-V',
    'VERBOSE=1'])

test(
  'pancake cram test',
  pancake_cram_script,
  args : [
    '--xunit-file=' + join_paths(meson.build_root(), 'pancake-cram-ccs.xml'),
    '--keep-tmpdir',
    '--verbose'] +
    pancake_test_cram_sources,
  env : [
    'BIN_DIR=' + join_paths(meson.build_root(), 'src'),
    'PROJECT_DIR=' + join_paths([meson.current_source_dir(), '..'])],
  timeout : 36000) # with '-O0 -g' tests can be *very* slow
