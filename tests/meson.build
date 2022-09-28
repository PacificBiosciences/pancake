pancake_test_cpp_sources = files([
  'unit/test_AdapterFinder.cpp',
  'unit/test_AlignerBatchCPU.cpp',
  'unit/test_AlignmentSeeded.cpp',
  'unit/test_AlignmentTools.cpp',
  'unit/test_DPChain.cpp',
  'unit/test_FastaSequenceCachedStore.cpp',
  'unit/test_FileIO.cpp',
  'unit/test_LCS.cpp',
  'unit/test_LIS.cpp',
  'unit/test_MapperBatchCPU.cpp',
  'unit/test_MapperCLR.cpp',
  'unit/test_MapperHiFi.cpp',
  'unit/test_MapperUtility.cpp',
  'unit/test_Minimizers.cpp',
  'unit/test_Overlap.cpp',
  'unit/test_Pancake.cpp',
  'unit/test_RangeTools.cpp',
  'unit/test_RunLengthEncoding.cpp',
  'unit/test_Secondary.cpp',
  'unit/test_SeedIndex.cpp',
  'unit/test_SeedDBReader.cpp',
  'unit/test_SeedDBReaderCached.cpp',
  'unit/test_SeedDBReaderCachedBlock.cpp',
  'unit/test_SeedDBReaderRawBlock.cpp',
  'unit/test_SeqDBReader.cpp',
  'unit/test_SeqDBReaderCached.cpp',
  'unit/test_SeqDBReaderCachedBlock.cpp',
  'unit/test_SesAlignBanded.cpp',
  'unit/test_SesDistanceBanded.cpp',
  'unit/test_Ses2AlignBanded.cpp',
  'unit/test_Ses2DistanceBanded.cpp',
  'unit/test_Twobit.cpp',
  'unit/test_Util.cpp',
  'unit/TestHelperUtils.cpp',
])

if opt_cuda_edelweiss
    pancake_test_cpp_sources += files([
      'unit/test_AlignerBatchGPUEdelweiss.cpp',
      'unit/test_MapperBatchGPUEdelweiss.cpp',
    ])
endif
if opt_cuda_genomeworks
    pancake_test_cpp_sources += files([
      'unit/test_AlignerBatchGPU.cpp',
      'unit/test_MapperBatchGPU.cpp',
    ])
endif
if opt_sse41
    pancake_test_cpp_sources += files([
        'unit/test_DPChainSimd.cpp',
    ])
endif

pancake_cram_tests = [
  'bugfixes',
  'hifi_ovl',
  'hifi_ovl_mapping',
  'hifi_ovl_multiple_hits_per_target',
  'map_clr',
  'robustness',
  'seqdb_dump',
  'seqdb_info',
  'seqdb_filter',
  'seqdb_writer',
  'seeddb_writer',
  'seqfetch',
]

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
  include_directories : [pancake_include_directories, pancake_third_party_include_directories, include_directories('include')],
  link_with : [pancake_lib],
  cpp_args : pancake_warning_flags + pancake_flags,
  install : false)

#########
# tests #
#########

test(
  'pancake formatting check',
  pancake_clang_formatter,
  args : [
    '--all'],
  workdir : meson.project_source_root())

test(
  'pancake gtest unittests',
  pancake_test,
  args : [
    '--gtest_output=xml:' + join_paths(meson.project_build_root(), 'pancake-gtest-unittests.xml')],
  timeout: 120)

foreach t : pancake_cram_tests
  test(
    'pancake cram test - ' + t,
    pancake_cram_script,
    args : [
      '--xunit-file=' + join_paths(meson.project_build_root(), 'pancake-cram-' + t.underscorify() + '-test.xml'),
      '--keep-tmpdir',
      '--verbose'] + files('cram/' + t + '.t'),
    env : {
      'BIN_DIR' : join_paths(meson.project_build_root(), 'src'),
      'PROJECT_DIR' : join_paths([meson.current_source_dir(), '..']),
      'TEST_DATA_DIR' : join_paths([meson.current_source_dir(), 'data']),
    },
    timeout : 36000) # with '-O0 -g' tests can be *very* slow
endforeach
