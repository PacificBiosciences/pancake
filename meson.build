project(
  'pancake',
  ['cpp', 'c'],
  version : '1.3.0',
  default_options : [
    'buildtype=release',
    'warning_level=3',
    'cpp_std=c++17',
    'c_std=c11',
    'b_ndebug=if-release'],
  license : 'BSD-3',
  meson_version : '>= 0.52.0')

############
# CXXFLAGS #
############

pancake_warning_flags = []
cpp = meson.get_compiler('cpp')
foreach cflag: [
  '-Wduplicated-cond',
  '-Wduplicated-branches',
  '-Wlogical-op',
  '-Wrestrict',
  '-Wnull-dereference',
  '-Wuseless-cast',
  '-Wdouble-promotion',
  '-Wshadow',
  '-fno-new-ttp-matching',
  '-Wformat=1']
    if cpp.has_argument(cflag)
      pancake_warning_flags += cflag
    endif
endforeach

pancake_flags = pancake_warning_flags + ['-DHAVE_KALLOC']

opt_tests = get_option('tests')
opt_cuda = get_option('gpu-cuda')

################
# dependencies #
################

# Threads.
pancake_thread_dep = dependency('threads', required : true)

# boost
pancake_boost_dep = dependency('boost', required : true, include_type : 'system')

# zlib
pancacke_zlib_dep = dependency('zlib', required : true, fallback : ['zlib', 'zlib_dep'])

# htslib
pancake_htslib_dep = dependency('htslib', required : true, version : '>=1.4', fallback : ['htslib', 'htslib_dep'])

# pbcopper
pancake_pbcopper_dep = dependency('pbcopper', required : true, fallback : ['pbcopper', 'pbcopper_dep'])

# pbbam
pancake_pbbam_dep = dependency('pbbam', fallback : ['pbbam', 'pbbam_dep'])

# clock_gettime on old glibc systems
pancake_rt_dep = cpp.find_library('rt', required : false)

# Collect all deps.
pancake_lib_deps = [
  pancake_thread_dep,
  pancake_boost_dep,
  pancacke_zlib_dep,
  pancake_htslib_dep,
  pancake_pbcopper_dep,
  pancake_pbbam_dep,
  pancake_rt_dep,
]

if opt_cuda
    add_languages('cuda')
    nvcc = meson.get_compiler('cuda')
    add_project_arguments('-DGW_ENABLE_CACHING_ALLOCATOR', language : ['cpp', 'cuda'])
    pancake_gw_dep = dependency('genomeworks', fallback : ['genomeworks', 'genomeworks_dep'])
    pancake_lib_deps += [
        pancake_gw_dep,
    ]
endif

###########
# headers #
###########

subdir('include')

#####################
# sources + library #
#####################

subdir('src')

##########
# extras #
##########

if not meson.is_subproject()
  # tests
  if opt_tests
    pancake_clang_formatter = find_program('tools/check-formatting')

    pancake_cram_script = find_program('cram', required : false)
    if not pancake_cram_script.found()
      warning('Using bundled cram script')
      pancake_cram_script = find_program('scripts/cram', required : true)
    endif

    subdir('tests')
  endif
endif

###################
# dependency info #
###################
if not meson.is_subproject()
  # need to add pbcopper into 'Requires:' field,
  # but Meson currently only allows this if it's a
  # 'pkgconfig-dependency object' and not a subproject
  pancake_requires = []
  foreach dep : pancake_lib_deps
    if dep.type_name() == 'pkgconfig'
      pancake_requires += [dep]
    endif
  endforeach

  import('pkgconfig').generate(
    pancake_lib,
    version : meson.project_version(),
    name : 'pancake',
    requires : pancake_requires,
    filebase : 'pancake',
    description : 'A fast overlapper for HiFi reads.')
endif

pancake_dep = declare_dependency(
  include_directories : pancake_include_directories,
  link_with : pancake_lib,
  dependencies: pancake_lib_deps,
  version : meson.project_version())
