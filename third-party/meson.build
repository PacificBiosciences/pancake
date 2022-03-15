####################
# third-party libs #
####################

if not meson.is_subproject()
  # flat_hash_map
  install_headers(
    files([
      'flat_hash_map/bytell_hash_map.hpp',
      'flat_hash_map/flat_hash_map.hpp',
      'flat_hash_map/unordered_map.hpp',
      ]),
      subdir : 'pancake/third-party/flat_hash_map')

  # intervaltree
  install_headers(
    files([
      'intervaltree/IntervalTree.h',
      ]),
      subdir : 'pancake/third-party/intervaltree')

  # istl
  install_headers(
    files([
      'istl/fenwick.hpp',
      'istl/lcs.hpp',
      'istl/lis.hpp',
      'istl/range_tools.hpp',
      ]),
      subdir : 'pancake/third-party/istl')

  # kxsort
  install_headers(
    files([
      'kxsort/kxsort.h',
      ]),
      subdir : 'pancake/third-party/kxsort')

  # sparsehash/internal
  install_headers(
    files([
      'sparsehash/internal/densehashtable.h',
      'sparsehash/internal/hashtable-common.h',
      'sparsehash/internal/libc_allocator_with_realloc.h',
      'sparsehash/internal/sparseconfig.h',
      'sparsehash/internal/sparsehashtable.h',
      ]),
      subdir : 'pancake/third-party/sparsehash/internal/')

  # sparsehash/sparsehash/internal
  install_headers(
    files([
      'sparsehash/sparsehash/internal/densehashtable.h',
      'sparsehash/sparsehash/internal/hashtable-common.h',
      'sparsehash/sparsehash/internal/libc_allocator_with_realloc.h',
      'sparsehash/sparsehash/internal/sparseconfig.h',
      'sparsehash/sparsehash/internal/sparsehashtable.h',
      ]),
      subdir : 'pancake/third-party/sparsehash/sparsehash/internal')

  # sparsehash/sparsehash
  install_headers(
    files([
      'sparsehash/sparsehash/template_util.h',
      'sparsehash/sparsehash/type_traits.h',
      ]),
      subdir : 'pancake/third-party/sparsehash/sparsehash')

  # sparsehash
  install_headers(
    files([
      'sparsehash/template_util.h',
      'sparsehash/type_traits.h',
      ]),
      subdir : 'pancake/third-party/sparsehash')
endif

pancake_third_party_include_directories = include_directories(['.', 'sparsehash'])
