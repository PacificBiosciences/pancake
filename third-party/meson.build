####################
# third-party libs #
####################

if not meson.is_subproject()
  # flat_hash_map
  install_headers(
    files([
      'pancake/third-party/flat_hash_map/bytell_hash_map.hpp',
      'pancake/third-party/flat_hash_map/flat_hash_map.hpp',
      'pancake/third-party/flat_hash_map/unordered_map.hpp',
      ]),
      subdir : 'pancake/third-party/flat_hash_map')

  # intervaltree
  install_headers(
    files([
      'pancake/third-party/intervaltree/IntervalTree.h',
      ]),
      subdir : 'pancake/third-party/intervaltree')

  # istl
  install_headers(
    files([
      'pancake/third-party/istl/fenwick.hpp',
      'pancake/third-party/istl/lcs.hpp',
      'pancake/third-party/istl/lis.hpp',
      'pancake/third-party/istl/range_tools.hpp',
      ]),
      subdir : 'pancake/third-party/istl')

  # kxsort
  install_headers(
    files([
      'pancake/third-party/kxsort/kxsort.h',
      ]),
      subdir : 'pancake/third-party/kxsort')

  # pdqsort
  install_headers(
    files([
      'pancake/third-party/pdqsort/pdqsort.h',
      ]),
      subdir : 'pancake/third-party/pdqsort')

  # sparsehash/internal
  install_headers(
    files([
      'pancake/third-party/sparsehash/internal/densehashtable.h',
      'pancake/third-party/sparsehash/internal/hashtable-common.h',
      'pancake/third-party/sparsehash/internal/libc_allocator_with_realloc.h',
      'pancake/third-party/sparsehash/internal/sparseconfig.h',
      'pancake/third-party/sparsehash/internal/sparsehashtable.h',
      ]),
      subdir : 'pancake/third-party/sparsehash/internal/')

  # sparsehash/sparsehash/internal
  install_headers(
    files([
      'pancake/third-party/sparsehash/sparsehash/internal/densehashtable.h',
      'pancake/third-party/sparsehash/sparsehash/internal/hashtable-common.h',
      'pancake/third-party/sparsehash/sparsehash/internal/libc_allocator_with_realloc.h',
      'pancake/third-party/sparsehash/sparsehash/internal/sparseconfig.h',
      'pancake/third-party/sparsehash/sparsehash/internal/sparsehashtable.h',
      ]),
      subdir : 'pancake/third-party/sparsehash/sparsehash/internal')

  # sparsehash/sparsehash
  install_headers(
    files([
      'pancake/third-party/sparsehash/sparsehash/template_util.h',
      'pancake/third-party/sparsehash/sparsehash/type_traits.h',
      ]),
      subdir : 'pancake/third-party/sparsehash/sparsehash')

  # sparsehash
  install_headers(
    files([
      'pancake/third-party/sparsehash/template_util.h',
      'pancake/third-party/sparsehash/type_traits.h',
      ]),
      subdir : 'pancake/third-party/sparsehash')
endif

pancake_third_party_include_directories = include_directories(['.', 'pancake/third-party/sparsehash'])
