###########
# headers #
###########

if not meson.is_subproject()

  #'pacbio/Version.h',

  install_headers(
    files([
      'lib/math.hpp',
      ]),
      subdir : 'lib/')

  # lib/flat_hash_map
  install_headers(
    files([
      'lib/flat_hash_map/bytell_hash_map.hpp',
      'lib/flat_hash_map/flat_hash_map.hpp',
      'lib/flat_hash_map/unordered_map.hpp',
      ]),
      subdir : 'lib/flat_hash_map')

  # lib/intervaltree
  install_headers(
    files([
      'lib/intervaltree/IntervalTree.h',
      ]),
      subdir : 'lib/intervaltree')

  # lib/istl
  install_headers(
    files([
      'lib/istl/lis.hpp',
      ]),
      subdir : 'lib/istl')

  # lib/kxsort
  install_headers(
    files([
      'lib/kxsort/kxsort.h',
      ]),
      subdir : 'lib/kxsort')

  # lib/sparsehash/internal
  install_headers(
    files([
      'lib/sparsehash/internal/densehashtable.h',
      'lib/sparsehash/internal/hashtable-common.h',
      'lib/sparsehash/internal/libc_allocator_with_realloc.h',
      'lib/sparsehash/internal/sparseconfig.h',
      'lib/sparsehash/internal/sparsehashtable.h',
      ]),
      subdir : 'lib/sparsehash/internal/')

  # lib/sparsehash/sparsehash/internal
  install_headers(
    files([
      'lib/sparsehash/sparsehash/internal/densehashtable.h',
      'lib/sparsehash/sparsehash/internal/hashtable-common.h',
      'lib/sparsehash/sparsehash/internal/libc_allocator_with_realloc.h',
      'lib/sparsehash/sparsehash/internal/sparseconfig.h',
      'lib/sparsehash/sparsehash/internal/sparsehashtable.h',
      ]),
      subdir : 'lib/sparsehash/sparsehash/internal')

  # lib/sparsehash/sparsehash
  install_headers(
    files([
      'lib/sparsehash/sparsehash/template_util.h',
      'lib/sparsehash/sparsehash/type_traits.h',
      ]),
      subdir : 'lib/sparsehash/sparsehash')

  # lib/sparsehash
  install_headers(
    files([
      'lib/sparsehash/template_util.h',
      'lib/sparsehash/type_traits.h',
      ]),
      subdir : 'lib/sparsehash')

  # pacbio/alignment
  install_headers(
    files([
      'pacbio/alignment/AlignmentTools.h',
      'pacbio/alignment/DiffCounts.h',
      'pacbio/alignment/SesAlignBanded.hpp',
      'pacbio/alignment/SesDistanceBanded.h',
      'pacbio/alignment/Ses2AlignBanded.hpp',
      'pacbio/alignment/Ses2DistanceBanded.hpp',
      'pacbio/alignment/SesOptions.h',
      'pacbio/alignment/SesResults.h',
      ]),
      subdir : 'pacbio/alignment')

  # pacbio/overlaphifi
  install_headers(
    files([
      'pacbio/overlaphifi/OverlapHifiSettings.h',
      ]),
      subdir : 'pacbio/overlaphifi')

  # pacbio/pancake
  install_headers(
    files([
      'pacbio/pancake/AlignerBase.h',
      'pacbio/pancake/AlignerBatchBase.h',
      'pacbio/pancake/AlignerBatchCPU.h',
      'pacbio/pancake/AlignerEdlib.h',
      'pacbio/pancake/AlignerFactory.h',
      'pacbio/pancake/AlignerKSW2.h',
      'pacbio/pancake/AlignerSES1.h',
      'pacbio/pancake/AlignerSES2.h',
      'pacbio/pancake/AlignmentParameters.h',
      'pacbio/pancake/AlignmentResult.h',
      'pacbio/pancake/AlignmentSeeded.h',
      'pacbio/pancake/CompressedSequence.h',
      'pacbio/pancake/ContiguousFilePart.h',
      'pacbio/pancake/DPChain.h',
      'pacbio/pancake/FastaSequenceCached.h',
      'pacbio/pancake/FastaSequenceCachedStore.h',
      'pacbio/pancake/FastaSequenceId.h',
      'pacbio/pancake/Lookups.h',
      'pacbio/pancake/MapperBase.h',
      'pacbio/pancake/MapperBatchBase.h',
      'pacbio/pancake/MapperBatchCPU.h',
      'pacbio/pancake/MapperBatchUtility.h',
      'pacbio/pancake/MapperCLR.h',
      'pacbio/pancake/MapperHiFi.h',
      'pacbio/pancake/MapperUtility.h',
      'pacbio/pancake/Minimizers.h',
      'pacbio/pancake/Overlap.h',
      'pacbio/pancake/OverlapWriterBase.h',
      'pacbio/pancake/OverlapWriterFactory.h',
      'pacbio/pancake/OverlapWriterFormat.h',
      'pacbio/pancake/OverlapWriterIPAOvl.h',
      'pacbio/pancake/OverlapWriterM4.h',
      'pacbio/pancake/OverlapWriterPAF.h',
      'pacbio/pancake/OverlapWriterSAM.h',
      'pacbio/pancake/Range.h',
      'pacbio/pancake/Secondary.h',
      'pacbio/pancake/Seed.h',
      'pacbio/pancake/SeedDBIndexCache.h',
      'pacbio/pancake/SeedDBParameters.h',
      'pacbio/pancake/SeedDBReader.h',
      'pacbio/pancake/SeedDBReaderCached.h',
      'pacbio/pancake/SeedDBReaderCachedBlock.h',
      'pacbio/pancake/SeedDBReaderRawBlock.h',
      'pacbio/pancake/SeedDBWriter.h',
      'pacbio/pancake/SeedHit.h',
      'pacbio/pancake/SeedHitWriter.h',
      'pacbio/pancake/SeedIndex.h',
      'pacbio/pancake/SeqDBIndexCache.h',
      'pacbio/pancake/SeqDBReader.h',
      'pacbio/pancake/SeqDBReaderBase.h',
      'pacbio/pancake/SeqDBReaderCached.h',
      'pacbio/pancake/SeqDBReaderCachedBlock.h',
      'pacbio/pancake/SeqDBWriter.h',
      'pacbio/pancake/SeqDBWriterBase.h',
      'pacbio/pancake/SequenceSeeds.h',
      'pacbio/pancake/SequenceSeedsCached.h',
      'pacbio/pancake/Twobit.h',
      ]),
      subdir : 'pacbio/pancake')

  # pacbio/alignment
  install_headers(
    files([
      'pacbio/util/CommonTypes.h',
      'pacbio/util/Conversion.h',
      'pacbio/util/FileIO.h',
      'pacbio/util/GenomicUnit.h',
      'pacbio/util/RunLengthEncoding.h',
      'pacbio/util/SeqLengthStats.h',
      'pacbio/util/TicToc.h',
      'pacbio/util/Util.h',
      ]),
      subdir : 'pacbio/util')

  if opt_cuda
      install_headers(
        files([
          'pacbio/pancake/AlignerBatchGPU.h',
          'pacbio/pancake/MapperBatchGPU.h',
          ]),
          subdir : 'pacbio/pancake')
  endif
endif

pancake_include_directories = include_directories(['.', 'pacbio', 'lib/sparsehash'])
