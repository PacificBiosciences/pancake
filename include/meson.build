###########
# headers #
###########

if not meson.is_subproject()
  # pancake
  install_headers(
    files([
      'pancake/AlignmentTools.hpp',
      'pancake/AlignerBase.hpp',
      'pancake/AlignerBatchBase.hpp',
      'pancake/AlignerBatchCPU.hpp',
      'pancake/AlignerEdlib.hpp',
      'pancake/AlignerFactory.hpp',
      'pancake/AlignerKSW2.hpp',
      'pancake/AlignerSES1.hpp',
      'pancake/AlignerSES2.hpp',
      'pancake/AlignmentParameters.hpp',
      'pancake/AlignmentResult.hpp',
      'pancake/AlignmentSeeded.hpp',
      'pancake/CompressedSequence.hpp',
      'pancake/ContiguousFilePart.hpp',
      'pancake/DiffCounts.hpp',
      'pancake/DPChain.hpp',
      'pancake/FastaSequenceCached.hpp',
      'pancake/FastaSequenceCachedStore.hpp',
      'pancake/FastaSequenceId.hpp',
      'pancake/Lookups.hpp',
      'pancake/MapperBase.hpp',
      'pancake/MapperBatchBase.hpp',
      'pancake/MapperBatchCPU.hpp',
      'pancake/MapperBatchUtility.hpp',
      'pancake/MapperCLR.hpp',
      'pancake/MapperHiFi.hpp',
      'pancake/MapperUtility.hpp',
      'pancake/Minimizers.hpp',
      'pancake/Overlap.hpp',
      'pancake/OverlapHifiSettings.hpp',
      'pancake/OverlapWriterBase.hpp',
      'pancake/OverlapWriterFactory.hpp',
      'pancake/OverlapWriterFormat.hpp',
      'pancake/OverlapWriterIPAOvl.hpp',
      'pancake/OverlapWriterM4.hpp',
      'pancake/OverlapWriterPAF.hpp',
      'pancake/OverlapWriterSAM.hpp',
      'pancake/Range.hpp',
      'pancake/Secondary.hpp',
      'pancake/Seed.hpp',
      'pancake/SeedDBIndexCache.hpp',
      'pancake/SeedDBParameters.hpp',
      'pancake/SeedDBReader.hpp',
      'pancake/SeedDBReaderCached.hpp',
      'pancake/SeedDBReaderCachedBlock.hpp',
      'pancake/SeedDBReaderRawBlock.hpp',
      'pancake/SeedDBWriter.hpp',
      'pancake/SeedHit.hpp',
      'pancake/SeedHitWriter.hpp',
      'pancake/SeedIndex.hpp',
      'pancake/SeqDBIndexCache.hpp',
      'pancake/SeqDBReader.hpp',
      'pancake/SeqDBReaderBase.hpp',
      'pancake/SeqDBReaderCached.hpp',
      'pancake/SeqDBReaderCachedBlock.hpp',
      'pancake/SeqDBWriter.hpp',
      'pancake/SeqDBWriterBase.hpp',
      'pancake/SequenceSeeds.hpp',
      'pancake/SequenceSeedsCached.hpp',
      'pancake/SesAlignBanded.hpp',
      'pancake/SesDistanceBanded.hpp',
      'pancake/Ses2AlignBanded.hpp',
      'pancake/Ses2DistanceBanded.hpp',
      'pancake/SesOptions.hpp',
      'pancake/SesResults.hpp',
      'pancake/Twobit.hpp',
      ]),
      subdir : 'pancake')

  # alignment
  install_headers(
    files([
      'pancake/util/CommonTypes.hpp',
      'pancake/util/DebugTools.hpp',
      'pancake/util/Conversion.hpp',
      'pancake/util/FileIO.hpp',
      'pancake/util/GenomicUnit.hpp',
      'pancake/util/Math.hpp',
      'pancake/util/RunLengthEncoding.hpp',
      'pancake/util/SeqLengthStats.hpp',
      'pancake/util/TicToc.hpp',
      'pancake/util/Util.hpp',
      ]),
      subdir : 'pancake/util')

  if opt_cuda
      install_headers(
        files([
          'pancake/AlignerBatchGPU.hpp',
          'pancake/AlignerBatchGPUEdelweiss.hpp',
          'pancake/MapperBatchGPU.hpp',
          'pancake/MapperBatchGPUEdelweiss.hpp',
          ]),
          subdir : 'pancake')
  endif
  if opt_sse41
      install_headers(
        files([
            'pancake/DPChainSimd.hpp',
          ]),
          subdir : 'pancake')
  endif
endif

pancake_include_directories = include_directories(['.'])
