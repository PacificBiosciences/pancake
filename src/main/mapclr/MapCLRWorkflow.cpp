// Authors: Ivan Sovic

#include "MapCLRWorkflow.hpp"
#include <pbcopper/logging/LogLevel.h>
#include <pbcopper/logging/Logging.h>
#include <pbcopper/parallel/FireAndForget.h>
#include <pbcopper/parallel/WorkQueue.h>
#include <pancake/MapperCLR.hpp>
#include <pancake/OverlapWriterFactory.hpp>
#include <pancake/Seed.hpp>
#include <pancake/SeedDBIndexCache.hpp>
#include <pancake/SeedDBReaderCachedBlock.hpp>
#include <pancake/SeedDBReaderRawBlock.hpp>
#include <pancake/SeedIndex.hpp>
#include <pancake/SeqDBIndexCache.hpp>
#include <pancake/SeqDBReaderCached.hpp>
#include <pancake/util/TicToc.hpp>
#include <sstream>
#include "MapCLRSettings.hpp"

namespace PacBio {
namespace Pancake {

namespace MapCLRCLI {

void Worker(const PacBio::Pancake::SeqDBReaderCachedBlock& targetSeqDBReader,
            const PacBio::Pancake::SeedIndex& index,
            const PacBio::Pancake::SeqDBReaderCachedBlock& querySeqDBReader,
            const PacBio::Pancake::SeedDBReaderCachedBlock& querySeedDBReader, MapperCLR& mapper,
            int64_t freqCutoff, int32_t start, int32_t end, std::vector<MapperBaseResult>& results)
{
    const int32_t numRecords = querySeqDBReader.records().size();
    if (start < 0 || end < 0 || start > end || start > numRecords || end > numRecords) {
        std::ostringstream oss;
        oss << "Invalid start/end indexes provided to the Worker. start = " << start
            << ", end = " << end << ", numRecords = " << numRecords;
        throw std::runtime_error(oss.str());
    }

    for (int32_t i = start; i < end; ++i) {
        const auto& querySeq = querySeqDBReader.records()[i];
        const auto& querySeeds = querySeedDBReader.GetSeedsForSequence(querySeq.Id());
        results[i] = mapper.MapAndAlignSingleQuery(targetSeqDBReader.records(), index, querySeq,
                                                   querySeeds, i, freqCutoff);
    }
}
}  // namespace MapCLRCLI

int MapCLRWorkflow::Runner(const PacBio::CLI_v2::Results& options)
{
    const MapCLRSettings settings{options};

    const std::string targetSeqDBFile = settings.TargetDBPrefix + ".seqdb";
    const std::string targetSeedDBFile = settings.TargetDBPrefix + ".seeddb";
    const std::string querySeqDBFile = settings.QueryDBPrefix + ".seqdb";
    const std::string querySeedDBFile = settings.QueryDBPrefix + ".seeddb";

    TicToc ttInit;
    PBLOG_INFO << "Loading the input DBs.";

    // Load the target DB caches.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> targetSeqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(targetSeqDBFile);
    PBLOG_INFO << "After loading target seq cache: " << ttInit.VerboseSecs(true);
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> targetSeedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(targetSeedDBFile);
    PBLOG_INFO << "After loading target seed cache: " << ttInit.VerboseSecs(true);
    PBLOG_INFO << "Target seed params: k = " << targetSeedDBCache->seedParams.KmerSize
               << ", w = " << targetSeedDBCache->seedParams.MinimizerWindow
               << ", s = " << targetSeedDBCache->seedParams.Spacing
               << ", hpc = " << targetSeedDBCache->seedParams.UseHPC
               << ", rc = " << targetSeedDBCache->seedParams.UseRC;
    if (targetSeedDBCache->seedParams.UseHPC) {
        throw std::runtime_error(
            "The target SeedDB was generated using the --use-hpc option which means that the "
            "sequence coordiante space is physically compressed. This is not supported by this "
            "tool.");
    }

    // Load the query DB caches.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> querySeqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(querySeqDBFile);
    PBLOG_INFO << "After loading query seq cache: " << ttInit.VerboseSecs(true);
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> querySeedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(querySeedDBFile);
    PBLOG_INFO << "After loading query seed cache: " << ttInit.VerboseSecs(true);
    PBLOG_INFO << "Query seed params: k = " << targetSeedDBCache->seedParams.KmerSize
               << ", w = " << targetSeedDBCache->seedParams.MinimizerWindow
               << ", s = " << targetSeedDBCache->seedParams.Spacing
               << ", hpc = " << targetSeedDBCache->seedParams.UseHPC
               << ", rc = " << targetSeedDBCache->seedParams.UseRC;
    if (querySeedDBCache->seedParams.UseHPC) {
        throw std::runtime_error(
            "The query SeedDB was generated using the --use-hpc option which means that the "
            "sequence coordiante space is physically compressed. This is not supported by this "
            "tool.");
    }

    // Create the target sequence and seed readers.
    PacBio::Pancake::SeqDBReaderCachedBlock targetSeqDBReader(targetSeqDBCache, false);
    targetSeqDBReader.LoadBlocks({settings.TargetBlockId});
    const PacBio::Pancake::SeedDBReaderRawBlock targetSeedDBReader(targetSeedDBCache);

    // Read the seeds for the target block. This will be moved into the SeedIndex below.
    std::vector<PacBio::Pancake::SeedRaw> targetSeeds =
        targetSeedDBReader.GetBlock(settings.TargetBlockId);

    ttInit.Stop();
    PBLOG_INFO << "Loaded the target index and seqs in " << ttInit.GetSecs() << " sec";
    PBLOG_INFO << "Target seqs: " << targetSeqDBReader.records().size();
    PBLOG_INFO << "Target seeds: " << targetSeeds.size();

    // Build the seed index.
    TicToc ttIndex;
    PacBio::Pancake::SeedIndex index(std::move(targetSeeds));
    ttIndex.Stop();
    PBLOG_INFO << "Built the seed index in " << ttIndex.GetSecs() << " sec";

    // Seed statistics, and computing the cutoff.
    TicToc ttSeedStats;
    int64_t freqMax = 0;
    int64_t freqCutoff = 0;
    double freqAvg = 0.0;
    double freqMedian = 0.0;
    index.ComputeFrequencyStats(settings.MapperSettings.map.freqPercentile, freqMax, freqAvg,
                                freqMedian, freqCutoff);
    ttSeedStats.Stop();
    PBLOG_INFO << "Computed the seed frequency statistics in " << ttSeedStats.GetSecs() << " sec";

    PBLOG_INFO << "Seed statistic: freqMax = " << freqMax << ", freqAvg = " << freqAvg
               << ", freqMedian = " << freqMedian << ", freqCutoff = " << freqCutoff;

    PBLOG_INFO << "Align: " << (settings.MapperSettings.align.align ? "true" : "false");

    PBLOG_INFO << "Beginning to map the sequences.";
    PBLOG_INFO << "Using " << settings.NumThreads << " threads.";

    // Create one mapper per thread.
    std::vector<MapperCLR> mappers;
    for (size_t i = 0; i < settings.NumThreads; ++i) {
        mappers.emplace_back(MapperCLR(settings.MapperSettings));
    }

    TicToc ttMap;

    // Create the overlap writer.
    auto writer = PacBio::Pancake::OverlapWriterFactory(settings.OutFormat, stdout,
                                                        settings.WriteIds, settings.WriteCigar);
    writer->WriteHeader(targetSeqDBReader);

    // Process all blocks.
    PacBio::Pancake::SeqDBReaderCachedBlock querySeqDBReader(querySeqDBCache, false);
    PacBio::Pancake::SeedDBReaderCachedBlock querySeedDBReader(querySeedDBCache);
    const int32_t endBlockId = (settings.QueryBlockEndId <= 0) ? querySeqDBCache->blockLines.size()
                                                               : settings.QueryBlockEndId;
    for (int32_t queryBlockId = settings.QueryBlockStartId; queryBlockId < endBlockId;
         queryBlockId += settings.CombineBlocks) {

        std::vector<int32_t> blocksToLoad;
        for (int32_t blockId = queryBlockId;
             blockId < std::min(endBlockId, (queryBlockId + settings.CombineBlocks)); ++blockId) {
            blocksToLoad.emplace_back(blockId);
        }
        if (blocksToLoad.empty()) {
            throw std::runtime_error("There are zero blocks to load!");
        }
        std::string blocksToLoadStr;
        {
            std::ostringstream oss;
            oss << "{" << blocksToLoad[0];
            for (size_t i = 1; i < blocksToLoad.size(); ++i) {
                oss << ", " << blocksToLoad[i];
            }
            oss << "}";
            blocksToLoadStr = oss.str();
        }

        PBLOG_INFO << "Loading the query blocks: " << blocksToLoadStr << ".";

        // Create the query readers for the current block.
        TicToc ttQueryLoad;
        querySeqDBReader.LoadBlocks(blocksToLoad);
        PBLOG_INFO << "Loaded the query SeqDB cache block after " << ttQueryLoad.GetSecs(true)
                   << " sec";

        querySeedDBReader.LoadBlock(blocksToLoad);
        ttQueryLoad.Stop();
        PBLOG_INFO << "Loaded the query SeedDB cache block after " << ttQueryLoad.GetSecs()
                   << " sec";

        PBLOG_INFO << "Loaded all query blocks in " << ttQueryLoad.GetSecs() << " sec";
        PBLOG_INFO << "About to map query blocks: " << blocksToLoadStr
                   << ": num_seqs = " << querySeqDBReader.records().size();

        // Map/align in parallel.
        {
            TicToc ttQueryBlockMapping;

            // Compute jobs per thread.
            const int32_t numRecords = querySeqDBReader.records().size();
            const std::vector<std::pair<int32_t, int32_t>> jobsPerThread =
                PacBio::Pancake::DistributeJobLoad<int32_t>(settings.NumThreads, numRecords);

            // Mapping.
            std::vector<MapperBaseResult> results(numRecords);
            PacBio::Parallel::FireAndForget faf(settings.NumThreads);
            for (size_t i = 0; i < jobsPerThread.size(); ++i) {
                const int32_t jobStart = jobsPerThread[i].first;
                const int32_t jobEnd = jobsPerThread[i].second;
                faf.ProduceWith(MapCLRCLI::Worker, std::cref(targetSeqDBReader), std::cref(index),
                                std::cref(querySeqDBReader), std::cref(querySeedDBReader),
                                mappers[i], freqCutoff, jobStart, jobEnd, std::ref(results));
            }
            faf.Finalize();

            // Write the results.
            for (size_t i = 0; i < querySeqDBReader.records().size(); ++i) {
                const auto& result = results[i];
                const auto& querySeq = querySeqDBReader.records()[i];
                for (const auto& chainedRegion : result.mappings) {
                    writer->Write(*chainedRegion->mapping, targetSeqDBReader, querySeq);
                }
            }

            ttQueryBlockMapping.Stop();

            PBLOG_INFO << "Mapped query block in " << ttQueryBlockMapping.GetSecs() << " sec";
        }
    }
    ttMap.Stop();
    PBLOG_INFO << "Mapped all query blocks in " << ttMap.GetSecs() << " sec";

    return EXIT_SUCCESS;
}

}  // namespace Pancake
}  // namespace PacBio
