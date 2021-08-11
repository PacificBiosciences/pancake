// Authors: Ivan Sovic

#include "OverlapHifiWorkflow.h"
#include <pacbio/overlaphifi/OverlapHifiSettings.h>
#include <pacbio/pancake/MapperHiFi.h>
#include <pacbio/pancake/OverlapWriterFactory.h>
#include <pacbio/pancake/Seed.h>
#include <pacbio/pancake/SeedDBIndexCache.h>
#include <pacbio/pancake/SeedDBReaderCachedBlock.h>
#include <pacbio/pancake/SeedDBReaderRawBlock.h>
#include <pacbio/pancake/SeedIndex.h>
#include <pacbio/pancake/SeqDBIndexCache.h>
#include <pacbio/pancake/SeqDBReaderCached.h>
#include <pacbio/util/TicToc.h>
#include <pbcopper/logging/LogLevel.h>
#include <pbcopper/logging/Logging.h>
#include <pbcopper/parallel/FireAndForget.h>
#include <pbcopper/parallel/WorkQueue.h>
#include <sstream>

namespace PacBio {
namespace Pancake {

void Worker(const PacBio::Pancake::SeqDBReaderCachedBlock& targetSeqDBReader,
            const PacBio::Pancake::SeedIndex& index,
            const PacBio::Pancake::SeqDBReaderCachedBlock& querySeqDBReader,
            const PacBio::Pancake::SeedDBReaderCachedBlock& querySeedDBReader,
            const OverlapHifiSettings& /*settings*/, const OverlapHiFi::Mapper& mapper,
            int64_t freqCutoff, bool generateFlippedOverlaps, int32_t start, int32_t end,
            std::vector<OverlapHiFi::MapperResult>& results)
{
    int32_t numRecords = querySeqDBReader.records().size();
    if (start < 0 || end < 0 || start > end || start > numRecords || end > numRecords) {
        std::ostringstream oss;
        oss << "Invalid start/end indexes provided to the Worker. start = " << start
            << ", end = " << end << ", numRecords = " << numRecords;
        throw std::runtime_error(oss.str());
    }

    // Map the reads.
    for (int32_t i = start; i < end; ++i) {
        const auto& querySeq = querySeqDBReader.records()[i];
        const auto& querySeeds = querySeedDBReader.GetSeedsForSequence(querySeq.Id());
        results[i] = mapper.Map(targetSeqDBReader.recordStore(), index, querySeq, querySeeds,
                                freqCutoff, generateFlippedOverlaps);
    }
}

int OverlapHifiWorkflow::Runner(const PacBio::CLI_v2::Results& options)
{
    OverlapHifiSettings settings{options};

    std::string targetSeqDBFile = settings.TargetDBPrefix + ".seqdb";
    std::string targetSeedDBFile = settings.TargetDBPrefix + ".seeddb";
    std::string querySeqDBFile = settings.QueryDBPrefix + ".seqdb";
    std::string querySeedDBFile = settings.QueryDBPrefix + ".seeddb";

    PBLOG_INFO << "Loading the input DBs.";
    TicToc ttInit;

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
    if (targetSeedDBCache->seedParams.UseHPC != settings.UseHPC) {
        throw std::runtime_error(
            "The --use-hpc option was either used to compute the target SeedDB but not specified "
            "for overlapping, or vice versa.");
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
    if (querySeedDBCache->seedParams.UseHPC != settings.UseHPC) {
        throw std::runtime_error(
            "The --use-hpc option was either used to compute the query SeedDB but not specified "
            "for overlapping, or vice versa.");
    }

    // Create the target readers.
    PacBio::Pancake::SeqDBReaderCachedBlock targetSeqDBReader(targetSeqDBCache, settings.UseHPC);
    targetSeqDBReader.LoadBlocks({settings.TargetBlockId});
    PacBio::Pancake::SeedDBReaderRawBlock targetSeedDBReader(targetSeedDBCache);

    // Read the seeds for the target block.
    std::vector<PacBio::Pancake::SeedDB::SeedRaw> targetSeeds =
        targetSeedDBReader.GetBlock(settings.TargetBlockId);

    ttInit.Stop();
    PBLOG_INFO << "Loaded the target index and seqs in " << ttInit.GetSecs() << " sec.";

    PBLOG_INFO << "Target seqs: " << targetSeqDBReader.records().size();
    PBLOG_INFO << "Target seeds: " << targetSeeds.size();

    // Build the seed index.
    TicToc ttIndex;
    PacBio::Pancake::SeedIndex index(std::move(targetSeeds));
    ttIndex.Stop();
    PBLOG_INFO << "Built the seed index in " << ttIndex.GetSecs() << " sec.";

    // Seed statistics, and computing the cutoff.
    TicToc ttSeedStats;
    int64_t freqMax = 0;
    int64_t freqCutoff = 0;
    double freqAvg = 0.0;
    double freqMedian = 0.0;
    index.ComputeFrequencyStats(settings.FreqPercentile, freqMax, freqAvg, freqMedian, freqCutoff);
    ttSeedStats.Stop();
    PBLOG_INFO << "Computed the seed frequency statistics in " << ttSeedStats.GetSecs() << " sec.";

    PBLOG_INFO << "Seed statistic: freqMax = " << freqMax << ", freqAvg = " << freqAvg
               << ", freqMedian = " << freqMedian << ", freqCutoff = " << freqCutoff;

    PBLOG_INFO << "Using traceback alignment: " << (settings.UseTraceback ? "on" : "off");

    PBLOG_INFO << "Beginning to map the sequences.";
    PBLOG_INFO << "Using " << settings.NumThreads << " threads.";
    std::vector<OverlapHiFi::Mapper> mappers;
    for (size_t i = 0; i < settings.NumThreads; ++i) {
        mappers.emplace_back(OverlapHiFi::Mapper(settings));
    }
    TicToc ttMap;

    auto writer = PacBio::Pancake::OverlapWriterFactory(settings.OutFormat, stdout,
                                                        settings.WriteIds, settings.WriteCigar);

    writer->WriteHeader(targetSeqDBReader);

    const int32_t endBlockId = (settings.QueryBlockEndId <= 0) ? querySeqDBCache->blockLines.size()
                                                               : settings.QueryBlockEndId;

    // Process all blocks.
    PacBio::Pancake::SeqDBReaderCachedBlock querySeqDBReader(querySeqDBCache, settings.UseHPC);
    PacBio::Pancake::SeedDBReaderCachedBlock querySeedDBReader(querySeedDBCache);
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
        // PacBio::Pancake::SeqDBReaderCached querySeqDBReader(querySeqDBCache, queryBlockId);
        querySeqDBReader.LoadBlocks(blocksToLoad);
        PBLOG_INFO << "Loaded the query SeqDB cache block after " << ttQueryLoad.GetSecs(true)
                   << " sec.";
        querySeedDBReader.LoadBlock(blocksToLoad);
        ttQueryLoad.Stop();
        PBLOG_INFO << "Loaded the query SeedDB cache block after " << ttQueryLoad.GetSecs()
                   << " sec.";
        PBLOG_INFO << "Loaded all query blocks in " << ttQueryLoad.GetSecs() << " sec.";

        std::ostringstream oss;
        oss << "About to map query blocks: " << blocksToLoadStr
            << ": num_seqs = " << querySeqDBReader.records().size();
        PBLOG_INFO << oss.str();

        // Parallel processing.
        {
            TicToc ttQueryBlockMapping;

            // Determine how many records should land in each thread, spread roughly evenly.
            const int32_t numRecords = static_cast<int32_t>(querySeqDBReader.records().size());
            const std::vector<std::pair<int32_t, int32_t>> jobsPerThread =
                PacBio::Pancake::DistributeJobLoad<int32_t>(settings.NumThreads, numRecords);

            // Storage for the results of the batch run.
            std::vector<OverlapHiFi::MapperResult> results(numRecords);

            // Run the mapping in parallel.
            PacBio::Parallel::FireAndForget faf(settings.NumThreads);
            for (size_t i = 0; i < jobsPerThread.size(); ++i) {
                const int32_t jobStart = jobsPerThread[i].first;
                const int32_t jobEnd = jobsPerThread[i].second;
                faf.ProduceWith(Worker, std::cref(targetSeqDBReader), std::cref(index),
                                std::cref(querySeqDBReader), std::cref(querySeedDBReader),
                                std::cref(settings), std::cref(mappers[i]), freqCutoff,
                                settings.WriteReverseOverlaps, jobStart, jobEnd, std::ref(results));
            }
            faf.Finalize();

            // Write the results.
            for (size_t i = 0; i < querySeqDBReader.records().size(); ++i) {
                const auto& result = results[i];
                const auto& querySeq = querySeqDBReader.records()[i];

                for (const auto& ovl : result.overlaps) {
                    writer->Write(*ovl, targetSeqDBReader, querySeq);
                }
            }

            ttQueryBlockMapping.Stop();

            PBLOG_INFO << "Mapped query block in " << ttQueryBlockMapping.GetSecs() << " sec.";
        }
    }
    ttMap.Stop();
    PBLOG_INFO << "Mapped all query blocks in " << ttMap.GetSecs() << " sec.";

    return EXIT_SUCCESS;
}

}  // namespace Pancake
}  // namespace PacBio
