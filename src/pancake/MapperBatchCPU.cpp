// Authors: Ivan Sovic

#include <pacbio/alignment/AlignmentTools.h>
#include <pacbio/pancake/AlignerBase.h>
#include <pacbio/pancake/MapperBatchCPU.h>
#include <pacbio/pancake/MapperBatchUtility.h>
#include <pacbio/pancake/OverlapWriterFactory.h>
#include <pacbio/util/TicToc.h>
#include <pbcopper/logging/Logging.h>
#include <pbcopper/parallel/FireAndForget.h>
#include <pbcopper/parallel/WorkQueue.h>
#include <pbcopper/utility/Stopwatch.h>
#include <algorithm>
#include <array>
#include <iostream>
#include <sstream>
#include <tuple>

namespace PacBio {
namespace Pancake {

MapperBatchCPU::MapperBatchCPU(const MapperCLRAlignSettings& alignSettings, int32_t numThreads)
    : MapperBatchCPU{alignSettings, nullptr}
{
    fafFallback_ = std::make_unique<Parallel::FireAndForget>(numThreads);
    faf_ = fafFallback_.get();
}

MapperBatchCPU::MapperBatchCPU(const MapperCLRAlignSettings& alignSettings,
                               Parallel::FireAndForget* faf)
    : alignSettings_{alignSettings}, faf_{faf}, fafFallback_{nullptr}
{
}

MapperBatchCPU::~MapperBatchCPU()
{
    if (fafFallback_) {
        fafFallback_->Finalize();
    }
}

std::vector<std::vector<MapperBaseResult>> MapperBatchCPU::MapAndAlign(
    const std::vector<MapperBatchChunk>& batchData)
{
    return MapAndAlignImpl_(batchData, alignSettings_, faf_);
}

std::vector<std::vector<MapperBaseResult>> MapperBatchCPU::MapAndAlignImpl_(
    const std::vector<MapperBatchChunk>& batchChunks, const MapperCLRAlignSettings& alignSettings,
    Parallel::FireAndForget* faf)
{
    PacBio::Utility::Stopwatch timer;
    // Determine how many records should land in each thread, spread roughly evenly.
    const int32_t numRecords = batchChunks.size();
    const int32_t numThreads = faf ? faf->NumThreads() : 1;
    const std::vector<std::pair<int32_t, int32_t>> jobsPerThread =
        PacBio::Pancake::DistributeJobLoad<int32_t>(numThreads, numRecords);

    // Run the mapping in parallel.
    std::vector<std::vector<MapperBaseResult>> results(batchChunks.size());
    const auto Submit = [&jobsPerThread, &batchChunks, &results](int32_t i) {
        const int32_t jobStart = jobsPerThread[i].first;
        const int32_t jobEnd = jobsPerThread[i].second;
        WorkerMapper_(batchChunks, jobStart, jobEnd, results);
    };

    Parallel::Dispatch(faf, jobsPerThread.size(), Submit);
    PBLOG_INFO << "CPU Mapping            : " << timer.ElapsedTime();
    timer.Reset();

    int64_t cpuTime = 0;
    if (alignSettings.align) {
        // Compute the reverse complements for alignment.
        std::vector<std::vector<std::string>> querySeqsRev =
            ComputeReverseComplements(batchChunks, results, faf);
        timer.Freeze();
        cpuTime += timer.ElapsedNanoseconds();
        PBLOG_INFO << "CPU RevComp            : " << timer.ElapsedTime();
        timer.Reset();

        // Prepare the sequences for alignment.
        std::vector<PairForBatchAlignment> partsGlobal;
        std::vector<PairForBatchAlignment> partsSemiglobal;
        std::vector<AlignmentStitchInfo> alnStitchInfo;
        int32_t longestSequenceForAln = 0;
        PrepareSequencesForBatchAlignment(batchChunks, querySeqsRev, results,
                                          alignSettings.selfHitPolicy, partsGlobal, partsSemiglobal,
                                          alnStitchInfo, longestSequenceForAln);
        PBLOG_TRACE << "partsGlobal.size() = " << partsGlobal.size();
        PBLOG_TRACE << "partsSemiglobal.size() = " << partsSemiglobal.size();
        timer.Freeze();
        cpuTime += timer.ElapsedNanoseconds();
        PBLOG_INFO << "CPU Prepare            : " << timer.ElapsedTime();
        timer.Reset();

        // Internal alignment on CPU.
        std::vector<AlignmentResult> internalAlns;
        int64_t prepareTime = 0;
        int64_t alignTime = 0;
        AlignPartsOnCpu(alignSettings.alignerTypeGlobal, alignSettings.alnParamsGlobal,
                        alignSettings.alignerTypeExt, alignSettings.alnParamsExt, partsGlobal, faf,
                        internalAlns, prepareTime, alignTime);
        PBLOG_TRACE << "internalAlns.size() = " << internalAlns.size();
        timer.Freeze();
        cpuTime += timer.ElapsedNanoseconds();
        PBLOG_INFO << "CPU Internal Prepare   : "
                   << PacBio::Utility::Stopwatch::PrettyPrintNanoseconds(prepareTime);
        PBLOG_INFO << "CPU Internal Alignment : "
                   << PacBio::Utility::Stopwatch::PrettyPrintNanoseconds(alignTime);
        timer.Reset();

        // Flank alignment on CPU.
        std::vector<AlignmentResult> flankAlns;
        prepareTime = 0;
        alignTime = 0;
        AlignPartsOnCpu(alignSettings.alignerTypeGlobal, alignSettings.alnParamsGlobal,
                        alignSettings.alignerTypeExt, alignSettings.alnParamsExt, partsSemiglobal,
                        faf, flankAlns, prepareTime, alignTime);
        PBLOG_TRACE << "flankAlns.size() = " << flankAlns.size();
        timer.Freeze();
        cpuTime += timer.ElapsedNanoseconds();
        PBLOG_INFO << "CPU Flanks Prepare     : "
                   << PacBio::Utility::Stopwatch::PrettyPrintNanoseconds(prepareTime);
        PBLOG_INFO << "CPU Flanks Alignment   : "
                   << PacBio::Utility::Stopwatch::PrettyPrintNanoseconds(alignTime);
        timer.Reset();

        StitchAlignmentsInParallel(results, batchChunks, querySeqsRev, internalAlns, flankAlns,
                                   alnStitchInfo, faf);
        timer.Freeze();
        cpuTime += timer.ElapsedNanoseconds();
        PBLOG_INFO << "CPU Stitch             : " << timer.ElapsedTime();
        timer.Reset();

        SetUnalignedAndMockedMappings(
            results, alignSettings.selfHitPolicy == MapperSelfHitPolicy::PERFECT_ALIGNMENT,
            alignSettings.alnParamsGlobal.matchScore);
        timer.Freeze();
        cpuTime += timer.ElapsedNanoseconds();
        PBLOG_INFO << "CPU Mock               : " << timer.ElapsedTime();
        timer.Reset();

        UpdateSecondaryAndFilter(results, faf, batchChunks);
        timer.Freeze();
        cpuTime += timer.ElapsedNanoseconds();
        PBLOG_INFO << "CPU Update             : " << timer.ElapsedTime();
        timer.Reset();
    }
    PBLOG_INFO << "CPU Time               : "
               << PacBio::Utility::Stopwatch::PrettyPrintNanoseconds(cpuTime);

    return results;
}

void MapperBatchCPU::WorkerMapper_(const std::vector<MapperBatchChunk>& batchChunks,
                                   int32_t startId, int32_t endId,
                                   std::vector<std::vector<MapperBaseResult>>& results)
{
    for (int32_t i = startId; i < endId; ++i) {
        const auto& chunk = batchChunks[i];

        // Create a copy of the settings so that we can turn off the alignment.
        MapperCLRSettings settingsCopy;
        settingsCopy.map = chunk.mapSettings;
        settingsCopy.align.align = false;

        // Create the mapper.
        MapperCLR mapper(settingsCopy);

        results[i] = mapper.MapAndAlign(chunk.targetSeqs, chunk.querySeqs);
    }
}

void UpdateSecondaryAndFilter(std::vector<std::vector<MapperBaseResult>>& mappingResults,
                              Parallel::FireAndForget* faf,
                              const std::vector<MapperBatchChunk>& batchChunks)
{

    const int32_t numThreads = faf ? faf->NumThreads() : 1;
    const int32_t numEntries = mappingResults.size();
    const std::vector<std::pair<int32_t, int32_t>> jobsPerThread =
        PacBio::Pancake::DistributeJobLoad<int32_t>(numThreads, numEntries);

    // Results are a vector for every chunk (one chunk is one ZMW).
    const auto Submit = [&](int32_t jobId) {
        const int32_t jobStart = jobsPerThread[jobId].first;
        const int32_t jobEnd = jobsPerThread[jobId].second;

        for (int32_t resultId = jobStart; resultId < jobEnd; ++resultId) {
            const auto& settings = batchChunks[resultId].mapSettings;
            auto& result = mappingResults[resultId];
            // One chunk can have multiple queries (subreads).
            for (size_t qId = 0; qId < result.size(); ++qId) {
                // Secondary/supplementary flagging.
                WrapFlagSecondaryAndSupplementary(result[qId].mappings,
                                                  settings.secondaryAllowedOverlapFractionQuery,
                                                  settings.secondaryAllowedOverlapFractionTarget,
                                                  settings.secondaryMinScoreFraction);

                const int32_t numPrimary =
                    CondenseMappings(result[qId].mappings, settings.bestNSecondary);

                // If this occurs, that means that a filtering stage removed the primary alignment for some reason.
                // This can happen in case self-hits are skipped in the overlapping use case (the self-hit is the
                // highest scoring one, and it would get marked for removal).
                // This reruns labeling to produce the next best primary.
                if (numPrimary == 0) {
                    WrapFlagSecondaryAndSupplementary(
                        result[qId].mappings, settings.secondaryAllowedOverlapFractionQuery,
                        settings.secondaryAllowedOverlapFractionTarget,
                        settings.secondaryMinScoreFraction);
                }
            }
        }
    };
    Parallel::Dispatch(faf, jobsPerThread.size(), Submit);
}

int32_t AlignPartsOnCpu(const AlignerType& alignerTypeGlobal,
                        const AlignmentParameters& alnParamsGlobal,
                        const AlignerType& alignerTypeExt, const AlignmentParameters& alnParamsExt,
                        const std::vector<PairForBatchAlignment>& parts, const int32_t numThreads,
                        std::vector<AlignmentResult>& retAlns, int64_t& prepareTime,
                        int64_t& alignTime)
{
    Parallel::FireAndForget faf(numThreads);
    const int32_t result =
        AlignPartsOnCpu(alignerTypeGlobal, alnParamsGlobal, alignerTypeExt, alnParamsExt, parts,
                        &faf, retAlns, prepareTime, alignTime);
    faf.Finalize();
    return result;
}
int32_t AlignPartsOnCpu(const AlignerType& alignerTypeGlobal,
                        const AlignmentParameters& alnParamsGlobal,
                        const AlignerType& alignerTypeExt, const AlignmentParameters& alnParamsExt,
                        const std::vector<PairForBatchAlignment>& parts,
                        Parallel::FireAndForget* faf, std::vector<AlignmentResult>& retAlns,
                        int64_t& prepareTime, int64_t& alignTime)
{
    PacBio::Utility::Stopwatch timer;
    retAlns.resize(parts.size());

    std::vector<size_t> partIds;

    AlignerBatchCPU aligner(faf, alignerTypeGlobal, alnParamsGlobal, alignerTypeExt, alnParamsExt);

    for (size_t i = 0; i < parts.size(); ++i) {
        const auto& part = parts[i];
        if (retAlns[i].valid) {
            continue;
        }
        partIds.emplace_back(i);
        if (part.regionType == RegionType::FRONT) {
            // Reverse the sequences for front flank alignment. No need to complement.
            std::string query(part.query, part.queryLen);
            std::reverse(query.begin(), query.end());
            std::string target(part.target, part.targetLen);
            std::reverse(target.begin(), target.end());
            aligner.AddSequencePair(query.c_str(), part.queryLen, target.c_str(), part.targetLen,
                                    part.regionType == RegionType::GLOBAL);
        } else {
            aligner.AddSequencePair(part.query, part.queryLen, part.target, part.targetLen,
                                    part.regionType == RegionType::GLOBAL);
        }
    }
    prepareTime += timer.ElapsedNanoseconds();
    timer.Reset();

    aligner.AlignAll();
    alignTime += timer.ElapsedNanoseconds();
    timer.Reset();

    const std::vector<AlignmentResult>& partInternalAlns = aligner.GetAlnResults();
    int32_t numNotValid = 0;
    for (size_t i = 0; i < partInternalAlns.size(); ++i) {
        const auto& aln = partInternalAlns[i];
        if (aln.valid == false) {
            ++numNotValid;
        }
        retAlns[partIds[i]] = std::move(partInternalAlns[i]);
    }
    prepareTime += timer.ElapsedNanoseconds();
    return numNotValid;
}

}  // namespace Pancake
}  // namespace PacBio
