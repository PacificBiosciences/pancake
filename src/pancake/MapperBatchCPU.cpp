// Authors: Ivan Sovic

#include <pancake/MapperBatchCPU.hpp>

#include <pancake/AlignerBase.hpp>
#include <pancake/AlignmentTools.hpp>
#include <pancake/MapperBatchUtility.hpp>
#include <pancake/MapperUtility.hpp>
#include <pancake/OverlapWriterFactory.hpp>
#include <pancake/util/TicToc.hpp>

#include <pbcopper/logging/Logging.h>
#include <pbcopper/parallel/FireAndForget.h>
#include <pbcopper/parallel/WorkQueue.h>

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
{}

MapperBatchCPU::~MapperBatchCPU()
{
    if (fafFallback_) {
        fafFallback_->Finalize();
    }
}

std::vector<std::vector<MapperBaseResult>> MapperBatchCPU::MapAndAlign(
    const std::vector<MapperBatchChunk>& batchData)
{
    try {
        return MapAndAlignImpl_(batchData, alignSettings_, faf_);
    } catch (const std::exception& e) {
        // Log, but do not fail. Important for clients of this class.
        // Return a vector of the size of the input, but with empty results for each query.
        PBLOG_DEBUG << "MapperBatchCPU generated an exception in " << std::string(__FUNCTION__)
                    << ". Message: " << e.what();
        std::vector<std::vector<MapperBaseResult>> results(batchData.size());
        for (size_t i = 0; i < batchData.size(); ++i) {
            results[i].resize(batchData[i].querySeqs.Size());
        }
        return results;
    }
}

std::vector<std::vector<MapperBaseResult>> MapperBatchCPU::MapAndAlignImpl_(
    const std::vector<MapperBatchChunk>& batchChunks, const MapperCLRAlignSettings& alignSettings,
    Parallel::FireAndForget* faf) const
{
    // Determine how many records should land in each thread, spread roughly evenly.
    const int32_t numRecords = batchChunks.size();
    const int32_t numThreads = faf ? faf->NumThreads() : 1;
    const std::vector<std::pair<int32_t, int32_t>> jobsPerThread =
        PacBio::Pancake::DistributeJobLoad<int32_t>(numThreads, numRecords);

    // Run the mapping in parallel.
    std::vector<std::vector<MapperBaseResult>> results(batchChunks.size());
    const auto Submit = [&jobsPerThread, &batchChunks, &results, this](int32_t i) {
        const int32_t jobStart = jobsPerThread[i].first;
        const int32_t jobEnd = jobsPerThread[i].second;
        this->WorkerMapper_(batchChunks, jobStart, jobEnd, results);
    };

    Parallel::Dispatch(faf, jobsPerThread.size(), Submit);

    if (alignSettings.align) {
        // Compute the reverse complements for alignment.
        std::vector<std::vector<FastaSequenceId>> querySeqsRev =
            ComputeQueryReverseComplements(batchChunks, results, true, faf);
        // Convert the reverse sequences to FastaSequenceCachedStore.
        std::vector<FastaSequenceCachedStore> querySeqsRevStore;
        for (const auto& chunkRevQueries : querySeqsRev) {
            querySeqsRevStore.emplace_back(FastaSequenceCachedStore(chunkRevQueries));
        }

        // Prepare the sequences for alignment.
        std::vector<PairForBatchAlignment> partsGlobal;
        std::vector<PairForBatchAlignment> partsSemiglobal;
        std::vector<AlignmentStitchInfo> alnStitchInfo;
        int32_t longestSequenceForAln = 0;
        PrepareSequencesForBatchAlignmentInParallel(
            faf, batchChunks, querySeqsRevStore, results, alignSettings.selfHitPolicy, false,
            partsGlobal, partsSemiglobal, alnStitchInfo, longestSequenceForAln);
        PBLOG_TRACE << "partsGlobal.size() = " << partsGlobal.size();
        PBLOG_TRACE << "partsSemiglobal.size() = " << partsSemiglobal.size();

        // Internal alignment on CPU.
        std::vector<AlignmentResult> internalAlns;
        AlignPartsOnCpu(alignSettings.alignerTypeGlobal, alignSettings.alnParamsGlobal,
                        alignSettings.alignerTypeExt, alignSettings.alnParamsExt, partsGlobal, faf,
                        internalAlns);
        PBLOG_TRACE << "internalAlns.size() = " << internalAlns.size();

        // Flank alignment on CPU.
        std::vector<AlignmentResult> flankAlns;
        AlignPartsOnCpu(alignSettings.alignerTypeGlobal, alignSettings.alnParamsGlobal,
                        alignSettings.alignerTypeExt, alignSettings.alnParamsExt, partsSemiglobal,
                        faf, flankAlns);
        PBLOG_TRACE << "flankAlns.size() = " << flankAlns.size();

        StitchAlignmentsInParallel(results, batchChunks, querySeqsRevStore, internalAlns, flankAlns,
                                   alnStitchInfo, faf);

        SetUnalignedAndMockedMappings(
            results, alignSettings.selfHitPolicy == MapperSelfHitPolicy::PERFECT_ALIGNMENT,
            alignSettings.alnParamsGlobal.matchScore);

        UpdateSecondaryAndFilter(results, faf, batchChunks);
    }

    return results;
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

int32_t AlignPartsOnCpu(const AlignerType alignerTypeGlobal,
                        const AlignmentParameters& alnParamsGlobal,
                        const AlignerType alignerTypeExt, const AlignmentParameters& alnParamsExt,
                        const std::vector<PairForBatchAlignment>& parts, const int32_t numThreads,
                        std::vector<AlignmentResult>& retAlns)
{
    Parallel::FireAndForget faf(numThreads);
    const int32_t result = AlignPartsOnCpu(alignerTypeGlobal, alnParamsGlobal, alignerTypeExt,
                                           alnParamsExt, parts, &faf, retAlns);
    faf.Finalize();
    return result;
}
int32_t AlignPartsOnCpu(const AlignerType alignerTypeGlobal,
                        const AlignmentParameters& alnParamsGlobal,
                        const AlignerType alignerTypeExt, const AlignmentParameters& alnParamsExt,
                        const std::vector<PairForBatchAlignment>& parts,
                        Parallel::FireAndForget* faf, std::vector<AlignmentResult>& retAlns)
{
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

            aligner.AddSequencePairForExtensionAlignment(query, target);
        } else {
            if (part.regionType == RegionType::GLOBAL) {
                aligner.AddSequencePairForGlobalAlignment(
                    std::string_view(part.query, part.queryLen),
                    std::string_view(part.target, part.targetLen));
            } else {
                aligner.AddSequencePairForExtensionAlignment(
                    std::string_view(part.query, part.queryLen),
                    std::string_view(part.target, part.targetLen));
            }
        }
    }
    aligner.AlignAll();

    const std::vector<AlignmentResult>& partInternalAlns = aligner.GetAlnResults();
    int32_t numNotValid = 0;
    for (size_t i = 0; i < partInternalAlns.size(); ++i) {
        const auto& aln = partInternalAlns[i];
        if (aln.valid == false) {
            ++numNotValid;
        }
        retAlns[partIds[i]] = std::move(partInternalAlns[i]);
    }
    return numNotValid;
}

}  // namespace Pancake
}  // namespace PacBio
