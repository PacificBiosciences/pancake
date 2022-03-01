// Copyright (c) 2019, Pacific Biosciences of California, Inc.
// All rights reserved.
// See LICENSE.txt.
//
// Contributions from NVIDIA are Copyright (c) 2021, NVIDIA Corporation.
// All rights reserved.
// SPDX-License-Identifier: BSD-3-Clause-Clear
//
// Authors: Ivan Sovic
#include <pacbio/pancake/AlignerBatchGPU.h>
#include <pacbio/pancake/MapperBatchGPU.h>

#include <pacbio/alignment/AlignmentTools.h>
#include <pacbio/pancake/AlignerBase.h>
#include <pacbio/pancake/AlignmentSeeded.h>
#include <pacbio/pancake/OverlapWriterFactory.h>
#include <pacbio/util/TicToc.h>
#include <pbcopper/logging/Logging.h>
#include <pbcopper/parallel/FireAndForget.h>
#include <pbcopper/parallel/WorkQueue.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <claraparabricks/genomeworks/utils/allocator.hpp>
#include <iostream>
#include <sstream>
#include <tuple>

namespace PacBio {
namespace Pancake {

MapperBatchGPU::MapperBatchGPU(const MapperCLRAlignSettings& alignSettings, int32_t numThreads,
                               int32_t gpuStartBandwidth, int32_t gpuMaxBandwidth,
                               uint32_t gpuDeviceId, int64_t gpuMemoryBytes,
                               bool alignRemainingOnCpu)
    : alignSettings_{alignSettings}
    , gpuStartBandwidth_(gpuStartBandwidth)
    , gpuMaxBandwidth_(gpuMaxBandwidth)
    , alignRemainingOnCpu_(alignRemainingOnCpu)
    , faf_{nullptr}
    , fafFallback_(nullptr)
    , aligner_{nullptr}
{
    fafFallback_ = std::make_unique<Parallel::FireAndForget>(numThreads);
    faf_ = fafFallback_.get();
    aligner_ = std::make_unique<AlignerBatchGPU>(alignSettings.alnParamsGlobal, gpuStartBandwidth,
                                                 gpuDeviceId, gpuMemoryBytes);
}

MapperBatchGPU::MapperBatchGPU(const MapperCLRAlignSettings& alignSettings,
                               Parallel::FireAndForget* faf, int32_t gpuStartBandwidth,
                               int32_t gpuMaxBandwidth, uint32_t gpuDeviceId,
                               int64_t gpuMemoryBytes, bool alignRemainingOnCpu)
    : alignSettings_{alignSettings}
    , gpuStartBandwidth_(gpuStartBandwidth)
    , gpuMaxBandwidth_(gpuMaxBandwidth)
    , alignRemainingOnCpu_(alignRemainingOnCpu)
    , faf_{faf}
    , fafFallback_(nullptr)
    , aligner_{std::make_unique<AlignerBatchGPU>(alignSettings.alnParamsGlobal, gpuStartBandwidth,
                                                 gpuDeviceId, gpuMemoryBytes)}
{}

MapperBatchGPU::~MapperBatchGPU()
{
    if (fafFallback_) {
        fafFallback_->Finalize();
    }
}

std::vector<std::vector<MapperBaseResult>> MapperBatchGPU::MapAndAlign(
    const std::vector<MapperBatchChunk>& batchData)
{
    assert(aligner_);

    try {
        return MapAndAlignImpl_(batchData, alignSettings_, alignRemainingOnCpu_, gpuStartBandwidth_,
                                gpuMaxBandwidth_, *aligner_, faf_);
    } catch (const std::exception& e) {
        // Log, but do not fail. Important for clients of this class.
        // Return a vector of the size of the input, but with empty subvectors.
        PBLOG_DEBUG << "MapperBatchGPU generated an exception in " << std::string(__FUNCTION__)
                    << ". Message: " << e.what();
        std::vector<std::vector<MapperBaseResult>> results(batchData.size());
        return results;
    }
}

std::vector<std::vector<MapperBaseResult>> MapperBatchGPU::MapAndAlignImpl_(
    const std::vector<MapperBatchChunk>& batchChunks, const MapperCLRAlignSettings& alignSettings,
    bool alignRemainingOnCpu, int32_t gpuStartBandwidth, int32_t gpuMaxBandwidth,
    AlignerBatchGPU& aligner, Parallel::FireAndForget* faf)
{
    const int32_t numThreads = faf ? faf->NumThreads() : 1;
    // Determine how many records should land in each thread, spread roughly evenly.
    const int32_t numRecords = batchChunks.size();
    const std::vector<std::pair<int32_t, int32_t>> jobsPerThread =
        PacBio::Pancake::DistributeJobLoad<int32_t>(numThreads, numRecords);

    // Run the mapping in parallel.
    std::vector<std::vector<MapperBaseResult>> results(batchChunks.size());

    const auto Submit = [&batchChunks, &jobsPerThread, &results](int32_t idx) {
        const int32_t jobStart = jobsPerThread[idx].first;
        const int32_t jobEnd = jobsPerThread[idx].second;
        WorkerMapper_(batchChunks, jobStart, jobEnd, results);
    };
    Parallel::Dispatch(faf, jobsPerThread.size(), Submit);

    // Align the mappings if required.
    if (alignSettings.align) {
        // Sanity check.
        if (gpuStartBandwidth <= 0) {
            throw std::runtime_error(
                "The startBandwidth needs to be a positive number. gpuStartBandwidth = " +
                std::to_string(gpuStartBandwidth));
        }

        // Compute the reverse complements for alignment.
        std::vector<std::vector<FastaSequenceId>> querySeqsRev =
            ComputeQueryReverseComplements(batchChunks, results, true, faf);
        // Convert the reverse sequences to FastaSequenceCachedStore.
        std::vector<FastaSequenceCachedStore> querySeqsRevStore;
        for (const auto& chunkRevQueries : querySeqsRev) {
            querySeqsRevStore.emplace_back(FastaSequenceCachedStore(chunkRevQueries));
        }

        PBLOG_TRACE << "Preparing parts for alignment.";

        // Prepare the sequences for alignment.
        std::vector<PairForBatchAlignment> partsGlobal;
        std::vector<PairForBatchAlignment> partsSemiglobal;
        std::vector<AlignmentStitchInfo> alnStitchInfo;
        int32_t longestSequenceForAln = 0;
        PrepareSequencesForBatchAlignment(batchChunks, querySeqsRevStore, results,
                                          alignSettings.selfHitPolicy, partsGlobal, partsSemiglobal,
                                          alnStitchInfo, longestSequenceForAln);
        PBLOG_TRACE << "partsGlobal.size() = " << partsGlobal.size();
        PBLOG_TRACE << "partsSemiglobal.size() = " << partsSemiglobal.size();

        // Global alignment on GPU. Try using different bandwidths,
        // increasing the bandwidth for the failed parts each iteration.
        std::vector<AlignmentResult> internalAlns;
        int32_t currentBandwidth = gpuStartBandwidth;
        int32_t numInternalNotValid = 0;
        const int32_t maxBandwidth =
            (gpuMaxBandwidth <= 0) ? (longestSequenceForAln * 2 + 1) : gpuMaxBandwidth;

        while (true) {
            PBLOG_TRACE << "Trying bandwidth: " << currentBandwidth;
            aligner.ResetMaxBandwidth(currentBandwidth);
            numInternalNotValid = AlignPartsOnGPU_(aligner, partsGlobal, internalAlns);
            if (numInternalNotValid == 0) {
                break;
            }
            if (currentBandwidth >= maxBandwidth) {
                break;
            }
            currentBandwidth *= 2;
            // Ensure that the last bandwidth tried is the maxBandwidth.
            if (currentBandwidth > maxBandwidth) {
                currentBandwidth = maxBandwidth;
            }
        }
        // Fallback to the CPU if there are any unaligned parts left.
        if (alignRemainingOnCpu && numInternalNotValid > 0) {
            PBLOG_TRACE << "Trying to align remaining parts on CPU.";
            const int32_t numNotValidInternal =
                AlignPartsOnCpu(alignSettings.alignerTypeGlobal, alignSettings.alnParamsGlobal,
                                alignSettings.alignerTypeExt, alignSettings.alnParamsExt,
                                partsGlobal, faf, internalAlns);
            PBLOG_TRACE << "Total not valid: " << numNotValidInternal << " / "
                        << internalAlns.size() << "\n";
        }
        PBLOG_TRACE << "internalAlns.size() = " << internalAlns.size();

        // Flank alignment on CPU.
        std::vector<AlignmentResult> flankAlns;
        const int32_t numNotValidFlanks =
            AlignPartsOnCpu(alignSettings.alignerTypeGlobal, alignSettings.alnParamsGlobal,
                            alignSettings.alignerTypeExt, alignSettings.alnParamsExt,
                            partsSemiglobal, faf, flankAlns);
        PBLOG_TRACE << "Total not valid: " << numNotValidFlanks << " / " << flankAlns.size()
                    << "\n";

        StitchAlignmentsInParallel(results, batchChunks, querySeqsRevStore, internalAlns, flankAlns,
                                   alnStitchInfo, faf);

        SetUnalignedAndMockedMappings(
            results, alignSettings.selfHitPolicy == MapperSelfHitPolicy::PERFECT_ALIGNMENT,
            alignSettings.alnParamsGlobal.matchScore);

        UpdateSecondaryAndFilter(results, faf, batchChunks);
    }

    return results;
}

void MapperBatchGPU::WorkerMapper_(const std::vector<MapperBatchChunk>& batchChunks,
                                   int32_t startId, int32_t endId,
                                   std::vector<std::vector<MapperBaseResult>>& results)
{

    // Create one mapper for this thread so we can reuse the memory.
    MapperCLRSettings settingsCopy;
    MapperCLR mapper(settingsCopy);

    for (int32_t i = startId; i < endId; ++i) {
        const auto& chunk = batchChunks[i];

        // Create a copy of the settings so that we can turn off the alignment.
        settingsCopy.map = chunk.mapSettings;
        settingsCopy.align.align = false;
        mapper.UpdateSettings(settingsCopy);

        try {
            results[i] = mapper.MapAndAlign(chunk.targetSeqs, chunk.querySeqs);
        } catch (const std::exception& e) {
            // Log, but do not fail, so that we do not kill all other alignments in the batch.
            PBLOG_DEBUG << "MapperBatchGPU caught an exception generated by MapperCLR in "
                        << std::string(__FUNCTION__) << ". Worker chunk i = " << i
                        << ". Message: " << e.what();
            results[i] = std::vector<PacBio::Pancake::MapperBaseResult>(chunk.querySeqs.Size());
        }
    }
}

int32_t MapperBatchGPU::AlignPartsOnGPU_(AlignerBatchGPU& aligner,
                                         const std::vector<PairForBatchAlignment>& parts,
                                         std::vector<AlignmentResult>& retInternalAlns)
{
    retInternalAlns.resize(parts.size());

    int32_t totalNumNotValid = 0;
    size_t partId = 0;
    while (partId < parts.size()) {
        aligner.Clear();

        std::vector<size_t> partIds;
        partIds.reserve(parts.size());

        PBLOG_TRACE << "Preparing sequences for GPU alignment.";
        for (; partId < parts.size(); ++partId) {
            const PairForBatchAlignment& part = parts[partId];

            if (retInternalAlns[partId].valid) {
                continue;
            }
            partIds.emplace_back(partId);

            StatusAddSequencePair rv;
            if (part.regionType == RegionType::FRONT) {
                // Reverse the sequences for front flank alignment. No need to complement.
                std::string query(part.query, part.queryLen);
                std::reverse(query.begin(), query.end());
                std::string target(part.target, part.targetLen);
                std::reverse(target.begin(), target.end());
                rv = aligner.AddSequencePairForGlobalAlignment(query.c_str(), part.queryLen,
                                                               target.c_str(), part.targetLen);
            } else {
                rv = aligner.AddSequencePairForGlobalAlignment(part.query, part.queryLen,
                                                               part.target, part.targetLen);
            }

            if (rv == StatusAddSequencePair::EXCEEDED_MAX_ALIGNMENTS) {
                break;

            } else if (rv != StatusAddSequencePair::OK) {
                throw std::runtime_error(
                    "Error occurred while trying to add sequences for batch alignment to "
                    "AlignerBatchGPU.");
            }
        }

        PBLOG_TRACE << "Aligning batch of " << aligner.BatchSize() << " sequence pairs.";
        aligner.AlignAll();

        std::vector<AlignmentResult>& partInternalAlns = aligner.GetAlnResults();

        int32_t numNotValid = 0;
        for (size_t i = 0; i < partInternalAlns.size(); ++i) {
            const auto& aln = partInternalAlns[i];
            if (aln.valid == false) {
                ++numNotValid;
            }
            // retInternalAlns.emplace_back(std::move(partInternalAlns[i]));
            retInternalAlns[partIds[i]] = std::move(partInternalAlns[i]);
        }
        totalNumNotValid += numNotValid;
    }
    PBLOG_TRACE << "Total not valid: " << totalNumNotValid << " / " << retInternalAlns.size()
                << "\n";
    return totalNumNotValid;
}

}  // namespace Pancake
}  // namespace PacBio
