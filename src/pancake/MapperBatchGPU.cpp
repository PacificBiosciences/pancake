// Copyright (c) 2019, Pacific Biosciences of California, Inc.
// All rights reserved.
// See LICENSE.txt.
//
// Contributions from NVIDIA are Copyright (c) 2021, NVIDIA Corporation.
// All rights reserved.
// SPDX-License-Identifier: BSD-3-Clause-Clear
//
// Authors: Ivan Sovic

#include <pancake/MapperBatchGPU.hpp>

#include <pancake/AlignerBase.hpp>
#include <pancake/AlignerBatchGPU.hpp>
#include <pancake/AlignmentSeeded.hpp>
#include <pancake/AlignmentTools.hpp>
#include <pancake/MapperBatchGPUUtility.hpp>
#include <pancake/OverlapWriterFactory.hpp>
#include <pancake/util/TicToc.hpp>

#include <pbcopper/logging/Logging.h>
#include <pbcopper/parallel/FireAndForget.h>
#include <pbcopper/parallel/WorkQueue.h>
#include <claraparabricks/genomeworks/utils/allocator.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <sstream>
#include <tuple>

namespace PacBio {
namespace Pancake {

MapperBatchGPU::MapperBatchGPU(const MapperCLRAlignSettings& alignSettings, int32_t numThreads,
                               int32_t gpuStartBandwidth, int32_t gpuMaxBandwidth,
                               uint32_t gpuDeviceId, int64_t gpuMemoryBytes,
                               const int32_t maxAllowedGapForGpu, bool alignRemainingOnCpu)
    : alignSettings_{alignSettings}
    , gpuStartBandwidth_(gpuStartBandwidth)
    , gpuMaxBandwidth_(gpuMaxBandwidth)
    , maxAllowedGapForGpu_(maxAllowedGapForGpu)
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
                               int64_t gpuMemoryBytes, int32_t maxAllowedGapForGpu,
                               bool alignRemainingOnCpu)
    : alignSettings_{alignSettings}
    , gpuStartBandwidth_(gpuStartBandwidth)
    , gpuMaxBandwidth_(gpuMaxBandwidth)
    , maxAllowedGapForGpu_(maxAllowedGapForGpu)
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
        return MapAndAlignImpl_(batchData, alignSettings_, maxAllowedGapForGpu_,
                                alignRemainingOnCpu_, gpuStartBandwidth_, gpuMaxBandwidth_,
                                *aligner_, faf_);
    } catch (const std::exception& e) {
        // Log, but do not fail. Important for clients of this class.
        // Return a vector of the size of the input, but with empty results for each query.
        PBLOG_DEBUG << "MapperBatchGPU generated an exception in " << std::string(__FUNCTION__)
                    << ". Message: " << e.what();
        std::vector<std::vector<MapperBaseResult>> results(batchData.size());
        for (size_t i = 0; i < batchData.size(); ++i) {
            results[i].resize(batchData[i].querySeqs.Size());
        }
        return results;
    }
}

std::vector<std::vector<MapperBaseResult>> MapperBatchGPU::MapAndAlignImpl_(
    const std::vector<MapperBatchChunk>& batchChunks, const MapperCLRAlignSettings& alignSettings,
    const int32_t maxAllowedGapForGpu, bool alignRemainingOnCpu, int32_t gpuStartBandwidth,
    int32_t gpuMaxBandwidth, AlignerBatchGPU& aligner, Parallel::FireAndForget* faf) const
{
    const int32_t numThreads = faf ? faf->NumThreads() : 1;
    // Determine how many records should land in each thread, spread roughly evenly.
    const int32_t numRecords = batchChunks.size();
    const std::vector<std::pair<int32_t, int32_t>> jobsPerThread =
        PacBio::Pancake::DistributeJobLoad<int32_t>(numThreads, numRecords);

    // Run the mapping in parallel.
    std::vector<std::vector<MapperBaseResult>> results(batchChunks.size());

    const auto Submit = [&batchChunks, &jobsPerThread, &results, this](int32_t idx) {
        const int32_t jobStart = jobsPerThread[idx].first;
        const int32_t jobEnd = jobsPerThread[idx].second;
        this->WorkerMapper_(batchChunks, jobStart, jobEnd, results);
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
        PrepareSequencesForBatchAlignmentInParallel(
            faf, batchChunks, querySeqsRevStore, results, alignSettings.selfHitPolicy, true,
            partsGlobal, partsSemiglobal, alnStitchInfo, longestSequenceForAln);
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
            numInternalNotValid =
                AlignPartsOnGPU(partsGlobal, aligner, maxAllowedGapForGpu, internalAlns);
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

}  // namespace Pancake
}  // namespace PacBio
