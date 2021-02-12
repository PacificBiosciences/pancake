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
#include <claraparabricks/genomeworks/utils/allocator.hpp>
#include <iostream>
#include <sstream>
#include <tuple>

namespace PacBio {
namespace Pancake {

MapperBatchGPU::MapperBatchGPU(const MapperCLRSettings& settings, int32_t numThreads,
                               int32_t gpuStartBandwidth, int32_t gpuMaxBandwidth,
                               uint32_t gpuDeviceId, double gpuMaxFreeMemoryFraction,
                               int64_t gpuMaxMemoryCap, bool alignRemainingOnCpu)
    : MapperBatchGPU(settings, nullptr, gpuStartBandwidth, gpuMaxBandwidth, gpuDeviceId,
                     gpuMaxFreeMemoryFraction, gpuMaxMemoryCap, alignRemainingOnCpu)
{
    fafFallback_ = std::make_unique<Parallel::FireAndForget>(numThreads);
    faf_ = fafFallback_.get();
}

MapperBatchGPU::MapperBatchGPU(const MapperCLRSettings& settings, Parallel::FireAndForget* faf,
                               int32_t gpuStartBandwidth, int32_t gpuMaxBandwidth,
                               uint32_t gpuDeviceId, double gpuMaxFreeMemoryFraction,
                               int64_t gpuMaxMemoryCap, bool alignRemainingOnCpu)
    : settings_{settings}
    , gpuStartBandwidth_(gpuStartBandwidth)
    , gpuMaxBandwidth_(gpuMaxBandwidth)
    , gpuDeviceId_(gpuDeviceId)
    , gpuMaxFreeMemoryFraction_(gpuMaxFreeMemoryFraction)
    , gpuMaxMemoryCap_(gpuMaxMemoryCap)
    , gpuStream_(0)
    , alignRemainingOnCpu_(alignRemainingOnCpu)
    , faf_{faf}
{
    const int64_t requestedMemory = std::min(
        gpuMaxMemoryCap_, PacBio::Pancake::ComputeMaxGPUMemory(1, gpuMaxFreeMemoryFraction_));

    GW_CU_CHECK_ERR(cudaSetDevice(gpuDeviceId_));
    GW_CU_CHECK_ERR(cudaStreamCreate(&gpuStream_));

    deviceAllocator_ =
        claraparabricks::genomeworks::create_default_device_allocator(requestedMemory, gpuStream_);
}

MapperBatchGPU::~MapperBatchGPU()
{
    if (gpuStream_ != 0) {
        GW_CU_CHECK_ERR(cudaStreamDestroy(gpuStream_));
    }
    if (fafFallback_) {
        fafFallback_->Finalize();
    }
}

std::vector<std::vector<MapperBaseResult>> MapperBatchGPU::MapAndAlign(
    const std::vector<MapperBatchChunk>& batchData)
{
    return MapAndAlignImpl_(batchData, settings_, alignRemainingOnCpu_, gpuStartBandwidth_,
                            gpuMaxBandwidth_, gpuDeviceId_, gpuStream_, deviceAllocator_, faf_);
}

std::vector<std::vector<MapperBaseResult>> MapperBatchGPU::MapAndAlignImpl_(
    const std::vector<MapperBatchChunk>& batchChunks, const MapperCLRSettings& settings,
    bool alignRemainingOnCpu, int32_t gpuStartBandwidth, int32_t gpuMaxBandwidth,
    uint32_t gpuDeviceId, cudaStream_t& gpuStream,
    claraparabricks::genomeworks::DefaultDeviceAllocator& deviceAllocator,
    Parallel::FireAndForget* faf)
{
    const int32_t numThreads = faf ? faf->NumThreads() : 1;
    // Determine how many records should land in each thread, spread roughly evenly.
    const int32_t numRecords = batchChunks.size();
    const std::vector<std::pair<int32_t, int32_t>> jobsPerThread =
        PacBio::Pancake::DistributeJobLoad<int32_t>(numThreads, numRecords);

    // Create a mapper for each thread.
    MapperCLRSettings settingsCopy = settings;
    settingsCopy.align = false;
    std::vector<std::unique_ptr<MapperCLR>> mappers;
    for (size_t i = 0; i < jobsPerThread.size(); ++i) {
        auto mapper = std::make_unique<MapperCLR>(settingsCopy);
        mappers.emplace_back(std::move(mapper));
    }

    // Run the mapping in parallel.
    std::vector<std::vector<MapperBaseResult>> results(batchChunks.size());

    const auto Submit = [&batchChunks, &mappers, &jobsPerThread, &results](int32_t idx) {
        const int32_t jobStart = jobsPerThread[idx].first;
        const int32_t jobEnd = jobsPerThread[idx].second;
        WorkerMapper_(batchChunks, jobStart, jobEnd, *mappers[idx], results);
    };
    Parallel::Dispatch(faf, jobsPerThread.size(), Submit);

    // Align the mappings if required.
    if (settings.align) {
        // Sanity check.
        if (gpuStartBandwidth <= 0) {
            throw std::runtime_error(
                "The startBandwidth needs to be a positive number. gpuStartBandwidth = " +
                std::to_string(gpuStartBandwidth));
        }

        // Compute the reverse complements for alignment.
        std::vector<std::vector<std::string>> querySeqsRev;
        for (const auto& chunk : batchChunks) {
            std::vector<std::string> revSeqs;
            for (const auto& query : chunk.querySeqs) {
                revSeqs.emplace_back(
                    PacBio::Pancake::ReverseComplement(query.c_str(), 0, query.size()));
            }
            querySeqsRev.emplace_back(std::move(revSeqs));
        }

        PBLOG_TRACE << "Preparing parts for alignment.";

        // Prepare the sequences for alignment.
        std::vector<PairForBatchAlignment> partsGlobal;
        std::vector<PairForBatchAlignment> partsSemiglobal;
        std::vector<AlignmentStitchInfo> alnStitchInfo;
        int32_t longestSequenceForAln = 0;
        PrepareSequencesForBatchAlignment(batchChunks, querySeqsRev, results, partsGlobal,
                                          partsSemiglobal, alnStitchInfo, longestSequenceForAln);
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
            numInternalNotValid =
                AlignPartsOnGPU_(currentBandwidth, gpuDeviceId, gpuStream, deviceAllocator,
                                 settings.alnParamsGlobal, partsGlobal, internalAlns);
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
            const int32_t numNotValidInternal = AlignPartsOnCpu(
                settings.alignerTypeGlobal, settings.alnParamsGlobal, settings.alignerTypeExt,
                settings.alnParamsExt, partsGlobal, numThreads, internalAlns);
            PBLOG_TRACE << "Total not valid: " << numNotValidInternal << " / "
                        << internalAlns.size() << "\n";
        }
        PBLOG_TRACE << "internalAlns.size() = " << internalAlns.size();

        // Flank alignment on CPU.
        std::vector<AlignmentResult> flankAlns;
        const int32_t numNotValidFlanks = AlignPartsOnCpu(
            settings.alignerTypeGlobal, settings.alnParamsGlobal, settings.alignerTypeExt,
            settings.alnParamsExt, partsSemiglobal, numThreads, flankAlns);
        PBLOG_TRACE << "Total not valid: " << numNotValidFlanks << " / " << flankAlns.size()
                    << "\n";

        StitchAlignments(results, batchChunks, querySeqsRev, internalAlns, flankAlns,
                         alnStitchInfo);

        UpdateSecondaryAndFilter(results, settings.secondaryAllowedOverlapFractionQuery,
                                 settings.secondaryAllowedOverlapFractionTarget,
                                 settings.secondaryMinScoreFraction, settings.bestNSecondary);
    }

    return results;
}

void MapperBatchGPU::WorkerMapper_(const std::vector<MapperBatchChunk>& batchChunks,
                                   int32_t startId, int32_t endId, MapperCLR& mapper,
                                   std::vector<std::vector<MapperBaseResult>>& results)
{
    for (int32_t i = startId; i < endId; ++i) {
        const auto& chunk = batchChunks[i];
        results[i] = mapper.MapAndAlign(chunk.targetSeqs, chunk.querySeqs);
    }
}

int32_t MapperBatchGPU::AlignPartsOnGPU_(
    int32_t gpuMaxBandwidth, uint32_t gpuDeviceId, cudaStream_t& gpuStream,
    claraparabricks::genomeworks::DefaultDeviceAllocator& deviceAllocator,
    const AlignmentParameters& alnParamsGlobal, const std::vector<PairForBatchAlignment>& parts,
    std::vector<AlignmentResult>& retInternalAlns)
{
    auto cudaAligner = claraparabricks::genomeworks::cudaaligner::create_aligner(
        claraparabricks::genomeworks::cudaaligner::AlignmentType::global_alignment, gpuMaxBandwidth,
        gpuStream, gpuDeviceId, deviceAllocator, -1);

    std::unique_ptr<AlignerBatchGPU> aligner =
        std::make_unique<AlignerBatchGPU>(alnParamsGlobal, std::move(cudaAligner));

    retInternalAlns.resize(parts.size());

    int32_t totalNumNotValid = 0;
    size_t partId = 0;
    while (partId < parts.size()) {
        aligner->Clear();

        std::vector<size_t> partIds;

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
                rv = aligner->AddSequencePair(query.c_str(), part.queryLen, target.c_str(),
                                              part.targetLen);
            } else {
                rv = aligner->AddSequencePair(part.query, part.queryLen, part.target,
                                              part.targetLen);
            }

            if (rv == StatusAddSequencePair::EXCEEDED_MAX_ALIGNMENTS) {
                break;

            } else if (rv != StatusAddSequencePair::OK) {
                throw std::runtime_error(
                    "Error occurred while trying to add sequences for batch alignment to "
                    "AlignerBatchGPU.");
            }
        }

        PBLOG_TRACE << "Aligning batch of " << aligner->BatchSize() << " sequence pairs.";
        aligner->AlignAll();

        const std::vector<AlignmentResult>& partInternalAlns = aligner->GetAlnResults();

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
