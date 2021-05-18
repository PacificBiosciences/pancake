// Authors: Ivan Sovic
#include <pacbio/pancake/AlignerBatchGPUGenomeWorksKSW2.h>
#include <pacbio/pancake/MapperBatchGPUGenomeWorksKSW2.h>

#include <pacbio/alignment/AlignmentTools.h>
#include <pacbio/pancake/AlignerBase.h>
#include <pacbio/pancake/AlignmentSeeded.h>
#include <pacbio/pancake/OverlapWriterFactory.h>
#include <pacbio/util/TicToc.h>
#include <pbcopper/logging/Logging.h>
#include <pbcopper/parallel/FireAndForget.h>
#include <pbcopper/parallel/WorkQueue.h>
#include <pbcopper/utility/Stopwatch.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <claraparabricks/genomeworks/utils/allocator.hpp>
#include <iostream>
#include <sstream>
#include <tuple>

namespace PacBio {
namespace Pancake {

MapperBatchGPUGenomeWorksKSW2::MapperBatchGPUGenomeWorksKSW2(
    const MapperCLRAlignSettings& alignSettings, int32_t numThreads, int32_t gpuStartBandwidth,
    int32_t gpuMaxBandwidth, uint32_t gpuDeviceId, int64_t gpuMemoryBytes, bool alignRemainingOnCpu)
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
    aligner_ = std::make_unique<AlignerBatchGPUGenomeWorksKSW2>(faf_, alignSettings.alnParamsGlobal,
                                                                gpuDeviceId, gpuMemoryBytes);
}

MapperBatchGPUGenomeWorksKSW2::MapperBatchGPUGenomeWorksKSW2(
    const MapperCLRAlignSettings& alignSettings, Parallel::FireAndForget* faf,
    int32_t gpuStartBandwidth, int32_t gpuMaxBandwidth, uint32_t gpuDeviceId,
    int64_t gpuMemoryBytes, bool alignRemainingOnCpu)
    : alignSettings_{alignSettings}
    , gpuStartBandwidth_(gpuStartBandwidth)
    , gpuMaxBandwidth_(gpuMaxBandwidth)
    , alignRemainingOnCpu_(alignRemainingOnCpu)
    , faf_{faf}
    , fafFallback_(nullptr)
    , aligner_{std::make_unique<AlignerBatchGPUGenomeWorksKSW2>(faf, alignSettings.alnParamsGlobal,
                                                                gpuDeviceId, gpuMemoryBytes)}
{
}

MapperBatchGPUGenomeWorksKSW2::~MapperBatchGPUGenomeWorksKSW2()
{
    if (fafFallback_) {
        fafFallback_->Finalize();
    }
}

std::vector<std::vector<MapperBaseResult>> MapperBatchGPUGenomeWorksKSW2::MapAndAlign(
    const std::vector<MapperBatchChunk>& batchData)
{
    assert(aligner_);
    return MapAndAlignImpl_(batchData, alignSettings_, alignRemainingOnCpu_, gpuStartBandwidth_,
                            gpuMaxBandwidth_, *aligner_, faf_);
}

std::vector<std::vector<MapperBaseResult>> MapperBatchGPUGenomeWorksKSW2::MapAndAlignImpl_(
    const std::vector<MapperBatchChunk>& batchChunks, const MapperCLRAlignSettings& alignSettings,
    bool alignRemainingOnCpu, int32_t gpuStartBandwidth, int32_t gpuMaxBandwidth,
    AlignerBatchGPUGenomeWorksKSW2& aligner, Parallel::FireAndForget* faf)
{
    PacBio::Utility::Stopwatch timer;
    int64_t cpuTime = 0;
    int64_t gpuTime = 0;

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
    timer.Freeze();
    cpuTime += timer.ElapsedNanoseconds();
    PBLOG_INFO << "CPU Mapping            : " << timer.ElapsedTime();

    // Align the mappings if required.
    if (alignSettings.align) {
        int64_t seedAlnCpuTime = 0;
        int64_t seedAlnGpuTime = 0;
        timer.Reset();
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
        timer.Freeze();
        seedAlnCpuTime += timer.ElapsedNanoseconds();
        PBLOG_INFO << "CPU RevComp            : " << timer.ElapsedTime();
        timer.Reset();

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
        timer.Freeze();
        cpuTime += timer.ElapsedNanoseconds();
        PBLOG_INFO << "CPU Prepare            : " << timer.ElapsedTime();
        timer.Reset();

        // Global alignment on GPU. Try using different bandwidths,
        // increasing the bandwidth for the failed parts each iteration.
        std::vector<AlignmentResult> internalAlns;
        int32_t currentBandwidth = gpuStartBandwidth;
        int32_t numInternalNotValid = 0;
        const int32_t maxBandwidth =
            (gpuMaxBandwidth <= 0) ? (longestSequenceForAln * 2 + 1) : gpuMaxBandwidth;

        while (true) {
            PBLOG_TRACE << "Trying bandwidth: " << currentBandwidth;
            // std::cerr << "Trying bandwidth: " << currentBandwidth << "\n";
            timer.Reset();
            aligner.ResetMaxBandwidth(currentBandwidth);
            timer.Freeze();
            cpuTime += timer.ElapsedNanoseconds();
            PBLOG_INFO << "CPU Reset BW           : " << timer.ElapsedTime();
            timer.Reset();
            numInternalNotValid = AlignPartsOnGPU_(aligner, partsGlobal, internalAlns,
                                                   seedAlnCpuTime, seedAlnGpuTime);
            // std::cerr << "numInternalNotValid = " << numInternalNotValid << "\n";
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
        cpuTime += seedAlnCpuTime;
        gpuTime += seedAlnGpuTime;
        PBLOG_INFO << "CPU Internal Prepare   : "
                   << PacBio::Utility::Stopwatch::PrettyPrintNanoseconds(seedAlnCpuTime);
        PBLOG_INFO << "GPU Internal Alignment : "
                   << PacBio::Utility::Stopwatch::PrettyPrintNanoseconds(seedAlnGpuTime);
        timer.Reset();

        // std::cerr << "alignRemainingOnCpu = " << alignRemainingOnCpu << "\n";
        // Fallback to the CPU if there are any unaligned parts left.
        int64_t prepareTime = 0;
        int64_t alignTime = 0;
        if (alignRemainingOnCpu && numInternalNotValid > 0) {
            PBLOG_TRACE << "Trying to align remaining parts on CPU.";
            const int32_t numNotValidInternal =
                AlignPartsOnCpu(alignSettings.alignerTypeGlobal, alignSettings.alnParamsGlobal,
                                alignSettings.alignerTypeExt, alignSettings.alnParamsExt,
                                partsGlobal, faf, internalAlns, prepareTime, alignTime);
            PBLOG_TRACE << "Total not valid: " << numNotValidInternal << " / "
                        << internalAlns.size() << "\n";
            // std::cerr << "Total not valid: " << numNotValidInternal << " / " << internalAlns.size()
            //          << "\n";
        }
        PBLOG_TRACE << "internalAlns.size() = " << internalAlns.size();

        timer.Freeze();
        cpuTime += timer.ElapsedNanoseconds();
        PBLOG_INFO << "CPU Fallback           : " << timer.ElapsedTime() << ' '
                   << numInternalNotValid;
        timer.Reset();

        // Flank alignment on CPU.
        std::vector<AlignmentResult> flankAlns;
        prepareTime = 0;
        alignTime = 0;
        const int32_t numNotValidFlanks =
            AlignPartsOnCpu(alignSettings.alignerTypeGlobal, alignSettings.alnParamsGlobal,
                            alignSettings.alignerTypeExt, alignSettings.alnParamsExt,
                            partsSemiglobal, faf, flankAlns, prepareTime, alignTime);
        PBLOG_TRACE << "Total not valid: " << numNotValidFlanks << " / " << flankAlns.size()
                    << "\n";
        timer.Freeze();
        cpuTime += timer.ElapsedNanoseconds();
        PBLOG_INFO << "CPU Flanks Prepare     : "
                   << PacBio::Utility::Stopwatch::PrettyPrintNanoseconds(prepareTime);
        PBLOG_INFO << "CPU Flanks Alignment   : "
                   << PacBio::Utility::Stopwatch::PrettyPrintNanoseconds(alignTime);
        timer.Reset();

        StitchAlignmentsInParallel(results, batchChunks, querySeqsRevStore, internalAlns, flankAlns,
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
    PBLOG_INFO << "GPU Time               : "
               << PacBio::Utility::Stopwatch::PrettyPrintNanoseconds(gpuTime);

    return results;
}

void MapperBatchGPUGenomeWorksKSW2::WorkerMapper_(
    const std::vector<MapperBatchChunk>& batchChunks, int32_t startId, int32_t endId,
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

int32_t MapperBatchGPUGenomeWorksKSW2::AlignPartsOnGPU_(
    AlignerBatchGPUGenomeWorksKSW2& aligner, const std::vector<PairForBatchAlignment>& parts,
    std::vector<AlignmentResult>& retInternalAlns, int64_t& cpuTime, int64_t& gpuTime)
{
    PacBio::Utility::Stopwatch timer;
    retInternalAlns.resize(parts.size());
    cpuTime += timer.ElapsedNanoseconds();

    int32_t totalNumNotValid = 0;
    size_t partId = 0;
    while (partId < parts.size()) {
        timer.Reset();
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
                    "AlignerBatchGPUGenomeWorksKSW2.");
            }
        }
        cpuTime += timer.ElapsedNanoseconds();

        PBLOG_TRACE << "Aligning batch of " << aligner.BatchSize() << " sequence pairs.";
        std::pair<int64_t, int64_t> seedTimings = aligner.AlignAll();
        cpuTime += seedTimings.first;
        gpuTime += seedTimings.second;
        timer.Reset();

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
        cpuTime += timer.ElapsedNanoseconds();
    }
    PBLOG_TRACE << "Total not valid: " << totalNumNotValid << " / " << retInternalAlns.size()
                << "\n";
    return totalNumNotValid;
}

}  // namespace Pancake
}  // namespace PacBio
