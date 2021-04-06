#include <pacbio/pancake/AlignerBatchFPGA.h>
#include <pacbio/pancake/MapperBatchFPGA.h>

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
#include <iostream>
#include <sstream>
#include <tuple>

namespace PacBio {
namespace Pancake {

MapperBatchFPGA::MapperBatchFPGA(const MapperCLRAlignSettings& alignSettings, int32_t numThreads,
                                 uint32_t fpgaDeviceId, bool alignRemainingOnCpu)
    : alignSettings_{alignSettings}
    , alignRemainingOnCpu_(alignRemainingOnCpu)
    , faf_{nullptr}
    , fafFallback_(nullptr)
    , aligner_{nullptr}
{
    fafFallback_ = std::make_unique<Parallel::FireAndForget>(numThreads);
    faf_ = fafFallback_.get();
    aligner_ = std::make_unique<AlignerBatchFPGA>(alignSettings.alnParamsGlobal, fpgaDeviceId);
}

MapperBatchFPGA::MapperBatchFPGA(const MapperCLRAlignSettings& alignSettings,
                                 Parallel::FireAndForget* faf, uint32_t fpgaDeviceId,
                                 bool alignRemainingOnCpu)
    : alignSettings_{alignSettings}
    , alignRemainingOnCpu_(alignRemainingOnCpu)
    , faf_{faf}
    , fafFallback_(nullptr)
    , aligner_{std::make_unique<AlignerBatchFPGA>(alignSettings.alnParamsGlobal, fpgaDeviceId)}
{}

MapperBatchFPGA::~MapperBatchFPGA()
{
    if (fafFallback_) {
        fafFallback_->Finalize();
    }
}

std::vector<std::vector<MapperBaseResult>> MapperBatchFPGA::MapAndAlign(
    const std::vector<MapperBatchChunk>& batchData)
{
    assert(aligner_);
    return MapAndAlignImpl_(batchData, alignSettings_, alignRemainingOnCpu_, *aligner_, faf_);
}

std::vector<std::vector<MapperBaseResult>> MapperBatchFPGA::MapAndAlignImpl_(
    const std::vector<MapperBatchChunk>& batchChunks, const MapperCLRAlignSettings& alignSettings,
    bool alignRemainingOnCpu, AlignerBatchFPGA& aligner, Parallel::FireAndForget* faf)
{
#ifdef PANCAKE_TIMINGS
    PacBio::Utility::Stopwatch timer;
#endif
    int64_t cpuTime = 0;
    int64_t fpgaTime = 0;

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
#ifdef PANCAKE_TIMINGS
    timer.Freeze();
    cpuTime += timer.ElapsedNanoseconds();
    PBLOG_INFO << "CPU Mapping            : " << timer.ElapsedTime();
#endif

    // Align the mappings if required.
    if (alignSettings.align) {
        int64_t seedAlnCpuTime = 0;
        int64_t seedAlnfpgaTime = 0;
#ifdef PANCAKE_TIMINGS
        timer.Reset();
#endif

        // Compute the reverse complements for alignment.
        std::vector<std::vector<FastaSequenceId>> querySeqsRev =
            ComputeQueryReverseComplements(batchChunks, results, true, faf);
        // Convert the reverse sequences to FastaSequenceCachedStore.
        std::vector<FastaSequenceCachedStore> querySeqsRevStore;
        for (const auto& chunkRevQueries : querySeqsRev) {
            querySeqsRevStore.emplace_back(FastaSequenceCachedStore(chunkRevQueries));
        }
#ifdef PANCAKE_TIMINGS
        timer.Freeze();
        seedAlnCpuTime += timer.ElapsedNanoseconds();
        PBLOG_INFO << "CPU RevComp            : " << timer.ElapsedTime();
        timer.Reset();
#endif

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
#ifdef PANCAKE_TIMINGS
        timer.Freeze();
        cpuTime += timer.ElapsedNanoseconds();
        PBLOG_INFO << "CPU Prepare            : " << timer.ElapsedTime();
        timer.Reset();
#endif

        // Global alignment on FPGA. Try using different bandwidths,
        // increasing the bandwidth for the failed parts each iteration.
        std::vector<AlignmentResult> internalAlns;
        int32_t numInternalNotValid = 0;
        numInternalNotValid = AlignPartsOnFPGA_(aligner, partsGlobal, internalAlns, seedAlnCpuTime,
                                                seedAlnfpgaTime, false);

#ifdef PANCAKE_TIMINGS
        cpuTime += seedAlnCpuTime;
        fpgaTime += seedAlnfpgaTime;
        PBLOG_INFO << "CPU Internal Prepare   : "
                   << PacBio::Utility::Stopwatch::PrettyPrintNanoseconds(seedAlnCpuTime);
        PBLOG_INFO << "FPGA Internal Alignment: "
                   << PacBio::Utility::Stopwatch::PrettyPrintNanoseconds(seedAlnfpgaTime);
        timer.Reset();
#endif

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
        }
        PBLOG_TRACE << "internalAlns.size() = " << internalAlns.size();
#ifdef PANCAKE_TIMINGS
        timer.Freeze();
        cpuTime += timer.ElapsedNanoseconds();
        PBLOG_INFO << "CPU Fallback           : " << timer.ElapsedTime() << ' '
                   << numInternalNotValid;
        timer.Reset();
#endif

        // Flank alignment on FPGA.
        std::vector<AlignmentResult> flankAlns;
        seedAlnCpuTime = 0;
        seedAlnfpgaTime = 0;
        // const int32_t numNotValidFlanks = AlignPartsOnFPGA_(aligner, partsSemiglobal, flankAlns,
        //                                                     seedAlnCpuTime, seedAlnfpgaTime, true);
        const int32_t numNotValidFlanks =
            AlignPartsOnCpu(alignSettings.alignerTypeGlobal, alignSettings.alnParamsGlobal,
                            alignSettings.alignerTypeExt, alignSettings.alnParamsExt,
                            partsSemiglobal, faf, flankAlns, prepareTime, alignTime);
        PBLOG_TRACE << "Total not valid: " << numNotValidFlanks << " / " << flankAlns.size()
                    << "\n";
#ifdef PANCAKE_TIMINGS
        cpuTime += seedAlnCpuTime;
        fpgaTime += seedAlnfpgaTime;
        PBLOG_INFO << "CPU Flanks Prepare     : "
                   << PacBio::Utility::Stopwatch::PrettyPrintNanoseconds(seedAlnCpuTime);
        PBLOG_INFO << "FPGA Flanks Alignment  : "
                   << PacBio::Utility::Stopwatch::PrettyPrintNanoseconds(seedAlnfpgaTime);
        timer.Reset();
#endif

        StitchAlignmentsInParallel(results, batchChunks, querySeqsRevStore, internalAlns, flankAlns,
                                   alnStitchInfo, faf);
#ifdef PANCAKE_TIMINGS
        timer.Freeze();
        cpuTime += timer.ElapsedNanoseconds();
        PBLOG_INFO << "CPU Stitch             : " << timer.ElapsedTime();
        timer.Reset();
#endif

        SetUnalignedAndMockedMappings(
            results, alignSettings.selfHitPolicy == MapperSelfHitPolicy::PERFECT_ALIGNMENT,
            alignSettings.alnParamsGlobal.matchScore);
#ifdef PANCAKE_TIMINGS
        timer.Freeze();
        cpuTime += timer.ElapsedNanoseconds();
        PBLOG_INFO << "CPU Mock               : " << timer.ElapsedTime();
        timer.Reset();
#endif

        UpdateSecondaryAndFilter(results, faf, batchChunks);
#ifdef PANCAKE_TIMINGS
        timer.Freeze();
        cpuTime += timer.ElapsedNanoseconds();
        PBLOG_INFO << "CPU Update             : " << timer.ElapsedTime();
        timer.Reset();
#endif
    }
#ifdef PANCAKE_TIMINGS
    PBLOG_INFO << "CPU Time               : "
               << PacBio::Utility::Stopwatch::PrettyPrintNanoseconds(cpuTime);
    PBLOG_INFO << "FPGA Time              : "
               << PacBio::Utility::Stopwatch::PrettyPrintNanoseconds(fpgaTime);
#endif

    return results;
}

void MapperBatchFPGA::WorkerMapper_(const std::vector<MapperBatchChunk>& batchChunks,
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

int32_t MapperBatchFPGA::AlignPartsOnFPGA_(AlignerBatchFPGA& aligner,
                                           const std::vector<PairForBatchAlignment>& parts,
                                           std::vector<AlignmentResult>& retInternalAlns,
                                           int64_t& cpuTime, int64_t& fpgaTime,
                                           const bool extension)
{
#ifdef PANCAKE_TIMINGS
    PacBio::Utility::Stopwatch timer;
#endif
    retInternalAlns.resize(parts.size());
#ifdef PANCAKE_TIMINGS
    cpuTime += timer.ElapsedNanoseconds();
#endif

    int32_t totalNumNotValid = 0;
    size_t partId = 0;
    while (partId < parts.size()) {
#ifdef PANCAKE_TIMINGS
        timer.Reset();
#endif
        aligner.Clear();

        std::vector<size_t> partIds;
        partIds.reserve(parts.size());

        PBLOG_TRACE << "Preparing sequences for FPGA alignment.";
        for (; partId < parts.size(); ++partId) {
            const PairForBatchAlignment& part = parts[partId];

            if (retInternalAlns[partId].valid) {
                continue;
            }
            partIds.emplace_back(partId);

            if (part.regionType == RegionType::FRONT) {
                // Reverse the sequences for front flank alignment. No need to complement.
                std::string query(part.query, part.queryLen);
                std::reverse(query.begin(), query.end());
                std::string target(part.target, part.targetLen);
                std::reverse(target.begin(), target.end());
                aligner.AddSequencePair(query.c_str(), part.queryLen, target.c_str(),
                                        part.targetLen);
            } else {
                aligner.AddSequencePair(part.query, part.queryLen, part.target, part.targetLen);
            }
        }
#ifdef PANCAKE_TIMINGS
        cpuTime += timer.ElapsedNanoseconds();
#endif

        PBLOG_TRACE << "Aligning batch of " << aligner.BatchSize() << " sequence pairs.";
        std::pair<int64_t, int64_t> seedTimings = aligner.AlignAll(extension);
#ifdef PANCAKE_TIMINGS
        cpuTime += seedTimings.first;
        fpgaTime += seedTimings.second;
        timer.Reset();
#endif

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
#ifdef PANCAKE_TIMINGS
        cpuTime += timer.ElapsedNanoseconds();
#endif
    }
    PBLOG_TRACE << "Total not valid: " << totalNumNotValid << " / " << retInternalAlns.size()
                << "\n";
    return totalNumNotValid;
}

}  // namespace Pancake
}  // namespace PacBio
