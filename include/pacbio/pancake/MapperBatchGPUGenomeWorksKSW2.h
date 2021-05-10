// Author: Ivan Sovic

#ifndef PANCAKE_MAPPER_BATCH_GPU_GENOMEWORKS_KSW2_H
#define PANCAKE_MAPPER_BATCH_GPU_GENOMEWORKS_KSW2_H

#include <pacbio/pancake/AlignerBatchGPUGenomeWorksKSW2.h>
#include <pacbio/pancake/FastaSequenceCached.h>
#include <pacbio/pancake/MapperBase.h>
#include <pacbio/pancake/MapperBatchBase.h>
#include <pacbio/pancake/MapperBatchCPU.h>
#include <pacbio/pancake/MapperBatchUtility.h>
#include <pacbio/pancake/MapperCLR.h>
#include <pbcopper/parallel/FireAndForget.h>
#include <claraparabricks/genomeworks/utils/cudautils.hpp>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace PacBio {
namespace Pancake {

class MapperBatchGPUGenomeWorksKSW2 : public MapperBatchBase
{
public:
    MapperBatchGPUGenomeWorksKSW2(const MapperCLRAlignSettings& alignSettings,
                                  Parallel::FireAndForget* faf, int32_t gpuStartBandwidth,
                                  int32_t gpuMaxBandwidth, uint32_t gpuDeviceId,
                                  int64_t gpuMemoryBytes, bool alignRemainingOnCpu);
    MapperBatchGPUGenomeWorksKSW2(const MapperCLRAlignSettings& alignSettings, int32_t numThreads,
                                  int32_t gpuStartBandwidth, int32_t gpuMaxBandwidth,
                                  uint32_t gpuDeviceId, int64_t gpuMemoryBytes,
                                  bool alignRemainingOnCpu);
    ~MapperBatchGPUGenomeWorksKSW2() override;

    std::vector<std::vector<MapperBaseResult>> MapAndAlign(
        const std::vector<MapperBatchChunk>& batchData) override;

private:
    MapperCLRAlignSettings alignSettings_;
    int32_t numThreads_;
    int32_t gpuStartBandwidth_;
    int32_t gpuMaxBandwidth_;
    bool alignRemainingOnCpu_;
    Parallel::FireAndForget* faf_;
    std::unique_ptr<Parallel::FireAndForget> fafFallback_;
    std::unique_ptr<AlignerBatchGPUGenomeWorksKSW2> aligner_;

    static std::vector<std::vector<MapperBaseResult>> MapAndAlignImpl_(
        const std::vector<MapperBatchChunk>& batchChunks,
        const MapperCLRAlignSettings& alignSettings, bool alignRemainingOnCpu,
        int32_t gpuStartBandwidth, int32_t gpuMaxBandwidth, AlignerBatchGPUGenomeWorksKSW2& aligner,
        Parallel::FireAndForget* faf);

    static void WorkerMapper_(const std::vector<MapperBatchChunk>& batchChunks, int32_t startId,
                              int32_t endId, std::vector<std::vector<MapperBaseResult>>& results);

    static int32_t AlignPartsOnGPU_(AlignerBatchGPUGenomeWorksKSW2& aligner,
                                    const std::vector<PairForBatchAlignment>& parts,
                                    std::vector<AlignmentResult>& retInternalAlns,
                                    int64_t& seedAlnCpuTime, int64_t& seedAlnGpuTime);
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_MAPPER_BATCH_GPU_GENOMEWORKS_KSW2_H
