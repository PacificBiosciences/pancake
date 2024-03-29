// Author: Ivan Sovic

#ifndef PANCAKE_MAPPER_BATCH_GPU_HPP
#define PANCAKE_MAPPER_BATCH_GPU_HPP

#include <pancake/AlignerBatchGPU.hpp>
#include <pancake/FastaSequenceCached.hpp>
#include <pancake/MapperBase.hpp>
#include <pancake/MapperBatchBase.hpp>
#include <pancake/MapperBatchCPU.hpp>
#include <pancake/MapperBatchUtility.hpp>
#include <pancake/MapperCLR.hpp>

#include <pbcopper/parallel/FireAndForget.h>
#include <claraparabricks/genomeworks/utils/cudautils.hpp>

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace PacBio {
namespace Pancake {

class MapperBatchGPU : public MapperBatchBase
{
public:
    MapperBatchGPU(const MapperCLRAlignSettings& alignSettings, Parallel::FireAndForget* faf,
                   int32_t gpuStartBandwidth, int32_t gpuMaxBandwidth, uint32_t gpuDeviceId,
                   int64_t gpuMemoryBytes, int32_t maxAllowedGapForGpu, bool alignRemainingOnCpu);
    MapperBatchGPU(const MapperCLRAlignSettings& alignSettings, int32_t numThreads,
                   int32_t gpuStartBandwidth, int32_t gpuMaxBandwidth, uint32_t gpuDeviceId,
                   int64_t gpuMemoryBytes, int32_t maxAllowedGapForGpu, bool alignRemainingOnCpu);
    ~MapperBatchGPU() override;

    std::vector<std::vector<MapperBaseResult>> MapAndAlign(
        const std::vector<MapperBatchChunk>& batchData) override;

private:
    MapperCLRAlignSettings alignSettings_;
    int32_t numThreads_;
    int32_t gpuStartBandwidth_;
    int32_t gpuMaxBandwidth_;
    int32_t maxAllowedGapForGpu_;
    bool alignRemainingOnCpu_;
    Parallel::FireAndForget* faf_;
    std::unique_ptr<Parallel::FireAndForget> fafFallback_;
    std::unique_ptr<AlignerBatchGPU> aligner_;

    std::vector<std::vector<MapperBaseResult>> MapAndAlignImpl_(
        const std::vector<MapperBatchChunk>& batchChunks,
        const MapperCLRAlignSettings& alignSettings, int32_t maxAllowedGapForGpu,
        bool alignRemainingOnCpu, int32_t gpuStartBandwidth, int32_t gpuMaxBandwidth,
        AlignerBatchGPU& aligner, Parallel::FireAndForget* faf) const;
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_MAPPER_BATCH_GPU_HPP
