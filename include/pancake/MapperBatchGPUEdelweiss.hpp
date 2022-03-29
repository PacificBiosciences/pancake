// Author: Ivan Sovic

#ifndef PANCAKE_MAPPER_BATCH_GPU_EDELWEISS_HPP
#define PANCAKE_MAPPER_BATCH_GPU_EDELWEISS_HPP

#include <pancake/AlignerBatchGPUEdelweiss.hpp>
#include <pancake/FastaSequenceCached.hpp>
#include <pancake/MapperBase.hpp>
#include <pancake/MapperBatchBase.hpp>
#include <pancake/MapperBatchCPU.hpp>
#include <pancake/MapperBatchUtility.hpp>
#include <pancake/MapperCLR.hpp>

#include <pbcopper/parallel/FireAndForget.h>

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace PacBio {
namespace Pancake {

class MapperBatchGPUEdelweiss : public MapperBatchBase
{
public:
    MapperBatchGPUEdelweiss(const MapperCLRAlignSettings& alignSettings,
                            Parallel::FireAndForget* faf, int32_t gpuStartBandwidth,
                            int32_t gpuMaxBandwidth, int32_t gpuDeviceId, int64_t gpuMemoryBytes,
                            int32_t maxAllowedGapForGpu, bool alignRemainingOnCpu);
    MapperBatchGPUEdelweiss(const MapperCLRAlignSettings& alignSettings, int32_t numThreads,
                            int32_t gpuStartBandwidth, int32_t gpuMaxBandwidth, int32_t gpuDeviceId,
                            int64_t gpuMemoryBytes, int32_t maxAllowedGapForGpu,
                            bool alignRemainingOnCpu);
    ~MapperBatchGPUEdelweiss() override;

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
    std::unique_ptr<AlignerBatchGPUEdelweiss> aligner_;

    static std::vector<std::vector<MapperBaseResult>> MapAndAlignImpl_(
        const std::vector<MapperBatchChunk>& batchChunks,
        const MapperCLRAlignSettings& alignSettings, int32_t maxAllowedGapForGpu,
        bool alignRemainingOnCpu, int32_t gpuStartBandwidth, int32_t gpuMaxBandwidth,
        AlignerBatchGPUEdelweiss& aligner, Parallel::FireAndForget* faf);

    static void WorkerMapper_(const std::vector<MapperBatchChunk>& batchChunks, int32_t startId,
                              int32_t endId, std::vector<std::vector<MapperBaseResult>>& results);
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_MAPPER_BATCH_GPU_EDELWEISS_HPP
