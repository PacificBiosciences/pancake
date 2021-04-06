#ifndef PANCAKE_MAPPER_BATCH_FPGA_H
#define PANCAKE_MAPPER_BATCH_FPGA_H

#include <pacbio/pancake/AlignerBatchFPGA.h>
#include <pacbio/pancake/FastaSequenceCached.h>
#include <pacbio/pancake/MapperBase.h>
#include <pacbio/pancake/MapperBatchBase.h>
#include <pacbio/pancake/MapperBatchCPU.h>
#include <pacbio/pancake/MapperBatchUtility.h>
#include <pacbio/pancake/MapperCLR.h>
#include <pbcopper/parallel/FireAndForget.h>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace PacBio {
namespace Pancake {

class MapperBatchFPGA : public MapperBatchBase
{
public:
    MapperBatchFPGA(const MapperCLRAlignSettings& alignSettings, Parallel::FireAndForget* faf,
                    uint32_t fpgaDeviceId, bool alignRemainingOnCpu);
    MapperBatchFPGA(const MapperCLRAlignSettings& alignSettings, int32_t numThreads,
                    uint32_t fpgaDeviceId, bool alignRemainingOnCpu);
    ~MapperBatchFPGA() override;

    std::vector<std::vector<MapperBaseResult>> MapAndAlign(
        const std::vector<MapperBatchChunk>& batchData) override;

private:
    MapperCLRAlignSettings alignSettings_;
    int32_t numThreads_;
    bool alignRemainingOnCpu_;
    Parallel::FireAndForget* faf_;
    std::unique_ptr<Parallel::FireAndForget> fafFallback_;
    std::unique_ptr<AlignerBatchFPGA> aligner_;

    static std::vector<std::vector<MapperBaseResult>> MapAndAlignImpl_(
        const std::vector<MapperBatchChunk>& batchChunks,
        const MapperCLRAlignSettings& alignSettings, bool alignRemainingOnCpu,
        AlignerBatchFPGA& aligner, Parallel::FireAndForget* faf);

    static void WorkerMapper_(const std::vector<MapperBatchChunk>& batchChunks, int32_t startId,
                              int32_t endId, std::vector<std::vector<MapperBaseResult>>& results);

    static int32_t AlignPartsOnFPGA_(AlignerBatchFPGA& aligner,
                                     const std::vector<PairForBatchAlignment>& parts,
                                     std::vector<AlignmentResult>& retInternalAlns,
                                     int64_t& seedAlnCpuTime, int64_t& seedAlnFpgaTime,
                                     const bool extension);
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_MAPPER_BATCH_FPGA_H
