// Author: Ivan Sovic

#ifndef PANCAKE_MAPPER_BATCH_CPU_HPP
#define PANCAKE_MAPPER_BATCH_CPU_HPP

#include <pancake/AlignerBatchCPU.hpp>
#include <pancake/FastaSequenceCached.hpp>
#include <pancake/MapperBase.hpp>
#include <pancake/MapperBatchBase.hpp>
#include <pancake/MapperBatchUtility.hpp>
#include <pancake/MapperCLR.hpp>

#include <pbcopper/parallel/FireAndForget.h>

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace PacBio {
namespace Pancake {

class MapperBatchCPU : public MapperBatchBase
{
public:
    MapperBatchCPU(const MapperCLRAlignSettings& alignSettings, int32_t numThreads);
    MapperBatchCPU(const MapperCLRAlignSettings& alignSettings, Parallel::FireAndForget* faf);

    ~MapperBatchCPU();

    std::vector<std::vector<MapperBaseResult>> MapAndAlign(
        const std::vector<MapperBatchChunk>& batchData) override;

private:
    MapperCLRAlignSettings alignSettings_;
    Parallel::FireAndForget* faf_;
    std::unique_ptr<Parallel::FireAndForget> fafFallback_;

    static std::vector<std::vector<MapperBaseResult>> MapAndAlignImpl_(
        const std::vector<MapperBatchChunk>& batchChunks,
        const MapperCLRAlignSettings& alignSettings, Parallel::FireAndForget* faf);

    static void WorkerMapper_(const std::vector<MapperBatchChunk>& batchChunks, int32_t startId,
                              int32_t endId, std::vector<std::vector<MapperBaseResult>>& results);
};

void UpdateSecondaryAndFilter(std::vector<std::vector<MapperBaseResult>>& mappingResults,
                              Parallel::FireAndForget* faf,
                              const std::vector<MapperBatchChunk>& batchChunks);

int32_t AlignPartsOnCpu(const AlignerType& alignerTypeGlobal,
                        const AlignmentParameters& alnParamsGlobal,
                        const AlignerType& alignerTypeExt, const AlignmentParameters& alnParamsExt,
                        const std::vector<PairForBatchAlignment>& parts, const int32_t numThreads,
                        std::vector<AlignmentResult>& retAlns);

int32_t AlignPartsOnCpu(const AlignerType& alignerTypeGlobal,
                        const AlignmentParameters& alnParamsGlobal,
                        const AlignerType& alignerTypeExt, const AlignmentParameters& alnParamsExt,
                        const std::vector<PairForBatchAlignment>& parts,
                        Parallel::FireAndForget* faf, std::vector<AlignmentResult>& retAlns);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_MAPPER_BATCH_CPU_HPP
