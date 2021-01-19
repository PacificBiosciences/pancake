// Author: Ivan Sovic

#ifndef PANCAKE_MAPPER_BATCH_H
#define PANCAKE_MAPPER_BATCH_H

#include <pacbio/pancake/AlignerBatchCPU.h>
#include <pacbio/pancake/FastaSequenceCached.h>
#include <pacbio/pancake/MapperBase.h>
#include <pacbio/pancake/MapperCLR.h>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace PacBio {
namespace Pancake {

struct MapperBatchChunk
{
    std::vector<FastaSequenceCached> targetSeqs;
    std::vector<FastaSequenceCached> querySeqs;
};

enum class BatchAlignerRegionType
{
    GLOBAL,
    SEMIGLOBAL,
    BOTH,
};

class MapperBatchCPU
{
public:
    MapperBatchCPU(const MapperCLRSettings& settings, int32_t numThreads);
    ~MapperBatchCPU();

    std::vector<std::vector<MapperBaseResult>> DummyMapAndAlign(
        const std::vector<MapperBatchChunk>& batchData);

    std::vector<std::vector<MapperBaseResult>> MapAndAlign(
        const std::vector<MapperBatchChunk>& batchData);

private:
    MapperCLRSettings settings_;
    int32_t numThreads_;

    static std::vector<std::vector<MapperBaseResult>> DummyMapAndAlignImpl_(
        const std::vector<MapperBatchChunk>& batchChunks, MapperCLRSettings settings,
        int32_t numThreads);

    static std::vector<std::vector<MapperBaseResult>> MapAndAlignImpl_(
        const std::vector<MapperBatchChunk>& batchChunks, MapperCLRSettings settings,
        int32_t numThreads);

    static void PrepareSequencesForAlignment_(
        AlignerBatchCPU& aligner, const std::vector<MapperBatchChunk>& batchChunks,
        const std::vector<std::vector<MapperBaseResult>>& mappingResults,
        const BatchAlignerRegionType& regionsToAdd);

    static void StitchAlignments_(const std::vector<MapperBatchChunk>& batchChunks,
                                  const std::vector<std::vector<MapperBaseResult>>& mappingResults,
                                  const std::vector<AlignmentResult>& internalAlns,
                                  const std::vector<AlignmentResult>& flankAlns);

    static void UpdateSecondaryAndFilter_(
        std::vector<std::vector<MapperBaseResult>>& mappingResults,
        double secondaryAllowedOverlapFractionQuery, double secondaryAllowedOverlapFractionTarget,
        double secondaryMinScoreFraction, int32_t bestNSecondary);

    static void WorkerMapper_(const std::vector<MapperBatchChunk>& batchChunks, int32_t startId,
                              int32_t endId, const std::unique_ptr<MapperCLR>& mapper,
                              std::vector<std::vector<MapperBaseResult>>& results);
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_MAPPER_BATCH_H
