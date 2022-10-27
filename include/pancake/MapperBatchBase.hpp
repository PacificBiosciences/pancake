#ifndef PANCAKE_MAPPER_BATCH_BASE_HPP
#define PANCAKE_MAPPER_BATCH_BASE_HPP

#include <pancake/MapperBase.hpp>
#include <pancake/MapperBatchUtility.hpp>

#include <vector>

namespace PacBio {
namespace Pancake {

class MapperBatchBase
{
public:
    virtual ~MapperBatchBase() = default;

    virtual std::vector<std::vector<MapperBaseResult>> MapAndAlign(
        const std::vector<MapperBatchChunk>& batchData) = 0;

protected:
    /**
     * @brief This is a worker function used to compute mappings in parallel, and meant to be used by a
     *          worker queue. It operates on a subset of batchChunks and requires that the resuls vector
     *          is preallocated in advance (to the length of batchChunks), so that the results can be stored.
     *
     * @param batchChunks Input batch chunks.
     * @param startId ID of the first chunk to process.
     * @param endId ID of the end chunk for processing.
     * @param results Return vector, needs to be preallocated before this function is called, otherwise no work will be performed.
     */
    virtual void WorkerMapper_(const std::vector<MapperBatchChunk>& batchChunks, int32_t startId,
                               int32_t endId,
                               std::vector<std::vector<MapperBaseResult>>& results) const;
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_MAPPER_BATCH_BASE_HPP
