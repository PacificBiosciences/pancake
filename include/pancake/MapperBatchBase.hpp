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
};
}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_MAPPER_BATCH_BASE_HPP
