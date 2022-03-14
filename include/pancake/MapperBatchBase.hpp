#ifndef PANCAKE_ABSTRACT_MAPPER_BATCH_H
#define PANCAKE_ABSTRACT_MAPPER_BATCH_H

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

#endif  // PANCAKE_ABSTRACT_MAPPER_BATCH_H
