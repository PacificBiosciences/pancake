#ifndef PANCAKE_ABSTRACT_MAPPER_BATCH_H
#define PANCAKE_ABSTRACT_MAPPER_BATCH_H

#include <pacbio/pancake/MapperBase.h>
#include <pacbio/pancake/MapperBatchUtility.h>
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
