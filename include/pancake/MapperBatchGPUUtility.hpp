// Author: Ivan Sovic

#ifndef PANCAKE_MAPPER_BATCH_GPU_UTILITY_HPP
#define PANCAKE_MAPPER_BATCH_GPU_UTILITY_HPP

#include <pancake/AlignerBatchBase.hpp>
#include <pancake/MapperBase.hpp>
#include <pancake/MapperBatchBase.hpp>

#include <vector>

namespace PacBio {
namespace Pancake {

int32_t AlignPartsOnGPU(const std::vector<PairForBatchAlignment>& parts, AlignerBatchBase& aligner,
                        std::vector<AlignmentResult>& retAlns);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_MAPPER_BATCH_GPU_UTILITY_HPP
