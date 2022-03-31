// Author: Ivan Sovic

#ifndef PANCAKE_MAPPER_BATCH_GPU_UTILITY_HPP
#define PANCAKE_MAPPER_BATCH_GPU_UTILITY_HPP

#include <pancake/AlignerBatchBase.hpp>
#include <pancake/MapperBase.hpp>
#include <pancake/MapperBatchBase.hpp>

#include <vector>

namespace PacBio {
namespace Pancake {

/**
 * @brief Takes a vector of sequence pairs for alignment and performs alignment of those pairs. In case
 *          there already is a valid alignment for that pair in the retAlns vector, alignment is not performed.
 *          Also, alignment is skipped for those sequence pairs which exceed the maxAllowedGap threshold.
 *
 * @param parts A vector of sequence pairs for alignment.
 * @param aligner Aligner to use on the sequence pairs.
 * @param maxAllowedGap Maximum allowed size difference between the query/target sequences. Value -1 deactivates this filter.
 * @param retAlns Return value via parameter, a vector of alignments, one for each input part. Relation is 1:1.
 * @return int32_t Number of alignments which were not valid at the end of this function. If everything was aligned, returns 0.
 */
int32_t AlignPartsOnGPU(const std::vector<PairForBatchAlignment>& parts, AlignerBatchBase& aligner,
                        int32_t maxAllowedGap, std::vector<AlignmentResult>& retAlns);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_MAPPER_BATCH_GPU_UTILITY_HPP
