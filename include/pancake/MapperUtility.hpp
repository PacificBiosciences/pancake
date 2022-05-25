// Author: Ivan Sovic

#ifndef PANCAKE_MAPPER_UTILITY_HPP
#define PANCAKE_MAPPER_UTILITY_HPP

#include <pancake/FastaSequenceCachedStore.hpp>
#include <pancake/Minimizers.hpp>
#include <pancake/Overlap.hpp>
#include <pancake/SeedHit.hpp>
#include <pancake/SeedIndex.hpp>

#include <cstdint>
#include <vector>

namespace PacBio {
namespace Pancake {

/*
* \brief Utility function which constructs an overlap from a given chain of seed hits.
* Overlap coordinates are determined based on the bounding box around the seed hits.
*/
OverlapPtr MakeOverlap(const std::vector<SeedHit>& sortedHits, int32_t queryId, int32_t queryLen,
                       const FastaSequenceCachedStore& targetSeqs, int32_t beginId, int32_t endId,
                       int32_t minTargetPosId, int32_t maxTargetPosId);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_MAPPER_UTILITY_HPP
