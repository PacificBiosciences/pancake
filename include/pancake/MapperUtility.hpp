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

/**
 * \brief Given several types of occurrence thresholds, computes the final occurrence threshold that can be used
 *          for collecting seed hits. Uses the following formula to compute the return value:
 *              cutoff = max(seedOccurrenceMin, min(seedOccurrenceMax, seedOccurrenceMemMax, percentileCutoff))
 * \param seedHitHistogram Histogram of seed hit occurrences, given as a vector of pairs: <seed occurrence, number of query seeds with this occurrence>.
 *                          Can be an empty vector. Used only when seedOccurrenceMaxMemory > 0. The vector should be sorted in the ascending order of seed occurrence.
 * \param seedOccurrenceMin Minimum value for the occurrence threshold. If the other cutoff values result in a smaller value than the minimum, then the return value is pinned to this.
 * \param seedOccurrenceMax Maximum allowed occurrence of a seed to keep for mapping. Value <= 0 turns off this threshold.
 * \param seedOccurrenceMaxMemory Maximum allowed memory to be consumed by collected seed hits. This is used to dynamically compute the maximum
 *                                  occurrence cutoff based on the seed hit histogram. Seeds are chosen in the sorted order by their occurrence
 *                                  until the memory threshold is reached. Value <= 0 turns off this threshold.
 * \param seedOccurrenceUserSpecified User-provided cutoff threshold, computed, for example, from the frequency percentile threshold
 *                                          (e.g. top 0.002% of most abundant seeds in the index should be skipped.).
 * \return Single value that consolidates the occurrence cutoff. If the cutoff is not applied (e.g. all or some provided values are zero), it returns std::numeric_limits<int64_t>::max().
*/
int64_t ComputeOccurrenceThreshold(const std::vector<std::pair<int64_t, int64_t>>& seedHitHistogram,
                                   int64_t seedOccurrenceMin, int64_t seedOccurrenceMax,
                                   int64_t seedOccurrenceMaxMemory,
                                   int64_t seedOccurrenceUserSpecified, bool debugVerbose);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_MAPPER_UTILITY_HPP
