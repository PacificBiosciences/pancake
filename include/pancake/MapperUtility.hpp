// Author: Ivan Sovic

#ifndef PANCAKE_MAPPER_UTILITY_HPP
#define PANCAKE_MAPPER_UTILITY_HPP

#include <pancake/AlignmentRegion.hpp>
#include <pancake/FastaSequenceCachedStore.hpp>
#include <pancake/MapperBase.hpp>
#include <pancake/Minimizers.hpp>
#include <pancake/Overlap.hpp>
#include <pancake/SeedHit.hpp>
#include <pancake/SeedIndex.hpp>

#include <pancake/third-party/pdqsort/pdqsort.h>

#include <cstdint>
#include <vector>

namespace PacBio {
namespace Pancake {

/**
 * \brief Utility function which constructs an overlap from a given chain of seed hits.
 * Overlap coordinates are determined based on the bounding box around the seed hits.
*/
OverlapPtr MakeOverlap(const std::vector<SeedHit>& sortedHits, int32_t queryId, int32_t queryLen,
                       const FastaSequenceCachedStore& targetSeqs, int32_t beginId, int32_t endId,
                       int32_t minTargetPosId, int32_t maxTargetPosId);

/**
 * \brief Removes the mappigns marked for removal (priority > 1) and keeps at most
 *          at most bestNSecondary secondary mappings (priority == 1).
*/
int32_t CondenseMappings(std::vector<std::unique_ptr<ChainedRegion>>& mappings,
                         int32_t bestNSecondary);

/**
 * \brief A helper function which creates a mocked self-mapping based on the query and target
 * IDs and lengths. The result is an unaligned mapping which spans full length of the
 * query and target.
 * Throws if the query and target lengths are different.
 * This function also does not initialize the Atype and Btype labels (it sets them to Unknown).
*/
std::unique_ptr<ChainedRegion> CreateMockedMapping(int32_t queryId, int32_t queryLen,
                                                   int32_t targetId, int32_t targetLen);

/**
 * @brief For a given input mapping of a query (A-read) to a target (B-read), create a
 *          mocked perfect alignment, but only if the lengths of A and B reads are the same.
 *          Otherwise, it throws.
 *
 * @param ovl Input mapping (overlap). Used only for Aid, Bid, Alen and Blen.
 * @param matchScore Alignment match score used to compute a mocked alignment score of the perfect alignment.
 * @return OverlapPtr New overlap containing the new perfect alignment.
*/
OverlapPtr CreateMockedAlignment(const OverlapPtr& ovl, int32_t matchScore);

/**
 * \brief Wraps the labelling of secondary and supplementary chained regions.
*/
void WrapFlagSecondaryAndSupplementary(
    std::vector<std::unique_ptr<ChainedRegion>>& allChainedRegions,
    double secondaryAllowedOverlapFractionQuery, double secondaryAllowedOverlapFractionTarget,
    double secondaryMinScoreFraction);

/**
 * \brief Filters symmetric and self seed hits, based on sequence IDs.
 * Be careful when using this function in case when the query and target DBs are not the same.
*/
std::vector<SeedHit> FilterSymmetricAndSelfHits(const std::vector<SeedHit>& hits, int32_t queryId,
                                                bool skipSelfHits, bool skipSymmetricOverlaps);

/**
 * \brief Merges the neighboring chains if they are not overlaping in neither the query nor
 * target coordinates.
 * The right of the two merged chained regions is set to nullptr.
*/
void LongMergeChains(std::vector<std::unique_ptr<ChainedRegion>>& chainedRegions,
                     int32_t maxBandwidth);

/**
 * @brief For a given mapped (but unaligned) segment of a query, this function extracts regions in between seed
 *          hits that will be used for alignment. Three types of regions are created: front, global (internal) and back.
 *          In its core, this function is a wrapper around ExtractAlignmentRegions but performs additional
 *          sanity checks and accepts a ChainedRegion as the input instead of a vector of seed hits.
 *
 * @param singleMapping Input mapped region.
 * @param minAlignmentSpan Minimum span of an alignment region, if possible. If seed hits are dense, then the alignment region
 *                          will cover multiple seed hits to achieve this span.
 * @param maxFlankExtensionDist Maximum distance for flank extension alignment, hard cap in case the remaining query flank
 *                              is very long, to limit time and memory.
 * @param flankExtensionFactor Factor to select the maximum distance for extension alignment in the sequence which has the longer flank.
 *                              Example: query has 1000 bp left of unmapped flank end, but target is a long reference of 100Mbp. This will
 *                              limit the maximum size of the target sequence for alignment to flankExtensionFactor * 1000bp.
 * @return std::vector<AlignmentRegion> Vector of extracted alignment regions.
*/
std::vector<AlignmentRegion> ExtractAlignmentRegions(const ChainedRegion& singleMapping,
                                                     int32_t minAlignmentSpan,
                                                     int32_t maxFlankExtensionDist,
                                                     double flankExtensionFactor);

/**
 * \brief Takes a chain of sorted seed hits and extracts regions ready for alignment.
 * The regions are extracted in three ways:
 *      - Front of the sequences (front flank) for semiglobal (extension) alignment.
 *      - Internal portions for global alignment. Each of those portions will be at least `minAlignmentSpan` in
 *          both query/target spans, unless the total span of the chain is not as long.
 *      - Back of the sequences (back flank) for semiglobal (extension) alignment.
 * Results of this function can then be input to AlignRegionsGeneric to produce alignments of each
 * of the regions.
 *
 * Inputs:
 * \param inSortedHits Chain of seed hits which should be monotonically increasing
 *          in both query and target coordinates, sorted by coordinates. The chain needs to have all
 *          hits on the same strand and on the same target ID. All of the seed hits need to
 *          match the strand specified via the `isRev` parameter.
 * \param qLen Total query length.
 * \param tLen Total target length.
 * \param isRev The alignment is on the reverse strand of the query sequence.
 * \param minAlignmentSpan Internal (global) alignment portions will be at least this wide in both
 *                          target and query coordinates, provided that the chain spans this far.
 * \param maxFlankExtensionDist Limit the maximum flank span for alignment to this value. If < 0, no limit is set.
 * \param flankExtensionFactor Heuristic factor to extend the potential target flank span for alignment (to allow for insertions).
 *                              Multiplies the query flank span to compute the potential target flank span. Value should be >= 1.0.
 * \returns Vector of alignment regions.
*/
std::vector<AlignmentRegion> ExtractAlignmentRegions(const std::vector<SeedHit>& inSortedHits,
                                                     int32_t qLen, int32_t tLen, bool isRev,
                                                     int32_t minAlignmentSpan,
                                                     int32_t maxFlankExtensionDist,
                                                     double flankExtensionFactor);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_MAPPER_UTILITY_HPP
