// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNMENT_SEEDED_HPP
#define PANCAKE_ALIGNMENT_SEEDED_HPP

#include <pancake/AlignerFactory.hpp>
#include <pancake/AlignmentRegion.hpp>
#include <pancake/DPChain.hpp>
#include <pancake/Overlap.hpp>
#include <pancake/Seed.hpp>

#include <pbcopper/data/Cigar.h>

#include <cstdint>
#include <memory>
#include <string_view>
#include <vector>

namespace PacBio {
namespace Pancake {

class AlignRegionsGenericResult
{
public:
    Data::Cigar cigar;
    int32_t offsetFrontQuery = 0;
    int32_t offsetBackQuery = 0;
    int32_t offsetFrontTarget = 0;
    int32_t offsetBackTarget = 0;
    int32_t score = 0;
};

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

AlignmentResult AlignSingleRegion(std::string_view targetSeq, std::string_view querySeqFwd,
                                  std::string_view querySeqRev, AlignerBasePtr& alignerGlobal,
                                  AlignerBasePtr& alignerExt, const AlignmentRegion& region);

AlignRegionsGenericResult AlignRegionsGeneric(std::string_view targetSeq,
                                              std::string_view querySeqFwd,
                                              std::string_view querySeqRev,
                                              const std::vector<AlignmentRegion>& regions,
                                              AlignerBasePtr& alignerGlobal,
                                              AlignerBasePtr& alignerExt);

OverlapPtr AlignmentSeeded(const OverlapPtr& ovl, const std::vector<AlignmentRegion>& alnRegions,
                           std::string_view targetSeq, std::string_view querySeqFwd,
                           std::string_view querySeqRev, AlignerBasePtr& alignerGlobal,
                           AlignerBasePtr& alignerExt);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ALIGNMENT_SEEDED_HPP
