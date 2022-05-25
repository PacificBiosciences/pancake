// Authors: Ivan Sovic

#include <pancake/MapperUtility.hpp>

#include <pbcopper/logging/Logging.h>

#include <iostream>
#include <sstream>

namespace PacBio {
namespace Pancake {

OverlapPtr MakeOverlap(const std::vector<SeedHit>& sortedHits, int32_t queryId, int32_t queryLen,
                       const FastaSequenceCachedStore& targetSeqs, int32_t beginId, int32_t endId,
                       int32_t minTargetPosId, int32_t maxTargetPosId)
{

    const auto& beginHit = sortedHits[minTargetPosId];
    const auto& endHit = sortedHits[maxTargetPosId];

    const int32_t targetId = beginHit.targetId;
    const int32_t numSeeds = endId - beginId;

    if (endHit.targetId != beginHit.targetId) {
        std::ostringstream oss;
        oss << "The targetId of the first and last seed does not match, in MakeOverlap. "
               "beginHit.targetId "
            << beginHit.targetId << ", endHit.targetId = " << endHit.targetId;
        throw std::runtime_error(oss.str());
    }

    const float score = numSeeds;
    const float identity = 0.0;
    const int32_t editDist = -1;

    const int32_t targetLen = targetSeqs.GetSequence(targetId).size();

    OverlapPtr ret =
        CreateOverlap(queryId, targetId, score, identity, beginHit.targetRev, beginHit.queryPos,
                      endHit.queryPos + endHit.querySpan, queryLen, false, beginHit.targetPos,
                      endHit.targetPos + endHit.targetSpan, targetLen, editDist, numSeeds,
                      OverlapType::Unknown, OverlapType::Unknown);

    ret->NormalizeStrand();

    return ret;
}

}  // namespace Pancake
}  // namespace PacBio
