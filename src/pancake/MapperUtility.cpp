// Authors: Ivan Sovic

#include <pacbio/pancake/MapperUtility.h>
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
        createOverlap(queryId, targetId, score, identity, beginHit.targetRev, beginHit.queryPos,
                      endHit.queryPos + endHit.querySpan, queryLen, false, beginHit.targetPos,
                      endHit.targetPos + endHit.targetSpan, targetLen, editDist, numSeeds,
                      OverlapType::Unknown, OverlapType::Unknown);

    ret->NormalizeStrand();

    return ret;
}

int64_t ComputeOccurrenceThreshold(const std::vector<std::pair<int64_t, int64_t>>& seedHitHistogram,
                                   const int64_t seedOccurrenceMin, const int64_t seedOccurrenceMax,
                                   const int64_t seedOccurrenceMaxMemory,
                                   const int64_t seedOccurrenceUserSpecified,
                                   const bool debugVerbose)
{
    int64_t occThresholdMemMax = std::numeric_limits<int64_t>::max();
    if (seedOccurrenceMaxMemory > 0) {
        const int64_t maxHitsToFit = std::ceil(static_cast<double>(seedOccurrenceMaxMemory) /
                                               static_cast<double>(sizeof(SeedHit)));
        int64_t totalHits = 0;
        for (size_t i = 0; i < seedHitHistogram.size(); ++i) {
            const int64_t currentBinSize = seedHitHistogram[i].first * seedHitHistogram[i].second;
            // std::cerr << "[i = " << i << " / " << seedHitHistogram.size() << "] hits = " << seedHitHistogram[i].first << ", numSeeds = "
            // << seedHitHistogram[i].second << ", currentBinSize = " << currentBinSize << ", totalHits = " << totalHits
            // << ", maxHitsToFit = " << maxHitsToFit << "\n";
            if ((totalHits + currentBinSize) > maxHitsToFit) {
                break;
            }
            totalHits += currentBinSize;
            occThresholdMemMax = seedHitHistogram[i].first + 1;
        }
    }

    const int64_t occThresholdMax =
        seedOccurrenceMax > 0 ? seedOccurrenceMax : std::numeric_limits<int64_t>::max();
    const int64_t occPercentileCutoff = seedOccurrenceUserSpecified > 0
                                            ? seedOccurrenceUserSpecified
                                            : std::numeric_limits<int64_t>::max();

    const int64_t occThreshold =
        std::max(seedOccurrenceMin,
                 std::min(std::min(occThresholdMax, occThresholdMemMax), occPercentileCutoff));

    if (debugVerbose) {
        std::cerr << "Seed hit occurrence threshold: occThreshold = " << occThreshold
                  << " (occUserSpec = " << seedOccurrenceUserSpecified
                  << ", occMin = " << seedOccurrenceMin << ", occMax = " << seedOccurrenceMax
                  << ", occMemMax = " << occThresholdMemMax
                  << ", maxMemoryInBytes = " << seedOccurrenceMaxMemory << ")"
                  << "\n";
    }

    return occThreshold;
}

}  // namespace Pancake
}  // namespace PacBio
