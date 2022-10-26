// Authors: Ivan Sovic

#include <pancake/MapperUtility.hpp>

#include <pancake/AlignmentSeeded.hpp>
#include <pancake/Secondary.hpp>

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

    // TODO: Remove tech debt regarding orientation of the A and B reads, and remove the following line.
    ret->NormalizeStrand();

    return ret;
}

int32_t CondenseMappings(std::vector<std::unique_ptr<ChainedRegion>>& mappings,
                         int32_t bestNSecondary)
{
    // Filter mappings.
    size_t numValid = 0;
    int32_t numPrimary = 0;
    int32_t numSelectedSecondary = 0;
    for (size_t i = 0; i < mappings.size(); ++i) {
        if (mappings[i] == nullptr) {
            continue;
        }
        const auto& region = mappings[i];
        if (region->mapping == nullptr) {
            continue;
        }
        if (region->priority > 1) {
            continue;
        }
        if (bestNSecondary >= 0 && region->priority == 1 &&
            numSelectedSecondary >= bestNSecondary) {
            continue;
        }
        if (region->priority == 0) {
            ++numPrimary;
        }

        if (bestNSecondary < 0 ||
            (region->priority == 1 && numSelectedSecondary < bestNSecondary)) {
            ++numSelectedSecondary;
        }
        if (i != numValid) {
            std::swap(mappings[i], mappings[numValid]);
        }
        ++numValid;
    }
    mappings.resize(numValid);
    return numPrimary;
}

std::unique_ptr<ChainedRegion> CreateMockedMapping(const int32_t queryId, const int32_t queryLen,
                                                   const int32_t targetId, const int32_t targetLen)
{
    if (queryLen != targetLen) {
        throw std::runtime_error(
            "Cannot mock a perfect mapping between the two sequences with same ID, the "
            "lengths are different. The sequences might be mislabelled. queryLen = " +
            std::to_string(queryLen) + ", targetLen = " + std::to_string(targetLen));
    }

    const int32_t numSeeds = queryLen;
    const int32_t alnScore = queryLen;
    const float identity = 0.0;
    const int32_t editDist = -1;
    const bool isRev = false;

    ChainedHits newChain{targetId,
                         false,
                         {
                             SeedHit(targetId, false, 0, 0, 0, 0, 0),
                             SeedHit(targetId, false, targetLen, queryLen, 0, 0, 0),
                         },
                         alnScore,
                         queryLen,
                         targetLen};

    OverlapPtr newOvl = CreateOverlap(queryId, targetId, alnScore, identity, isRev, 0, queryLen,
                                      queryLen, false, 0, targetLen, targetLen, editDist, numSeeds,
                                      OverlapType::Unknown, OverlapType::Unknown);

    auto newChainedRegion = std::make_unique<ChainedRegion>();
    newChainedRegion->chain = std::move(newChain);
    newChainedRegion->regionsForAln = {};  // This will be generated later.
    newChainedRegion->mapping = std::move(newOvl);
    newChainedRegion->priority = 0;
    newChainedRegion->isSupplementary = false;
    newChainedRegion->isMockedMapping = true;
    newChainedRegion->isMockedAlignment = false;

    return newChainedRegion;
}

OverlapPtr CreateMockedAlignment(const OverlapPtr& ovl, const int32_t matchScore)
{

    if (ovl->Alen != ovl->Blen) {
        std::ostringstream oss;
        oss << "Cannot mock a perfect alignment between the two sequences with same ID, the "
               "lengths are different. The sequences might be mislabelled. Overlap: "
            << *ovl;
        throw std::runtime_error(oss.str());
    }

    const int32_t score = ovl->Alen * matchScore;
    const float identity = 1.0;
    const int32_t editDist = 0;
    const int32_t numSeeds = ovl->Alen;
    OverlapPtr newOvl = CreateOverlap(
        ovl->Aid, ovl->Bid, score, identity, false, 0, ovl->Alen, ovl->Alen, false, 0, ovl->Blen,
        ovl->Blen, editDist, numSeeds, OverlapType::Unknown, OverlapType::Unknown,
        // OverlapType::Contains, OverlapType::Contains,
        Data::Cigar(std::to_string(ovl->Alen) + "="), "", "", false, false, false);
    return newOvl;
}

void WrapFlagSecondaryAndSupplementary(
    std::vector<std::unique_ptr<ChainedRegion>>& allChainedRegions,
    double secondaryAllowedOverlapFractionQuery, double secondaryAllowedOverlapFractionTarget,
    double secondaryMinScoreFraction)
{
    /*
     * Edits the objects in place.
    */
    // Copy the overlaps so we can satisfy the FlagSecondaryAndSupplementary API.
    std::vector<OverlapPtr> tmpOverlaps;
    std::vector<int32_t> inputOverlapToTmpOverlap(allChainedRegions.size(), -1);
    for (size_t i = 0; i < allChainedRegions.size(); ++i) {
        if (allChainedRegions[i] == nullptr || allChainedRegions[i]->mapping == nullptr) {
            continue;
        }
        inputOverlapToTmpOverlap[i] = tmpOverlaps.size();
        tmpOverlaps.emplace_back(CreateOverlap(allChainedRegions[i]->mapping));
    }
    // Flag the secondary and supplementary overlaps.
    std::vector<OverlapPriority> overlapPriorities = FlagSecondaryAndSupplementary(
        tmpOverlaps, secondaryAllowedOverlapFractionQuery, secondaryAllowedOverlapFractionTarget,
        secondaryMinScoreFraction);

    // Set the flags.
    for (size_t i = 0; i < allChainedRegions.size(); ++i) {
        if (allChainedRegions[i] == nullptr || allChainedRegions[i]->mapping == nullptr) {
            continue;
        }
        const int32_t tmpOverlapId = inputOverlapToTmpOverlap[i];
        if (tmpOverlapId < 0) {
            continue;
        }
        auto& cr = allChainedRegions[i];
        cr->mapping->IsSecondary = (overlapPriorities[tmpOverlapId].priority > 0);
        cr->mapping->IsSupplementary = overlapPriorities[tmpOverlapId].isSupplementary;
        cr->priority = overlapPriorities[tmpOverlapId].priority;
        cr->isSupplementary = overlapPriorities[tmpOverlapId].isSupplementary;
    }
}

std::vector<SeedHit> FilterSymmetricAndSelfHits(const std::vector<SeedHit>& hits,
                                                const int32_t queryId, const bool skipSelfHits,
                                                const bool skipSymmetricOverlaps)
{
    std::vector<SeedHit> newHits(hits.size());
    size_t newHitId = 0;
    for (const SeedHit& hit : hits) {
        if ((skipSelfHits && hit.targetId == queryId) ||
            (skipSymmetricOverlaps && hit.targetId > queryId)) {
            continue;
        }
        newHits[newHitId] = hit;
        ++newHitId;
    }
    newHits.resize(newHitId);
    return newHits;
}

void LongMergeChains(std::vector<std::unique_ptr<ChainedRegion>>& chainedRegions,
                     const int32_t maxBandwidth)
{
    if (maxBandwidth < 0) {
        throw std::runtime_error(
            "(LongMergeChains_) maxBandwidth cannot be negative. maxBandwidth = " +
            std::to_string(maxBandwidth));
    }

    if (chainedRegions.empty()) {
        return;
    }

    std::vector<int32_t> candidates;
    for (size_t i = 0; i < chainedRegions.size(); ++i) {
        if (chainedRegions[i] == nullptr) {
            continue;
        }
        auto& cr = chainedRegions[i];
        if (cr->priority < 0 || cr->priority > 1) {
            continue;
        }
        candidates.emplace_back(i);
    }

    pdqsort(candidates.begin(), candidates.end(), [&](const auto& a, const auto& b) {
        const auto& ar = chainedRegions[a];
        const auto& br = chainedRegions[b];
        return std::tuple(ar->mapping->Bid, ar->mapping->Brev, ar->mapping->Astart,
                          ar->mapping->Bstart) < std::tuple(br->mapping->Bid, br->mapping->Brev,
                                                            br->mapping->Astart,
                                                            br->mapping->Bstart);
    });

    std::unordered_set<int32_t> doSort;

    int32_t lastId = candidates.at(0);
    for (int32_t i = 1; i < static_cast<int32_t>(candidates.size()); ++i) {
        const int32_t currId = candidates[i];
        auto& last = chainedRegions[lastId];
        const auto& curr = chainedRegions[currId];

        // Skip if there is an overlap (or if the target is wrong), to preserve colinearity.
        if (curr->mapping->Bid != last->mapping->Bid ||
            curr->mapping->Brev != last->mapping->Brev ||
            curr->mapping->Astart <= last->mapping->Aend ||
            curr->mapping->Bstart <= last->mapping->Bend) {
            lastId = currId;
            continue;
        }

        const int32_t gap = std::abs((curr->mapping->Bstart - last->mapping->Bend) -
                                     (curr->mapping->Astart - last->mapping->Aend));

        if (gap > maxBandwidth) {
            lastId = currId;
            continue;
        }

        last->chain.hits.insert(last->chain.hits.end(), curr->chain.hits.begin(),
                                curr->chain.hits.end());
        last->chain.score += curr->chain.score;
        last->chain.coveredBasesQuery += curr->chain.coveredBasesQuery;
        last->chain.coveredBasesTarget += curr->chain.coveredBasesTarget;

        last->mapping->Aend = curr->mapping->Aend;
        last->mapping->Bend = curr->mapping->Bend;
        last->mapping->Score += curr->mapping->Score;
        last->mapping->NumSeeds += curr->mapping->NumSeeds;

        last->isMockedMapping |= curr->isMockedMapping;
        last->isMockedAlignment |= curr->isMockedAlignment;

        chainedRegions[currId] = nullptr;

        doSort.emplace(lastId);
    }

    // Sort only the extended chains.
    for (const auto& i : doSort) {
        if (chainedRegions[i] == nullptr) {
            continue;
        }
        pdqsort(chainedRegions[i]->chain.hits.begin(), chainedRegions[i]->chain.hits.end());
    }

    // Remove the merged ones.
    std::vector<std::unique_ptr<ChainedRegion>> ret;
    for (size_t i = 0; i < chainedRegions.size(); ++i) {
        if (chainedRegions[i] == nullptr) {
            continue;
        }
        if (chainedRegions[i]->mapping == nullptr) {
            continue;
        }
        ret.emplace_back(std::move(chainedRegions[i]));
    }

    std::swap(chainedRegions, ret);
}

std::vector<AlignmentRegion> CollectAlignmentRegions(const ChainedRegion& singleMapping,
                                                     int32_t minAlignmentSpan,
                                                     int32_t maxFlankExtensionDist,
                                                     double flankExtensionFactor)
{
    if (singleMapping.mapping == nullptr) {
        return {};
    }

    const auto& ovl = singleMapping.mapping;
    const auto& chain = singleMapping.chain;
    const auto& sortedHits = chain.hits;

    // Sanity checks.
    if (ovl->Arev) {
        throw std::runtime_error("(CollectAlignmentRegions) The ovl->Arev should always be false!");
    }
    const int32_t Astart = ovl->Brev ? (ovl->Alen - ovl->Aend) : ovl->Astart;
    const int32_t Aend = ovl->Brev ? (ovl->Alen - ovl->Astart) : ovl->Aend;
    const int32_t Bstart = ovl->Brev ? (ovl->Blen - ovl->Bend) : ovl->Bstart;
    const int32_t Bend = ovl->Brev ? (ovl->Blen - ovl->Bstart) : ovl->Bend;

    if (Astart != sortedHits.front().queryPos ||
        Aend != (sortedHits.back().queryPos + sortedHits.back().querySpan) ||
        Bstart != sortedHits.front().targetPos ||
        Bend != (sortedHits.back().targetPos + sortedHits.back().targetSpan)) {
        std::ostringstream oss;
        oss << "(" << __FUNCTION__
            << ") Provided overlap coordinates do not match the first/last seed hit!"
            << " ovl: " << *ovl << "; sortedHits.front() = " << sortedHits.front()
            << "; sortedHits.back() = " << sortedHits.back() << "; Astart = " << Astart
            << ", Aend = " << Aend << ", Bstart = " << Bstart << ", Bend = " << Bend;
        throw std::runtime_error(oss.str());
    }
    if (ovl->Brev != sortedHits.front().targetRev || ovl->Brev != sortedHits.back().targetRev) {
        std::ostringstream oss;
        oss << "(" << __FUNCTION__
            << ") Strand of the provided overlap does not match the first/last "
               "seed hit."
            << " ovl->Brev = " << (ovl->Brev ? "true" : "false")
            << ", sortedHits.front().targetRev = "
            << (sortedHits.front().targetRev ? "true" : "false")
            << ", sortedHits.back().targetRev = "
            << (sortedHits.back().targetRev ? "true" : "false");
        throw std::runtime_error(oss.str());
    }

    // Extract the regions.
    std::vector<AlignmentRegion> regions =
        ExtractAlignmentRegions(chain.hits, ovl->Alen, ovl->Blen, ovl->Brev, minAlignmentSpan,
                                maxFlankExtensionDist, flankExtensionFactor);

    return regions;
}

}  // namespace Pancake
}  // namespace PacBio
