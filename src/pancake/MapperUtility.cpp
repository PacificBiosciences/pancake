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

std::vector<AlignmentRegion> ExtractAlignmentRegions(const ChainedRegion& singleMapping,
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
        throw std::runtime_error("(ExtractAlignmentRegions) The ovl->Arev should always be false!");
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

std::vector<AlignmentRegion> ExtractAlignmentRegions(const std::vector<SeedHit>& inSortedHits,
                                                     const int32_t qLen, const int32_t tLen,
                                                     const bool isRev,
                                                     const int32_t minAlignmentSpan,
                                                     int32_t maxFlankExtensionDist,
                                                     const double flankExtensionFactor)
{
    if (qLen < 0) {
        throw std::runtime_error("Invalid function parameter in ExtractAlignmentRegions! qLen = " +
                                 std::to_string(qLen) + ", should be >= 0.");
    }
    if (tLen < 0) {
        throw std::runtime_error("Invalid function parameter in ExtractAlignmentRegions! tLen = " +
                                 std::to_string(tLen) + ", should be >= 0.");
    }
    if (flankExtensionFactor < 1.0) {
        throw std::runtime_error(
            "Invalid function parameter in ExtractAlignmentRegions! flankExtensionFactor = " +
            std::to_string(flankExtensionFactor) + ", should be >= 1.0.");
    }

    if (inSortedHits.empty()) {
        return {};
    }

    if (maxFlankExtensionDist < 0) {
        maxFlankExtensionDist = std::max(qLen, tLen);
    }

    const std::vector<SeedHit>& hits = inSortedHits;

    std::vector<AlignmentRegion> ret;
    int32_t numRegions = 0;
    const int32_t globalAlnQueryStart = hits.front().queryPos;
    const int32_t globalAlnTargetStart = hits.front().targetPos;
    const int32_t globalAlnQueryEnd = hits.back().queryPos;
    const int32_t globalAlnTargetEnd = hits.back().targetPos;

    // Extract the front.
    if (globalAlnQueryStart > 0 && globalAlnTargetStart > 0) {
        // Determine the maximum flank we want to align, in both target and query.
        const int32_t projTStart =
            std::max(0.0, globalAlnTargetStart - flankExtensionFactor * globalAlnQueryStart);
        const int32_t projQStart = 0;
        const int32_t qExtLen = std::min(globalAlnQueryStart - projQStart, maxFlankExtensionDist);
        const int32_t tExtLen = std::min(globalAlnTargetStart - projTStart, maxFlankExtensionDist);

        // Find the frame of the sequences to compare.
        const int32_t qStart = globalAlnQueryStart - qExtLen;
        const int32_t tStart = globalAlnTargetStart - tExtLen;
        const int32_t qSpan = globalAlnQueryStart - qStart;
        const int32_t tSpan = globalAlnTargetStart - tStart;

        AlignmentRegion region;
        region.qStart = qStart;
        region.qSpan = qSpan;
        region.tStart = tStart;
        region.tSpan = tSpan;
        region.queryRev = isRev;
        region.type = RegionType::FRONT;
        region.regionId = numRegions;
        ++numRegions;
        ret.emplace_back(region);
    }

    if (!hits.empty()) {
        const auto& h1 = hits.front();
        // Sanity check that the first hit matches the strand specified via function call.
        if (h1.targetRev != isRev) {
            std::ostringstream oss;
            oss << "[" << __FUNCTION__
                << "] Hit strand does not match the strand specified as the "
                   "parameter to this function. First hit: "
                << h1;
            throw std::runtime_error(oss.str());
        }
        // Sanity check that the seed hit itself is valid.
        if ((h1.queryPos + h1.querySpan) > qLen || (h1.targetPos + h1.targetSpan) > tLen) {
            std::ostringstream oss;
            oss << "[" << __FUNCTION__
                << "] Seed hit coordinates/span are not valid, they span "
                   "out of the bounds of query or target. First hit: "
                << h1;
            throw std::runtime_error(oss.str());
        }
    }

    // Align between seeds.
    const int32_t nHits = hits.size();
    std::vector<Data::Cigar> cigarChunks;
    int32_t startId = 0;
    for (int32_t i = 1; i < nHits; ++i) {
        // Shorthands.
        const auto& h1 = hits[startId];
        const auto& h2 = hits[i];
        const auto& hPrev = hits[i - 1];

        // Compute the maximum indel gap between the two neighboring seed hits.
        const int32_t distQuery = h2.queryPos - hPrev.queryPos;
        const int32_t distTarget = h2.targetPos - hPrev.targetPos;
        const int32_t currGap = std::abs(distTarget - distQuery);

        // Sanity check that the strands are valid.
        if (h2.targetRev != hPrev.targetRev || h2.targetId != hPrev.targetId) {
            std::ostringstream oss;
            oss << "[" << __FUNCTION__
                << "] Hit target ID or strand is not consistent in the chain. Previous hit: "
                << hPrev << ". Current hit (i = " << i << "): " << h2;
            throw std::runtime_error(oss.str());
        }
        // Sanity check that the seed hit has the same strand as the specified via function arguments.
        if (h1.targetRev != isRev) {
            std::ostringstream oss;
            oss << "[" << __FUNCTION__
                << "] Hit strand does not match the strand specified as the "
                   "parameter to this function. Current hit (i = "
                << i << "): " << h2;
            throw std::runtime_error(oss.str());
        }
        // Sanity check that the coordinates are valid.
        if (h2.targetPos < hPrev.targetPos || h2.queryPos < hPrev.queryPos) {
            std::ostringstream oss;
            oss << "[" << __FUNCTION__
                << "] The chain of seed hits is not monotonically "
                   "increasing in terms of coordinates. Previous hit: "
                << hPrev << ". Current hit (i = " << i << "): " << h2;
            throw std::runtime_error(oss.str());
        }
        // Sanity check that the seed hit itself is valid.
        if ((h2.queryPos + h2.querySpan) > qLen || (h2.targetPos + h2.targetSpan) > tLen) {
            std::ostringstream oss;
            oss << "[" << __FUNCTION__
                << "] Seed hit coordinates/span are not valid, they span "
                   "out of the bounds of query or target. Current hit (i = "
                << i << "): " << h2;
            throw std::runtime_error(oss.str());
        }

        // Compute the new region.
        AlignmentRegion region;
        region.qStart = h1.queryPos;
        region.qSpan = h2.queryPos - h1.queryPos;
        region.tStart = h1.targetPos;
        region.tSpan = h2.targetPos - h1.targetPos;
        region.type = RegionType::GLOBAL;
        region.queryRev = isRev;
        region.regionId = numRegions;
        region.maxGap = std::max(region.maxGap, currGap);

        // Sanity check that the spans are valid.
        if (region.qSpan < 0 || region.tSpan < 0) {
            std::ostringstream oss;
            oss << "Region span not valid, in ExtractAlignmentRegions! qStart = " << region.qStart
                << ", qSpan = " << region.qSpan << ", tStart = " << region.tStart
                << ", tSpan = " << region.tSpan << "\n";
            throw std::runtime_error(oss.str());
        }

        // Skip short alignment portions.
        // if ((i + 1) < nHits && h2.CheckFlagLongJoin() == false &&
        if ((i + 1) < nHits &&
            (region.qSpan < minAlignmentSpan || region.tSpan < minAlignmentSpan)) {
            continue;
        }

        // Update the start ID for the next iteration.
        startId = i;
        ++numRegions;

        // Add the new region.
        ret.emplace_back(region);
    }

    // Back chunk.
    if (globalAlnQueryEnd < qLen && globalAlnTargetEnd < tLen) {
        // Determine the maximum flank we want to align, in both target and query.
        const int32_t qFlankLen = qLen - globalAlnQueryEnd;
        const int32_t projTEnd = std::min(
            tLen, static_cast<int32_t>(globalAlnTargetEnd + flankExtensionFactor * qFlankLen));
        const int32_t projQEnd = qLen;
        const int32_t qExtLen = std::min(projQEnd - globalAlnQueryEnd, maxFlankExtensionDist);
        const int32_t tExtLen = std::min(projTEnd - globalAlnTargetEnd, maxFlankExtensionDist);

        // Compute the new region.
        AlignmentRegion region;
        region.qStart = globalAlnQueryEnd;
        region.qSpan = qExtLen;
        region.tStart = globalAlnTargetEnd;
        region.tSpan = tExtLen;
        region.type = RegionType::BACK;
        region.queryRev = isRev;
        region.regionId = numRegions;
        ++numRegions;
        ret.emplace_back(region);
    }

    return ret;
}

}  // namespace Pancake
}  // namespace PacBio
