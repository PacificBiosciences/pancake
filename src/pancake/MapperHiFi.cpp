// Authors: Ivan Sovic

#include <lib/kxsort/kxsort.h>
#include <pacbio/alignment/AlignmentTools.h>
#include <pacbio/alignment/DiffCounts.h>
#include <pacbio/alignment/SesDistanceBanded.h>
#include <pacbio/pancake/MapperHiFi.h>
#include <pacbio/pancake/Minimizers.h>
#include <pacbio/pancake/OverlapWriterBase.h>
#include <pacbio/pancake/Secondary.h>
#include <pacbio/pancake/SeedHitWriter.h>
#include <pacbio/util/RunLengthEncoding.h>
#include <pacbio/util/TicToc.h>
#include <pacbio/util/Util.h>
#include <pbcopper/logging/Logging.h>
#include <pbcopper/third-party/edlib.h>
#include <algorithm>
#include <iostream>
#include <lib/istl/lis.hpp>
#include <pacbio/alignment/Ses2AlignBanded.hpp>
#include <pacbio/alignment/Ses2DistanceBanded.hpp>
#include <pacbio/alignment/SesAlignBanded.hpp>
#include <sstream>

namespace PacBio {
namespace Pancake {
namespace OverlapHiFi {

// #define PANCAKE_DEBUG
// #define PANCAKE_DEBUG_ALN

auto AlignWithTraceback(const char* query, size_t queryLen, const char* target, size_t targetLen,
                        int32_t maxDiffs, int32_t bandwidth,
                        std::shared_ptr<Alignment::SESScratchSpace> ss = nullptr)
{
    return Alignment::SES2AlignBanded<Alignment::SESAlignMode::Semiglobal,
                                      Alignment::SESTrimmingMode::Disabled,
                                      Alignment::SESTracebackMode::Enabled>(
        query, queryLen, target, targetLen, maxDiffs, bandwidth, ss);
}

auto AlignNoTraceback(const char* query, size_t queryLen, const char* target, size_t targetLen,
                      int32_t maxDiffs, int32_t bandwidth,
                      std::shared_ptr<Alignment::SESScratchSpace> ss = nullptr)
{
    return Alignment::SES2AlignBanded<Alignment::SESAlignMode::Semiglobal,
                                      Alignment::SESTrimmingMode::Disabled,
                                      Alignment::SESTracebackMode::Disabled>(
        query, queryLen, target, targetLen, maxDiffs, bandwidth, ss);
}

static const int32_t MIN_DIFFS_CAP = 10;
static const int32_t MIN_BANDWIDTH_CAP = 10;
// static const int32_t MASK_DEGREE = 3;

MapperResult Mapper::Map(const PacBio::Pancake::FastaSequenceCachedStore& targetSeqs,
                         const PacBio::Pancake::SeedIndex& index,
                         const PacBio::Pancake::FastaSequenceCached& querySeq,
                         const PacBio::Pancake::SequenceSeedsCached& querySeeds, int64_t freqCutoff,
                         bool generateFlippedOverlap) const
{
    return MapSingleQuery_(targetSeqs, index, querySeq, querySeeds, freqCutoff,
                           generateFlippedOverlap);
}

MapperResult Mapper::MapSingleQuery_(const PacBio::Pancake::FastaSequenceCachedStore& targetSeqs,
                                     const PacBio::Pancake::SeedIndex& index,
                                     const PacBio::Pancake::FastaSequenceCached& querySeq,
                                     const PacBio::Pancake::SequenceSeedsCached& querySeeds,
                                     int64_t freqCutoff, bool generateFlippedOverlap) const
{
#ifdef PANCAKE_DEBUG
    PBLOG_INFO << "Mapping query ID = " << querySeq.Id() << ", header = " << querySeq.Name();
#endif

    if (querySeq.Size() < settings_.MinQueryLen) {
        return {};
    }

    TicToc ttCollectHits;
    std::vector<SeedHit> hits;
    index.CollectHits(querySeeds.Seeds(), querySeeds.Size(), querySeq.Size(), hits, freqCutoff);
    ttCollectHits.Stop();

    TicToc ttSortHits;
    std::sort(hits.begin(), hits.end(), [](const auto& a, const auto& b) {
        return PackSeedHitWithDiagonalToTuple(a) < PackSeedHitWithDiagonalToTuple(b);
    });
    ttSortHits.Stop();

    // PBLOG_INFO << "Hits: " << hits.size();

    TicToc ttChain;
    auto overlaps =
        FormAnchors2_(hits, querySeq, targetSeqs, index, settings_.ChainBandwidth,
                      settings_.MinNumSeeds, settings_.MinChainSpan, index.GetMinSeedSpan() * 3,
                      settings_.SkipSelfHits, settings_.SkipSymmetricOverlaps);
    ttChain.Stop();
#ifdef PANCAKE_DEBUG
    PBLOG_INFO << "Formed diagonal anchors: " << overlaps.size();
#endif

    // Filter out multiple hits per query-target pair (e.g. tandem repeats) by
    // taking only the longest overlap chain.
    TicToc ttFilterTandem;
    if (settings_.OneHitPerTarget) {
        overlaps = FilterTandemOverlaps_(overlaps);
    }
    ttFilterTandem.Stop();
#ifdef PANCAKE_DEBUG
    PBLOG_INFO << "Overlaps after tandem filtering: " << overlaps.size();
#endif

    const std::string reverseQuerySeq =
        PacBio::Pancake::ReverseComplement(querySeq.Bases(), querySeq.Size(), 0, querySeq.Size());

    TicToc ttAlign;
    overlaps = AlignOverlaps_(
        targetSeqs, querySeq, reverseQuerySeq, overlaps, settings_.AlignmentBandwidth,
        settings_.AlignmentMaxD, settings_.UseTraceback, settings_.NoSNPsInIdentity,
        settings_.NoIndelsInIdentity, settings_.MaskHomopolymers, settings_.MaskSimpleRepeats,
        settings_.MaskHomopolymerSNPs, settings_.MaskHomopolymersArbitrary, settings_.TrimAlignment,
        settings_.TrimWindowSize, settings_.TrimWindowMatchFraction, settings_.TrimToFirstMatch,
        sesScratch_);
    ttAlign.Stop();

    TicToc ttMarkSecondary;
    if (settings_.MarkSecondary) {
        // Flag the secondary and supplementary overlaps.
        // Overlaps don't have to be sorted, the maximum is found internally.
        std::vector<OverlapPriority> overlapPriorities = FlagSecondaryAndSupplementary(
            overlaps, settings_.SecondaryAllowedOverlapFraction,
            settings_.SecondaryAllowedOverlapFraction, settings_.SecondaryMinScoreFraction);

        // Generate a new, filtered list of overlaps. The non-primary, non-secondary and non-supplementary
        // alignments are filtered out.
        std::vector<OverlapPtr> newOverlaps;
        for (size_t i = 0; i < overlaps.size(); ++i) {
            if (overlapPriorities[i].priority > 1) {
                continue;
            }
            auto& ovl = overlaps[i];
            ovl->IsSupplementary = overlapPriorities[i].isSupplementary;
            ovl->IsSecondary = (overlapPriorities[i].priority > 0);
            newOverlaps.emplace_back(std::move(overlaps[i]));
        }
        std::swap(overlaps, newOverlaps);
    }
    ttMarkSecondary.Stop();

    // Filter the overlaps.
    TicToc ttFilter;
    overlaps = FilterOverlaps_(
        overlaps, settings_.MinNumSeeds, settings_.MinIdentity, settings_.MinMappedLength,
        settings_.MinQueryLen, settings_.MinTargetLen, settings_.ChainBandwidth,
        settings_.AllowedDovetailDist, settings_.AllowedHeuristicExtendDist, settings_.BestN);
    ttFilter.Stop();

    // Generating flipped overlaps.
    TicToc ttFlip;
    if (generateFlippedOverlap) {
        std::vector<OverlapPtr> flippedOverlaps = GenerateFlippedOverlaps_(
            targetSeqs, querySeq, reverseQuerySeq, overlaps, settings_.NoSNPsInIdentity,
            settings_.NoIndelsInIdentity, settings_.MaskHomopolymers, settings_.MaskSimpleRepeats,
            settings_.MaskHomopolymerSNPs, settings_.MaskHomopolymersArbitrary);
        for (size_t i = 0; i < flippedOverlaps.size(); ++i) {
            overlaps.emplace_back(std::move(flippedOverlaps[i]));
        }
    }
    ttFlip.Stop();

#ifdef PANCAKE_DEBUG
    PBLOG_INFO << "Overlaps after alignment: " << overlaps.size();
    for (const auto& ovl : overlaps) {
        PBLOG_INFO << *ovl;
    }
#endif

#ifdef PANCAKE_DEBUG
    PBLOG_INFO << "Overlaps after filtering: " << overlaps.size();
    for (const auto& ovl : overlaps) {
        OverlapWriterBase::PrintOverlapAsM4(stderr, *ovl, querySeq.Name(),
                                            targetSeqs.GetSequence(ovl->Bid).Name(), false, false);
    }
    PBLOG_INFO << "Num anchors: " << overlaps.size();
    PBLOG_INFO << "Collected " << hits.size() << " hits.";
    PBLOG_INFO << "Time - collecting hits: " << ttCollectHits.GetMillisecs() << " ms";
    PBLOG_INFO << "Time - sorting: " << ttSortHits.GetMillisecs() << " ms";
    PBLOG_INFO << "Time - chaining: " << ttChain.GetMillisecs() << " ms";
    PBLOG_INFO << "Time - tandem filter: " << ttFilterTandem.GetMillisecs() << " ms";
    PBLOG_INFO << "Time - alignment: " << ttAlign.GetMillisecs() << " ms";
    PBLOG_INFO << "Time - filter: " << ttFilter.GetMillisecs() << " ms";
    DebugWriteSeedHits_("temp/debug/mapper-0-seed_hits.csv", hits, 30, querySeq.Name(),
                        querySeq.Size(), "target", 0);
#endif

    MapperResult result;
    std::swap(result.overlaps, overlaps);
    return result;
}

OverlapPtr Mapper::MakeOverlap_(const std::vector<SeedHit>& sortedHits,
                                const PacBio::Pancake::FastaSequenceCached& querySeq,
                                const PacBio::Pancake::FastaSequenceCachedStore& targetSeqs,
                                int32_t numSeeds, int32_t minTargetPosId, int32_t maxTargetPosId)
{

    const auto& beginHit = sortedHits[minTargetPosId];
    const auto& endHit = sortedHits[maxTargetPosId];

    const int32_t targetId = beginHit.targetId;

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
        createOverlap(querySeq.Id(), targetId, score, identity, beginHit.targetRev,
                      beginHit.queryPos, endHit.queryPos + endHit.querySpan, querySeq.Size(), false,
                      beginHit.targetPos, endHit.targetPos + endHit.targetSpan, targetLen, editDist,
                      numSeeds, OverlapType::Unknown, OverlapType::Unknown);

    ret->NormalizeStrand();

    return ret;
}

std::vector<OverlapPtr> Mapper::FormAnchors_(
    const std::vector<SeedHit>& sortedHits, const PacBio::Pancake::FastaSequenceCached& querySeq,
    const PacBio::Pancake::FastaSequenceCachedStore& targetSeqs, int32_t chainBandwidth,
    int32_t minNumSeeds, int32_t minChainSpan, bool skipSelfHits, bool skipSymmetricOverlaps)
{
#ifdef PANCAKE_DEBUG
    std::cerr << "[Function: " << __FUNCTION__ << "]\n";
#endif

    if (sortedHits.empty()) {
        return {};
    }

    std::vector<OverlapPtr> overlaps;

    const int32_t numHits = static_cast<int32_t>(sortedHits.size());
    int32_t beginId = 0;
    int32_t beginDiag = sortedHits[beginId].Diagonal();

    // This is a combination of <targetPos, queryPos>, intended for simple comparison
    // without defining a custom complicated comparison operator.
    uint64_t minTargetQueryPosCombo = (static_cast<uint64_t>(sortedHits[beginId].targetPos) << 32) |
                                      (static_cast<uint64_t>(sortedHits[beginId].queryPos));
    uint64_t maxTargetQueryPosCombo = minTargetQueryPosCombo;
    int32_t minPosId = 0;
    int32_t maxPosId = 0;

    for (int32_t i = 0; i < numHits; ++i) {
        const auto& prevHit = sortedHits[beginId];
        const auto& currHit = sortedHits[i];
        const int32_t currDiag = currHit.Diagonal();
        const int32_t diagDiff = abs(currDiag - beginDiag);
        const uint64_t targetQueryPosCombo =
            (static_cast<uint64_t>(sortedHits[i].targetPos) << 32) |
            (static_cast<uint64_t>(sortedHits[i].queryPos));

#ifdef PANCAKE_DEBUG
        std::cerr << "[hit " << i << "] tid = " << currHit.targetId
                  << ", trev = " << currHit.targetRev << ", tpos = " << currHit.targetPos
                  << ", qpos = " << currHit.queryPos << ", flag = " << currHit.flags
                  << ", diag = " << currDiag << ", beginDiag = " << beginDiag
                  << ", diagDiff = " << diagDiff << "; minPosId = " << minPosId
                  << ", maxPosId = " << maxPosId << "\n";
#endif

        if (currHit.targetId != prevHit.targetId || currHit.targetRev != prevHit.targetRev ||
            diagDiff > chainBandwidth) {
            auto ovl =
                MakeOverlap_(sortedHits, querySeq, targetSeqs, i - beginId, minPosId, maxPosId);
            beginId = i;
            beginDiag = currDiag;

#ifdef PANCAKE_DEBUG
            std::cerr << "ovl->NumSeeds = " << ovl->NumSeeds << " (" << minNumSeeds
                      << "), minChainSpan = " << minChainSpan << ", ovl->ASpan() = " << ovl->ASpan()
                      << ", ovl->BSpan() = " << ovl->BSpan() << ", skipSelfHits = " << skipSelfHits
                      << ", ovl->Aid = " << ovl->Aid << ", ovl->Bid = " << ovl->Bid
                      << ", skipSymmetricOverlaps = " << skipSymmetricOverlaps << "\n";
            std::cerr << *ovl << "\n";
#endif

            // Add a new overlap.
            if (ovl->NumSeeds >= minNumSeeds && ovl->ASpan() > minChainSpan &&
                ovl->BSpan() > minChainSpan &&
                (skipSelfHits == false || (skipSelfHits && ovl->Bid != ovl->Aid)) &&
                (skipSymmetricOverlaps == false ||
                 (skipSymmetricOverlaps && ovl->Bid < ovl->Aid))) {
                overlaps.emplace_back(std::move(ovl));
            }
            minPosId = maxPosId = i;
            minTargetQueryPosCombo = maxTargetQueryPosCombo = targetQueryPosCombo;
        }

        // Track the minimum and maximum target positions for each diagonal.
        if (targetQueryPosCombo < minTargetQueryPosCombo) {
            minPosId = i;
            minTargetQueryPosCombo = targetQueryPosCombo;
        }
        if (targetQueryPosCombo > maxTargetQueryPosCombo) {
            maxPosId = i;
            maxTargetQueryPosCombo = targetQueryPosCombo;
        }
    }

    if ((numHits - beginId) > 0) {
        auto ovl =
            MakeOverlap_(sortedHits, querySeq, targetSeqs, numHits - beginId, minPosId, maxPosId);

#ifdef PANCAKE_DEBUG
        std::cerr << "ovl->NumSeeds = " << ovl->NumSeeds << " (" << minNumSeeds
                  << "), minChainSpan = " << minChainSpan << ", ovl->ASpan() = " << ovl->ASpan()
                  << ", ovl->BSpan() = " << ovl->BSpan() << ", skipSelfHits = " << skipSelfHits
                  << ", ovl->Aid = " << ovl->Aid << ", ovl->Bid = " << ovl->Bid
                  << ", skipSymmetricOverlaps = " << skipSymmetricOverlaps << "\n";
#endif

        // Add a new overlap.
        if (ovl->NumSeeds >= minNumSeeds && ovl->ASpan() > minChainSpan &&
            ovl->BSpan() > minChainSpan &&
            (skipSelfHits == false || (skipSelfHits && ovl->Bid != ovl->Aid)) &&
            (skipSymmetricOverlaps == false || (skipSymmetricOverlaps && ovl->Bid < ovl->Aid))) {

            overlaps.emplace_back(std::move(ovl));
        }
    }

    return overlaps;
}

void RefineBadEnds(const std::vector<SeedHit>& chainedHits, int32_t beginId, int32_t endId,
                   int32_t kmerSize, int32_t bandwidth, int32_t minMatch, int32_t& retFirst,
                   int32_t& retLast)
{
    const int32_t n = chainedHits.size();
    if (beginId < 0 || beginId >= n || endId < 0 || endId > n || beginId >= endId) {
        std::ostringstream oss;
        oss << "Invalid beginId or endId in a call t o RefineBadEnds. beginId = " << beginId
            << ", endId = " << endId << ", chainedHits.size() = " << n;
        throw std::runtime_error(oss.str());
    }

    retFirst = beginId;
    retLast = endId;

    if (chainedHits.size() < 3) {
        return;
    }

    int32_t coveredBasesQuery = 0;
    int32_t coveredBasesTarget = 0;
    CalcHitCoverage(chainedHits, beginId, endId, coveredBasesQuery, coveredBasesTarget);
    const int32_t minCoveredBases = std::min(coveredBasesQuery, coveredBasesTarget);

    // Front.
    {
        int32_t numMatches = kmerSize;
        int32_t totalSpan = kmerSize;
        for (int32_t i = (beginId + 1); i < (endId - 1); ++i) {
            if (chainedHits[i].CheckFlagLongJoin()) {
                break;
            }
            const int32_t qDist = chainedHits[i].queryPos - chainedHits[i - 1].queryPos;
            const int32_t tDist = chainedHits[i].targetPos - chainedHits[i - 1].targetPos;
            const int32_t minDist = std::min(qDist, tDist);
            const int32_t maxDist = std::max(qDist, tDist);
            const int32_t gap = maxDist - minDist;
            if (gap > (totalSpan >> 1)) {
                retFirst = i;
            }
            totalSpan += minDist;
            numMatches += std::min(minDist, kmerSize);
            if (totalSpan >= (bandwidth << 1) ||
                (numMatches >= minMatch && numMatches >= bandwidth) ||
                numMatches >= (minCoveredBases >> 1)) {
                break;
            }
        }
    }

    // Back.
    {
        int32_t numMatches = kmerSize;
        int32_t totalSpan = kmerSize;
        for (int32_t i = endId - 2; i > beginId; --i) {
            if (chainedHits[i + 1].CheckFlagLongJoin()) {
                break;
            }
            const int32_t qDist = chainedHits[i + 1].queryPos - chainedHits[i].queryPos;
            const int32_t tDist = chainedHits[i + 1].targetPos - chainedHits[i].targetPos;
            const int32_t minDist = std::min(qDist, tDist);
            const int32_t maxDist = std::max(qDist, tDist);
            const int32_t gap = maxDist - minDist;
            if (gap > (totalSpan >> 1)) {
                retLast = i + 1;
            }
            totalSpan += minDist;
            numMatches += std::min(minDist, kmerSize);
            if (totalSpan >= (bandwidth << 1) ||
                (numMatches >= minMatch && numMatches >= bandwidth) ||
                numMatches >= (minCoveredBases >> 1)) {
                break;
            }
        }
    }
}

std::vector<OverlapPtr> Mapper::FormAnchors2_(
    const std::vector<SeedHit>& sortedHits, const PacBio::Pancake::FastaSequenceCached& querySeq,
    const PacBio::Pancake::FastaSequenceCachedStore& targetSeqs,
    const PacBio::Pancake::SeedIndex& index, int32_t chainBandwidth, int32_t minNumSeeds,
    int32_t minChainSpan, int32_t minMatch, bool skipSelfHits, bool skipSymmetricOverlaps)
{
#ifdef PANCAKE_DEBUG
    std::cerr << "[Function: " << __FUNCTION__ << "]\n";
#endif

    if (sortedHits.empty()) {
        return {};
    }

    auto WrapMakeOverlap = [](
        const std::vector<SeedHit>& _sortedHits, const int32_t beginId, const int32_t endId,
        const PacBio::Pancake::FastaSequenceCached& _querySeq,
        const PacBio::Pancake::FastaSequenceCachedStore& _targetSeqs, int32_t _chainBandwidth,
        int32_t kmerSize, int32_t _minMatch) -> OverlapPtr {
        if (endId <= beginId) {
            return nullptr;
        }

        static std::function<bool(const SeedHit& a, const SeedHit& b)> ComparisonLIS = [](
            const SeedHit& a, const SeedHit& b) {
            // This needs to always return the upper-left element as the smaller one.
            return a.targetPos < b.targetPos && a.queryPos < b.queryPos;
        };

        // Extract the subset so we can sort it.
        std::vector<SeedHit> groupHits(_sortedHits.begin() + beginId, _sortedHits.begin() + endId);

        std::sort(groupHits.begin(), groupHits.end(), [](const SeedHit& a, const SeedHit& b) {
            return std::pair(a.targetPos, a.queryPos) < std::pair(b.targetPos, b.queryPos);
        });

        // Longest Increasing Subsequence of the diagonal bin.
        std::vector<PacBio::Pancake::SeedHit> lisHits =
            istl::LIS(groupHits, 0, groupHits.size(), ComparisonLIS);

        int32_t finalFirst = 0;
        int32_t finalLast = 0;
        RefineBadEnds(lisHits, 0, lisHits.size(), kmerSize, _chainBandwidth, _minMatch, finalFirst,
                      finalLast);

#ifdef PANCAKE_DEBUG
        std::cerr << "LIS hits for a group, beginId = " << beginId << ", endId = " << endId << "\n";
        for (size_t i = 0; i < lisHits.size(); ++i) {
            std::cerr << "  [lis hit " << i << "] " << lisHits[i] << "\n";
        }
        std::cerr << "RefineBadEnds: finalFirst = " << finalFirst << ", finalLast = " << finalLast
                  << "\n";
#endif

        // Make the overlap.
        auto ovl = MakeOverlap_(lisHits, _querySeq, _targetSeqs, (finalLast - finalFirst),
                                finalFirst, finalLast - 1);

#ifdef PANCAKE_DEBUG
        if (ovl->NumSeeds > 200) {
            WriteSeedHits("temp-debug/hits-q" + std::to_string(_querySeq.Id()) + "-1-group_" +
                              std::to_string(beginId) + "_" + std::to_string(endId) +
                              "-all_hits.csv",
                          _sortedHits, beginId, endId, beginId, "query-" + std::to_string(ovl->Aid),
                          ovl->Alen, "target-" + std::to_string(ovl->Bid), ovl->Blen, false);
            WriteSeedHits("temp-debug/hits-q" + std::to_string(_querySeq.Id()) + "-2-group_" +
                              std::to_string(beginId) + "_" + std::to_string(endId) + "-lis.csv",
                          lisHits, 0, lisHits.size(), beginId, "query-" + std::to_string(ovl->Aid),
                          ovl->Alen, "target-" + std::to_string(ovl->Bid), ovl->Blen, false);
        }
#endif

        return ovl;
    };

    std::vector<OverlapPtr> overlaps;

    const int32_t numHits = static_cast<int32_t>(sortedHits.size());
    int32_t beginId = 0;
    int32_t beginDiag = sortedHits[beginId].Diagonal();

    for (int32_t i = 0; i < numHits; ++i) {
        const auto& prevHit = sortedHits[beginId];
        const auto& currHit = sortedHits[i];
        const int32_t currDiag = currHit.Diagonal();
        const int32_t diagDiff = abs(currDiag - beginDiag);

#ifdef PANCAKE_DEBUG
        std::cerr << "[hit " << i << "] tid = " << currHit.targetId
                  << ", trev = " << currHit.targetRev << ", tpos = " << currHit.targetPos
                  << ", qpos = " << currHit.queryPos << ", flag = " << currHit.flags
                  << ", diag = " << currDiag << ", beginDiag = " << beginDiag
                  << ", diagDiff = " << diagDiff << "\n";
#endif

        if (currHit.targetId != prevHit.targetId || currHit.targetRev != prevHit.targetRev ||
            diagDiff > chainBandwidth) {

            auto ovl = WrapMakeOverlap(sortedHits, beginId, i, querySeq, targetSeqs, chainBandwidth,
                                       index.GetMinSeedSpan(), minMatch);
            beginId = i;
            beginDiag = currDiag;

#ifdef PANCAKE_DEBUG
            std::cerr << "ovl->NumSeeds = " << ovl->NumSeeds << " (" << minNumSeeds
                      << "), minChainSpan = " << minChainSpan << ", ovl->ASpan() = " << ovl->ASpan()
                      << ", ovl->BSpan() = " << ovl->BSpan() << ", skipSelfHits = " << skipSelfHits
                      << ", ovl->Aid = " << ovl->Aid << ", ovl->Bid = " << ovl->Bid
                      << ", skipSymmetricOverlaps = " << skipSymmetricOverlaps << "\n";
            std::cerr << *ovl << "\n";
#endif

            // Add a new overlap.
            if (ovl->NumSeeds >= minNumSeeds && ovl->ASpan() > minChainSpan &&
                ovl->BSpan() > minChainSpan &&
                (skipSelfHits == false || (skipSelfHits && ovl->Bid != ovl->Aid)) &&
                (skipSymmetricOverlaps == false ||
                 (skipSymmetricOverlaps && ovl->Bid < ovl->Aid))) {
                overlaps.emplace_back(std::move(ovl));
            }
        }
    }

    if ((numHits - beginId) > 0) {
        auto ovl = WrapMakeOverlap(sortedHits, beginId, numHits, querySeq, targetSeqs,
                                   chainBandwidth, index.GetMinSeedSpan(), minMatch);

#ifdef PANCAKE_DEBUG
        std::cerr << "ovl->NumSeeds = " << ovl->NumSeeds << " (" << minNumSeeds
                  << "), minChainSpan = " << minChainSpan << ", ovl->ASpan() = " << ovl->ASpan()
                  << ", ovl->BSpan() = " << ovl->BSpan() << ", skipSelfHits = " << skipSelfHits
                  << ", ovl->Aid = " << ovl->Aid << ", ovl->Bid = " << ovl->Bid
                  << ", skipSymmetricOverlaps = " << skipSymmetricOverlaps << "\n";
#endif

        // Add a new overlap.
        if (ovl->NumSeeds >= minNumSeeds && ovl->ASpan() > minChainSpan &&
            ovl->BSpan() > minChainSpan &&
            (skipSelfHits == false || (skipSelfHits && ovl->Bid != ovl->Aid)) &&
            (skipSymmetricOverlaps == false || (skipSymmetricOverlaps && ovl->Bid < ovl->Aid))) {

            overlaps.emplace_back(std::move(ovl));
        }
    }

    return overlaps;
}

std::vector<OverlapPtr> Mapper::FilterOverlaps_(const std::vector<OverlapPtr>& overlaps,
                                                int32_t minNumSeeds, float minIdentity,
                                                int32_t minMappedSpan, int32_t minQueryLen,
                                                int32_t minTargetLen, int32_t diagonalBandwidth,
                                                int32_t allowedDovetailDist,
                                                int32_t allowedExtendDist, int32_t bestN)
{
    std::vector<OverlapPtr> newOverlaps;
    for (const auto& ovl : overlaps) {
        if (100 * ovl->Identity < minIdentity || ovl->ASpan() < minMappedSpan ||
            ovl->BSpan() < minMappedSpan || ovl->NumSeeds < minNumSeeds ||
            ovl->Alen < minQueryLen || ovl->Blen < minTargetLen) {
            continue;
        }
        auto newOvl = createOverlap(ovl);
        newOvl->Atype =
            DetermineOverlapType(ovl->Arev, ovl->AstartFwd(), ovl->AendFwd(), ovl->Alen, ovl->Brev,
                                 ovl->BstartFwd(), ovl->BendFwd(), ovl->Blen, allowedDovetailDist);
        // Arev and Brev are intentionally out of place here! The A-read's orientation should always
        // be FWD, so the B-read is the one that determines the direction.
        newOvl->Btype =
            DetermineOverlapType(ovl->Arev, ovl->BstartFwd(), ovl->BendFwd(), ovl->Blen, ovl->Brev,
                                 ovl->AstartFwd(), ovl->AendFwd(), ovl->Alen, allowedDovetailDist);
        HeuristicExtendOverlapFlanks(newOvl, allowedExtendDist);
        newOverlaps.emplace_back(std::move(newOvl));
    }

    // Sort by diagonal, to filter the duplicate overlaps.
    std::stable_sort(newOverlaps.begin(), newOverlaps.end(), [](const auto& a, const auto& b) {
        return std::make_tuple(a->Bid, a->Brev, (a->Bstart - a->Astart)) <
               std::make_tuple(b->Bid, b->Brev, (b->Bstart - b->Astart));
    });
    for (size_t i = 0; i < newOverlaps.size(); ++i) {
        if (newOverlaps[i] == nullptr) {
            continue;
        }
        for (size_t j = (i + 1); j < newOverlaps.size(); ++j) {
            if (newOverlaps[j] == nullptr) {
                continue;
            }
            // Stop the loop if we reached a different target or orientation.
            if (newOverlaps[j]->Bid != newOverlaps[i]->Bid ||
                newOverlaps[j]->Brev != newOverlaps[i]->Brev) {
                break;
            }
            // Break if the diagonal is too far away from the current overlap.
            const int32_t diagI = newOverlaps[i]->Bstart - newOverlaps[i]->Astart;
            const int32_t diagJ = newOverlaps[j]->Bstart - newOverlaps[j]->Astart;
            if (std::abs(diagI - diagJ) > diagonalBandwidth) {
                break;
            }
            // Two overlaps are within the bandwidth. Remove the one with lower score.
            // Score is negative, as per legacy Falcon convention.
            if (newOverlaps[i]->Score > newOverlaps[j]->Score) {
                newOverlaps[i] = nullptr;
                break;
            } else {
                newOverlaps[j] = nullptr;
            }
        }
    }

    // Collect remaining overlaps.
    std::vector<OverlapPtr> ret;
    for (size_t i = 0; i < newOverlaps.size(); ++i) {
        if (newOverlaps[i] == nullptr) {
            continue;
        }
        ret.emplace_back(std::move(newOverlaps[i]));
    }
    // Score is negative, as per legacy Falcon convention.
    std::stable_sort(ret.begin(), ret.end(),
                     [](const auto& a, const auto& b) { return a->Score < b->Score; });
    // Keep best N.
    int32_t nToKeep = (bestN > 0) ? std::min(bestN, static_cast<int32_t>(ret.size())) : ret.size();
    ret.resize(nToKeep);

    return ret;
}

std::vector<OverlapPtr> Mapper::FilterTandemOverlaps_(const std::vector<OverlapPtr>& overlaps)
{
    if (overlaps.empty()) {
        return {};
    }

    // Make an internal copy for sorting.
    std::vector<OverlapPtr> overlaps_copy;
    for (const auto& ovl : overlaps) {
        overlaps_copy.emplace_back(createOverlap(ovl));
    }

    // Sort by length.
    std::sort(overlaps_copy.begin(), overlaps_copy.end(), [](const auto& a, const auto& b) {
        return a->Bid < b->Bid ||
               (a->Bid == b->Bid &&
                std::max(a->ASpan(), a->BSpan()) > std::max(b->ASpan(), b->BSpan()));
    });

    // Accumulate the results.
    std::vector<OverlapPtr> ret;
    ret.emplace_back(createOverlap(overlaps_copy.front()));
    for (const auto& ovl : overlaps_copy) {
        if (ovl->Bid == ret.back()->Bid) {
            continue;
        }
        ret.emplace_back(createOverlap(ovl));
    }

    return ret;
}

std::vector<OverlapPtr> Mapper::AlignOverlaps_(
    const PacBio::Pancake::FastaSequenceCachedStore& targetSeqs,
    const PacBio::Pancake::FastaSequenceCached& querySeq, const std::string reverseQuerySeq,
    const std::vector<OverlapPtr>& overlaps, double alignBandwidth, double alignMaxDiff,
    bool useTraceback, bool noSNPs, bool noIndels, bool maskHomopolymers, bool maskSimpleRepeats,
    bool maskHomopolymerSNPs, bool maskHomopolymersArbitrary, bool trimAlignment,
    int32_t trimWindowSize, double trimMatchFraction, bool trimToFirstMatch,
    std::shared_ptr<PacBio::Pancake::Alignment::SESScratchSpace> sesScratch)
{
    std::vector<OverlapPtr> ret;

    for (size_t i = 0; i < overlaps.size(); ++i) {
#ifdef PANCAKE_DEBUG_ALN
        PBLOG_INFO << "Aligning overlap: [" << i << "] " << *overlaps[i];
#endif
        const auto& targetSeq = targetSeqs.GetSequence(overlaps[i]->Bid);
        OverlapPtr newOverlap = AlignOverlap_(
            targetSeq, querySeq, reverseQuerySeq, overlaps[i], alignBandwidth, alignMaxDiff,
            useTraceback, noSNPs, noIndels, maskHomopolymers, maskSimpleRepeats,
            maskHomopolymerSNPs, maskHomopolymersArbitrary, trimAlignment, trimWindowSize,
            trimMatchFraction, trimToFirstMatch, sesScratch);
        if (newOverlap != nullptr) {
            ret.emplace_back(std::move(newOverlap));
#ifdef PANCAKE_DEBUG_ALN
            PBLOG_INFO << "After alignment: " << *overlaps[i];
#endif
        }

#ifdef PANCAKE_DEBUG_ALN
        PBLOG_INFO << "\n";
#endif
    }

    return ret;
}

std::vector<OverlapPtr> Mapper::GenerateFlippedOverlaps_(
    const PacBio::Pancake::FastaSequenceCachedStore& targetSeqs,
    const PacBio::Pancake::FastaSequenceCached& querySeq, const std::string reverseQuerySeq,
    const std::vector<OverlapPtr>& overlaps, bool noSNPs, bool noIndels, bool maskHomopolymers,
    bool maskSimpleRepeats, bool maskHomopolymerSNPs, bool maskHomopolymersArbitrary)
{
    std::vector<OverlapPtr> ret;

    for (size_t i = 0; i < overlaps.size(); ++i) {
        const auto& targetSeq = targetSeqs.GetSequence(overlaps[i]->Bid);
        const auto& ovl = overlaps[i];

        OverlapPtr newOverlapFlipped = CreateFlippedOverlap(ovl);
        NormalizeAndExtractVariantsInPlace_(newOverlapFlipped, targetSeq, querySeq, reverseQuerySeq,
                                            noSNPs, noIndels, maskHomopolymers, maskSimpleRepeats,
                                            maskHomopolymerSNPs, maskHomopolymersArbitrary);

        if (newOverlapFlipped == nullptr) {
            std::ostringstream oss;
            oss << "Problem generating the flipped overlap, it's nullptr! Before flipping: "
                << *ovl;
            throw std::runtime_error(oss.str());
        }

        ret.emplace_back(std::move(newOverlapFlipped));
    }

    return ret;
}

std::string Mapper::FetchTargetSubsequence_(const PacBio::Pancake::FastaSequenceCached& targetSeq,
                                            int32_t seqStart, int32_t seqEnd, bool revCmp)
{
    return FetchTargetSubsequence_(targetSeq.Bases(), targetSeq.Size(), seqStart, seqEnd, revCmp);
}

std::string Mapper::FetchTargetSubsequence_(const char* seq, int32_t seqLen, int32_t seqStart,
                                            int32_t seqEnd, bool revCmp)
{
    if (seqEnd == seqStart) {
        return {};
    }
    if (seqStart < 0 || seqEnd < 0 || seqStart > seqLen || seqEnd > seqLen || seqEnd < seqStart) {
        std::ostringstream oss;
        oss << "Invalid seqStart or seqEnd in a call to FetchTargetSubsequence_. seqStart = "
            << seqStart << ", seqEnd = " << seqEnd << ", seqLen = " << seqLen
            << ", revCmp = " << revCmp << ".";
        std::cerr << oss.str() << "\n";
        throw std::runtime_error(oss.str());
    }
    seqEnd = (seqEnd == 0) ? seqLen : seqEnd;
    std::string ret(seq + seqStart, seqEnd - seqStart);
    if (revCmp) {
        ret = PacBio::Pancake::ReverseComplement(ret, 0, ret.size());
    }
    return ret;
}

OverlapPtr Mapper::AlignOverlap_(
    const PacBio::Pancake::FastaSequenceCached& targetSeq,
    const PacBio::Pancake::FastaSequenceCached& querySeq, const std::string reverseQuerySeq,
    const OverlapPtr& ovl, double alignBandwidth, double alignMaxDiff, bool useTraceback,
    bool noSNPs, bool noIndels, bool maskHomopolymers, bool maskSimpleRepeats,
    bool maskHomopolymerSNPs, bool maskHomopolymersArbitrary, bool trimAlignment,
    int32_t trimWindowSize, double trimMatchFraction, bool trimToFirstMatch,
    std::shared_ptr<PacBio::Pancake::Alignment::SESScratchSpace> sesScratch)
{

    if (ovl == nullptr) {
        return nullptr;
    }

#ifdef PANCAKE_DEBUG_ALN
    PBLOG_INFO << "Initial: " << *ovl;
#endif

    OverlapPtr ret = createOverlap(ovl);
    PacBio::Pancake::Alignment::SesResults sesResultRight;
    PacBio::Pancake::Alignment::SesResults sesResultLeft;

    ///////////////////////////
    /// Align forward pass. ///
    ///////////////////////////
    {
        const int32_t qStart = ovl->Astart;
        const int32_t qEnd = ovl->Alen;
        const int32_t qSpan = qEnd - qStart;
        const int32_t tStartFwd = ovl->Brev ? (ovl->Blen - ovl->Bend) : ovl->Bstart;
        const int32_t tEndFwd = ovl->Brev ? (ovl->Blen - ovl->Bstart) : ovl->Bend;
        std::string tseq;
        if (ovl->Brev) {
            // Extract reverse complemented target sequence.
            // The reverse complement begins at the last mapped position (tEndFwd),
            // and ends at the the tStartFwd reduced by an allowed overhang.
            // The allowed overhang is designed to absorb the entire unmapped end
            // of the query, and at most 2x that length in the target. (The 2xhang
            // is just to be safe, because no proper alignment should be 2x longer
            // in one sequence than in the other).
            int32_t minHangLen = std::min(ovl->Alen - ovl->Aend, tStartFwd);
            int32_t extractBegin = std::max(0, tStartFwd - minHangLen * 2);
            int32_t extractEnd = tEndFwd;
            tseq = FetchTargetSubsequence_(targetSeq, extractBegin, extractEnd, ovl->Brev);
        } else {
            // Take the sequence starting from the start position, and reaching
            // until the end of the query (or target, which ever is the shorter).
            // Extract 2x larger overhang to be safe - no proper alignment should be
            // 2x longer in one sequence than the other.
            int32_t minHangLen = std::min(ovl->Blen - tEndFwd, ovl->Alen - ovl->Aend);
            int32_t extractBegin = tStartFwd;
            int32_t extractEnd = std::min(ovl->Blen, tEndFwd + minHangLen * 2);
            tseq = FetchTargetSubsequence_(targetSeq, extractBegin, extractEnd, ovl->Brev);
        }
        const int32_t tSpan = tseq.size();
        const int32_t dMax =
            std::max(MIN_DIFFS_CAP, static_cast<int32_t>(ovl->Alen * alignMaxDiff));
        const int32_t bandwidth =
            std::max(MIN_BANDWIDTH_CAP,
                     static_cast<int32_t>(std::min(ovl->Blen, ovl->Alen) * alignBandwidth));

        if (useTraceback) {
            sesResultRight = AlignWithTraceback(querySeq.Bases() + qStart, qSpan, tseq.c_str(),
                                                tSpan, dMax, bandwidth, sesScratch);
        } else {
            sesResultRight = AlignNoTraceback(querySeq.Bases() + qStart, qSpan, tseq.c_str(), tSpan,
                                              dMax, bandwidth, sesScratch);
        }

        ret->Aend = sesResultRight.lastQueryPos;
        ret->Bend = sesResultRight.lastTargetPos;
        ret->Aend += ovl->Astart;
        ret->Bend += ovl->Bstart;
        ret->EditDistance = sesResultRight.numDiffs;
        ret->Score = -(std::min(ret->ASpan(), ret->BSpan()) - ret->EditDistance);

#ifdef PANCAKE_DEBUG_ALN
        PBLOG_INFO << "dMax = " << dMax << ", bandwidth = " << bandwidth;
        PBLOG_INFO << "Right: diffs = " << sesResultRight.numDiffs;
        PBLOG_INFO << "After right: " << *ret;
#endif
    }

    ///////////////////////////
    /// Align reverse pass. ///
    ///////////////////////////
    {
        // Reverse query coordinates.
        const int32_t qStart = ret->Alen - ret->Astart;
        const int32_t qEnd = ret->Alen;
        const int32_t qSpan = qEnd - qStart;
        const int32_t tStartFwd = ret->Brev ? (ret->Blen - ret->Bend) : ret->Bstart;
        const int32_t tEndFwd = ret->Brev ? (ret->Blen - ret->Bstart) : ret->Bend;
        std::string tseq;
        if (ovl->Brev) {
            int32_t minHangLen = std::min(ovl->Blen - tEndFwd, qStart);
            int32_t extractBegin = tEndFwd;
            int32_t extractEnd = std::min(ret->Blen, tEndFwd + minHangLen * 2);
            tseq = FetchTargetSubsequence_(targetSeq, extractBegin, extractEnd, !ret->Brev);
        } else {
            int32_t minHangLen = std::min(ovl->Astart, tStartFwd);
            int32_t extractBegin = std::max(0, tStartFwd - minHangLen * 2);
            int32_t extractEnd = tStartFwd;
            tseq = FetchTargetSubsequence_(targetSeq, extractBegin, extractEnd, !ret->Brev);
        }
        const int32_t tSpan = tseq.size();
        const int32_t dMax = std::max(MIN_DIFFS_CAP, static_cast<int32_t>(ovl->Alen * alignMaxDiff -
                                                                          sesResultRight.numDiffs));
        const int32_t bandwidth =
            std::max(MIN_BANDWIDTH_CAP,
                     static_cast<int32_t>(std::min(ovl->Blen, ovl->Alen) * alignBandwidth));

        if (useTraceback) {
            sesResultLeft = AlignWithTraceback(reverseQuerySeq.c_str() + qStart, qSpan,
                                               tseq.c_str(), tSpan, dMax, bandwidth, sesScratch);
        } else {
            sesResultLeft = AlignNoTraceback(reverseQuerySeq.c_str() + qStart, qSpan, tseq.c_str(),
                                             tSpan, dMax, bandwidth, sesScratch);
        }

        ret->Astart = ovl->Astart - sesResultLeft.lastQueryPos;
        ret->Bstart = ovl->Bstart - sesResultLeft.lastTargetPos;
        std::reverse(sesResultLeft.cigar.begin(), sesResultLeft.cigar.end());
    }

    PacBio::Pancake::Alignment::DiffCounts diffs =
        sesResultRight.diffCounts + sesResultLeft.diffCounts;

    // Compute edit distance, identity and score.
    ret->EditDistance = -1;
    if (useTraceback) {
        diffs.Identity(noSNPs, noIndels, ret->Identity, ret->EditDistance);
        ret->Score = -diffs.numEq;
    } else {
        ret->EditDistance = sesResultLeft.numDiffs + sesResultRight.numDiffs;
        ret->Score = -(std::min(ret->ASpan(), ret->BSpan()) - ret->EditDistance);
        const float span = std::max(ret->ASpan(), ret->BSpan());
        ret->Identity =
            ((span != 0.0f) ? ((span - static_cast<float>(ret->EditDistance)) / span) : -0.0f);
    }

    // Merge the CIGAR strings.
    ret->Cigar = std::move(sesResultLeft.cigar);
    if (sesResultRight.cigar.size() > 0) {
        AppendToCigar(ret->Cigar, sesResultRight.cigar.front().Type(),
                      sesResultRight.cigar.front().Length());
        ret->Cigar.insert(ret->Cigar.end(), sesResultRight.cigar.begin() + 1,
                          sesResultRight.cigar.end());
    }

    if (trimAlignment && ret->Cigar.size() > 0) {
        PacBio::BAM::Cigar newCigar;
        TrimmingInfo trimInfo;
        TrimCigar(ret->Cigar, trimWindowSize, std::max(1.0, trimWindowSize * trimMatchFraction),
                  trimToFirstMatch, newCigar, trimInfo);
        ret->Astart += trimInfo.queryFront;
        ret->Bstart += trimInfo.targetFront;
        ret->Aend -= trimInfo.queryBack;
        ret->Bend -= trimInfo.targetBack;
        std::swap(ret->Cigar, newCigar);
    }

    NormalizeAndExtractVariantsInPlace_(ret, targetSeq, querySeq, reverseQuerySeq, noSNPs, noIndels,
                                        maskHomopolymers, maskSimpleRepeats, maskHomopolymerSNPs,
                                        maskHomopolymersArbitrary);

#ifdef PANCAKE_DEBUG_ALN
    PBLOG_INFO << "Final: " << *ret;
#endif

    return ret;
}

void Mapper::NormalizeAndExtractVariantsInPlace_(
    OverlapPtr& ovl, const PacBio::Pancake::FastaSequenceCached& targetSeq,
    const PacBio::Pancake::FastaSequenceCached& querySeq, const std::string /*reverseQuerySeq*/,
    bool noSNPs, bool noIndels, bool maskHomopolymers, bool maskSimpleRepeats,
    bool maskHomopolymerSNPs, bool maskHomopolymersArbitrary)
{
    // Extract the variant strings.
    if (ovl->Cigar.empty()) {
        return;
    }

    const char* Aseq = querySeq.Bases();
    const char* Bseq = targetSeq.Bases();
    int32_t Blen = targetSeq.Size();

    if (ovl->IsFlipped) {
        Aseq = targetSeq.Bases();
        Bseq = querySeq.Bases();
        Blen = querySeq.Size();
    }

    // / This works, but is slower because it requires to copy the target sequence every time.
    // / Since we already have the reversed query, we can simply reverse the CIGAR string and provide
    // / the reversed query, as below.
    PacBio::Pancake::Alignment::DiffCounts diffsPerBase;
    PacBio::Pancake::Alignment::DiffCounts diffsPerEvent;

    auto tseq = FetchTargetSubsequence_(Bseq, Blen, ovl->BstartFwd(), ovl->BendFwd(), ovl->Brev);

    const char* querySub = Aseq + ovl->Astart;
    int32_t querySubLen = ovl->ASpan();
    const char* targetSub = tseq.c_str();
    int32_t targetSubLen = tseq.size();
    auto& cigar = ovl->Cigar;

    // If the B-read is reversed, then we need to reverse and left-align the CIGAR string
    // to be consistent with the other alignments which were in forward.
    // The reason is that in the rest of the workflow up to this point we treated the
    // query as either fwd or reverse, but finally we actually take query as FWD in the output.
    // Normalization needs to happen in this context, so that all left-aligned indels line up.
    //    if (ovl->Brev) {
    // //        auto tseq =
    // //            FetchTargetSubsequence_(Bseq, Blen, ovl->BstartFwd(), ovl->BendFwd(), ovl->Brev);
    // //        cigar = NormalizeCigar(querySub, querySubLen, tseq.c_str(), targetSubLen, cigar);
    //        cigar = NormalizeCigar(querySub, querySubLen, targetSub, targetSubLen, cigar);
    //    } else {
    //        cigar = NormalizeCigar(querySub, querySubLen, targetSub, targetSubLen, cigar);
    //    }

    cigar = NormalizeCigar(querySub, querySubLen, targetSub, targetSubLen, cigar);

    ExtractVariantString(querySub, querySubLen, targetSub, targetSubLen, cigar, maskHomopolymers,
                         maskSimpleRepeats, maskHomopolymerSNPs, maskHomopolymersArbitrary,
                         ovl->Avars, ovl->Bvars, diffsPerBase, diffsPerEvent);

    const auto& diffs = diffsPerBase;
    diffs.Identity(noSNPs, noIndels, ovl->Identity, ovl->EditDistance);
    ovl->Score = -diffs.numEq;
}

// void NormalizeAndExtractVariantsInPlaceDeprecated_(
//     OverlapPtr& ovl, const PacBio::Pancake::FastaSequenceCached& targetSeq,
//     const PacBio::Pancake::FastaSequenceCached& querySeq, const std::string reverseQuerySeq,
//     bool noSNPs, bool noIndels, bool maskHomopolymers, bool maskSimpleRepeats)
// {
//     // Extract the variant strings.
//     if (ovl->Cigar.empty()) {
//         return;
//     }

//     /// This works, but is slower because it requires to copy the target sequence every time.
//     /// Since we already have the reversed query, we can simply reverse the CIGAR string and provide
//     /// the reversed query, as below.
//     // auto tseq = FetchTargetSubsequence_(targetSeq, ovl->BstartFwd(), ovl->BendFwd(), ovl->Brev);
//     // ovl->VariantString = ExtractVariantString(querySeq.Bases() + ovl->Astart, ovl->ASpan(), tseq.c_str(), tseq.size(),
//     //               ovl->Cigar, false, false);

//     PacBio::Pancake::Alignment::DiffCounts diffsPerBase;
//     PacBio::Pancake::Alignment::DiffCounts diffsPerEvent;

//     if (ovl->Brev) {
//         // Get the correct subsequences as C-strings.
//         const char* querySub = reverseQuerySeq.c_str() + (ovl->Alen - ovl->Aend);
//         int32_t querySubLen = ovl->ASpan();
//         const char* targetSub = targetSeq.Bases() + ovl->BstartFwd();
//         int32_t targetSubLen = ovl->BSpan();
//         // Reverse the CIGAR (it was generated earlier by reversing the target).
//         auto tempCigar = ovl->Cigar;
//         std::reverse(tempCigar.begin(), tempCigar.end());

//         // Normalize the CIGAR operations.
//         tempCigar = NormalizeCigar(querySub, querySubLen, targetSub, targetSubLen, tempCigar);
//         ovl->Cigar = tempCigar;
//         std::reverse(ovl->Cigar.begin(), ovl->Cigar.end());

//         // Extract variants for the updated CIGAR.
//         ExtractVariantString(querySub, querySubLen, targetSub, targetSubLen, tempCigar,
//                              maskHomopolymers, maskSimpleRepeats, ovl->Avars, ovl->Bvars,
//                              diffsPerBase, diffsPerEvent);

//         ovl->Avars = Pancake::ReverseComplement(ovl->Avars, 0, ovl->Avars.size());
//         ovl->Bvars = Pancake::ReverseComplement(ovl->Bvars, 0, ovl->Bvars.size());
//     } else {
//         // Get the correct subsequences as C-strings.
//         const char* querySub = querySeq.Bases() + ovl->Astart;
//         int32_t querySubLen = ovl->ASpan();
//         const char* targetSub = targetSeq.Bases() + ovl->BstartFwd();
//         int32_t targetSubLen = ovl->BSpan();
//         auto& tempCigar = ovl->Cigar;
//         // Normalize the CIGAR operations.
//         tempCigar = NormalizeCigar(querySub, querySubLen, targetSub, targetSubLen, tempCigar);
//         // Extract variants for the updated CIGAR.
//         ExtractVariantString(querySub, querySubLen, targetSub, targetSubLen, tempCigar,
//                              maskHomopolymers, maskSimpleRepeats, ovl->Avars, ovl->Bvars,
//                              diffsPerBase, diffsPerEvent);
//     }

//     const auto& diffs = diffsPerBase;
//     diffs.Identity(noSNPs, noIndels, ovl->Identity, ovl->EditDistance);
//     ovl->Score = -diffs.numEq;
// }

void Mapper::DebugWriteSeedHits_(const std::string& outPath, const std::vector<SeedHit>& hits,
                                 int32_t seedLen, const std::string& queryName, int64_t queryLen,
                                 const std::string& targetName, int64_t targetLen)
{
    std::ofstream ofs(outPath);
    // Simply walk away if the file cannot be open.
    // Avoid writing debug output if a specific path is not available.
    if (ofs.is_open() == false) {
        return;
    }
    ofs << queryName << "\t0\t" << queryLen << "\t" << targetName << "\t0\t" << targetLen << "\t0.0"
        << std::endl;
    for (size_t i = 0; i < hits.size(); ++i) {
        int32_t clusterId = hits[i].targetId * 2 + (hits[i].targetRev ? 1 : 0);
        ofs << hits[i].queryPos << "\t" << hits[i].targetPos << "\t"
            << "\t" << clusterId << std::endl;
        ofs << hits[i].queryPos + seedLen << "\t" << hits[i].targetPos + seedLen << "\t"
            << clusterId << std::endl;
    }
}

std::vector<MapperResult> MapHiFi(const std::vector<std::string>& targetSeqs,
                                  const std::vector<std::string>& querySeqs,
                                  const PacBio::Pancake::SeedDB::SeedDBParameters& seedParams,
                                  const OverlapHifiSettings& settings)
{
    PacBio::Pancake::FastaSequenceCachedStore targetSeqsCached;
    for (int32_t i = 0; i < static_cast<int32_t>(targetSeqs.size()); ++i) {
        targetSeqsCached.AddRecord(
            FastaSequenceCached(std::to_string(i), targetSeqs[i].c_str(), targetSeqs[i].size(), i));
    }

    PacBio::Pancake::FastaSequenceCachedStore querySeqsCached;
    for (int32_t i = 0; i < static_cast<int32_t>(querySeqs.size()); ++i) {
        querySeqsCached.AddRecord(
            FastaSequenceCached(std::to_string(i), querySeqs[i].c_str(), querySeqs[i].size(), i));
    }

    return MapHiFi(targetSeqsCached, querySeqsCached, seedParams, settings);
}

std::vector<MapperResult> MapHiFi(const FastaSequenceCachedStore& targetSeqs,
                                  const FastaSequenceCachedStore& querySeqs,
                                  const PacBio::Pancake::SeedDB::SeedDBParameters& seedParams,
                                  const OverlapHifiSettings& settings)
{
    const bool generateFlippedOverlap = settings.WriteReverseOverlaps;

    // Construct the target sequence index.
    std::vector<PacBio::Pancake::Int128t> seeds;
    PacBio::Pancake::SeedDB::GenerateMinimizers(seeds, targetSeqs.records(), seedParams.KmerSize,
                                                seedParams.MinimizerWindow, seedParams.Spacing,
                                                seedParams.UseRC, seedParams.UseHPCForSeedsOnly);
    std::unique_ptr<SeedIndex> seedIndex = std::make_unique<SeedIndex>(std::move(seeds));

    // Seed statistics, and computing the cutoff.
    TicToc ttSeedStats;
    int64_t freqMax = 0;
    int64_t freqCutoff = 0;
    double freqAvg = 0.0;
    double freqMedian = 0.0;
    seedIndex->ComputeFrequencyStats(settings.FreqPercentile, freqMax, freqAvg, freqMedian,
                                     freqCutoff);
    ttSeedStats.Stop();
#ifdef PANCAKE_MAP_HIFI_DEBUG
    PBLOG_INFO << "Computed the seed frequency statistics in " << ttSeedStats.GetSecs() << " sec.";
#endif

    // Create the overlapper.
    OverlapHiFi::Mapper mapper(settings);

    // Run mapping for each query.
    std::vector<MapperResult> results;
    for (int32_t i = 0; i < static_cast<int32_t>(querySeqs.records().size()); ++i) {
        const auto& query = querySeqs.records()[i];
        int32_t queryId = query.Id();

#ifdef PANCAKE_MAP_HIFI_DEBUG
        PBLOG_INFO << "[i = " << i << "] New query: queryId = " << queryId
                   << ", length = " << query.size();
#endif

        // Compute the query seeds.
        std::vector<PacBio::Pancake::Int128t> querySeeds;
        int32_t seqLen = query.size();
        const uint8_t* seq = reinterpret_cast<const uint8_t*>(query.data());
        int rv = SeedDB::GenerateMinimizers(
            querySeeds, seq, seqLen, 0, queryId, seedParams.KmerSize, seedParams.MinimizerWindow,
            seedParams.Spacing, seedParams.UseRC, seedParams.UseHPCForSeedsOnly);
        if (rv) {
            throw std::runtime_error("Generating minimizers failed for the query sequence i = " +
                                     std::to_string(i) + ", id = " + std::to_string(queryId));
        }
        SequenceSeedsCached querySeedsCached("", &querySeeds[0], querySeeds.size(), queryId);

        // Map/align.
        MapperResult queryResults = mapper.Map(targetSeqs, *seedIndex, query, querySeedsCached,
                                               freqCutoff, generateFlippedOverlap);

        // for (const auto& m : queryResults.mappings) {
        //     m->mapping->Aid = queryId;
        // }

        // Collect.
        results.emplace_back(std::move(queryResults));

#ifdef PANCAKE_MAP_HIFI_DEBUG
        PBLOG_INFO << "\n\n\n";
#endif
    }
    return results;
}

}  // namespace OverlapHiFi
}  // namespace Pancake
}  // namespace PacBio
