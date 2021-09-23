/*
 * DPChain.cpp
 *
 *  Created on: Dec 16, 2017
 *      Author: Ivan Sovic
 *
 * Originally implemented in the Raptor graph-based mapper.
 */

#include <pacbio/pancake/DPChain.h>
#include <iostream>
#include <lib/math.hpp>
#include <sstream>

// #define DPCHAIN_DEBUG

namespace PacBio {
namespace Pancake {

constexpr int32_t PlusInf = std::numeric_limits<int32_t>::max() - 10000;  // Leave a margin.

/*
 * \brief Performs the forward pass of the DP chaining procedure (no backtrack included).
 *          The DP results are returned via parameters: dp, pred and chainId.
*/
int32_t ChainHitsForward(const SeedHit* hits, const int32_t hitsSize, const int32_t chainMaxSkip,
                         const int32_t chainMaxPredecessors, const int32_t seedJoinDist,
                         const int32_t diagMargin, std::vector<int32_t>& dp,
                         std::vector<int32_t>& pred, std::vector<int32_t>& chainId)
{
    dp.clear();
    pred.clear();
    chainId.clear();

    // Zeroth element will be the "Null" state.
    // For each node, it's chain ID is the same as of it's predecessor.
    dp.resize(hitsSize + 1, 0);
    pred.resize(hitsSize + 1, 0);
    chainId.resize(hitsSize + 1, -1);

    int32_t numChains = 0;

    double avgQuerySpan = 0.0;
    for (int32_t i = 0; i < hitsSize; i++) {
        avgQuerySpan += static_cast<double>(hits[i].querySpan);
    }
    avgQuerySpan = (hitsSize > 0) ? avgQuerySpan / static_cast<double>(hitsSize) : 0.0;

    const double linFactor = 0.01 * avgQuerySpan;

    for (int32_t i = 1; i < (hitsSize + 1); i++) {
        // The [i - 1] is correct here. This is because the "i" loop starts at 1,
        // just to make sure there are at least 2 points.
        const int32_t x_i_start = hits[i - 1].queryPos;
        const int32_t y_i_start = hits[i - 1].targetPos;
        const int32_t t_id_i = hits[i - 1].targetId;
        const bool t_rev_i = hits[i - 1].targetRev;

        // Add the initial gap open penalty.
        const int32_t x_i_score = hits[i - 1].querySpan;

        int32_t newDpVal = x_i_score;
        int32_t newDpPred = 0;
        int32_t newDpChain = numChains;
        int32_t numSkippedPredecessors = 0;
        int32_t numProcessed = 0;

        const int32_t minJ =
            (chainMaxPredecessors <= 0) ? 0 : std::max(0, (i - 1 - chainMaxPredecessors));

        for (int32_t j = (i - 1); j > minJ; j--) {
            // The [j - 1] is correct here. The outter loop looks at hits[i-1], and this one has to go
            // one before that.
            const int32_t x_j_start = hits[j - 1].queryPos;
            const int32_t y_j_start = hits[j - 1].targetPos;
            const int32_t t_id_j = hits[j - 1].targetId;
            const bool t_rev_j = hits[j - 1].targetRev;
            const int32_t x_j_span = hits[j - 1].querySpan;

            const int32_t distX = x_i_start - x_j_start;  // If < 0, it's not a predecessor.
            const int32_t distY = y_i_start - y_j_start;

            const int32_t gapDist = (distX < distY) ? (distY - distX) : (distX - distY);

            if (t_id_j != t_id_i || t_rev_j != t_rev_i) {
                break;
            }
            if (distY > seedJoinDist) {
                break;
            }

            if (x_i_start <= x_j_start || y_i_start <= y_j_start) {
                continue;
            }
            if (gapDist > diagMargin) {
                continue;
            }
            if (distX > seedJoinDist) {
                continue;
            }

            numProcessed += 1;

            const int32_t linPart = (gapDist * linFactor);
            const int32_t logPart = ((gapDist == 0) ? 0 : raptor::utility::ilog2_32(gapDist));
            const int32_t edge_score = linPart + (logPart / 2.0);
            const int32_t x_j_score = std::min(x_j_span, std::min<int32_t>(abs(distX), abs(distY)));
            const int32_t score_ij = dp[j] + x_j_score - edge_score;

            if (score_ij >= newDpVal) {
                newDpPred = j;
                newDpVal = score_ij;
                newDpChain = chainId[j];

                // This is the main difference to how I previously calculated the scan_depth.
                numSkippedPredecessors -= 1;
                numSkippedPredecessors = std::max(0, numSkippedPredecessors);

            } else {
                numSkippedPredecessors += 1;
                if (numSkippedPredecessors > chainMaxSkip) {
                    break;
                }
            }
        }

        dp[i] = newDpVal;
        pred[i] = newDpPred;
        chainId[i] = newDpChain;
        if (newDpChain == numChains) {
            numChains += 1;
        }
    }

    return numChains;
}

std::vector<ChainedHits> ChainHits(const SeedHit* hits, int32_t hitsSize, int32_t chainMaxSkip,
                                   int32_t chainMaxPredecessors, int32_t seedJoinDist,
                                   int32_t diagMargin, int32_t minNumSeeds, int32_t minCovBases,
                                   int32_t minDPScore)
{
    /**
     * Hits need to be sorted in this order of priority:
     *      target_id, target_rev, target_pos, query_pos
    */

    std::vector<ChainedHits> chains;
    const int32_t n_hits = hitsSize;

    if (n_hits == 0) {
        return chains;
    }

    if (chainMaxSkip <= 0) {
        return chains;
    }

    std::vector<int32_t> dp;
    std::vector<int32_t> pred;
    std::vector<int32_t> chainId;

    const int32_t numChains = ChainHitsForward(hits, hitsSize, chainMaxSkip, chainMaxPredecessors,
                                               seedJoinDist, diagMargin, dp, pred, chainId);

    ////////////////////
    /// Backtrack.   ///
    ////////////////////
    // Find the maximum of every chain for backtracking.
    std::vector<int32_t> chainMaxima(numChains, -PlusInf);
    for (int32_t i = 1; i < (n_hits + 1); i++) {
        if (chainMaxima[chainId[i]] == -PlusInf || dp[i] >= dp[chainMaxima[chainId[i]]]) {
            chainMaxima[chainId[i]] = i;
        }
    }

    for (int32_t i = 0; i < static_cast<int32_t>(chainMaxima.size()); i++) {
        // Trace back from the maxima.
        int32_t node_id = chainMaxima[i];
        const int32_t score = dp[node_id];

        if (score < minDPScore) {
            continue;
        }

        std::vector<int32_t> nodes;
        while (node_id > 0) {
            nodes.emplace_back(node_id - 1);  // The "- 1" is because of the DP offset.
            node_id = pred[node_id];
        }
        // Reverse the backtracked nodes.
        std::reverse(nodes.begin(), nodes.end());

        // Skip if needed.
        if (nodes.empty() || static_cast<int32_t>(nodes.size()) < minNumSeeds) {
            continue;
        }

        /////////////////////////
        /// Create the chain. ///
        /////////////////////////
        ChainedHits chain;
        const int32_t currTargetId = hits[nodes.front()].targetId;
        const bool currTargetRev = hits[nodes.front()].targetRev;
        if (chain.targetId == -1 || chain.targetId != currTargetId ||
            chain.targetRev != currTargetRev) {
            chain = ChainedHits(currTargetId, currTargetRev);
        }

        for (const auto& node : nodes) {
            chain.hits.emplace_back(hits[node]);
        }

        // Penalize the distance from the end of the query.
        // Otherwise, shorted chains near the beginning would
        // prevail longer ones in some cases.
        // int32_t chain_dist_to_end = qseq.get_sequence_length() - chain->hits().back().QueryPos();
        // int32_t chain_score = score - chain_dist_to_end * params->chain_penalty_gap;
        // chain->score(chain_score);
        chain.score = score;

        CalcHitCoverage(chain.hits, 0, chain.hits.size(), chain.coveredBasesQuery,
                        chain.coveredBasesTarget);

        // int32_t qspan = chain.hits.back().queryPos - chain.hits.front().queryPos;
        // double frac = (qspan == 0) ? 0 : ((double)chain.coveredBasesQuery) / ((double)qspan);

        // Add the new chain.
        if (chain.coveredBasesQuery >= minCovBases && chain.coveredBasesTarget >= minCovBases) {
            chains.emplace_back(std::move(chain));
        }
        /////////////////////////
    }

#ifdef DEBUG_DP_VERBOSE_
    printf("The DP:\n");
    for (int32_t i = 0; i < dp.size(); i++) {
        printf("[%d] dp[i] = %d, pred[i] = %d, chainId[i] = %d\n", i, dp[i], pred[i], chainId[i]);
    }
#endif

    return chains;
}

inline int32_t ComputeGap(const SeedHit& hitStart, const SeedHit& hitEnd)
{
    // return (hitEnd.targetPos - hitStart.targetPos) - (hitEnd.queryPos - hitStart.queryPos);
    return (hitEnd.queryPos - hitStart.queryPos) - (hitEnd.targetPos - hitStart.targetPos);
}

double ComputeChainDivergence(const std::vector<SeedHit>& hits)
{
    /*
     * Measures the sum of absolute gaps between hits.
     * The lower the divergence, the closer the hits are to a diagonal.
    */
    double divergence = 0.0;
    if (hits.empty()) {
        return divergence;
    }

    for (int32_t i = 1; i < static_cast<int32_t>(hits.size()); ++i) {
        const int32_t gap = ComputeGap(hits[i - 1], hits[i]);
        divergence += std::abs(static_cast<double>(gap));
    }

    return divergence;
}

std::vector<int32_t> CollectLongGaps(const std::vector<SeedHit>& hits, int32_t minGap)
{
    /*
        * Finds hits which have a gap > minGap from the previous hit.
    */
    std::vector<int32_t> breakpoints;
    for (int32_t i = 1; i < static_cast<int32_t>(hits.size()); ++i) {
        const int32_t gap = ComputeGap(hits[i - 1], hits[i]);
        if (std::abs(gap) > minGap) {
            breakpoints.emplace_back(i);
        }
    }
    return breakpoints;
}

ChainedHits RefineChainedHits(const ChainedHits& chain, int32_t minGap, int32_t diffThreshold,
                              int32_t maxForwardSeedDist, int32_t maxForwardSeedCount)
{
    /*
     * \param minGap Minimum gap distsance between two seeds to mark it as a breakpoint.
     * \param maxForwardSeedDist Stop the inner for loop if two seeds are more than this apart.
     * \param maxForwardSeedCount Heuristic value to limit the number of succesive seeds to be checked in a for loop.
    */

    const int8_t FLAG_MASK_IGNORE = (1 << 0);
    const auto& hits = chain.hits;
    const std::vector<int32_t> breakpoints = CollectLongGaps(hits, minGap);
    const int32_t nBreakpoints = breakpoints.size();

    int32_t maxVal = 0, maxStart = -1, maxEnd = -1;
    std::vector<int8_t> flags(hits.size(), 0);

    for (int32_t bpId = 0; bpId <= nBreakpoints; ++bpId) {
        // Either filter out a maximum-gap span, or break.
        if (bpId == nBreakpoints || bpId >= maxEnd) {
            if (maxEnd > 0) {
                for (int32_t j = breakpoints[maxStart]; j < breakpoints[maxEnd]; ++j) {
                    flags[j] |= FLAG_MASK_IGNORE;
                }
            }
            maxVal = 0;
            maxStart = -1;
            maxEnd = -1;
            if (bpId == nBreakpoints) {
                break;
            }
        }

        const int32_t seedId = breakpoints[bpId];
        const int32_t gap = ComputeGap(hits[seedId - 1], hits[seedId]);
        int32_t nIns = (gap > 0) ? gap : 0;
        int32_t nDel = (gap > 0) ? 0 : -gap;
        const int32_t qStart = hits[seedId - 1].queryPos;
        const int32_t tStart = hits[seedId - 1].targetPos;
        int32_t maxDiff = 0;
        int32_t maxDiffBpId = -1;

        // Go through next breakpoints and check how concordant they are.
        for (int32_t nextBpId = (bpId + 1);
             nextBpId < nBreakpoints && nextBpId < (bpId + maxForwardSeedCount); ++nextBpId) {
            const int32_t nextSeedId = breakpoints[nextBpId];
            // Check that the next breakpoint is not too far away.
            if ((hits[nextSeedId].queryPos - qStart) > maxForwardSeedDist ||
                (hits[nextSeedId].targetPos - tStart) > maxForwardSeedDist) {
                break;
            }
            const int32_t nextGap = ComputeGap(hits[nextSeedId - 1], hits[nextSeedId]);
            if (nextGap > 0) {
                nIns += nextGap;
            } else {
                nDel += (-nextGap);
            }
            const int32_t diff = nIns + nDel - std::abs(nIns - nDel);
            // const int32_t diff = std::abs(nIns - nDel);
            if (diff > maxDiff) {
                maxDiff = diff;
                maxDiffBpId = nextBpId;
            }
        }
        if (maxDiff > diffThreshold && maxDiff > maxVal) {
            maxVal = maxDiff;
            maxStart = bpId;
            maxEnd = maxDiffBpId;
        }
    }

    // Construct a set of filtered hits.
    ChainedHits ret;
    for (size_t i = 0; i < hits.size(); ++i) {
        if (flags[i] != 0) {
            continue;
        }
        ret.hits.emplace_back(hits[i]);
    }
    ret.score = chain.score;
    ret.targetId = chain.targetId;
    ret.targetRev = chain.targetRev;
    CalcHitCoverage(ret.hits, 0, ret.hits.size(), ret.coveredBasesQuery, ret.coveredBasesTarget);

    return ret;
}

ChainedHits RefineChainedHits2(const ChainedHits& chain, int32_t minGap, int32_t maxForwardSeedDist)
{
    /*
    * This filters more extreme outliers.
    * For every breakpoint, we check the succeeding breakpoints and compute a value
    * m (approximate number of matches computed as the min(target_dist, query_dist) between the two breakpoints,
    * and the sum of gaps (gap1 which is the leading gap into the first breakpoint, and gap2 which is the
    * gap preceding the next breakpoing).
    * If the "number of matches" is lower than the sum of gaps, then this is a candidate for filtering.
    */

    const int8_t FLAG_MASK_IGNORE = (1 << 0);
    const int8_t FLAG_MASK_LONG_JOIN = (1 << 1);

    const auto& hits = chain.hits;
    const std::vector<int32_t> breakpoints = CollectLongGaps(hits, minGap);
    const int32_t nBreakpoints = breakpoints.size();
    if (nBreakpoints == 0) {
        return chain;
    }

    std::vector<int8_t> flags(hits.size(), 0);
    for (int32_t bpId = 0; bpId < nBreakpoints;) {
        const int32_t seedId = breakpoints[bpId];
        int32_t gap1 = std::abs(ComputeGap(hits[seedId - 1], hits[seedId]));
        int32_t qEnd1 = hits[seedId].queryPos + static_cast<int32_t>(hits[seedId].querySpan);
        int32_t tEnd1 = hits[seedId].targetPos + static_cast<int32_t>(hits[seedId].targetSpan);

        int32_t nextBpId = 0;
        for (nextBpId = (bpId + 1); nextBpId < nBreakpoints; ++nextBpId) {
            const int32_t nextSeedId = breakpoints[nextBpId];
            // Check that the next breakpoint is not too far away.
            if ((hits[nextSeedId].queryPos - qEnd1) > maxForwardSeedDist ||
                (hits[nextSeedId].targetPos - tEnd1) > maxForwardSeedDist) {
                break;
            }
            // Compute the difference in (qspan-tspan).
            const int32_t gap2 = std::abs(ComputeGap(hits[nextSeedId - 1], hits[nextSeedId]));
            // Start of the previous gap.
            const int32_t qStart2 = hits[nextSeedId - 1].queryPos;
            const int32_t tStart2 = hits[nextSeedId - 1].targetPos;
            // Compute approximate number of matches from the first breakpoint to here.
            const int32_t m = std::min(tStart2 - tEnd1, qStart2 - qEnd1);
            if (m > (gap1 + gap2)) {
                // If there are more matches than the total gap, stop going forward.
                break;
            }
            qEnd1 = hits[nextSeedId].queryPos + static_cast<int32_t>(hits[nextSeedId].querySpan);
            tEnd1 = hits[nextSeedId].targetPos + static_cast<int32_t>(hits[nextSeedId].targetSpan);
            gap1 = gap2;
        }
        if (nextBpId > (bpId + 1)) {
            for (int32_t j = breakpoints[bpId]; j < breakpoints[nextBpId - 1]; ++j) {
                flags[j] |= FLAG_MASK_IGNORE;
            }
            flags[breakpoints[nextBpId - 1]] |= FLAG_MASK_LONG_JOIN;
        }
        bpId = nextBpId;
    }

    // Construct a set of filtered hits.
    ChainedHits ret;
    for (size_t i = 0; i < hits.size(); ++i) {
        if ((flags[i] & FLAG_MASK_IGNORE) != 0) {
            continue;
        }
        ret.hits.emplace_back(hits[i]);
        if ((flags[i] & FLAG_MASK_LONG_JOIN) != 0) {
            ret.hits.back().SetFlagLongJoin();
        }
    }
    ret.score = chain.score;
    ret.targetId = chain.targetId;
    ret.targetRev = chain.targetRev;
    CalcHitCoverage(ret.hits, 0, ret.hits.size(), ret.coveredBasesQuery, ret.coveredBasesTarget);

    return ret;
}

ChainedHits RefineBadEnds(const ChainedHits& chain, int32_t bandwidth, int32_t minMatch)
{
    if (chain.hits.size() < 3) {
        return chain;
    }
    // const int32_t minCoveredBases = std::min(chain.coveredBasesQuery, chain.coveredBasesTarget);
    const int32_t minCoveredBases = chain.coveredBasesQuery;

    int32_t start = 0;
    {
        int32_t numMatches = chain.hits[0].querySpan;
        int32_t totalSpan = chain.hits[0].querySpan;
        for (int32_t i = 1; i < (static_cast<int32_t>(chain.hits.size()) - 1); ++i) {
            if (chain.hits[i].CheckFlagLongJoin()) {
                break;
            }
            const int32_t qDist = chain.hits[i].queryPos - chain.hits[i - 1].queryPos;
            const int32_t tDist = chain.hits[i].targetPos - chain.hits[i - 1].targetPos;
            const int32_t minDist = std::min(qDist, tDist);
            const int32_t maxDist = std::max(qDist, tDist);
            const int32_t gap = maxDist - minDist;
            if (gap > (totalSpan >> 1)) {
                start = i;
            }
            totalSpan += minDist;
            const int32_t qSpan = chain.hits[i].querySpan;
            numMatches += std::min(minDist, qSpan);
            if (totalSpan >= (bandwidth << 1) ||
                (numMatches >= minMatch && numMatches >= bandwidth) ||
                numMatches >= (minCoveredBases >> 1)) {
                break;
            }
        }
    }

    int32_t end = chain.hits.size();
    {
        int32_t numMatches = chain.hits.back().querySpan;
        int32_t totalSpan = chain.hits.back().querySpan;
        for (int32_t i = (static_cast<int32_t>(chain.hits.size()) - 2); i > start; --i) {
            if (chain.hits[i + 1].CheckFlagLongJoin()) {
                break;
            }
            const int32_t qDist = chain.hits[i + 1].queryPos - chain.hits[i].queryPos;
            const int32_t tDist = chain.hits[i + 1].targetPos - chain.hits[i].targetPos;
            const int32_t minDist = std::min(qDist, tDist);
            const int32_t maxDist = std::max(qDist, tDist);
            const int32_t gap = maxDist - minDist;
            if (gap > (totalSpan >> 1)) {
                end = i + 1;
            }
            totalSpan += minDist;
            const int32_t qSpan = chain.hits[i + 1].querySpan;
            numMatches += std::min(minDist, qSpan);
            if (totalSpan >= (bandwidth << 1) ||
                (numMatches >= minMatch && numMatches >= bandwidth) ||
                numMatches >= (minCoveredBases >> 1)) {
                break;
            }
        }
    }

    // Construct a set of filtered hits.
    ChainedHits ret;
    ret.hits.insert(ret.hits.end(), chain.hits.begin() + start, chain.hits.begin() + end);
    ret.score = chain.score;
    ret.targetId = chain.targetId;
    ret.targetRev = chain.targetRev;
    CalcHitCoverage(ret.hits, 0, ret.hits.size(), ret.coveredBasesQuery, ret.coveredBasesTarget);
    return ret;
}

std::vector<Range> GroupByTargetAndStrand(const std::vector<SeedHit>& sortedHits)
{
    if (sortedHits.empty()) {
        return {};
    }
    const int32_t numHits = static_cast<int32_t>(sortedHits.size());
    int32_t beginId = 0;
    std::vector<Range> groups;
    for (int32_t i = 0; i < numHits; ++i) {
        const auto& prevHit = sortedHits[beginId];
        const auto& currHit = sortedHits[i];
        if (currHit.targetId != prevHit.targetId || currHit.targetRev != prevHit.targetRev) {
            groups.emplace_back(Range{beginId, i});
            beginId = i;
        }
    }
    if ((numHits - beginId) > 0) {
        groups.emplace_back(Range{beginId, numHits});
    }
    return groups;
}

std::vector<Range> DiagonalGroup(const std::vector<SeedHit>& sortedHits, int32_t chainBandwidth,
                                 bool overlappingWindows)
{
    /*
     * Groups seed hits by:
     *      - Target ID.
     *      - Target strand.
     *      - Diagonal, groupping a current hit with previous by comparing the diagonal with the first diagonal in the current range.
     *      - Optionally, hits that could fall within neighboring windows will be groupped in both. This doesn't have to happen, diagonals can have clean cuts.
     *
     * For example, if there are 2 diagonals, one at 1000bp and the other at 1700bp (with chainBandwidth = 500 and overlappingWindows = true), this will result
     * in two ranges.
     * If there are 3 diagonals: (1) 1000bp, (2) 1300bp, (3) 1700bp; then there will be two ranges, and the middle diagonal will be included in both.
    */

    if (sortedHits.empty()) {
        return {};
    }

    std::vector<Range> groups;

    const int32_t numHits = static_cast<int32_t>(sortedHits.size());
    int32_t beginId = 0;
    int32_t beginDiag = sortedHits[beginId].Diagonal();

    // This is a combination of <targetPos, queryPos>, intended for simple comparison
    // without defining a custom comparison operator.
    uint64_t minTargetQueryPosCombo = (static_cast<uint64_t>(sortedHits[beginId].targetPos) << 32) |
                                      (static_cast<uint64_t>(sortedHits[beginId].queryPos));
    uint64_t maxTargetQueryPosCombo = minTargetQueryPosCombo;

    int32_t firstInBandwidth = 0;

    for (int32_t i = 0; i < numHits; ++i) {
        const auto& beginHit = sortedHits[beginId];
        const auto& currHit = sortedHits[i];
        const int32_t currDiag = currHit.Diagonal();
        const int32_t diagDiff = abs(currDiag - beginDiag);
        const uint64_t targetQueryPosCombo =
            (static_cast<uint64_t>(sortedHits[i].targetPos) << 32) |
            (static_cast<uint64_t>(sortedHits[i].queryPos));

        if (currHit.targetId != beginHit.targetId || currHit.targetRev != beginHit.targetRev ||
            diagDiff > chainBandwidth) {

            if (overlappingWindows) {
                groups.emplace_back(Range{firstInBandwidth, i});
            } else {
                groups.emplace_back(Range{beginId, i});
            }

            beginId = i;
            beginDiag = currDiag;

            minTargetQueryPosCombo = maxTargetQueryPosCombo = targetQueryPosCombo;

            // Find the earliest hit which is within the bandwidth window from the current hit.
            if (overlappingWindows) {
                for (; firstInBandwidth < i; ++firstInBandwidth) {
                    const auto& firstHit = sortedHits[firstInBandwidth];
                    const int32_t firstDiag = firstHit.Diagonal();
                    const int32_t diagDiffToFirst = abs(firstDiag - beginDiag);
                    if (currHit.targetId != firstHit.targetId ||
                        currHit.targetRev != firstHit.targetRev ||
                        diagDiffToFirst > chainBandwidth) {
                        continue;
                    }
                    break;
                }
            }
        }

        // Track the minimum and maximum target positions for each diagonal.
        if (targetQueryPosCombo < minTargetQueryPosCombo) {
            minTargetQueryPosCombo = targetQueryPosCombo;
        }
        if (targetQueryPosCombo > maxTargetQueryPosCombo) {
            maxTargetQueryPosCombo = targetQueryPosCombo;
        }
    }

    if ((numHits - beginId) > 0) {
        groups.emplace_back(Range{beginId, numHits});
    }

    return groups;
}

}  // namespace Pancake
}  // namespace PacBio
