/*
 * DPChain.h
 *
 *  Created on: Dec 16, 2017
 *      Author: Ivan Sovic
 *
 * Originally implemented in the Raptor graph-based mapper.
 */

#ifndef PANCAKE_DP_CHAIN_HPP
#define PANCAKE_DP_CHAIN_HPP

#include <pancake/Range.hpp>
#include <pancake/SeedHit.hpp>

#include <cstdint>
#include <memory>
#include <span>
#include <vector>
#ifdef PANCAKE_USE_SSE41
#include <emmintrin.h>
#endif

namespace PacBio {
namespace Pancake {

struct ChainingScratchSpace
{
    std::vector<int32_t> dp;
    std::vector<int32_t> pred;
    std::vector<int32_t> chainId;

#ifdef PANCAKE_USE_SSE41
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"
    std::vector<__m128i> dpSimd;
    std::vector<__m128i> predSimd;

    std::vector<__m128i> qp;             // Query pos.
    std::vector<__m128i> tp;             // Target pos.
    std::vector<__m128i> qs;             // Query span.
    std::vector<__m128i> tid;            // Target ID and strand.
    std::vector<__m128i> vectorIndices;  // Target ID and strand.

#pragma GCC diagnostic pop
#endif
};

struct ChainedHits
{
    int32_t targetId = -1;
    bool targetRev = false;
    std::vector<SeedHit> hits;
    int32_t score = 0;
    int32_t coveredBasesQuery = 0;
    int32_t coveredBasesTarget = 0;

    ChainedHits() = default;
    ChainedHits(int32_t _targetId, bool _targetRev) : targetId(_targetId), targetRev(_targetRev) {}
    ChainedHits(int32_t _targetId, bool _targetRev, const std::vector<SeedHit>& _hits,
                int32_t _score, int32_t _coveredBasesQuery, int32_t _coveredBasesTarget)
        : targetId(_targetId)
        , targetRev(_targetRev)
        , hits(_hits)
        , score(_score)
        , coveredBasesQuery(_coveredBasesQuery)
        , coveredBasesTarget(_coveredBasesTarget)
    {}
};

inline bool operator==(const ChainedHits& lhs, const ChainedHits& rhs)
{
    return std::tie(lhs.targetId, lhs.targetRev, lhs.hits, lhs.score, lhs.coveredBasesQuery,
                    lhs.coveredBasesTarget) == std::tie(rhs.targetId, rhs.targetRev, rhs.hits,
                                                        rhs.score, rhs.coveredBasesQuery,
                                                        rhs.coveredBasesTarget);
}

inline std::ostream& operator<<(std::ostream& os, const ChainedHits& b)
{
    os << "targetId = " << b.targetId << ", targetRev = " << (b.targetRev ? "true" : "false")
       << ", score = " << b.score << ", covBasesQuery = " << b.coveredBasesQuery
       << ", covBasesTarget = " << b.coveredBasesTarget << ", hits:\n";
    for (size_t i = 0; i < b.hits.size(); ++i) {
        os << "  [hit " << i << "] " << b.hits[i] << "\n";
    }
    return os;
}

inline int32_t SisdCompareLte(const int32_t a, const int32_t b)
{
    /**
     * Returns 0xFFFFFFFF if a <= b, 0x00000000 otherwise.
    */
    return (a - b - 1) >> 31;  // a <= b
}
inline int32_t SisdCompareLt(const int32_t a, const int32_t b)
{
    /**
     * Returns 0xFFFFFFFF if a < b, 0x00000000 otherwise.
    */
    return (a - b) >> 31;  // a < b
}
inline int32_t SisdCompareGte(const int32_t a, const int32_t b)
{
    /**
     * Returns 0xFFFFFFFF if a >= b, 0x00000000 otherwise.
    */
    return (b - a - 1) >> 31;  // a >= b
}
inline int32_t SisdCompareGt(const int32_t a, const int32_t b)
{
    /**
     * Returns 0xFFFFFFFF if a > b, 0x00000000 otherwise.
    */
    return (b - a) >> 31;  // a > b
}

inline uint32_t ilog2_32_clz_special_zero(const int32_t v)
{
    /**
     * \brief This is a modified integer log2 function, it returns 0 when v == 0.
     *          Of course, a proper mathematical log2 would return -inf, but this
     *          special case is required for chaining.
    */
    const int32_t vNonZero = std::max(v, 1);
    return (31 - __builtin_clz(vNonZero));
}

int32_t ChainHitsForward(std::span<const SeedHit> hits, int32_t chainMaxSkip,
                         int32_t chainMaxPredecessors, int32_t seedJoinDist, int32_t diagMargin,
                         std::vector<int32_t>& dp, std::vector<int32_t>& pred,
                         std::vector<int32_t>& chainId);

int32_t ChainHitsForwardFastSisd(std::span<const SeedHit> hits, int32_t chainMaxSkip,
                                 int32_t chainMaxPredecessors, int32_t seedJoinDist,
                                 int32_t diagMargin, std::vector<int32_t>& dp,
                                 std::vector<int32_t>& pred, std::vector<int32_t>& chainId);

std::vector<ChainedHits> ChainHitsBacktrack(std::span<const SeedHit> hits,
                                            std::span<const int32_t> dp,
                                            std::span<const int32_t> pred,
                                            std::span<const int32_t> chainId, int32_t numChains,
                                            int32_t minNumSeeds, int32_t minCovBases,
                                            int32_t minDPScore);

std::vector<ChainedHits> ChainHitsSisd(std::span<const SeedHit> hits, int32_t chainMaxSkip,
                                       int32_t chainMaxPredecessors, int32_t seedJoinDist,
                                       int32_t diagMargin, int32_t minNumSeeds, int32_t minCovBases,
                                       int32_t minDPScore, double& timeChaining,
                                       double& timeBacktrack,
                                       std::shared_ptr<ChainingScratchSpace> ss = nullptr);

double ComputeChainDivergence(const std::vector<SeedHit>& hits);

/**
 * @brief Removes stretches of chained hits which cover a region of a lot of local indel diffs.
 *          It begins by determining a list of potential breakpoints - seed hits which follow
 *          gaps longer than minGap.
 *          For every breakpoint Bi, it looks at the next maxForwardSeedCount breakpoints Bj (or until
 *          query or target distance between Bi and Bj is > maxForwardSeedDist) and computes
 *          the estimated gap as the "diff" (min(numIns, numDel), where numIns and numDel are the sums of approx.
 *          gaps from Bi to Bj). It tracks the maximum diff from Bi to any of Bj and logs j_max.
 *          If maxDiff > diffThreshold, then stretch of hits between Bi and Bj_max are marked for removal.
 *          In the next Bi iteration the i is simply incremented by 1, and it can overlap a previously flagged
 *          stretch for removal. Of all the overlapping ranges, only the one with the highest diff will be removed.
 *
 * @param chain Input chained seed hits.
 * @param minGap Minimum gap distaance between two seeds to mark the second one as a breakpoint.
 * @param diffThreshold Maximum allowed approximate gap threshold used to mark a stretch of seed hits for removal.
 * @param maxForwardSeedDist Stop looking at next seed hits if they are this much away in either query or target.
 * @param maxForwardSeedCount Heuristic value to limit the number of succesive seeds.
 * @return ChainedHits. Filtered seed hits, chained.
 */
ChainedHits RefineChainedHits(const ChainedHits& chain, int32_t minGap, int32_t diffThreshold,
                              int32_t maxForwardSeedDist, int32_t maxForwardSeedCount);

/**
 * @brief This filters more extreme outliers.
 *          For every breakpoint B[i], we iterate through the remaining breakpoints B[j] (j > i) and compute:
 *          (1) approximate number of matches `m` (computed as the min(target_dist, query_dist) between the two breakpoints B[j-1] and B[j]), and
 *          (2) `gap1 + gap2` where gap1 is the gap preceding breakpoint B[j-1], and gap2 the gap preceding breakpoint B[j]).
 *          The first `j` for which `m > (gap1 + gap2)` is satisfied (i.e. the diagonal slide is much longer than the gaps in between) is where iteration stops
 *          (lets call this value of `j` as `last_j`).
 *          All seed hits between `i` and `last_j` are filtered out because there are more indel gap jumps than actual matches/mismatches.
 *
 * @param chain Input chained seed hits
 * @param minGap Threshold to select potential breakpoints.
 * @param maxForwardSeedDist Stop looking at next seed hits if they are this much away in either query or target.
 * @return ChainedHits. Filtered seed hits, chained.
 */
ChainedHits RefineChainedHits2(const ChainedHits& chain, int32_t minGap,
                               int32_t maxForwardSeedDist);

/**
 * @brief Has two iterations: one from the left end and one from the right end. They are symmetrical.
 *          In one iteration starting from the front end, if the gap between the current and the previous
 *          seed hit is too large (>totalSpan/2), then all seed hits before the current one will be filtered out.
 *          Here, gap = abs(queryDist - targetDist), and totalSpan = sum(min(queryDist-targetDist)) for all seed hits.
 *          Iteration stops if:
 *              1. totalSpan >= 2 * maxDist (i.e. we covered a large enough portion), or
 *              2. numMatches >= minMatch && numMatches >= maxDist, where numMatches is the sum of bases covered by seed hits.
 *                  (i.e. we covered more than enough match bases). Or,
 *              3. numMatches >= chain.coveredBasesQuery / 2, where chain.coveredBasesQuery is the total number of covered bases
 *                  by seed hits in the query sequence.
 *
 * @param chain Input chain of seed hits.
 * @param maxDist Threshold for maximum distance from the first seed hit
 * @param minMatch Used to break the iteration, but only if number of matches >= maxDist.
 * @return ChainedHits. Filtered seed hits, chained.
 */
ChainedHits RefineBadEnds(const ChainedHits& chain, int32_t maxAllowedDist, int32_t minMatch);

std::vector<Range> GroupByTargetAndStrand(const std::vector<SeedHit>& sortedHits);

std::vector<Range> DiagonalGroup(const std::vector<SeedHit>& sortedHits, int32_t chainBandwidth,
                                 bool overlappingWindows);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_DP_CHAIN_HPP
