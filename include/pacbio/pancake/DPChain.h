/*
 * DPChain.h
 *
 *  Created on: Dec 16, 2017
 *      Author: Ivan Sovic
 *
 * Originally implemented in the Raptor graph-based mapper.
 */

#ifndef PANCAKE_DP_CHAIN_H
#define PANCAKE_DP_CHAIN_H

#include <pacbio/pancake/Range.h>
#include <pacbio/pancake/SeedHit.h>
#include <cstdint>
#include <vector>

namespace PacBio {
namespace Pancake {

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
    {
    }
};
inline bool operator==(const ChainedHits& lhs, const ChainedHits& rhs)
{
    return lhs.targetId == rhs.targetId && lhs.targetRev == rhs.targetRev && lhs.hits == rhs.hits &&
           lhs.score == rhs.score && lhs.coveredBasesQuery == rhs.coveredBasesQuery &&
           lhs.coveredBasesTarget == rhs.coveredBasesTarget;
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

inline int32_t SisdCompareLte(int32_t a, int32_t b)
{
    /**
     * Returns 0xFFFFFFFF if a <= b, 0x00000000 otherwise.
    */
    return (a - b - 1) >> 31;  // a <= b
}
inline int32_t SisdCompareLt(int32_t a, int32_t b)
{
    /**
     * Returns 0xFFFFFFFF if a < b, 0x00000000 otherwise.
    */
    return (a - b) >> 31;  // a < b
}
inline int32_t SisdCompareGte(int32_t a, int32_t b)
{
    /**
     * Returns 0xFFFFFFFF if a >= b, 0x00000000 otherwise.
    */
    return (b - a - 1) >> 31;  // a >= b
}
inline int32_t SisdCompareGt(int32_t a, int32_t b)
{
    /**
     * Returns 0xFFFFFFFF if a > b, 0x00000000 otherwise.
    */
    return (b - a) >> 31;  // a > b
}

inline uint32_t ilog2_32_clz_special_zero(int32_t v)
{
    /**
     * \brief This is a modified integer log2 function, it returns 0 when v == 0.
     *          Of course, a proper mathematical log2 would return -inf, but this
     *          special case is required for chaining.
    */
    // If v <= 0, then lt == 0x0, otherwise lt == 0xFFFFFFFF.
    const int32_t lt = (0 - v) >> 31;
    return (31 - __builtin_clz(v)) & lt;
}

int32_t ChainHitsForward(const SeedHit* hits, int32_t hitsSize, int32_t chainMaxSkip,
                         int32_t chainMaxPredecessors, int32_t seedJoinDist, int32_t diagMargin,
                         std::vector<int32_t>& dp, std::vector<int32_t>& pred,
                         std::vector<int32_t>& chainId);

int32_t ChainHitsForwardFastSisd(const SeedHit* hits, const int32_t hitsSize,
                                 const int32_t chainMaxSkip, const int32_t chainMaxPredecessors,
                                 const int32_t seedJoinDist, const int32_t diagMargin,
                                 std::vector<int32_t>& dp, std::vector<int32_t>& pred,
                                 std::vector<int32_t>& chainId);

std::vector<ChainedHits> ChainHitsBacktrack(const SeedHit* hits, int32_t hitsSize,
                                            const int32_t* dp, const int32_t* pred,
                                            const int32_t* chainId, const int32_t numChains,
                                            int32_t minNumSeeds, int32_t minCovBases,
                                            int32_t minDPScore);

std::vector<ChainedHits> ChainHits(const SeedHit* hits, int32_t hits_size, int32_t chain_max_skip,
                                   int32_t chain_max_predecessors, int32_t seed_join_dist,
                                   int32_t diag_margin, int32_t min_num_seeds,
                                   int32_t min_cov_bases, int32_t min_dp_score,
                                   double& timeChaining, double& timeBacktrack);

double ComputeChainDivergence(const std::vector<SeedHit>& hits);

ChainedHits RefineChainedHits(const ChainedHits& chain, int32_t minGap, int32_t diffThreshold,
                              int32_t maxForwardSeedDist, int32_t maxForwardSeedCount);

ChainedHits RefineChainedHits2(const ChainedHits& chain, int32_t minGap,
                               int32_t maxForwardSeedDist);

ChainedHits RefineBadEnds(const ChainedHits& chain, int32_t bandwidth, int32_t minMatch);

std::vector<Range> GroupByTargetAndStrand(const std::vector<SeedHit>& sortedHits);

std::vector<Range> DiagonalGroup(const std::vector<SeedHit>& sortedHits, int32_t chainBandwidth,
                                 bool overlappingWindows);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_DP_CHAIN_H