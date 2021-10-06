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

#include <emmintrin.h>
#include <pacbio/pancake/Range.h>
#include <pacbio/pancake/SeedHit.h>
#include <cstdint>
#include <memory>
#include <vector>

namespace PacBio {
namespace Pancake {

struct ChainingScratchSpace
{
    std::vector<int32_t> dp;
    std::vector<int32_t> pred;
    std::vector<int32_t> chainId;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"
    std::vector<__m128i> dpSimd;
    std::vector<__m128i> predSimd;
#pragma GCC diagnostic pop
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
    ChainedHits(const int32_t _targetId, const bool _targetRev)
        : targetId(_targetId), targetRev(_targetRev)
    {
    }
    ChainedHits(const int32_t _targetId, const bool _targetRev, const std::vector<SeedHit>& _hits,
                const int32_t _score, const int32_t _coveredBasesQuery,
                const int32_t _coveredBasesTarget)
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

int32_t ChainHitsForward(const SeedHit* hits, const int32_t hitsSize, const int32_t chainMaxSkip,
                         const int32_t chainMaxPredecessors, const int32_t seedJoinDist,
                         const int32_t diagMargin, std::vector<int32_t>& dp,
                         std::vector<int32_t>& pred, std::vector<int32_t>& chainId);

int32_t ChainHitsForwardFastSisd(const SeedHit* hits, const int32_t hitsSize,
                                 const int32_t chainMaxSkip, const int32_t chainMaxPredecessors,
                                 const int32_t seedJoinDist, const int32_t diagMargin,
                                 std::vector<int32_t>& dp, std::vector<int32_t>& pred,
                                 std::vector<int32_t>& chainId);

std::vector<ChainedHits> ChainHitsBacktrack(const SeedHit* hits, const int32_t hitsSize,
                                            const int32_t* dp, const int32_t* pred,
                                            const int32_t* chainId, const int32_t numChains,
                                            const int32_t minNumSeeds, const int32_t minCovBases,
                                            const int32_t minDPScore);

std::vector<ChainedHits> ChainHits(const SeedHit* hits, const int32_t hitsSize,
                                   const int32_t chainMaxSkip, const int32_t chainMaxPredecessors,
                                   const int32_t seedJoinDist, const int32_t diagMargin,
                                   const int32_t minNumSeeds, const int32_t minCovBases,
                                   const int32_t minDPScore, double& timeChaining,
                                   double& timeBacktrack,
                                   std::shared_ptr<ChainingScratchSpace> ss = nullptr);

double ComputeChainDivergence(const std::vector<SeedHit>& hits);

ChainedHits RefineChainedHits(const ChainedHits& chain, const int32_t minGap,
                              const int32_t diffThreshold, const int32_t maxForwardSeedDist,
                              const int32_t maxForwardSeedCount);

ChainedHits RefineChainedHits2(const ChainedHits& chain, const int32_t minGap,
                               const int32_t maxForwardSeedDist);

ChainedHits RefineBadEnds(const ChainedHits& chain, const int32_t bandwidth,
                          const int32_t minMatch);

std::vector<Range> GroupByTargetAndStrand(const std::vector<SeedHit>& sortedHits);

std::vector<Range> DiagonalGroup(const std::vector<SeedHit>& sortedHits,
                                 const int32_t chainBandwidth, const bool overlappingWindows);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_DP_CHAIN_H
