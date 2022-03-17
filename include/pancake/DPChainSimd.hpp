/*
 * DPChain.h
 *
 *  Created on: Sep 28, 2021
 *      Author: Ivan Sovic
 *
 */

#ifndef PANCAKE_DP_CHAIN_SIMD_HPP
#define PANCAKE_DP_CHAIN_SIMD_HPP

#include <pancake/DPChain.hpp>
#include <pancake/Range.hpp>
#include <pancake/SeedHit.hpp>

#include <cstdint>
#include <memory>
#include <vector>

namespace PacBio {
namespace Pancake {

int32_t ChainHitsForwardFastSimd(const SeedHit* hits, const int32_t hitsSize,
                                 const int32_t chainMaxSkip, const int32_t chainMaxPredecessors,
                                 const int32_t seedJoinDist, const int32_t diagMargin,
                                 std::vector<int32_t>& retDp, std::vector<int32_t>& retPred,
                                 std::vector<int32_t>& retChainId);

std::vector<ChainedHits> ChainHitsSimd(const SeedHit* hits, const int32_t hitsSize,
                                       const int32_t chainMaxSkip,
                                       const int32_t chainMaxPredecessors,
                                       const int32_t seedJoinDist, const int32_t diagMargin,
                                       const int32_t minNumSeeds, const int32_t minCovBases,
                                       const int32_t minDPScore, double& timeChaining,
                                       double& timeBacktrack,
                                       std::shared_ptr<ChainingScratchSpace> ss = nullptr);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_DP_CHAIN_SIMD_HPP
