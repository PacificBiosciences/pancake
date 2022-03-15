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
#include <span>
#include <vector>

namespace PacBio {
namespace Pancake {

int32_t ChainHitsForwardFastSimd(std::span<const SeedHit> hits, int32_t chainMaxSkip,
                                 int32_t chainMaxPredecessors, int32_t seedJoinDist,
                                 int32_t diagMargin, std::vector<int32_t>& retDp,
                                 std::vector<int32_t>& retPred, std::vector<int32_t>& retChainId);

std::vector<ChainedHits> ChainHitsSimd(std::span<const SeedHit> hits, int32_t chainMaxSkip,
                                       int32_t chainMaxPredecessors, int32_t seedJoinDist,
                                       int32_t diagMargin, int32_t minNumSeeds, int32_t minCovBases,
                                       int32_t minDPScore, double& timeChaining,
                                       double& timeBacktrack,
                                       std::shared_ptr<ChainingScratchSpace> ss = nullptr);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_DP_CHAIN_SIMD_HPP
