/*
 * DPChain.cpp
 *
 *  Created on: Sep 30, 2021
 *      Author: Ivan Sovic
 *
 * SSE-vectorized dynamic programming based seed hit chaining.
 */

#include <pacbio/pancake/DPChain.h>
#include <pacbio/util/TicToc.h>
#include <array>
#include <iostream>
#include <lib/math.hpp>
#include <sstream>

#include <emmintrin.h>
#include <smmintrin.h>
#include <tmmintrin.h>
#include <xmmintrin.h>

// #define PANCAKE_DPCHAIN_SIMD_DEBUG

#define PANCAKE_DPCHAIN_SIMD_SKIP_NONMONOTONIC_COORDS

namespace PacBio {
namespace Pancake {

int32_t ChainHitsForwardFastSimd(const SeedHit* hits, const int32_t hitsSize,
                                 const int32_t chainMaxSkip, const int32_t chainMaxPredecessors,
                                 const int32_t seedJoinDist, const int32_t diagMargin,
                                 std::vector<__m128i>& dp, std::vector<__m128i>& pred,
                                 std::vector<int32_t>& chainId)
{
    constexpr int32_t VECTOR_SIZE = 128;
    constexpr int32_t REGISTER_SIZE = 32;
    constexpr uint32_t NUM_REGISTERS = VECTOR_SIZE / REGISTER_SIZE;
    constexpr int32_t MAX_INT32 = std::numeric_limits<int32_t>::max();

    dp.clear();
    pred.clear();
    chainId.clear();

    if (hitsSize == 0) {
        return 0;
    }

    // Compute the linear factor for the DP.
    const int32_t finalChainMaxPredecessors =
        (chainMaxPredecessors > 0) ? chainMaxPredecessors : std::numeric_limits<int32_t>::max();
    float avgQuerySpan = 0.0f;
    for (int32_t i = 0; i < hitsSize; ++i) {
        avgQuerySpan += static_cast<float>(hits[i].querySpan);
    }
    avgQuerySpan = (hitsSize > 0) ? avgQuerySpan / static_cast<float>(hitsSize) : 0.0f;
    const float linFactor = 0.01f * avgQuerySpan;
    // const double avgQuerySpan = hits[0].querySpan;

    // Allocate the vectors.
    const int32_t dpPaddedSize = std::ceil(static_cast<float>(hitsSize) / NUM_REGISTERS);
    dp.resize(dpPaddedSize);
    pred.resize(dpPaddedSize);
    chainId.resize(hitsSize, -1);
    std::vector<__m128i> qp(dpPaddedSize);  // Query pos.
    std::vector<__m128i> tp(dpPaddedSize);  // Target pos.
    std::vector<__m128i> qs(dpPaddedSize);  // Query span.

    // Used to allow direct access to data, instead through operator[].
    __m128i* dpPtr = dp.data();
    __m128i* qpPtr = qp.data();
    __m128i* tpPtr = tp.data();
    __m128i* qsPtr = qs.data();

    // Helper pointers for direct data access, without using operator[].
    int32_t* dpInt32 = reinterpret_cast<int32_t*>(dpPtr);
    int32_t* predInt32 = reinterpret_cast<int32_t*>(pred.data());
    int32_t* chainIdInt32 = chainId.data();
    int32_t* qpInt32 = reinterpret_cast<int32_t*>(qpPtr);
    int32_t* tpInt32 = reinterpret_cast<int32_t*>(tpPtr);
    int32_t* qsInt32 = reinterpret_cast<int32_t*>(qsPtr);

    // Initialize the values needed for DP computation.
    for (int32_t i = 0; i < hitsSize; ++i) {
        dpInt32[i] = 0;
        predInt32[i] = -1;
        qpInt32[i] = hits[i].queryPos;
        tpInt32[i] = hits[i].targetPos;
        qsInt32[i] = hits[i].querySpan;
    }
    for (uint32_t i = hitsSize; i < dpPaddedSize * NUM_REGISTERS; ++i) {
        dpInt32[i] = 0;
        predInt32[i] = -1;
        qpInt32[i] = MAX_INT32;
        tpInt32[i] = MAX_INT32;
        qsInt32[i] = MAX_INT32;
    }

    // The distDiag values will have to be used to compute the log2.
    __m128i distDiag = _mm_set1_epi32(0);
    int32_t* distDiagPtr_0 = reinterpret_cast<int32_t*>(&distDiag);
    int32_t* distDiagPtr_1 = distDiagPtr_0 + 1;
    int32_t* distDiagPtr_2 = distDiagPtr_0 + 2;
    int32_t* distDiagPtr_3 = distDiagPtr_0 + 3;

    // Compute the addresses here, instead of the inner loop. Needed to compute numSkippedPredecessors.
    __m128i skipDiff = _mm_set1_epi32(0);
    int32_t* skipDiffPtr_0 = reinterpret_cast<int32_t*>(&skipDiff);
    int32_t* skipDiffPtr_1 = skipDiffPtr_0 + 1;
    int32_t* skipDiffPtr_2 = skipDiffPtr_0 + 2;
    int32_t* skipDiffPtr_3 = skipDiffPtr_0 + 3;

    // Constants for score computation.
    const __m128 linFactorVec = _mm_set1_ps(linFactor);
    const __m128i diagMarginVec =
        _mm_set_epi32(diagMargin - 1, diagMargin - 1, diagMargin - 1, diagMargin - 1);
    const __m128i seedJoinDistVec =
        _mm_set_epi32(seedJoinDist - 1, seedJoinDist - 1, seedJoinDist - 1, seedJoinDist - 1);

    const __m128i M128_MASK_1BIT = _mm_set_epi32(0, 0, 0, 1);
    const __m128i M128_EPI32_ALL_POSITIVE_1 = _mm_set_epi32(1, 1, 1, 1);
    const __m128i M128_EPI32_ALL_NEGATIVE_1 = _mm_set_epi32(-1, -1, -1, -1);
    // const __m128i M128_MASK_FULL = _mm_set1_epi32(0xFFFFFFFF);
    // const __m128i M128_MASK_ZERO = _mm_set1_epi32(0);
    // const __m128i M128_MASK_31 = _mm_set1_epi32(31);

    // clang-format off
    // Compute addresses before the loops. Needed to compute the maximum score.
    __m128i bestDpScore = _mm_set1_epi32(0);
    std::array<int32_t*, 4> bestDpScorePtr = {
        reinterpret_cast<int32_t *>(&bestDpScore),
        reinterpret_cast<int32_t *>(&bestDpScore) + 1,
        reinterpret_cast<int32_t *>(&bestDpScore) + 2,
        reinterpret_cast<int32_t *>(&bestDpScore) + 3,
    };

    __m128i bestDpPred = _mm_set1_epi32(-1);
    std::array<int32_t*, 4> bestDpPredPtr = {
        reinterpret_cast<int32_t *>(&bestDpPred),
        reinterpret_cast<int32_t *>(&bestDpPred) + 1,
        reinterpret_cast<int32_t *>(&bestDpPred) + 2,
        reinterpret_cast<int32_t *>(&bestDpPred) + 3,
    };
// clang-format on

#ifdef PANCAKE_DPCHAIN_SIMD_DEBUG
    auto PrintVectorInt32 = [](std::ostream& os, __m128i vals) {
        const int32_t* valsInt32 = reinterpret_cast<int32_t*>(&vals);
        os << valsInt32[0] << ", " << valsInt32[1] << ", " << valsInt32[2] << ", " << valsInt32[3];
    };
    auto PrintVectorFloat = [](std::ostream& os, __m128 vals) {
        const float* valsFloat = reinterpret_cast<float*>(&vals);
        os << valsFloat[0] << ", " << valsFloat[1] << ", " << valsFloat[2] << ", " << valsFloat[3];
    };
#endif

    int32_t numChains = 0;

    for (int32_t i = 0, minJ = 0; i < hitsSize; ++i) {
        const auto& hi = hits[i];

        // Move the furthest allowed predecessor.
        while ((minJ < i) && (hits[i].targetId != hits[minJ].targetId ||
                              hits[i].targetRev != hits[minJ].targetRev ||
                              (hits[i].targetPos - hits[minJ].targetPos) > seedJoinDist ||
                              (i - minJ) > finalChainMaxPredecessors)) {
            ++minJ;
        }

        // Current positions.
        const __m128i qpi = _mm_set1_epi32(hi.queryPos);
        const __m128i tpi = _mm_set1_epi32(hi.targetPos);
        const __m128i qpiMinusOne = _mm_sub_epi32(qpi, M128_EPI32_ALL_POSITIVE_1);
        const __m128i tpiMinusOne = _mm_sub_epi32(tpi, M128_EPI32_ALL_POSITIVE_1);

        bestDpScore = _mm_set1_epi32(hi.querySpan);
        bestDpPred = _mm_set1_epi32(-1);

#ifdef PANCAKE_DPCHAIN_SIMD_SKIP_NONMONOTONIC_COORDS
        const int32_t currChainMaxSkip = chainMaxSkip;
#elif
        // The position "i" can begin in any of the NUM_REGISTERS registers in the current SIMD vector. Any
        // register which begins at position "i" or later in that vector will not be valid. The maxSkip needs to be
        // increased to alow for that, or the heuristic will terminate early.
        const int32_t currChainMaxSkip = chainMaxSkip + (NUM_REGISTERS - (i % NUM_REGISTERS));
#endif

        // Loop "j" boundaries.
        const int32_t startJ4 = i / NUM_REGISTERS;
        const int32_t minJ4 = minJ / NUM_REGISTERS;

#ifdef PANCAKE_DPCHAIN_SIMD_DEBUG
        std::cerr << "[i = " << i << "] startJ4 = " << startJ4 << ", minJ4 = " << minJ4
                  << ", minJ = " << minJ << ", currChainMaxSkip = " << currChainMaxSkip
                  << "; hit = {" << hi << "}"
                  << "\n";
#endif

        // clang-format off
        for (int32_t j = startJ4, numSkippedPredecessors = 0; (j >= minJ4) && numSkippedPredecessors <= currChainMaxSkip; --j) {

#ifdef PANCAKE_DPCHAIN_SIMD_DEBUG
            std::cerr << "    [i = " << i << ", i&0x03 = " << (i & 0x03)
                      << ", startJ4 = " << startJ4 << ", j = " << j
                      << "] numSkippedPredecessors = " << numSkippedPredecessors << "\n";
#endif

            // Compute X and Y distsances, and the diagonal distance.
            const __m128i dpj = _mm_load_si128(&dpPtr[j]);
            const __m128i qsj = _mm_load_si128(&qsPtr[j]);
            const __m128i qpj = _mm_load_si128(&qpPtr[j]);
            const __m128i tpj = _mm_load_si128(&tpPtr[j]);
            const __m128i distQuery = _mm_sub_epi32(qpi, qpj);
            const __m128i distTarget = _mm_sub_epi32(tpi, tpj);
            distDiag = _mm_abs_epi32(_mm_sub_epi32(distTarget, distQuery));

            // Check any of the boundary criteria.
            const __m128i c0 = _mm_cmplt_epi32(qpiMinusOne, qpj);
            const __m128i c1 = _mm_cmplt_epi32(tpiMinusOne, tpj);
            const __m128i c2 = _mm_cmplt_epi32(diagMarginVec, distDiag);
            const __m128i c3 = _mm_cmplt_epi32(seedJoinDistVec, distQuery);
            const __m128i c = _mm_or_si128(c0, _mm_or_si128(c1, _mm_or_si128(c2, c3)));

            // Compute the linear part of the score.
            const __m128 distDiagFloat = _mm_cvtepi32_ps(distDiag);
            const __m128 linPartFloat = _mm_mul_ps(linFactorVec, distDiagFloat);
            // IMPORTANT: _mm_cvtps_epi32 rounds to the closest int, and not down.
            const __m128i linPart = _mm_cvtps_epi32(_mm_floor_ps(linPartFloat));

            // Compute the log part of the score. Not using SIMD because there are no
            // SSE alternatives to the __builtin_clz.
            const int32_t logPart0 = ilog2_32_clz_special_zero(*distDiagPtr_0);
            const int32_t logPart1 = ilog2_32_clz_special_zero(*distDiagPtr_1);
            const int32_t logPart2 = ilog2_32_clz_special_zero(*distDiagPtr_2);
            const int32_t logPart3 = ilog2_32_clz_special_zero(*distDiagPtr_3);
            __m128i logPart = _mm_set_epi32(logPart3, logPart2, logPart1, logPart0);
            logPart = _mm_srl_epi32(logPart, M128_MASK_1BIT);

            // // Version 2: Use SSE instructions for simple math operations inside ilog2_32_clz_special_zero.
            // const int32_t logPart0 = __builtin_clz(*distDiagPtr_0);
            // const int32_t logPart1 = __builtin_clz(*distDiagPtr_1);
            // const int32_t logPart2 = __builtin_clz(*distDiagPtr_2);
            // const int32_t logPart3 = __builtin_clz(*distDiagPtr_3);
            // __m128i logPart = _mm_set_epi32(logPart3, logPart2, logPart1, logPart0);
            // logPart = _mm_sub_epi32(M128_MASK_31, logPart);
            // logPart = _mm_and_si128(logPart, _mm_cmpgt_epi32(distDiag, M128_MASK_ZERO));
            // logPart = _mm_srl_epi32(logPart, M128_MASK_1BIT);

            // // Version 3: Reuse the distDiag variable.
            // const __m128i tmp = _mm_cmpgt_epi32(distDiag, M128_MASK_ZERO);
            // distDiagPtr_0[0] = __builtin_clz(distDiagPtr_0[0]);
            // distDiagPtr_0[1] = __builtin_clz(distDiagPtr_0[1]);
            // distDiagPtr_0[2] = __builtin_clz(distDiagPtr_0[2]);
            // distDiagPtr_0[3] = __builtin_clz(distDiagPtr_0[3]);
            // __m128i logPart = _mm_sub_epi32(M128_MASK_31, distDiag);
            // logPart = _mm_and_si128(logPart, tmp);
            // logPart = _mm_srl_epi32(logPart, M128_MASK_1BIT);

            // Compute the new score. Invalidate the values which are out of bounds (defined by c).
            const __m128i edgeScore = _mm_add_epi32(linPart, logPart);
            const __m128i matchScore = _mm_min_epi32(qsj, _mm_min_epi32(distQuery, distTarget));
            const __m128i score = _mm_or_si128(_mm_add_epi32(dpj, _mm_sub_epi32(matchScore, edgeScore)), c);

            // Pick the best score.
            const __m128i jVec = _mm_set1_epi32(j);
            const __m128i isBetter = _mm_cmpgt_epi32(score, bestDpScore);
            bestDpScore = _mm_blendv_epi8(bestDpScore, score, isBetter);
            bestDpPred = _mm_blendv_epi8(bestDpPred, jVec, isBetter);

            // Horizontal add to update the numSkippedPredecessors, and limit lower value to 0.
            skipDiff = _mm_blendv_epi8(M128_EPI32_ALL_POSITIVE_1, M128_EPI32_ALL_NEGATIVE_1, isBetter);
#ifdef PANCAKE_DPCHAIN_SIMD_SKIP_NONMONOTONIC_COORDS
            // NOTE: The following line would re-enable the "continue" behaviour (where coordinates out of order
            // would not be counted in numSkippedPredecessors:
            skipDiff = _mm_andnot_si128(c, skipDiff);
#endif
            numSkippedPredecessors += *skipDiffPtr_0;
            numSkippedPredecessors = std::max(numSkippedPredecessors, 0);
            numSkippedPredecessors += *skipDiffPtr_1;
            numSkippedPredecessors = std::max(numSkippedPredecessors, 0);
            numSkippedPredecessors += *skipDiffPtr_2;
            numSkippedPredecessors = std::max(numSkippedPredecessors, 0);
            numSkippedPredecessors += *skipDiffPtr_3;
            numSkippedPredecessors = std::max(numSkippedPredecessors, 0);

#ifdef PANCAKE_DPCHAIN_SIMD_DEBUG
            {
                std::cerr << "        range = [" << (j * NUM_REGISTERS) << ", " << ((j + 1) * NUM_REGISTERS) << "]\n";
                std::cerr << "        score = [";
                PrintVectorInt32(std::cerr, score);
                std::cerr << "], isBetter = [";
                PrintVectorInt32(std::cerr, isBetter);
                std::cerr << "], c = [";
                PrintVectorInt32(std::cerr, c);
                std::cerr << "], bestDpScore = [";
                PrintVectorInt32(std::cerr, bestDpScore);
                std::cerr << "], bestDpPred = [";
                PrintVectorInt32(std::cerr, bestDpPred);
                std::cerr << "]";
                std::cerr << "\n";
                std::cerr << "        logPart = [";
                PrintVectorInt32(std::cerr, logPart);
                std::cerr << "], linPart [";
                PrintVectorInt32(std::cerr, linPart);
                std::cerr << "], matchScore [";
                PrintVectorInt32(std::cerr, matchScore);
                std::cerr << "]\n";
                std::cerr << "        distDiagFloat = [";
                PrintVectorFloat(std::cerr, distDiagFloat);
                std::cerr << "], linPartFloat [";
                PrintVectorFloat(std::cerr, linPartFloat);
                std::cerr << "]\n";
                std::cerr << "        numSkippedPredecessors = " << numSkippedPredecessors
                          << ", currChainMaxSkip = " << currChainMaxSkip << "\n";
                std::cerr << "\n";
            }
#endif
        }
        // clang-format on

        // Find the maximum.
        dpInt32[i] = hi.querySpan;
        predInt32[i] = -1;
        for (size_t j = 0; j < bestDpScorePtr.size(); ++j) {
            if (*bestDpScorePtr[j] > dpInt32[i]) {
                dpInt32[i] = *bestDpScorePtr[j];
                predInt32[i] = (*bestDpPredPtr[j]) * NUM_REGISTERS + j;
            }
        }
        chainIdInt32[i] = (predInt32[i] >= 0) ? chainIdInt32[predInt32[i]] : numChains++;

#ifdef PANCAKE_DPCHAIN_SIMD_DEBUG
        std::cerr << "    - dp[i] = " << dpInt32[i] << ", pred[i] = " << predInt32[i]
                  << ", chainId[i] = " << chainIdInt32[i] << "\n";
        std::cerr << "    - DP: [";
        for (size_t j = 0; j < dp.size(); ++j) {
            PrintVectorInt32(std::cerr, dp[j]);
            if ((j + 1) < dp.size()) {
                std::cerr << ", ";
            }
        }
        std::cerr << "]\n";
        std::cerr << "    - PR: [";
        for (size_t j = 0; j < pred.size(); ++j) {
            PrintVectorInt32(std::cerr, pred[j]);
            if ((j + 1) < pred.size()) {
                std::cerr << ", ";
            }
        }
        std::cerr << "]\n";
        std::cerr << "\n";
#endif
    }

    return numChains;
}

std::vector<ChainedHits> ChainHitsSimd(
    const SeedHit* hits, const int32_t hitsSize, const int32_t chainMaxSkip,
    const int32_t chainMaxPredecessors, const int32_t seedJoinDist, const int32_t diagMargin,
    const int32_t minNumSeeds, const int32_t minCovBases, const int32_t minDPScore,
    double& timeChaining, double& timeBacktrack, std::shared_ptr<ChainingScratchSpace> ss)
{
/**
     * Hits need to be sorted in this order of priority:
     *      target_id, target_rev, target_pos, query_pos
    */

#ifdef PANCAKE_ENABLE_TIMINGS
    TicToc ttPartial;
#endif

    timeChaining = 0.0;
    timeBacktrack = 0.0;

    if (hitsSize == 0) {
        return {};
    }

    if (chainMaxSkip <= 0) {
        return {};
    }

    if (ss == nullptr) {
        ss = std::make_shared<ChainingScratchSpace>();
    }

#ifdef PANCAKE_DPCHAIN_SIMD_DEBUG
    {
        std::vector<int32_t> dp1;
        std::vector<int32_t> pred1;
        std::vector<int32_t> chainId1;

        std::cerr << "========= TEST NON-SIMD ==========\n";
        const int32_t numChains =
            ChainHitsForwardFastSisd(hits, hitsSize, chainMaxSkip, chainMaxPredecessors,
                                     seedJoinDist, diagMargin, dp1, pred1, chainId1);
        std::cerr << "==================================\n";
    }
#endif

    std::vector<__m128i>& dp = ss->dpSimd;
    std::vector<__m128i>& pred = ss->predSimd;
    std::vector<int32_t>& chainId = ss->chainId;

    const int32_t numChains =
        ChainHitsForwardFastSimd(hits, hitsSize, chainMaxSkip, chainMaxPredecessors, seedJoinDist,
                                 diagMargin, dp, pred, chainId);

#ifdef PANCAKE_ENABLE_TIMINGS
    ttPartial.Stop();
    timeChaining = ttPartial.GetMicrosecs();
    ttPartial.Start();
#endif

    ////////////////////
    /// Backtrack.   ///
    ////////////////////
    const int32_t* dpInt32 = reinterpret_cast<int32_t*>(&dp[0]);
    const int32_t* predInt32 = reinterpret_cast<int32_t*>(&pred[0]);
    const int32_t* chainIdInt32 = &chainId[0];

    std::vector<ChainedHits> chains =
        ChainHitsBacktrack(hits, hitsSize, dpInt32, predInt32, chainIdInt32, numChains, minNumSeeds,
                           minCovBases, minDPScore);

#ifdef DEBUG_DP_VERBOSE_
    printf("The DP:\n");
    for (int32_t i = 0; i < dp.size(); i++) {
        printf("[%d] dp[i] = %d, pred[i] = %d\n", i, dp[i], pred[i]);
    }
#endif

#ifdef PANCAKE_ENABLE_TIMINGS
    ttPartial.Stop();
    timeBacktrack = ttPartial.GetMicrosecs();
    ttPartial.Start();
#endif

    return chains;
}

}  // namespace Pancake
}  // namespace PacBio
