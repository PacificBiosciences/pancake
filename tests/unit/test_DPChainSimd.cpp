#include <gtest/gtest.h>

#include <cstdint>
#include <iostream>
#include <pancake/DPChain.hpp>
#include <pancake/DPChainSimd.hpp>
#include <vector>

using namespace PacBio::Pancake;

TEST(DPChainSimd, ChainHits_ArrayOfTests)
{
    struct TestData
    {
        const std::string testName;
        const int32_t chainMaxSkip = 25;
        const int32_t chainMaxPredecessors = 500;
        const int32_t chainMaxGap = 10000;
        const int32_t chainBandwidth = 500;
        const int32_t minNumSeeds = 3;
        const int32_t minCovBases = 0;
        const int32_t minDpScore = 40;
        const std::vector<PacBio::Pancake::SeedHit> seedHits;
        std::vector<ChainedHits> expected;
    };

    // clang-format off
    std::vector<TestData> allTests = {
        {
            "Empty input",
            25, 500, 10000, 500, 3, 0, 40,
            // Seed hits.
            {
            },
            // Results.
            {
            },
        },

        //////////////////////////////
        /// Single seed hit tests. ///
        //////////////////////////////
        {
            "Single seed hit. No chains reported because BOTH minimum DP score and minimum number of seeds thresholds remove it.",
            25, 500, 10000, 500, 3, 0, 40,
            // Seed hits.
            {
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                {0, 0, 50, 50, 15, 15, 0},
            },
            // Results.
            {
                // targetId, targetRev, hits, score, coveredQueryBases, coveredTargetBases
            },
        },
        {
            "Single seed hit. No chains reported because minimum number of seeds thresholds removes it.",
            25, 500, 10000, 500, 3, 0, 0,
            // Seed hits.
            {
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                {0, 0, 50, 50, 15, 15, 0},
            },
            // Results.
            {
                // targetId, targetRev, hits, score, coveredQueryBases, coveredTargetBases
            },
        },
        {
            "Single seed hit. No chains reported because minimum DP score removes it.",
            25, 500, 10000, 500, 0, 0, 40,
            // Seed hits.
            {
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                {0, 0, 50, 50, 15, 15, 0},
            },
            // Results.
            {
                // targetId, targetRev, hits, score, coveredQueryBases, coveredTargetBases
            },
        },
        {
            "Single seed hit. A chain is created because there is no lower threshold for minimum number of seeds and DP chaining score.",
            25, 500, 10000, 500, 0, 0, 0,
            // Seed hits.
            {
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                {0, 0, 50, 50, 15, 15, 0},
            },
            // Results.
            {
                // targetId, targetRev, hits, score, coveredQueryBases, coveredTargetBases
                ChainedHits(0, 0, {
                    {0, 0, 50, 50, 15, 15, 0},
                }, 15, 15, 15),
            },
        },

        /////////////////////////////
        /// Simple linear chains. ///
        /////////////////////////////
        {
            "Simple single linear chain with coordinates perfectly lined up.",
            25, 500, 10000, 500, 3, 0, 40,
            // Seed hits.
            {
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                {0, 0, 0, 0, 15, 15, 0},
                {0, 0, 50, 50, 15, 15, 0},
                {0, 0, 100, 100, 15, 15, 0},
                {0, 0, 150, 150, 15, 15, 0},
                {0, 0, 200, 200, 15, 15, 0},
            },
            // Results.
            {
                // targetId, targetRev, hits, score, coveredQueryBases, coveredTargetBases
                ChainedHits(0, 0, {
                    {0, 0, 0, 0, 15, 15, 0},
                    {0, 0, 50, 50, 15, 15, 0},
                    {0, 0, 100, 100, 15, 15, 0},
                    {0, 0, 150, 150, 15, 15, 0},
                    {0, 0, 200, 200, 15, 15, 0},
                }, 75, 75, 75),
            },
        },
        {
            "Simple linear chain, with indels between seed hits.",
            25, 500, 10000, 500, 3, 0, 40,
            // Seed hits.
            {
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                {0, 0, 0, 0, 15, 15, 0},
                {0, 0, 50, 60, 15, 15, 0},
                {0, 0, 55, 67, 15, 15, 0},
                {0, 0, 100, 125, 15, 15, 0},
                {0, 0, 150, 180, 15, 15, 0},
                {0, 0, 200, 230, 15, 15, 0},
            },
            // Results.
            {
                // targetId, targetRev, hits, score, coveredQueryBases, coveredTargetBases
                ChainedHits(0, 0, {
                    {0, 0, 0, 0, 15, 15, 0},
                    {0, 0, 50, 60, 15, 15, 0},
                    {0, 0, 55, 67, 15, 15, 0},
                    {0, 0, 100, 125, 15, 15, 0},
                    {0, 0, 150, 180, 15, 15, 0},
                    {0, 0, 200, 230, 15, 15, 0},
                }, 75, 82, 80), // Note that here we have query and then target coverage in this order, while a SeedHit has target position first, then query position.
            },
        },
        {
            "Two chains because of a large gap.",
            25, 500, 10000, 500, 3, 0, 40,
            // Seed hits.
            {
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                {0, 0, 0, 0, 15, 15, 0},
                {0, 0, 50, 50, 15, 15, 0},
                {0, 0, 100, 100, 15, 15, 0},
                {0, 0, 150, 150, 15, 15, 0},
                {0, 0, 200, 200, 15, 15, 0},

                {0, 0, 20000, 20000, 15, 15, 0},
                {0, 0, 20050, 20050, 15, 15, 0},
                {0, 0, 20100, 20100, 15, 15, 0},
                {0, 0, 20150, 20150, 15, 15, 0},
                {0, 0, 20200, 20200, 15, 15, 0},
            },
            // Results.
            {
                // targetId, targetRev, hits, score, coveredQueryBases, coveredTargetBases
                ChainedHits(0, 0, {
                    {0, 0, 0, 0, 15, 15, 0},
                    {0, 0, 50, 50, 15, 15, 0},
                    {0, 0, 100, 100, 15, 15, 0},
                    {0, 0, 150, 150, 15, 15, 0},
                    {0, 0, 200, 200, 15, 15, 0},
                }, 75, 75, 75),
                ChainedHits(0, 0, {
                    {0, 0, 20000, 20000, 15, 15, 0},
                    {0, 0, 20050, 20050, 15, 15, 0},
                    {0, 0, 20100, 20100, 15, 15, 0},
                    {0, 0, 20150, 20150, 15, 15, 0},
                    {0, 0, 20200, 20200, 15, 15, 0},
                }, 75, 75, 75),
            },
        },
        {
            "Two chains because of a diagonal skew.",
            25, 500, 10000, 500, 3, 0, 40,
            // Seed hits.
            {
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                {0, 0, 0, 0, 15, 15, 0},
                {0, 0, 50, 50, 15, 15, 0},
                {0, 0, 100, 100, 15, 15, 0},
                {0, 0, 150, 150, 15, 15, 0},
                {0, 0, 200, 200, 15, 15, 0},

                {0, 0, 500, 0, 15, 15, 0},
                {0, 0, 550, 50, 15, 15, 0},
                {0, 0, 600, 100, 15, 15, 0},
                {0, 0, 650, 150, 15, 15, 0},
                {0, 0, 700, 200, 15, 15, 0},
            },
            // Results.
            {
                // targetId, targetRev, hits, score, coveredQueryBases, coveredTargetBases
                ChainedHits(0, 0, {
                    {0, 0, 0, 0, 15, 15, 0},
                    {0, 0, 50, 50, 15, 15, 0},
                    {0, 0, 100, 100, 15, 15, 0},
                    {0, 0, 150, 150, 15, 15, 0},
                    {0, 0, 200, 200, 15, 15, 0},
                }, 75, 75, 75),
                ChainedHits(0, 0, {
                    {0, 0, 500, 0, 15, 15, 0},
                    {0, 0, 550, 50, 15, 15, 0},
                    {0, 0, 600, 100, 15, 15, 0},
                    {0, 0, 650, 150, 15, 15, 0},
                    {0, 0, 700, 200, 15, 15, 0},
                }, 75, 75, 75),
            },
        },
        {
            "Two chains because of a diagonal skew in query coordinates.",
            25, 500, 10000, 500, 3, 0, 40,
            // Seed hits.
            {
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                {0, 0, 0, 0, 15, 15, 0},
                {0, 0, 50, 50, 15, 15, 0},
                {0, 0, 100, 100, 15, 15, 0},
                {0, 0, 150, 150, 15, 15, 0},
                {0, 0, 200, 200, 15, 15, 0},

                {0, 0, 0, 500, 15, 15, 0},
                {0, 0, 50, 550, 15, 15, 0},
                {0, 0, 100, 600, 15, 15, 0},
                {0, 0, 150, 650, 15, 15, 0},
                {0, 0, 200, 700, 15, 15, 0},
            },
            // Results.
            {
                // targetId, targetRev, hits, score, coveredQueryBases, coveredTargetBases
                ChainedHits(0, 0, {
                    {0, 0, 0, 0, 15, 15, 0},
                    {0, 0, 50, 50, 15, 15, 0},
                    {0, 0, 100, 100, 15, 15, 0},
                    {0, 0, 150, 150, 15, 15, 0},
                    {0, 0, 200, 200, 15, 15, 0},
                }, 75, 75, 75),
                ChainedHits(0, 0, {
                    {0, 0, 0, 500, 15, 15, 0},
                    {0, 0, 50, 550, 15, 15, 0},
                    {0, 0, 100, 600, 15, 15, 0},
                    {0, 0, 150, 650, 15, 15, 0},
                    {0, 0, 200, 700, 15, 15, 0},
                }, 75, 75, 75),
            },
        },
        {
            "Two chains because hits are on different strands.",
            25, 500, 10000, 500, 3, 0, 40,
            // Seed hits.
            {
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                {0, 0, 0, 0, 15, 15, 0},
                {0, 0, 50, 50, 15, 15, 0},
                {0, 0, 100, 100, 15, 15, 0},
                {0, 0, 150, 150, 15, 15, 0},
                {0, 0, 200, 200, 15, 15, 0},

                {0, 1, 0, 0, 15, 15, 0},
                {0, 1, 50, 50, 15, 15, 0},
                {0, 1, 100, 100, 15, 15, 0},
                {0, 1, 150, 150, 15, 15, 0},
                {0, 1, 200, 200, 15, 15, 0},
            },
            // Results.
            {
                // targetId, targetRev, hits, score, coveredQueryBases, coveredTargetBases
                ChainedHits(0, 0, {
                    {0, 0, 0, 0, 15, 15, 0},
                    {0, 0, 50, 50, 15, 15, 0},
                    {0, 0, 100, 100, 15, 15, 0},
                    {0, 0, 150, 150, 15, 15, 0},
                    {0, 0, 200, 200, 15, 15, 0},
                }, 75, 75, 75),
                ChainedHits(0, 1, {
                    {0, 1, 0, 0, 15, 15, 0},
                    {0, 1, 50, 50, 15, 15, 0},
                    {0, 1, 100, 100, 15, 15, 0},
                    {0, 1, 150, 150, 15, 15, 0},
                    {0, 1, 200, 200, 15, 15, 0},
                }, 75, 75, 75),
            },
        },
        {
            "Two chains because hits are on different targets.",
            25, 500, 10000, 500, 3, 0, 40,
            // Seed hits.
            {
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                {0, 0, 0, 0, 15, 15, 0},
                {0, 0, 50, 50, 15, 15, 0},
                {0, 0, 100, 100, 15, 15, 0},
                {0, 0, 150, 150, 15, 15, 0},
                {0, 0, 200, 200, 15, 15, 0},

                {123, 0, 0, 0, 15, 15, 0},
                {123, 0, 50, 50, 15, 15, 0},
                {123, 0, 100, 100, 15, 15, 0},
                {123, 0, 150, 150, 15, 15, 0},
                {123, 0, 200, 200, 15, 15, 0},
            },
            // Results.
            {
                // targetId, targetRev, hits, score, coveredQueryBases, coveredTargetBases
                ChainedHits(0, 0, {
                    {0, 0, 0, 0, 15, 15, 0},
                    {0, 0, 50, 50, 15, 15, 0},
                    {0, 0, 100, 100, 15, 15, 0},
                    {0, 0, 150, 150, 15, 15, 0},
                    {0, 0, 200, 200, 15, 15, 0},
                }, 75, 75, 75),
                ChainedHits(123, 0, {
                    {123, 0, 0, 0, 15, 15, 0},
                    {123, 0, 50, 50, 15, 15, 0},
                    {123, 0, 100, 100, 15, 15, 0},
                    {123, 0, 150, 150, 15, 15, 0},
                    {123, 0, 200, 200, 15, 15, 0},
                }, 75, 75, 75),
            },
        },
        {
            "No chains because minimum number of covered bases is not satisfied.",
            25, 500, 10000, 500, 3, 100, 40,
            // Seed hits.
            {
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                {0, 0, 0, 0, 15, 15, 0},
                {0, 0, 50, 50, 15, 15, 0},
                {0, 0, 100, 100, 15, 15, 0},
                {0, 0, 150, 150, 15, 15, 0},
                {0, 0, 200, 200, 15, 15, 0},

                {123, 0, 0, 0, 15, 15, 0},
                {123, 0, 50, 50, 15, 15, 0},
                {123, 0, 100, 100, 15, 15, 0},
                {123, 0, 150, 150, 15, 15, 0},
                {123, 0, 200, 200, 15, 15, 0},
            },
            // Results.
            {
                // targetId, targetRev, hits, score, coveredQueryBases, coveredTargetBases
            },
        },

        ///////////////////////////////
        /// Low-complexity regions. ///
        ///////////////////////////////
        // {    TODO: This test is commented for now to log a potential issue with DP chaining when integer scores are used! Minimap2 likely has the same issue.
        //      In this case, low complexity regions may cause a false positive seed hit to be picked.
        //      When int32_t is used for scores, then this test picks the coordinate (50, 47) instead of (50, 50) for the low-complexity seed hit.
        //
        //     "Low complexity seed hit. With integer scores in the DP chaining function, this can result in a wrong seed hit picked for chaining (the one at (50, 47) instead of (50, 50). "
        //     "The chainMaxSkip does not play a role here, because all low-complexity predecessors to the seed at position (100, 100) will have a better score (up to (50, 50).",
        //     5, 500, 10000, 500, 3, 0, 40,
        //     // Seed hits.
        //     {
        //         // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
        //         {0, 0, 0, 0, 15, 15, 0},    // 0

        //         {0, 0, 50, 42, 15, 15, 0},  // 1
        //         {0, 0, 50, 43, 15, 15, 0},  // 2
        //         {0, 0, 50, 44, 15, 15, 0},  // 3
        //         {0, 0, 50, 45, 15, 15, 0},  // 4
        //         {0, 0, 50, 46, 15, 15, 0},  // 5
        //         {0, 0, 50, 47, 15, 15, 0},  // 6
        //         {0, 0, 50, 48, 15, 15, 0},  // 7
        //         {0, 0, 50, 49, 15, 15, 0},  // 8
        //         {0, 0, 50, 50, 15, 15, 0},  // -> This is the correct minimizer to fit in the linear chain, but it's drowned in neighboring hits (8 above and 8 below),
        //         {0, 0, 50, 51, 15, 15, 0},  // and none of those hits provide a better DP score. Because chainMaxSkip == 5 (<= 8), chaining stops at this seed.
        //         {0, 0, 50, 52, 15, 15, 0},  // 11
        //         {0, 0, 50, 53, 15, 15, 0},  // 12
        //         {0, 0, 50, 54, 15, 15, 0},  // 13
        //         {0, 0, 50, 55, 15, 15, 0},  // 14
        //         {0, 0, 50, 56, 15, 15, 0},  // 15
        //         {0, 0, 50, 57, 15, 15, 0},  // 16
        //         {0, 0, 50, 58, 15, 15, 0},  // 17

        //         {0, 0, 100, 100, 15, 15, 0},    // 18
        //         {0, 0, 150, 150, 15, 15, 0},    // 19
        //         {0, 0, 200, 200, 15, 15, 0},    // 20
        //     },
        //     // Results.
        //     {
        //         // targetId, targetRev, hits, score, coveredQueryBases, coveredTargetBases
        //         ChainedHits(0, 0, {
        //             {0, 0, 0, 0, 15, 15, 0},
        //             {0, 0, 50, 50, 15, 15, 0},
        //             {0, 0, 100, 100, 15, 15, 0},
        //             {0, 0, 150, 150, 15, 15, 0},
        //             {0, 0, 200, 200, 15, 15, 0},
        //         }, 75, 75, 75),
        //     },
        // },

///////////////
        {
            "Low complexity seed hit. Parameter chainMaxSkip is set lower than the copy number of the kmer (chainMaxSkip = 5, while the on-diagonal kmer is 8th at target position 5000)."
            "It was important to bump up the values of the chainMaxGap and chainBandwidth so that the 'continue' statements would not sidetrack the num_skipped_predecessors count."
            "To actually pick up both chains, the minNumSeeds and minDpScore had to be reduced.",
            5, 500, 100000, 50000, 2, 0, 30,
            // Seed hits.
            {
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                {0, 0, 4950, 4950, 15, 15, 0},  // 0    0
                {0, 0, 4965, 4965, 15, 15, 0},  // 1    0

                {0, 0, 5000, 4200, 15, 15, 0},  // 2    0
                {0, 0, 5000, 4300, 15, 15, 0},  // 3    0
                {0, 0, 5000, 4400, 15, 15, 0},  // 4    1
                {0, 0, 5000, 4500, 15, 15, 0},  // 5    1
                {0, 0, 5000, 4600, 15, 15, 0},  // 6    1
                {0, 0, 5000, 4700, 15, 15, 0},  // 7    1
                {0, 0, 5000, 4800, 15, 15, 0},  // 8    2
                {0, 0, 5000, 4900, 15, 15, 0},  // 9    2
                {0, 0, 5000, 5000, 15, 15, 0},  // 10   2   -> This is the correct minimizer to fit in the linear chain, but it's drowned in neighboring hits (8 above and 8 below),
                {0, 0, 5000, 5100, 15, 15, 0},  // 11   2       and none of those hits provide a better DP score. Because chainMaxSkip == 5 (<= 8), chaining stops at this seed.
                {0, 0, 5000, 5200, 15, 15, 0},  // 12   3   -> SIMD specific: This is the last one checked from seed 19, because the distance is a multiple of 4.
                {0, 0, 5000, 5300, 15, 15, 0},  // 13   3
                {0, 0, 5000, 5400, 15, 15, 0},  // 14   3
                {0, 0, 5000, 5500, 15, 15, 0},  // 15   3
                {0, 0, 5000, 5600, 15, 15, 0},  // 16   4
                {0, 0, 5000, 5700, 15, 15, 0},  // 17   4
                {0, 0, 5000, 5800, 15, 15, 0},  // 18   4

                {0, 0, 15100, 15100, 15, 15, 0},    // 19   4
                {0, 0, 15150, 15150, 15, 15, 0},    // 20   5
                {0, 0, 15200, 15200, 15, 15, 0},    // 21   5
            },
            // Results.
            {
                // targetId, targetRev, hits, score, coveredQueryBases, coveredTargetBases
                // NOTE: Faster DP chaining, which does not "continue" when a predecessor with smaller query or target coordinates is reached, would
                // would produce the following chain instead of the one below. Because there is no "continue", even predecessors with smaller
                // coordinates will be counted in numSkippedPredecessors, and the inner loop would end earlier.
                // This also applies if the 'c1' heuristic is not used, i.e. inner loop does not check if targetPos[i] < targetPos[j]. Then, all the
                // seed hits on the identical target position contribute to the numSkippedPredecessors.
                ChainedHits(0, 0, {
                    {0, 0, 4950, 4950, 15, 15, 0},
                    {0, 0, 4965, 4965, 15, 15, 0},
                }, 30, 30, 30),
                // // NOTE: The original DP chaining (which had continue statements in the inner loop) would produce the following chain instead of
                // // the one above. This is because the seed hits with non-valid coordinates (query/target coord is greater than the next seed hit's)
                // // are not processed ('continue') instead of counted in numSkippedPredecessors.
                // ChainedHits(0, 0, {
                //     {0, 0, 4950, 4950, 15, 15, 0},
                //     {0, 0, 4965, 4965, 15, 15, 0},
                //     {0, 0, 5000, 5000, 15, 15, 0},
                // }, 45, 45, 45),
                ChainedHits(0, 0, {
                    {0, 0, 15100, 15100, 15, 15, 0},
                    {0, 0, 15150, 15150, 15, 15, 0},
                    {0, 0, 15200, 15200, 15, 15, 0},
                }, 45, 45, 45),
            },
        },

        {
            "Low complexity seed hit. Parameter chainMaxSkip is set higher than the multiplicity of the repetitive kmer. Chaining should run fine.",
            25, 500, 100000, 50000, 2, 0, 30,
            // Seed hits.
            {
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                {0, 0, 4950, 4950, 15, 15, 0},  // 0
                {0, 0, 4965, 4965, 15, 15, 0},  // 1

                {0, 0, 5000, 4200, 15, 15, 0},  // 2
                {0, 0, 5000, 4300, 15, 15, 0},  // 3
                {0, 0, 5000, 4400, 15, 15, 0},  // 4
                {0, 0, 5000, 4500, 15, 15, 0},  // 5
                {0, 0, 5000, 4600, 15, 15, 0},  // 6
                {0, 0, 5000, 4700, 15, 15, 0},  // 7
                {0, 0, 5000, 4800, 15, 15, 0},  // 8
                {0, 0, 5000, 4900, 15, 15, 0},  // 9
                {0, 0, 5000, 5000, 15, 15, 0},  // -> This is the correct minimizer to fit in the linear chain, but it's drowned in neighboring hits (8 above and 8 below),
                {0, 0, 5000, 5100, 15, 15, 0},  // and none of those hits provide a better DP score. Because chainMaxSkip == 5 (<= 8), chaining stops at this seed.
                {0, 0, 5000, 5200, 15, 15, 0},  // 12
                {0, 0, 5000, 5300, 15, 15, 0},  // 13
                {0, 0, 5000, 5400, 15, 15, 0},  // 14
                {0, 0, 5000, 5500, 15, 15, 0},  // 15
                {0, 0, 5000, 5600, 15, 15, 0},  // 16
                {0, 0, 5000, 5700, 15, 15, 0},  // 17
                {0, 0, 5000, 5800, 15, 15, 0},  // 18

                {0, 0, 15100, 15100, 15, 15, 0},    // 19
                {0, 0, 15150, 15150, 15, 15, 0},    // 20
                {0, 0, 15200, 15200, 15, 15, 0},    // 21
            },
            // Results.
            {
                // targetId, targetRev, hits, score, coveredQueryBases, coveredTargetBases
                ChainedHits(0, 0, {
                    {0, 0, 4950, 4950, 15, 15, 0},
                    {0, 0, 4965, 4965, 15, 15, 0},
                    {0, 0, 5000, 5000, 15, 15, 0},
                    {0, 0, 15100, 15100, 15, 15, 0},
                    {0, 0, 15150, 15150, 15, 15, 0},
                    {0, 0, 15200, 15200, 15, 15, 0},
                }, 90, 90, 90),
            },
        },
        {
            "Low complexity seed hit. Test chainMaxPredecessors.",
            25, 5, 100000, 50000, 2, 0, 30,
            // Seed hits.
            {
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                {0, 0, 0, 0, 15, 15, 0},    // 0    // SIMD vector: 0
                {0, 0, 15, 15, 15, 15, 0},  // 1    // SIMD vector: 0

                {0, 0, 50, 42, 15, 15, 0},  // 2    // SIMD vector: 0
                {0, 0, 50, 43, 15, 15, 0},  // 3    // SIMD vector: 0
                {0, 0, 50, 44, 15, 15, 0},  // 4    // SIMD vector: 1
                {0, 0, 50, 45, 15, 15, 0},  // 5    // SIMD vector: 1
                {0, 0, 50, 46, 15, 15, 0},  // 6    // SIMD vector: 1   -> In SISD, this one would be taken as the end of the first chain, because the inner for loop will stop traversing after chainMaxPredecessors = 5 predecessors (so the next seed will never exit the low complexity region.)
                {0, 0, 50, 47, 15, 15, 0},  // 7    // SIMD vector: 1
                {0, 0, 50, 48, 15, 15, 0},  // 8    // SIMD vector: 2   -> In SIMD, this one is the end of the first chain, because the chainMaxPredecessors is rounded to a value divisible by 4. (Vector (8 - 5) / 4 = 0)
                {0, 0, 50, 49, 15, 15, 0},  // 9    // SIMD vector: 2
                {0, 0, 50, 50, 15, 15, 0},  // 10   // SIMD vector: 2   -> This is the correct minimizer to fit in the linear chain, but it's drowned in neighboring hits (8 above and 8 below),
                {0, 0, 50, 51, 15, 15, 0},  // 11   // SIMD vector: 2       and none of those hits provide a better DP score. Because chainMaxSkip == 5 (<= 8), chaining stops at this seed.
                {0, 0, 50, 52, 15, 15, 0},  // 12   // SIMD vector: 3
                {0, 0, 50, 53, 15, 15, 0},  // 13   // SIMD vector: 3   -> In SIMD, this is the beginning of the second chain because that is the beginning of the furthest reaching vector in the "j" loop. (Vector (19 - 5) / 4 = 3). Technically, seed hit 12 would be the start of the chain, but both 12 and 13 have the identical score of 30. The closest position to 'i' is preferred, so position 13 is chosen.
                {0, 0, 50, 54, 15, 15, 0},  // 14   // SIMD vector: 3   -> In SISD, this one would be the beginning of the second chain. This is the furthest reachable because of the chainMaxPredecessor heuristic set to 5.
                {0, 0, 50, 55, 15, 15, 0},  // 15   // SIMD vector: 3
                {0, 0, 50, 56, 15, 15, 0},  // 16   // SIMD vector: 4
                {0, 0, 50, 57, 15, 15, 0},  // 17   // SIMD vector: 4
                {0, 0, 50, 58, 15, 15, 0},  // 18   // SIMD vector: 4

                {0, 0, 100, 100, 15, 15, 0},  // 19 // SIMD vector: 4
                {0, 0, 150, 150, 15, 15, 0},  // 20 // SIMD vector: 5
                {0, 0, 200, 200, 15, 15, 0},  // 21 // SIMD vector: 5
            },
            // Results.
            {
                /////
                /// These results are the expected behaviour when the inner DP loop checks for targetPos[i] < targetPos[j] and skips those seeds with the 'continue' heuristic.
                /////
                // // targetId, targetRev, hits, score, coveredQueryBases, coveredTargetBases
                // ChainedHits(0, 0, {
                //     {0, 0, 0, 0, 15, 15, 0},
                //     {0, 0, 15, 15, 15, 15, 0},
                //     {0, 0, 50, 48, 15, 15, 0},
                // }, 45, 45, 45),
                // ChainedHits(0, 0, {
                //     {0, 0, 50, 53, 15, 15, 0},
                //     {0, 0, 100, 100, 15, 15, 0},
                //     {0, 0, 150, 150, 15, 15, 0},
                //     {0, 0, 200, 200, 15, 15, 0},
                // }, 60, 60, 60),

                /////
                /// These results are the expected behaviour when the inner DP loop DOES NOT check for targetPos[i] < targetPos[j]
                /// and DOES NOT skip those seeds with the 'continue' heuristic. This may bridge through some low complexity regions, but only if
                /// they are not too deep, because otherwise the chainMaxSkip heuristic will kick in.
                /////
                // targetId, targetRev, hits, score, coveredQueryBases, coveredTargetBases
                ChainedHits(0, 0, {
                    {0, 0, 0, 0, 15, 15, 0},
                    {0, 0, 15, 15, 15, 15, 0},

                    {0, 0, 50, 47, 15, 15, 0},
                    {0, 0, 50, 48, 15, 15, 0},
                    {0, 0, 50, 49, 15, 15, 0},
                    {0, 0, 50, 50, 15, 15, 0},
                    {0, 0, 50, 51, 15, 15, 0},
                    {0, 0, 50, 52, 15, 15, 0},
                    {0, 0, 50, 53, 15, 15, 0},

                    {0, 0, 100, 100, 15, 15, 0},
                    {0, 0, 150, 150, 15, 15, 0},
                    {0, 0, 200, 200, 15, 15, 0},
                }, 90, 96, 90),

            },
        },

        //////////////////////////////////////////
        /// Edge cases                         ///
        //////////////////////////////////////////
        {
            "Edge case in the FIRST vector. Two chains on different strands, but coordinates only appear to be monotonically rising. This tests if the chaining"
            "method will take care of the edge case in the first vector (where minJ is), where seed hits may be mixed.",
            25, 500, 10000, 500, 3, 0, 40,
            // Seed hits.
            {
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                {0, 0, 0, 0, 15, 15, 0},
                {0, 0, 50, 50, 15, 15, 0},
                {0, 0, 100, 100, 15, 15, 0},
                {0, 0, 150, 150, 15, 15, 0},
                {0, 0, 200, 200, 15, 15, 0},

                {0, 1, 220, 220, 15, 15, 0},
                {0, 1, 250, 250, 15, 15, 0},
                {0, 1, 300, 300, 15, 15, 0},
                {0, 1, 350, 350, 15, 15, 0},
                {0, 1, 400, 400, 15, 15, 0},
            },
            // Results.
            {
                // targetId, targetRev, hits, score, coveredQueryBases, coveredTargetBases
                ChainedHits(0, 0, {
                    {0, 0, 0, 0, 15, 15, 0},
                    {0, 0, 50, 50, 15, 15, 0},
                    {0, 0, 100, 100, 15, 15, 0},
                    {0, 0, 150, 150, 15, 15, 0},
                    {0, 0, 200, 200, 15, 15, 0},
                }, 75, 75, 75),
                ChainedHits(0, 1, {
                    {0, 1, 220, 220, 15, 15, 0},
                    {0, 1, 250, 250, 15, 15, 0},
                    {0, 1, 300, 300, 15, 15, 0},
                    {0, 1, 350, 350, 15, 15, 0},
                    {0, 1, 400, 400, 15, 15, 0},
                }, 75, 75, 75),
            },
        },
        {
            "Edge case in the LAST vector of the second chain. Three chains on different target IDs/strands. This tests if the chaining"
            "method will take care of the edge case in the last vector (concretely, of the second chain here), where the starting pivot"
            "is located ('i' coordinate), where seed hits may be mixed.",
            25, 500, 10000, 500, 3, 0, 40,
            // Seed hits.
            {
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                {0, 0, 0, 0, 15, 15, 0},        // 0    0
                {0, 0, 50, 50, 15, 15, 0},      // 1    0
                {0, 0, 100, 100, 15, 15, 0},    // 2    0
                {0, 0, 150, 150, 15, 15, 0},    // 3    0
                {0, 0, 200, 200, 15, 15, 0},    // 4    1

                {0, 1, 220, 220, 15, 15, 0},    // 5    1
                {0, 1, 250, 250, 15, 15, 0},    // 6    1
                {0, 1, 300, 300, 15, 15, 0},    // 7    1
                {0, 1, 350, 350, 15, 15, 0},    // 8    2
                {0, 1, 400, 400, 15, 15, 0},    // 9    2

                {1, 1, 370, 370, 15, 15, 0},    // 10   2
                {1, 1, 380, 380, 15, 15, 0},    // 11   2
                {1, 1, 390, 390, 15, 15, 0},    // 12   3
                {1, 1, 395, 395, 15, 15, 0},    // 13   3
            },
            // Results.
            {
                // targetId, targetRev, hits, score, coveredQueryBases, coveredTargetBases
                ChainedHits(0, 0, {
                    {0, 0, 0, 0, 15, 15, 0},
                    {0, 0, 50, 50, 15, 15, 0},
                    {0, 0, 100, 100, 15, 15, 0},
                    {0, 0, 150, 150, 15, 15, 0},
                    {0, 0, 200, 200, 15, 15, 0},
                }, 75, 75, 75),
                ChainedHits(0, 1, {
                    {0, 1, 220, 220, 15, 15, 0},
                    {0, 1, 250, 250, 15, 15, 0},
                    {0, 1, 300, 300, 15, 15, 0},
                    {0, 1, 350, 350, 15, 15, 0},
                    {0, 1, 400, 400, 15, 15, 0},
                }, 75, 75, 75),
                ChainedHits(1, 1, {
                    {1, 1, 370, 370, 15, 15, 0},
                    {1, 1, 380, 380, 15, 15, 0},
                    {1, 1, 390, 390, 15, 15, 0},
                    {1, 1, 395, 395, 15, 15, 0},
                }, 40, 40, 40),
            },
        },
        {
            "Edge case which can happen if the chaining DP loop compares 'newScore >= bestScore' instead of 'newScore > bestScore'."
            "If there are multiple predecessors with the identical score, then potentially the most distant one from the current position"
            "will be selected because of the sorting order. In reality, the closest one should be chosen so that the"
            "chain produces the highest possible coverage. This may be especially important for HiFi reads."
            "This can happen when there are multiple seed hits which are <= seedSpan apart."
            "So, if the window size is <= seedSpan, then this can cause potential loss of seeds in chaining."
            "The fix is trivial - use '>' instead of '>=' to track the maximum score.",
            25, 500, 10000, 500, 3, 0, 40,
            // Seed hits.
            {
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                {0, 1, 220, 220, 15, 15, 0},    // 0
                {0, 1, 250, 250, 15, 15, 0},    // 1
                {0, 1, 300, 300, 15, 15, 0},    // 2
                {0, 1, 350, 350, 15, 15, 0},    // 3
                {0, 1, 400, 400, 15, 15, 0},    // 4

                {1, 1, 370, 370, 15, 15, 0},    // 5
                {1, 1, 380, 380, 15, 15, 0},    // 6
                {1, 1, 390, 390, 15, 15, 0},    // 7
                {1, 1, 395, 395, 15, 15, 0},    // 8    -> Identical score for transitions between 8->7 and 8->6. Backtracking can pick 8->6 because of the sort order, but it would be better to have 8->7->6.
            },
            // Results.
            {
                // targetId, targetRev, hits, score, coveredQueryBases, coveredTargetBases
                ChainedHits(0, 1, {
                    {0, 1, 220, 220, 15, 15, 0},    // 0
                    {0, 1, 250, 250, 15, 15, 0},    // 1
                    {0, 1, 300, 300, 15, 15, 0},    // 2
                    {0, 1, 350, 350, 15, 15, 0},    // 3
                    {0, 1, 400, 400, 15, 15, 0},    // 4
                }, 75, 75, 75),
                ChainedHits(1, 1, {
                    {1, 1, 370, 370, 15, 15, 0},    // 5
                    {1, 1, 380, 380, 15, 15, 0},    // 6
                    {1, 1, 390, 390, 15, 15, 0},    // 7
                    {1, 1, 395, 395, 15, 15, 0},    // 8
                }, 40, 40, 40),
            },
        },
        {
            "Edge case. When several seed hits are processed at once (i.e. SIMD), the last vector (where pivot 'i' for the current node is located) will contain "
            "nodes which are not predecessors of 'i' (at most 3 extra successor nodes for SSE4)."
            "There are two cases:"
            "1) In case when queryPos[i+1] > queryPos[i] or targetPos[i+1] > targetPos[i], the explicit boundary condition checking will invalidate this case. This is FINE."
            "2) In case when queryPos[i+1] < queryPos[i] or targetPos[i+1] < targetPos[i], then THIS CAN CAUSE THE ISSUE COVERED BY THIS TEST. In this case, the"
            "    dp[i+1] can turn out to be > dp[i], but only if querySpan[i+1] > querySpan[i] (e.g. HPC seeds). The chaining code could pick a successor as the"
            "    'predecessor' which can crash the program. Also, this can most likely occur only when the analyzed vector is the first one, because otherwise"
            "    the actual predecessors dp[j] for j < i will likely have a better score."
            "    This needs to be handled as a special case.",
            25, 500, 10000, 500, 1, 0, 0,
            // Seed hits.
            {
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                {0, 1, 400, 400, 15, 15, 0},    // 0

                {1, 1, 270, 270, 17, 17, 0},    // 1
                {1, 1, 280, 280, 17, 17, 0},    // 2
                {1, 1, 290, 290, 17, 17, 0},    // 3
                {1, 1, 295, 295, 17, 17, 0},    // 4
            },
            // Results.
            {
                // targetId, targetRev, hits, score, coveredQueryBases, coveredTargetBases
                ChainedHits(0, 1, {
                    {0, 1, 400, 400, 15, 15, 0},
                }, 15, 15, 15),
                ChainedHits(1, 1, {
                    {1, 1, 270, 270, 17, 17, 0},
                    {1, 1, 280, 280, 17, 17, 0},
                    {1, 1, 290, 290, 17, 17, 0},
                    {1, 1, 295, 295, 17, 17, 0},
                }, 42, 42, 42),
            },
        },
    };

                // {targetId, false, 0, 0, span, span, 0},
    // clang-format on

    for (const auto& data : allTests) {
        // Name the test.
        SCOPED_TRACE(data.testName);
        std::cerr << "testName = " << data.testName << "\n";

        // Run unit under test.
        double timeChaining = 0.0;
        double timeBacktrack = 0.0;
        std::vector<ChainedHits> results =
            ChainHitsSimd(data.seedHits, data.chainMaxSkip, data.chainMaxPredecessors,
                          data.chainMaxGap, data.chainBandwidth, data.minNumSeeds, data.minCovBases,
                          data.minDpScore, timeChaining, timeBacktrack);
        // Sort to make sure the chains are really the same.
        std::sort(results.begin(), results.end(), [](const auto& a, const auto& b) {
            return std::tuple(a.targetId, a.targetRev, a.score,
                              (a.hits.empty() ? 0 : a.hits.front().targetPos),
                              (a.hits.empty() ? 0 : a.hits.front().queryPos)) >
                   std::tuple(b.targetId, b.targetRev, b.score,
                              (b.hits.empty() ? 0 : b.hits.front().targetPos),
                              (b.hits.empty() ? 0 : b.hits.front().queryPos));
        });

        // Sort the expected data to make sure the chains are really the same.
        auto expectedCopy = data.expected;
        std::sort(expectedCopy.begin(), expectedCopy.end(), [](const auto& a, const auto& b) {
            return std::tuple(a.targetId, a.targetRev, a.score,
                              (a.hits.empty() ? 0 : a.hits.front().targetPos),
                              (a.hits.empty() ? 0 : a.hits.front().queryPos)) >
                   std::tuple(b.targetId, b.targetRev, b.score,
                              (b.hits.empty() ? 0 : b.hits.front().targetPos),
                              (b.hits.empty() ? 0 : b.hits.front().queryPos));
        });

        // Evaluate.
        EXPECT_EQ(expectedCopy, results);
    }
}
