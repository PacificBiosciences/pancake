#include <gtest/gtest.h>

#include <pacbio/pancake/DPChain.h>
#include <cstdint>
#include <iostream>
#include <vector>

using namespace PacBio::Pancake;

TEST(DPChain, GroupByTargetAndStrand_ArrayOfTests)
{
    struct TestData
    {
        std::string testName;
        std::vector<PacBio::Pancake::SeedHit> sortedHits;
        bool expectedThrow = false;
        std::vector<PacBio::Pancake::Range> expectedResult;
    };

    // clang-format off
    std::vector<TestData> allTests{
        TestData{
            "Empty input", {}, false, {}
        },

        TestData{
            "Single point",
            {
                SeedHit(0, false, 1000, 1000, 15, 15, 0),
            },
            false,
            {
                {0, 1},
            }
        },

        TestData{
            "Two points, different target",
            {
                SeedHit(0, false, 1000, 1000, 15, 15, 0),
                SeedHit(1, false, 2000, 2000, 15, 15, 0),
            },
            false,
            {
                {0, 1},
                {1, 2},
            }
        },

        TestData{
            "Two points, same target, different strand",
            {
                SeedHit(0, false, 1000, 1000, 15, 15, 0),
                SeedHit(0, true, 2000, 2000, 15, 15, 0),
            },
            false,
            {
                {0, 1},
                {1, 2},
            }
        },

        TestData{
            "Multiple points, same target and strand",
            {
                SeedHit(0, false, 1000, 1000, 15, 15, 0),
                SeedHit(0, false, 2000, 2000, 15, 15, 0),
                SeedHit(0, false, 3000, 3000, 15, 15, 0),
                SeedHit(0, false, 4000, 4000, 15, 15, 0),
                SeedHit(0, false, 5000, 5000, 15, 15, 0),
            },
            false,
            {
                {0, 5},
            }
        },

        TestData{
            "Multiple points, three targets and strands",
            {
                SeedHit(0, false, 1000, 1000, 15, 15, 0),
                SeedHit(0, false, 2000, 2000, 15, 15, 0),
                SeedHit(0, false, 3000, 3000, 15, 15, 0),
                SeedHit(1, false, 4000, 4000, 15, 15, 0),
                SeedHit(1, false, 5000, 5000, 15, 15, 0),
                SeedHit(1, true, 6000, 6000, 15, 15, 0),
                SeedHit(2, true, 7000, 7000, 15, 15, 0),
                SeedHit(2, true, 8000, 8000, 15, 15, 0),
                SeedHit(3, true, 9000, 9000, 15, 15, 0),
            },
            false,
            {
                {0, 3},
                {3, 5},
                {5, 6},
                {6, 8},
                {8, 9},
            }
        },

    };
    // clang-format on

    for (const auto& data : allTests) {
        // Name the test.
        SCOPED_TRACE(data.testName);

        if (data.expectedThrow) {
            EXPECT_THROW({ GroupByTargetAndStrand(data.sortedHits); }, std::runtime_error);

        } else {
            // Run the unit under test.
            std::vector<PacBio::Pancake::Range> result = GroupByTargetAndStrand(data.sortedHits);

            // Evaluate.
            EXPECT_EQ(data.expectedResult, result);
        }
    }
}

TEST(DPChain, DiagonalGroup_ArrayOfTests)
{
    struct TestData
    {
        std::string testName;
        std::vector<PacBio::Pancake::SeedHit> sortedHits;
        int32_t chainBandwidth = 0;
        bool overlappingWindows = 0;
        bool expectedThrow = false;
        std::vector<PacBio::Pancake::Range> expectedResult;
    };

    // clang-format off
    std::vector<TestData> allTests{
        TestData{
            "Empty input", {}, 500, true, false, {}
        },

        TestData{
            "Single point",
            // sortedHits
            {
                SeedHit(0, false, 1000, 1000, 15, 15, 0),
            },
            // chainBandwidth, overlappingWindows
            500, false,
            // expectedThrow
            false,
            // expectedResult
            {
                {0, 1},
            }
        },

        TestData{
            "Multiple points, same target and strand",
            // sortedHits
            {
                SeedHit(0, false, 1000, 1000, 15, 15, 0),
                SeedHit(0, false, 2000, 2000, 15, 15, 0),
                SeedHit(0, false, 3000, 3000, 15, 15, 0),
                SeedHit(0, false, 4000, 4000, 15, 15, 0),
                SeedHit(0, false, 5000, 5000, 15, 15, 0),
            },
            // chainBandwidth, overlappingWindows
            500, false,
            // expectedThrow
            false,
            // expectedResult
            {
                {0, 5},
            }
        },

        TestData{
            "Multiple points, three targets and strands",
            // sortedHits
            {
                SeedHit(0, false, 1000, 1000, 15, 15, 0),
                SeedHit(0, false, 2000, 2000, 15, 15, 0),
                SeedHit(0, false, 3000, 3000, 15, 15, 0),
                SeedHit(1, false, 4000, 4000, 15, 15, 0),
                SeedHit(1, false, 5000, 5000, 15, 15, 0),
                SeedHit(1, true, 6000, 6000, 15, 15, 0),
                SeedHit(2, true, 7000, 7000, 15, 15, 0),
                SeedHit(2, true, 8000, 8000, 15, 15, 0),
                SeedHit(3, true, 9000, 9000, 15, 15, 0),
            },
            // chainBandwidth, overlappingWindows
            500, false,
            // expectedThrow
            false,
            // expectedResult
            {
                {0, 3},
                {3, 5},
                {5, 6},
                {6, 8},
                {8, 9},
            }
        },

        TestData{
            "Multiple points, same target and strand, chain bandwidth broken",
            // sortedHits
            {
                SeedHit(0, false, 1000, 1000, 15, 15, 0),
                SeedHit(0, false, 2000, 2000, 15, 15, 0),
                SeedHit(0, false, 3000, 3000, 15, 15, 0),
                SeedHit(0, false, 4700, 4000, 15, 15, 0),
                SeedHit(0, false, 5700, 5000, 15, 15, 0),
            },
            // chainBandwidth, overlappingWindows
            500, false,
            // expectedThrow
            false,
            // expectedResult
            {
                {0, 3},
                {3, 5},
            }
        },

        TestData{
            "Multiple points, multiple targets and strands, chain bandwidths broken",
            // sortedHits
            {
                SeedHit(0, false, 1000, 1000, 15, 15, 0),
                SeedHit(0, false, 2000, 2000, 15, 15, 0),
                SeedHit(0, false, 3000, 3000, 15, 15, 0),
                SeedHit(0, false, 4700, 4000, 15, 15, 0),
                SeedHit(0, false, 5700, 5000, 15, 15, 0),

                SeedHit(1, false, 1000, 1000, 15, 15, 0),
                SeedHit(1, false, 2000, 2000, 15, 15, 0),
                SeedHit(1, false, 3000, 3000, 15, 15, 0),
                SeedHit(1, false, 4700, 4000, 15, 15, 0),
                SeedHit(1, false, 5700, 5700, 15, 15, 0),
            },
            // chainBandwidth, overlappingWindows
            500, false,
            // expectedThrow
            false,
            // expectedResult
            {
                {0, 3},
                {3, 5},
                {5, 8},
                {8, 9},
                {9, 10},
            }
        },

        TestData{
            "Overlapping windows. Some hits are close to both the first and second range, and are included in both.",
            // sortedHits
            {
                // First diagonal range: from 0bp to 400bp.
                SeedHit(0, false, 1000, 1000, 15, 15, 0),
                SeedHit(0, false, 2000, 2000, 15, 15, 0),
                SeedHit(0, false, 3000, 3000, 15, 15, 0),
                SeedHit(0, false, 4000, 4000, 15, 15, 0),
                // These two seed hits should be included in the first
                // and the second range.
                SeedHit(0, false, 4300, 4000, 15, 15, 0),
                SeedHit(0, false, 4400, 4000, 15, 15, 0),

                // Second diagonal range: 700bp.
                SeedHit(0, false, 4700, 4000, 15, 15, 0),
                SeedHit(0, false, 5700, 5000, 15, 15, 0),
                SeedHit(0, false, 6700, 6000, 15, 15, 0),

                // Third diagonal range: 0bp again.
                SeedHit(0, false, 7700, 7700, 15, 15, 0),
                SeedHit(0, false, 8700, 8700, 15, 15, 0),
                SeedHit(0, false, 9700, 9700, 15, 15, 0),
            },
            // chainBandwidth, overlappingWindows
            500, true,
            // expectedThrow
            false,
            // expectedResult
            {
                {0, 6},
                {4, 9},
                {9, 12},
            }
        },
    };
    // clang-format on

    for (const auto& data : allTests) {
        // Name the test.
        SCOPED_TRACE(data.testName);

        if (data.expectedThrow) {
            EXPECT_THROW(
                { DiagonalGroup(data.sortedHits, data.chainBandwidth, data.overlappingWindows); },
                std::runtime_error);

        } else {
            // Run the unit under test.
            std::vector<PacBio::Pancake::Range> result =
                DiagonalGroup(data.sortedHits, data.chainBandwidth, data.overlappingWindows);

            // Evaluate.
            EXPECT_EQ(data.expectedResult, result);
        }
    }
}

TEST(DPChain, ChainHits_ArrayOfTests)
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
                }, 70, 82, 80), // Note that here we have query and then target coverage in this order, while a SeedHit has target position first, then query position.
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
        {
            "Low complexity seed hit. With integer scores in the DP chaining function, this can result in a wrong seed hit picked for chaining (the one at (50, 47) instead of (50, 50). "
            "The chainMaxSkip does not play a role here, because all low-complexity predecessors to the seed at position (100, 100) will have a better score (up to (50, 50).",
            5, 500, 10000, 500, 3, 0, 40,
            // Seed hits.
            {
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                {0, 0, 0, 0, 15, 15, 0},    // 0

                {0, 0, 50, 42, 15, 15, 0},  // 1
                {0, 0, 50, 43, 15, 15, 0},  // 2
                {0, 0, 50, 44, 15, 15, 0},  // 3
                {0, 0, 50, 45, 15, 15, 0},  // 4
                {0, 0, 50, 46, 15, 15, 0},  // 5
                {0, 0, 50, 47, 15, 15, 0},  // 6
                {0, 0, 50, 48, 15, 15, 0},  // 7
                {0, 0, 50, 49, 15, 15, 0},  // 8
                {0, 0, 50, 50, 15, 15, 0},  // -> This is the correct minimizer to fit in the linear chain, but it's drowned in neighboring hits (8 above and 8 below),
                {0, 0, 50, 51, 15, 15, 0},  // and none of those hits provide a better DP score. Because chainMaxSkip == 5 (<= 8), chaining stops at this seed.
                {0, 0, 50, 52, 15, 15, 0},  // 11
                {0, 0, 50, 53, 15, 15, 0},  // 12
                {0, 0, 50, 54, 15, 15, 0},  // 13
                {0, 0, 50, 55, 15, 15, 0},  // 14
                {0, 0, 50, 56, 15, 15, 0},  // 15
                {0, 0, 50, 57, 15, 15, 0},  // 16
                {0, 0, 50, 58, 15, 15, 0},  // 17

                {0, 0, 100, 100, 15, 15, 0},    // 18
                {0, 0, 150, 150, 15, 15, 0},    // 19
                {0, 0, 200, 200, 15, 15, 0},    // 20
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
            "Low complexity seed hit. Parameter chainMaxSkip is set lower than the copy number of the kmer (chainMaxSkip = 5, while the on-diagonal kmer is 8th at target position 5000)."
            "It was important to bump up the values of the chainMaxGap and chainBandwidth so that the 'continue' statements would not sidetrack the num_skipped_predecessors count."
            "To actually pick up both chains, the minNumSeeds and minDpScore had to be reduced.",
            5, 500, 100000, 50000, 2, 0, 30,
            // Seed hits.
            {
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                {0, 0, 4950, 4950, 15, 15, 0},    // 0

                {0, 0, 5000, 4200, 15, 15, 0},  // 1
                {0, 0, 5000, 4300, 15, 15, 0},  // 2
                {0, 0, 5000, 4400, 15, 15, 0},  // 3
                {0, 0, 5000, 4500, 15, 15, 0},  // 4
                {0, 0, 5000, 4600, 15, 15, 0},  // 5
                {0, 0, 5000, 4700, 15, 15, 0},  // 6
                {0, 0, 5000, 4800, 15, 15, 0},  // 7
                {0, 0, 5000, 4900, 15, 15, 0},  // 8
                {0, 0, 5000, 5000, 15, 15, 0},  // -> This is the correct minimizer to fit in the linear chain, but it's drowned in neighboring hits (8 above and 8 below),
                {0, 0, 5000, 5100, 15, 15, 0},  // and none of those hits provide a better DP score. Because chainMaxSkip == 5 (<= 8), chaining stops at this seed.
                {0, 0, 5000, 5200, 15, 15, 0},  // 11
                {0, 0, 5000, 5300, 15, 15, 0},  // 12
                {0, 0, 5000, 5400, 15, 15, 0},  // 13
                {0, 0, 5000, 5500, 15, 15, 0},  // 14
                {0, 0, 5000, 5600, 15, 15, 0},  // 15
                {0, 0, 5000, 5700, 15, 15, 0},  // 16
                {0, 0, 5000, 5800, 15, 15, 0},  // 17

                {0, 0, 15100, 15100, 15, 15, 0},    // 18
                {0, 0, 15150, 15150, 15, 15, 0},    // 19
                {0, 0, 15200, 15200, 15, 15, 0},    // 20
            },
            // Results.
            {
                // targetId, targetRev, hits, score, coveredQueryBases, coveredTargetBases
                ChainedHits(0, 0, {
                    {0, 0, 4950, 4950, 15, 15, 0},
                    {0, 0, 5000, 5000, 15, 15, 0},
                }, 30, 30, 30),
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
                {0, 0, 4950, 4950, 15, 15, 0},    // 0

                {0, 0, 5000, 4200, 15, 15, 0},  // 1
                {0, 0, 5000, 4300, 15, 15, 0},  // 2
                {0, 0, 5000, 4400, 15, 15, 0},  // 3
                {0, 0, 5000, 4500, 15, 15, 0},  // 4
                {0, 0, 5000, 4600, 15, 15, 0},  // 5
                {0, 0, 5000, 4700, 15, 15, 0},  // 6
                {0, 0, 5000, 4800, 15, 15, 0},  // 7
                {0, 0, 5000, 4900, 15, 15, 0},  // 8
                {0, 0, 5000, 5000, 15, 15, 0},  // -> This is the correct minimizer to fit in the linear chain, but it's drowned in neighboring hits (8 above and 8 below),
                {0, 0, 5000, 5100, 15, 15, 0},  // and none of those hits provide a better DP score. Because chainMaxSkip == 5 (<= 8), chaining stops at this seed.
                {0, 0, 5000, 5200, 15, 15, 0},  // 11
                {0, 0, 5000, 5300, 15, 15, 0},  // 12
                {0, 0, 5000, 5400, 15, 15, 0},  // 13
                {0, 0, 5000, 5500, 15, 15, 0},  // 14
                {0, 0, 5000, 5600, 15, 15, 0},  // 15
                {0, 0, 5000, 5700, 15, 15, 0},  // 16
                {0, 0, 5000, 5800, 15, 15, 0},  // 17

                {0, 0, 15100, 15100, 15, 15, 0},    // 18
                {0, 0, 15150, 15150, 15, 15, 0},    // 19
                {0, 0, 15200, 15200, 15, 15, 0},    // 20
            },
            // Results.
            {
                // targetId, targetRev, hits, score, coveredQueryBases, coveredTargetBases
                ChainedHits(0, 0, {
                    {0, 0, 4950, 4950, 15, 15, 0},
                    {0, 0, 5000, 5000, 15, 15, 0},
                    {0, 0, 15100, 15100, 15, 15, 0},
                    {0, 0, 15150, 15150, 15, 15, 0},
                    {0, 0, 15200, 15200, 15, 15, 0},
                }, 75, 75, 75),
            },
        },
        {
            "Low complexity seed hit. Test chainMaxPredecessors.",
            25, 5, 100000, 50000, 2, 0, 30,
            // Seed hits.
            {
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                {0, 0, 0, 0, 15, 15, 0},    // 0
                {0, 0, 15, 15, 15, 15, 0},    // 1

                {0, 0, 50, 42, 15, 15, 0},  // 2
                {0, 0, 50, 43, 15, 15, 0},  // 3
                {0, 0, 50, 44, 15, 15, 0},  // 4
                {0, 0, 50, 45, 15, 15, 0},  // 5
                {0, 0, 50, 46, 15, 15, 0},  // 6    -> This one will be taken as the end of the first chain, because the inner for loop will stop traversing after chainMaxPredecessors = 5 predecessors (so the next seed will never exit the low complexity region.)
                {0, 0, 50, 47, 15, 15, 0},  // 7
                {0, 0, 50, 48, 15, 15, 0},  // 8
                {0, 0, 50, 49, 15, 15, 0},  // 9
                {0, 0, 50, 50, 15, 15, 0},  // -> This is the correct minimizer to fit in the linear chain, but it's drowned in neighboring hits (8 above and 8 below),
                {0, 0, 50, 51, 15, 15, 0},  // and none of those hits provide a better DP score. Because chainMaxSkip == 5 (<= 8), chaining stops at this seed.
                {0, 0, 50, 52, 15, 15, 0},  // 12
                {0, 0, 50, 53, 15, 15, 0},  // 13
                {0, 0, 50, 54, 15, 15, 0},  // 14   -> This one will be the beginning of the second chain. This is the furthest reachable because of the chainMaxPredecessor heuristic set to 5.
                {0, 0, 50, 55, 15, 15, 0},  // 15
                {0, 0, 50, 56, 15, 15, 0},  // 16
                {0, 0, 50, 57, 15, 15, 0},  // 17
                {0, 0, 50, 58, 15, 15, 0},  // 18

                {0, 0, 100, 100, 15, 15, 0},    // 18
                {0, 0, 150, 150, 15, 15, 0},    // 19
                {0, 0, 200, 200, 15, 15, 0},    // 20
            },
            // Results.
            {
                // targetId, targetRev, hits, score, coveredQueryBases, coveredTargetBases
                ChainedHits(0, 0, {
                    {0, 0, 0, 0, 15, 15, 0},
                    {0, 0, 15, 15, 15, 15, 0},
                    {0, 0, 50, 46, 15, 15, 0},
                }, 44, 45, 45),
                ChainedHits(0, 0, {
                    {0, 0, 50, 54, 15, 15, 0},
                    {0, 0, 100, 100, 15, 15, 0},
                    {0, 0, 150, 150, 15, 15, 0},
                    {0, 0, 200, 200, 15, 15, 0},
                }, 59, 60, 60),
            },
        },

    };

                // {targetId, false, 0, 0, span, span, 0},
    // clang-format on

    for (const auto& data : allTests) {
        // Name the test.
        SCOPED_TRACE(data.testName);

        // Run unit under test.
        const std::vector<ChainedHits> results =
            ChainHits(data.seedHits.data(), data.seedHits.size(), data.chainMaxSkip,
                      data.chainMaxPredecessors, data.chainMaxGap, data.chainBandwidth,
                      data.minNumSeeds, data.minCovBases, data.minDpScore);

        // Evaluate.
        EXPECT_EQ(data.expected, results);
    }
}
