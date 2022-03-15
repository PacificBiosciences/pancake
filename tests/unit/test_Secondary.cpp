#include <gtest/gtest.h>
#include <pancake/Secondary.hpp>
#include <string>
#include <tuple>
#include <vector>

TEST(FlagSecondaryAndSupplementary, SmallTest1)
{
    // clang-format off
    std::vector<std::tuple<std::string, std::vector<std::string>, double, double, double,
                           std::vector<PacBio::Pancake::OverlapPriority>>>
        testData = {
            {   "Test 0 - empty input",
                {},
                0.50, 0.50,
                0.80,
                {}
            },
            {   "Test 1 - single primary alignment.",
                {
                    "0 1 -1000 100.0000 0 0 10000 10000 0 0 10000 4641652",
                },
                0.50, 0.50,
                0.80,
                {
                    {0, false},
                }
            },
            {   "Test 2 - simple secondary, full overlap between mappings in both target and query.",
                {
                    "0 1 -1000 100.0000 0 0 10000 10000 0 0 10000 4641652",
                    "0 2 -900 100.0000 0 0 10000 10000 0 0 10000 4641652",
                },
                0.50, 0.50,
                0.80,
                {
                    {0, false},
                    {1, false},
                }
            },
            {   "Test 3 - simple supplementary, no overlap between mappings in either target or query.",
                {
                    "0 1 -1000 100.0000 0 0 5000 10000 0 0 5000 4641652",
                    "0 2 -900 100.0000 0 5000 10000 10000 0 50000 55000 4641652",
                },
                0.50, 0.50,
                0.80,
                {
                    {0, false},
                    {0, true},
                }
            },
            {   "Test 4 - secondary check, partial overlap > minimum fraction.",
                {
                    "0 1 -1000 100.0000 0 0 10000 10000 0 0 10000 4641652",
                    "0 2 -900 100.0000 0 4000 10000 10000 0 4000 10000 4641652",
                },
                0.50, 0.50,
                0.80,
                {
                    {0, false},
                    {1, false},
                }
            },
            {   "Test 5 - partial overlap of length larger than minimum allowed fraction because the second mapping is fully covered by the overlap. The other overlap should be secondary.",
                {
                    "0 1 -1000 100.0000 0 0 10000 10000 0 0 10000 4641652",
                    "0 2 -900 100.0000 0 6000 10000 10000 0 6000 10000 4641652",
                },
                0.50, 0.50,
                0.80,
                {
                    {0, false},
                    {1, false},
                }
            },
            {   "Test 6 - partial overlap smaller than min allowed fraction for secondary. The other overlap should be supplementary.",
                {
                    "0 1 -1000 100.0000 0 0 7000 10000 0 0 7000 4641652",
                    "0 2 -900 100.0000 0 5000 10000 10000 0 5000 10000 4641652",
                },
                0.50, 0.50,
                0.80,
                {
                    {0, false},
                    {0, true},
                }
            },
            {   "Test 7 - filtering by alignment score. The second overlap has too small alignment score (< 0.80 * primary), so it's flagged out.",
                {
                    "0 1 -1000 100.0000 0 0 10000 10000 0 0 10000 4641652",
                    "0 2 -700 100.0000 0 6000 10000 10000 0 6000 10000 4641652",
                },
                0.50, 0.50,
                0.80,
                {
                    {0, false},
                    {2, false},
                }
            },
            {   "Test 8 - multiple mappings, some overlap, some don't. The output should have one primary, two supplementaries, one filtered and the rest are secondary. Also, the primary alignment is second on the list (not sorted).",
                {
                    "0 1 -1000 100.0000 0 0 1000 10000 0 0 1000 4641652",               // Supplementary.
                    "0 1 -8000 100.0000 0 1000 8000 10000 0 1000 8000 4641652",         // Primary.
                    "0 1 -1000 100.0000 0 8000 10000 10000 0 8000 10000 4641652",       // Supplementary.
                    "0 1 -7000 100.0000 0 3000 10000 10000 0 23000 30000 4641652",      // Secondary.
                    "0 1 -6000 100.0000 0 1000 7000 10000 0 21000 27000 4641652",       // Filtered, score < 0.80 * primary.
                    "0 1 -7000 100.0000 0 1000 8000 10000 0 23000 28000 4641652",       // Secondary.
                },
                0.50, 0.50,
                0.80,
                {
                    {0, true},
                    {0, false},
                    {0, true},
                    {1, false},
                    {2, false},
                    {1, false},
                }
            },
            {   "Test 9 - test secondary flagging when two mappings overlap only in target. E.g. read with missing adapter. Also, the second mapping is to the reverse complement of the target.",
                {
                    "0 1 -1000 100.0000 0 0 5000 10000 0 0 5000 10000",
                    "0 1 -900 100.0000 0 5000 8000 10000 0 2000 5000 10000",
                },
                0.50, 0.50,
                0.80,
                {
                    {0, false},
                    {1, false},
                }
            },
            {   "Test 10 - test secondary flagging when two mappings overlap only in target, but in reverse strands. Note: overlaps are in the FWD strand, as by standard.",
                {
                    "0 1 -1000 100.0000 0 0 5000 10000 0 0 5000 10000",
                    "0 1 -900 100.0000 0 5000 8000 10000 1 2000 5000 10000",
                },
                0.50, 0.50,
                0.80,
                {
                    {0, false},
                    {1, false},
                }
            },

            {   "Test 11a - Overlap in query is allowed, so this should be a supplementary mapping.",
                {
                    "0 1 -1000 100.0000 0 0 5000 10000 0 0 5000 10000",
                    "0 1 -900 100.0000 0 4000 8000 10000 0 4000 8000 10000",
                },
                0.50, 0.50,
                0.80,
                {
                    {0, false},
                    {0, true},
                }
            },
            {   "Test 11b - Overlap in query is not allowed, so this should be a secondary mapping.",
                {
                    "0 1 -1000 100.0000 0 0 5000 10000 0 0 5000 10000",
                    "0 1 -900 100.0000 0 4000 8000 10000 0 4000 8000 10000",
                },
                0.00, 0.50,
                0.80,
                {
                    {0, false},
                    {1, false},
                }
            },
            {   "Test 11c - Overlap in target is not allowed, so this should be a secondary mapping.",
                {
                    "0 1 -1000 100.0000 0 0 5000 10000 0 0 5000 10000",
                    "0 1 -900 100.0000 0 4000 8000 10000 0 4000 8000 10000",
                },
                0.50, 0.00,
                0.80,
                {
                    {0, false},
                    {1, false},
                }
            },
            {   "Test 11d - Overlap in query is not allowed but the two mappings are adjacent in query. They overlap in target, but that's allowed.",
                {
                    "0 1 -1000 100.0000 0 0 5000 10000 0 0 5000 10000",
                    "0 1 -900 100.0000 0 5000 8000 10000 0 4000 8000 10000",
                },
                0.00, 0.50,
                0.80,
                {
                    {0, false},
                    {0, true},
                }
            },
            {   "Test 11d - Strict no overlap allowed, and the second mapping is adjacent in both query and target.",
                {
                    "0 1 -1000 100.0000 0 0 5000 10000 0 0 5000 10000",
                    "0 1 -900 100.0000 0 5000 8000 10000 0 5000 8000 10000",
                },
                0.00, 0.00,
                0.80,
                {
                    {0, false},
                    {0, true},
                }
            },

    };
    // clang-format on

    for (const auto& data : testData) {
        const auto& testName = std::get<0>(data);
        const auto& inOverlapsStr = std::get<1>(data);
        const auto& allowedOverlapFractionQuery = std::get<2>(data);
        const auto& allowedOverlapFractionTarget = std::get<3>(data);
        const auto& minSecondaryScoreFraction = std::get<4>(data);
        const auto& expected = std::get<5>(data);

        SCOPED_TRACE(testName);

        // Convert the input overlap strings to the proper Overlap objects.
        std::vector<PacBio::Pancake::OverlapPtr> inOverlaps;
        for (const auto& ovlStr : inOverlapsStr) {
            auto ovl = PacBio::Pancake::ParseM4OverlapFromString(ovlStr);
            inOverlaps.emplace_back(std::move(ovl));
        }

        std::vector<PacBio::Pancake::OverlapPriority> results =
            PacBio::Pancake::FlagSecondaryAndSupplementary(inOverlaps, allowedOverlapFractionQuery,
                                                           allowedOverlapFractionTarget,
                                                           minSecondaryScoreFraction);

        // Evaluate.
        EXPECT_EQ(expected, results);
    }
}
