#include <gtest/gtest.h>

#include <pacbio/pancake/AlignmentSeeded.h>
#include <pacbio/pancake/Overlap.h>
#include <pacbio/pancake/OverlapWriterBase.h>

#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

using namespace PacBio::Pancake;

TEST(AlignmentSeeded, ExtractAlignmentRegions_ArrayOfTests)
{
    struct TestData
    {
        std::string testName;
        std::vector<PacBio::Pancake::SeedHit> sortedHits;
        int32_t queryLen = 0;
        int32_t targetLen = 0;
        bool isRev = false;
        int32_t minAlignmentSpan = 200;
        int32_t maxFlankExtensionDist = 5000;
        double flankExtensionFactor = 1.3;
        bool expectedThrow = false;
        std::vector<PacBio::Pancake::AlignmentRegion> expectedRegions;
    };

    // clang-format off
    std::vector<TestData> allTests{
        TestData{
            "Empty input", {}, 0, 0, false, 200, 5000, 1.3, false, {},
        },

        TestData{
            "Test parameters - invalid query length, should throw", {}, -1, 0, false, 200, 5000, 1.3, true, {},
        },
        TestData{
            "Test parameters - invalid target length, should throw", {}, 0, -1, false, 200, 5000, 1.3, true, {},
        },
        TestData{
            "Test parameters - invalid flankExtensionFactor, value is below 1.0", {}, 0, 0, false, 200, 5000, 0.8, true, {},
        },
        TestData{
            "Test parameters - invalid flankExtensionFactor, value is below 0.0", {}, 0, 0, false, 200, 5000, -1.0, true, {},
        },

        TestData{
            "Single seed - global aln is performed from seed to seed, so the final global coords start and end at the same point.",
            {   // sortedHits
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                SeedHit(0, false, 5, 5, 5, 5, 0),
            },
            // queryLen, targetLen, isRev, minAlignmentSpan, maxFlankExtensionDist, flankExtensionFactor
            20, 20, false, 200, 5000, 1.3,
            // expectedThrow
            false,
            // expectedRegions
            {
                {0, 5, 0, 5, false, RegionType::FRONT, 0},
                {5, 15, 5, 15, false, RegionType::BACK, 1},
            },
        },

        TestData{
            "Two seeds",
            {   // sortedHits
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                SeedHit(0, false, 5, 5, 5, 5, 0),
                SeedHit(0, false, 10, 10, 5, 5, 0),
            },
            // queryLen, targetLen, isRev, minAlignmentSpan, maxFlankExtensionDist, flankExtensionFactor
            20, 20, false, 200, 5000, 1.3,
            // expectedThrow
            false,
            // expectedRegions
            {
                // qStart, qSpan, tStart, tSpan, queryRev, regionType, regionId
                {0, 5, 0, 5, false, RegionType::FRONT, 0},
                {5, 5, 5, 5, false, RegionType::GLOBAL, 1},
                {10, 10, 10, 10, false, RegionType::BACK, 2},
            },
        },

        TestData{
            "Three seeds",
            {   // sortedHits
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                SeedHit(0, false, 5, 5, 15, 15, 0),
                SeedHit(0, false, 400, 400, 15, 15, 0),
                SeedHit(0, false, 900, 900, 15, 15, 0),
            },
            // queryLen, targetLen, isRev, minAlignmentSpan, maxFlankExtensionDist, flankExtensionFactor
            1000, 1000, false, 200, 5000, 1.3,
            // expectedThrow
            false,
            // expectedRegions
            {
                // qStart, qSpan, tStart, tSpan, queryRev, regionType, regionId
                {0, 5, 0, 5, false, RegionType::FRONT, 0},
                {5, 395, 5, 395, false, RegionType::GLOBAL, 1},
                {400, 500, 400, 500, false, RegionType::GLOBAL, 2},
                {900, 100, 900, 100, false, RegionType::BACK, 3},
            },
        },

        TestData{
            "Three seeds but middle one is skipped because distance from previous is too short",
            {   // sortedHits
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                SeedHit(0, false, 5, 5, 15, 15, 0),
                SeedHit(0, false, 100, 100, 15, 15, 0),
                SeedHit(0, false, 900, 900, 15, 15, 0),
            },
            // queryLen, targetLen, isRev, minAlignmentSpan, maxFlankExtensionDist, flankExtensionFactor
            1000, 1000, false, 200, 5000, 1.3,
            // expectedThrow
            false,
            // expectedRegions
            {
                // qStart, qSpan, tStart, tSpan, queryRev, regionType, regionId
                {0, 5, 0, 5, false, RegionType::FRONT, 0},
                {5, 895, 5, 895, false, RegionType::GLOBAL, 1},
                {900, 100, 900, 100, false, RegionType::BACK, 2},
            },
        },

        TestData{
            "Three seeds, same as previous, where the middle one is close to the previous hit. Here, we allow any distance for alignment.",
            {   // sortedHits
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                SeedHit(0, false, 5, 5, 15, 15, 0),
                SeedHit(0, false, 100, 100, 15, 15, 0),
                SeedHit(0, false, 900, 900, 15, 15, 0),
            },
            // queryLen, targetLen, isRev, minAlignmentSpan, maxFlankExtensionDist, flankExtensionFactor
            1000, 1000, false, 0, 5000, 1.3,
            // expectedThrow
            false,
            // expectedRegions
            {
                // qStart, qSpan, tStart, tSpan, queryRev, regionType, regionId
                {0, 5, 0, 5, false, RegionType::FRONT, 0},
                {5, 95, 5, 95, false, RegionType::GLOBAL, 1},
                {100, 800, 100, 800, false, RegionType::GLOBAL, 2},
                {900, 100, 900, 100, false, RegionType::BACK, 3},
            },
        },

        TestData{
            "Left flank is too far to be aligned fully, so it's limited to 5000bp",
            {   // sortedHits
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                SeedHit(0, false, 12000, 12000, 15, 15, 0),
                SeedHit(0, false, 14000, 14000, 15, 15, 0),
                SeedHit(0, false, 15000, 15000, 15, 15, 0),
            },
            // queryLen, targetLen, isRev, minAlignmentSpan, maxFlankExtensionDist, flankExtensionFactor
            17000, 17000, false, 200, 5000, 1.3,
            // expectedThrow
            false,
            // expectedRegions
            {
                // qStart, qSpan, tStart, tSpan, queryRev, regionType, regionId
                {7000, 5000, 7000, 5000, false, RegionType::FRONT, 0},
                {12000, 2000, 12000, 2000, false, RegionType::GLOBAL, 1},
                {14000, 1000, 14000, 1000, false, RegionType::GLOBAL, 2},
                {15000, 2000, 15000, 2000, false, RegionType::BACK, 3},
            },
        },

        TestData{
            "Right and left flanks are far from the edge, but we allow that with the maxFlankExtensionDist parameter",
            {   // sortedHits
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                SeedHit(0, false, 12000, 12000, 15, 15, 0),
                SeedHit(0, false, 14000, 14000, 15, 15, 0),
                SeedHit(0, false, 15000, 15000, 15, 15, 0),
            },
            // queryLen, targetLen, isRev, minAlignmentSpan, maxFlankExtensionDist, flankExtensionFactor
            25000, 25000, false, 200, -1, 1.3,
            // expectedThrow
            false,
            // expectedRegions
            {
                // qStart, qSpan, tStart, tSpan, queryRev, regionType, regionId
                {0, 12000, 0, 12000, false, RegionType::FRONT, 0},
                {12000, 2000, 12000, 2000, false, RegionType::GLOBAL, 1},
                {14000, 1000, 14000, 1000, false, RegionType::GLOBAL, 2},
                {15000, 10000, 15000, 10000, false, RegionType::BACK, 3},
            },
        },

        TestData{
            "Hits are not sorted, should throw",
            {   // sortedHits
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                SeedHit(0, false, 12000, 12000, 15, 15, 0),
                SeedHit(0, false, 14000, 14000, 15, 15, 0),
                SeedHit(0, false, 13000, 13000, 15, 15, 0),
                SeedHit(0, false, 15000, 15000, 15, 15, 0),
            },
            // queryLen, targetLen, isRev, minAlignmentSpan, maxFlankExtensionDist, flankExtensionFactor
            25000, 25000, false, 200, 5000, 1.3,
            // expectedThrow
            true,
            // expectedRegions
            {
            },
        },

        TestData{
            "Query is internal to target, and flank extension in target should cover 1.5x the flank length in query to allow for potential indels to be aligned",
            {   // sortedHits
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                SeedHit(0, false, 12000, 2000, 15, 15, 0),
                SeedHit(0, false, 14000, 4000, 15, 15, 0),
                SeedHit(0, false, 15000, 5000, 15, 15, 0),
            },
            // queryLen, targetLen, isRev, minAlignmentSpan, maxFlankExtensionDist, flankExtensionFactor
            7000, 25000, false, 200, 5000, 1.5,
            // expectedThrow
            false,
            // expectedRegions
            {
                // qStart, qSpan, tStart, tSpan, queryRev, regionType, regionId
                {0, 2000, 9000, 3000, false, RegionType::FRONT, 0},
                {2000, 2000, 12000, 2000, false, RegionType::GLOBAL, 1},
                {4000, 1000, 14000, 1000, false, RegionType::GLOBAL, 2},
                {5000, 2000, 15000, 3000, false, RegionType::BACK, 3},
            },
        },

        TestData{
            "Query is internal to target, but no extra flank extension in target is allowed",
            {   // sortedHits
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                SeedHit(0, false, 12000, 2000, 15, 15, 0),
                SeedHit(0, false, 14000, 4000, 15, 15, 0),
                SeedHit(0, false, 15000, 5000, 15, 15, 0),
            },
            // queryLen, targetLen, isRev, minAlignmentSpan, maxFlankExtensionDist, flankExtensionFactor
            7000, 25000, false, 200, 5000, 1.0,
            // expectedThrow
            false,
            // expectedRegions
            {
                // qStart, qSpan, tStart, tSpan, queryRev, regionType, regionId
                {0, 2000, 10000, 2000, false, RegionType::FRONT, 0},
                {2000, 2000, 12000, 2000, false, RegionType::GLOBAL, 1},
                {4000, 1000, 14000, 1000, false, RegionType::GLOBAL, 2},
                {5000, 2000, 15000, 2000, false, RegionType::BACK, 3},
            },
        },

        TestData{
            "Three seeds, reverse complement strand. Input coordinates are expected to be: (1) target coordinates always forward, "
            "(2) query coordinates always in the strand of alignment.",
            {   // sortedHits
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                SeedHit(0, true, 5, 5, 15, 15, 0),
                SeedHit(0, true, 400, 400, 15, 15, 0),
                SeedHit(0, true, 900, 900, 15, 15, 0),
            },
            // queryLen, targetLen, isRev, minAlignmentSpan, maxFlankExtensionDist, flankExtensionFactor
            1000, 1000, true, 200, 5000, 1.3,
            // expectedThrow
            false,
            // expectedRegions
            {
                // qStart, qSpan, tStart, tSpan, queryRev, regionType, regionId
                {0, 5, 0, 5, true, RegionType::FRONT, 0},
                {5, 395, 5, 395, true, RegionType::GLOBAL, 1},
                {400, 500, 400, 500, true, RegionType::GLOBAL, 2},
                {900, 100, 900, 100, true, RegionType::BACK, 3},
            },
        },
        TestData{
            "A mix of forward and reverse seed hits. This should not be possible here.",
            {   // sortedHits
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                SeedHit(0, false, 5, 5, 15, 15, 0),
                SeedHit(0, true, 400, 400, 15, 15, 0),
                SeedHit(0, true, 900, 900, 15, 15, 0),
            },
            // queryLen, targetLen, isRev, minAlignmentSpan, maxFlankExtensionDist, flankExtensionFactor
            1000, 1000, true, 200, 5000, 1.3,
            // expectedThrow
            true,
            // expectedRegions
            {
            },
        },
        TestData{
            "A mix of forward and reverse seed hits, but an internal hit is of the opposite strand. This should not be possible here.",
            {   // sortedHits
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                SeedHit(0, true, 5, 5, 15, 15, 0),
                SeedHit(0, false, 400, 400, 15, 15, 0),
                SeedHit(0, true, 900, 900, 15, 15, 0),
            },
            // queryLen, targetLen, isRev, minAlignmentSpan, maxFlankExtensionDist, flankExtensionFactor
            1000, 1000, true, 200, 5000, 1.3,
            // expectedThrow
            true,
            // expectedRegions
            {
            },
        },
        TestData{
            "Seed hit strand doesn't match the strand specified via the function parameters.",
            {   // sortedHits
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                SeedHit(0, true, 5, 5, 15, 15, 0),
                SeedHit(0, true, 400, 400, 15, 15, 0),
                SeedHit(0, true, 900, 900, 15, 15, 0),
            },
            // queryLen, targetLen, isRev, minAlignmentSpan, maxFlankExtensionDist, flankExtensionFactor
            1000, 1000, false, 200, 5000, 1.3,
            // expectedThrow
            true,
            // expectedRegions
            {
            },
        },
        TestData{
            "Chain of seed hits is not monotonically increasing.",
            {   // sortedHits
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                SeedHit(0, true, 5, 5, 15, 15, 0),
                SeedHit(0, true, 900, 900, 15, 15, 0),
                SeedHit(0, true, 400, 400, 15, 15, 0),
            },
            // queryLen, targetLen, isRev, minAlignmentSpan, maxFlankExtensionDist, flankExtensionFactor
            1000, 1000, true, 200, 5000, 1.3,
            // expectedThrow
            true,
            // expectedRegions
            {
            },
        },
        TestData{
            "Seed hit + seed span are out of bounds of the query/target lengths. (995 + 15 > 1000)",
            {   // sortedHits
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                SeedHit(0, true, 100, 100, 15, 15, 0),
                SeedHit(0, true, 600, 600, 15, 15, 0),
                SeedHit(0, true, 995, 995, 15, 15, 0),
            },
            // queryLen, targetLen, isRev, minAlignmentSpan, maxFlankExtensionDist, flankExtensionFactor
            1000, 1000, true, 200, 5000, 1.3,
            // expectedThrow
            true,
            // expectedRegions
            {
            },
        },
    };
    // clang-format on

    for (const auto& data : allTests) {
        // Name the test.
        SCOPED_TRACE(data.testName);

        if (data.expectedThrow) {
            EXPECT_THROW(
                {
                    std::vector<AlignmentRegion> result = ExtractAlignmentRegions(
                        data.sortedHits, data.queryLen, data.targetLen, data.isRev,
                        data.minAlignmentSpan, data.maxFlankExtensionDist,
                        data.flankExtensionFactor);
                },
                std::runtime_error);

        } else {
            // Run the unit under test.
            std::vector<AlignmentRegion> result = ExtractAlignmentRegions(
                data.sortedHits, data.queryLen, data.targetLen, data.isRev, data.minAlignmentSpan,
                data.maxFlankExtensionDist, data.flankExtensionFactor);

            // std::cerr << "Test name: " << data.testName << "\n";
            // std::cerr << "Results:\n";
            // for (size_t i = 0; i < result.size(); ++i) {
            //     std::cerr << "[i = " << i << "] " << result[i] << "\n";
            // }
            // std::cerr << "Expected:\n";
            // for (size_t i = 0; i < data.expectedRegions.size(); ++i) {
            //     std::cerr << "[i = " << i << "] " << data.expectedRegions[i] << "\n";
            // }
            // std::cerr << "\n";

            // Evaluate.
            EXPECT_EQ(data.expectedRegions, result);
        }
    }
}

TEST(AlignmentSeeded, AlignSingleRegion_ArrayOfTests)
{
    struct TestData
    {
        std::string testName;
        std::string querySeq;
        std::string targetSeq;
        AlignmentRegion region;
        bool expectedThrow = false;
        PacBio::Data::Cigar expectedCigar;
        int32_t expectedLastQueryPos = -1;   // Last aligned query position within this region.
        int32_t expectedLastTargetPos = -1;  // Last aligned target position within this region.
        bool expectedValid = false;
    };

    // clang-format off
    std::vector<TestData> allTests{
        TestData{
            "Empty input", "", "", {}, false, PacBio::Data::Cigar(), 0, 0, false
        },

        // TestData{
        //     "Invalid input, NULL query",
        //     NULL,
        //     "AAAATCCCCCTGTTTGGGGG",
        //     {0, 5, 0, 5, false, RegionType::FRONT, 0},
        //     true, PacBio::Data::Cigar(), 0, 0, false
        // },

        // TestData{
        //     "Invalid input, NULL target",
        //     "AAAAACCCCCTTTTTGGGGG",
        //     NULL,
        //     {0, 5, 0, 5, false, RegionType::FRONT, 0},
        //     true, PacBio::Data::Cigar(), 0, 0, false
        // },

        TestData{
            "Invalid input, invalid region",
            "AAAAACCCCCTTTTTGGGGG",
            "AAAATCCCCCTGTTTGGGGG",
            {0, -1, 0, -1, false, RegionType::FRONT, 0},
            true, PacBio::Data::Cigar(), 0, 0, false
        },

        TestData{
            "Aligning the front",
            "AAAAACCCCCTTTTTGGGGG",
            "AAAATCCCCCTGTTTGGGGG",
            {0, 5, 0, 5, false, RegionType::FRONT, 0},
            false, PacBio::Data::Cigar("4=1X"), 5, 5, true
        },

        TestData{
            "Aligning the middle, globally",
            "AAAAACCCCCTTTTTGGGGG",
            "AAAATCCCCCTGTTTGGGGG",
            {5, 10, 5, 10, false, RegionType::GLOBAL, 0},
            false, PacBio::Data::Cigar("6=1X3="), 10, 10, true
        },

        TestData{
            "Aligning the back",
            "AAAAACCCCCTTTTTGGGGG",
            "AAAATCCCCCTGTTTGGGGG",
            {10, 10, 10, 10, false, RegionType::BACK, 0},
            false, PacBio::Data::Cigar("1=1X8="), 10, 10, true
        },

        TestData{
            "Aligning the middle, globally, reverse complemented.",
            "CCCCCAAAAAGGGGGTTTTT",
            "AAAATCCCCCTGTTTGGGGG",
            {5, 10, 5, 10, true, RegionType::GLOBAL, 0},
            false, PacBio::Data::Cigar("6=1X3="), 10, 10, true
        },

        TestData{
            "Aligning the front, reversed",
            "CCCCCAAAAAGGGGGTTTTT",
            "AAAATCCCCCTGTTTGGGGG",
            {0, 5, 0, 5, true, RegionType::FRONT, 0},
            false, PacBio::Data::Cigar("4=1X"), 5, 5, true
        },

        TestData{
            "Aligning the back, reversed",
            "CCCCCAAAAAGGGGGTTTTT",
            "AAAATCCCCCTGTTTGGGGG",
            {10, 10, 10, 10, true, RegionType::BACK, 0},
            false, PacBio::Data::Cigar("1=1X8="), 10, 10, true
        },

    };
    // clang-format on

    // Create any aligner, that's not important for testing of the unit under test.
    AlignerType alignerTypeGlobal = AlignerType::KSW2;
    AlignmentParameters alnParamsGlobal;
    AlignerType alignerTypeExt = AlignerType::KSW2;
    AlignmentParameters alnParamsExt;
    auto alignerGlobal = AlignerFactory(alignerTypeGlobal, alnParamsGlobal);
    auto alignerExt = AlignerFactory(alignerTypeExt, alnParamsExt);

    for (const auto& data : allTests) {
        // Name the test.
        SCOPED_TRACE(data.testName);

        const std::string querySeqRev =
            PacBio::Pancake::ReverseComplement(data.querySeq, 0, data.querySeq.size());

        if (data.expectedThrow) {
            EXPECT_THROW(
                {
                    AlignSingleRegion(data.targetSeq.c_str(), data.targetSeq.size(),
                                      data.querySeq.c_str(), querySeqRev.c_str(),
                                      data.querySeq.size(), alignerGlobal, alignerExt, data.region);
                },
                std::runtime_error);

        } else {
            // Run the unit under test.
            AlignmentResult result = AlignSingleRegion(
                data.targetSeq.c_str(), data.targetSeq.size(), data.querySeq.c_str(),
                querySeqRev.c_str(), data.querySeq.size(), alignerGlobal, alignerExt, data.region);

            // std::cerr << "Test name: " << data.testName << "\n";
            // std::cerr << "result.cigar = " << result.cigar.ToStdString() << "\n";
            // std::cerr << "result.lastQueryPos = " << result.lastQueryPos << "\n";
            // std::cerr << "result.lastTargetPos = " << result.lastTargetPos << "\n";
            // std::cerr << "result.valid = " << result.valid << "\n";
            // std::cerr << "result = " << result << "\n";

            // std::cerr << "data.expectedCigar = " << data.expectedCigar.ToStdString() << "\n";
            // std::cerr << "data.expectedLastQueryPos = " << data.expectedLastQueryPos << "\n";
            // std::cerr << "data.expectedLastTargetPos = " << data.expectedLastTargetPos << "\n";
            // std::cerr << "data.expectedValid = " << data.expectedValid << "\n";
            // std::cerr << "\n";

            // Evaluate.
            EXPECT_EQ(data.expectedCigar, result.cigar);
            EXPECT_EQ(data.expectedLastQueryPos, result.lastQueryPos);
            EXPECT_EQ(data.expectedLastTargetPos, result.lastTargetPos);
            EXPECT_EQ(data.expectedValid, result.valid);
        }
    }
}

TEST(AlignmentSeeded, AlignmentSeeded_ArrayOfTests)
{
    struct TestData
    {
        std::string testName;
        std::string querySeq;
        std::string targetSeq;
        PacBio::Pancake::Overlap ovl;
        std::vector<PacBio::Pancake::SeedHit> sortedHits;
        int32_t minAlignmentSpan = 200;
        int32_t maxFlankExtensionDist = 5000;
        double flankExtensionFactor = 1.3;
        bool expectedThrow = false;
        PacBio::Pancake::Overlap expectedOvl;
    };

    // clang-format off
    std::vector<TestData> allTests{
        TestData{
            "Empty input, should throw", "", "", PacBio::Pancake::Overlap(), {}, 200, 5000, 1.3, true, PacBio::Pancake::Overlap()
        },

        TestData{
            "Alignment forward, the internal hit will be skipped.",
            "AAAAACCCCCTTTTTGGGGG",
            "AAAATCCCCCTGTTTGGGGG",

            PacBio::Pancake::Overlap(0, 0, 0.0, 0.0, false, 5, 16, 20, false, 5, 16, 20),
            {   // sortedHits
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                SeedHit(0, false, 5, 5, 5, 5, 0),
                SeedHit(0, false, 10, 10, 5, 5, 0),
                SeedHit(0, false, 16, 16, 4, 4, 0),
            },
            // minAlignmentSpan, maxFlankExtensionDist, flankExtensionFactor
            200, 5000, 1.3,
            // expectedThrow
            false,
            // expectedOvl
            PacBio::Pancake::Overlap(
                0, 0, 28, 0.90, false, 0, 20, 20, false, 0, 20, 20,
                2, 0, OverlapType::Unknown, OverlapType::Unknown,           // editDist, numSeeds, aType, bType
                PacBio::Data::Cigar("4=1X6=1X8="),                          // cigar
                "", "", false, false, false                                 // aVars, bVars, isFlipped, isSupplementary, isSecondary
            ),
        },

        TestData{
            "Alignment forward, using the internal hit too.",
            "AAAAACCCCCTTTTTGGGGG",
            "AAAATCCCCCTGTTTGGGGG",

            PacBio::Pancake::Overlap(0, 0, 0.0, 0.0, false, 5, 16, 20, false, 5, 16, 20),
            {   // sortedHits
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                SeedHit(0, false, 5, 5, 5, 5, 0),
                SeedHit(0, false, 10, 10, 5, 5, 0),
                SeedHit(0, false, 16, 16, 4, 4, 0),
            },
            // minAlignmentSpan, maxFlankExtensionDist, flankExtensionFactor
            0, 5000, 1.3,
            // expectedThrow
            false,
            // expectedOvl
            PacBio::Pancake::Overlap(
                0, 0, 28, 0.90, false, 0, 20, 20, false, 0, 20, 20,
                2, 0, OverlapType::Unknown, OverlapType::Unknown,           // editDist, numSeeds, aType, bType
                PacBio::Data::Cigar("4=1X6=1X8="),                          // cigar
                "", "", false, false, false                                 // aVars, bVars, isFlipped, isSupplementary, isSecondary
            ),
        },

        TestData{
            "Alignment reverse, the internal hit will be skipped.",
            "CCCCCAAAAAGGGGGTTTTT",
            "AAAATCCCCCTGTTTGGGGG",

            PacBio::Pancake::Overlap(0, 0, 0.0, 0.0, false, 5, 16, 20, true, 5, 16, 20),
            {   // sortedHits
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                SeedHit(0, true, 5, 5, 5, 5, 0),
                SeedHit(0, true, 10, 10, 5, 5, 0),
                SeedHit(0, true, 16, 16, 4, 4, 0),
            },
            // minAlignmentSpan, maxFlankExtensionDist, flankExtensionFactor
            200, 5000, 1.3,
            // expectedThrow
            false,
            // expectedOvl
            PacBio::Pancake::Overlap(
                0, 0, 28, 0.90, false, 0, 20, 20, true, 0, 20, 20,
                2, 0, OverlapType::Unknown, OverlapType::Unknown,           // editDist, numSeeds, aType, bType
                PacBio::Data::Cigar("8=1X6=1X4="),                          // cigar
                "", "", false, false, false                                 // aVars, bVars, isFlipped, isSupplementary, isSecondary
            ),
        },

        TestData{
            "Alignment reverse, using the internal hit too.",
            "CCCCCAAAAAGGGGGTTTTT",
            "AAAATCCCCCTGTTTGGGGG",

            PacBio::Pancake::Overlap(0, 0, 0.0, 0.0, false, 5, 16, 20, true, 5, 16, 20),
            {   // sortedHits
                // targetId, targetRev, targetPos, queryPos, targetSpan, querySpan, flags
                SeedHit(0, true, 5, 5, 5, 5, 0),
                SeedHit(0, true, 10, 10, 5, 5, 0),
                SeedHit(0, true, 16, 16, 4, 4, 0),
            },
            // minAlignmentSpan, maxFlankExtensionDist, flankExtensionFactor
            0, 5000, 1.3,
            // expectedThrow
            false,
            // expectedOvl
            PacBio::Pancake::Overlap(
                0, 0, 28, 0.90, false, 0, 20, 20, true, 0, 20, 20,
                2, 0, OverlapType::Unknown, OverlapType::Unknown,           // editDist, numSeeds, aType, bType
                PacBio::Data::Cigar("8=1X6=1X4="),                          // cigar
                "", "", false, false, false                                 // aVars, bVars, isFlipped, isSupplementary, isSecondary
            ),
        },

        TestData{
            "Wrong overlap compared to the seed hits. Should throw.",
            "AAAAACCCCCTTTTTGGGGG",
            "AAAATCCCCCTGTTTGGGGG",

            PacBio::Pancake::Overlap(),
            {   // sortedHits
            },
            // minAlignmentSpan, maxFlankExtensionDist, flankExtensionFactor
            200, 5000, 1.3,
            // expectedThrow
            true,
            // expectedOvl
            PacBio::Pancake::Overlap()
        },

    };
    // clang-format on

    // Create any aligner, that's not important for testing of the unit under test.
    AlignerType alignerTypeGlobal = AlignerType::KSW2;
    AlignmentParameters alnParamsGlobal;
    AlignerType alignerTypeExt = AlignerType::KSW2;
    AlignmentParameters alnParamsExt;
    auto alignerGlobal = AlignerFactory(alignerTypeGlobal, alnParamsGlobal);
    auto alignerExt = AlignerFactory(alignerTypeExt, alnParamsExt);

    for (const auto& data : allTests) {
        // Name the test.
        SCOPED_TRACE(data.testName);

        const std::string querySeqRev =
            PacBio::Pancake::ReverseComplement(data.querySeq, 0, data.querySeq.size());
        auto ovl = createOverlap(data.ovl);

        std::vector<AlignmentRegion> regions = ExtractAlignmentRegions(
            data.sortedHits, data.querySeq.size(), data.targetSeq.size(), data.ovl.Brev,
            data.minAlignmentSpan, data.maxFlankExtensionDist, data.flankExtensionFactor);

        // std::cerr << "Test name: " << data.testName << "\n";
        // std::cerr << "Input overlap: "<< OverlapWriterBase::PrintOverlapAsM4(data.ovl, "", "", true, true) << "\n";

        if (data.expectedThrow) {
            EXPECT_THROW(
                {
                    AlignmentSeeded(ovl, regions, data.targetSeq.c_str(), data.targetSeq.size(),
                                    data.querySeq.c_str(), querySeqRev.c_str(),
                                    data.querySeq.size(), alignerGlobal, alignerExt);
                },
                std::runtime_error);

        } else {

            // Run the unit under test.
            OverlapPtr result = AlignmentSeeded(
                ovl, regions, data.targetSeq.c_str(), data.targetSeq.size(), data.querySeq.c_str(),
                querySeqRev.c_str(), data.querySeq.size(), alignerGlobal, alignerExt);

            // std::cerr << "Expected overlap:\n"
            //           << OverlapWriterBase::PrintOverlapAsM4(data.expectedOvl, "", "", true, true)
            //           << "\n";
            // std::cerr << "Resulting overlap:\n"
            //           << OverlapWriterBase::PrintOverlapAsM4(result, "", "", true, true) << "\n";
            // std::cerr << "\n";

            // Evaluate.
            EXPECT_EQ(data.expectedOvl, *result);
        }
    }
}