#include <PancakeTestData.h>
#include <gtest/gtest.h>
#include <iostream>
#include <memory>
#include <pancake/MapperUtility.hpp>
#include "TestHelperUtils.hpp"

TEST(MapperUtility, ComputeOccurrenceThreshold_ArrayOfTests)
{
    // clang-format off
    struct TestData
    {
        const std::string testName;
        const std::vector<std::pair<int64_t, int64_t>> seedHitHistogram;
        const int64_t seedOccurrenceMin = 0;
        const int64_t seedOccurrenceMax = 0;
        const int64_t seedOccurrenceMaxMemory = 0;
        const int64_t seedOccurrencePercentileCutoff = 0;
        const int64_t expected = 0;
    };

    // const int32_t span = 15;
    // const int32_t targetId = 0;
    const std::vector<TestData> testData = {
        {   // Test 1.
            "Empty input",
            // seedHitHistogram
            {},
            // seedOccurrenceMin
            0,
            // seedOccurrenceMax
            0,
            // seedOccurrenceMaxMemory
            0,
            // seedOccurrencePercentileCutoff
            0,
            // Expected.
            std::numeric_limits<int64_t>::max()
        },

        {   // Test 2.
            "No thresholds. Histogram from a real low-complexity dataset: Target = 'm64004_191221_060015/1441/134152_156483', query = 'm64004_191221_060015/1441/156528_178747'",
            // seedHitHistogram
            {
                {0, 341}, {1, 76}, {2, 48}, {3, 9}, {29, 33}, {30, 33}, {31, 64},
                {32, 36}, {38, 170}, {39, 61}, {40, 76}, {41, 101}, {42, 113},
                {43, 117}, {79, 84}, {80, 127}, {82, 114}, {8899, 6518},
            },
            // seedOccurrenceMin
            0,
            // seedOccurrenceMax
            0,
            // seedOccurrenceMaxMemory
            0,
            // seedOccurrencePercentileCutoff
            0,
            // Expected.
            std::numeric_limits<int64_t>::max()
        },

        {   // Test 3.
            "Maximum cutoff of 1000. Histogram from a real low-complexity dataset: Target = 'm64004_191221_060015/1441/134152_156483', query = 'm64004_191221_060015/1441/156528_178747'",
            // seedHitHistogram
            {
                {0, 341}, {1, 76}, {2, 48}, {3, 9}, {29, 33}, {30, 33}, {31, 64},
                {32, 36}, {38, 170}, {39, 61}, {40, 76}, {41, 101}, {42, 113},
                {43, 117}, {79, 84}, {80, 127}, {82, 114}, {8899, 6518},
            },
            // seedOccurrenceMin
            0,
            // seedOccurrenceMax
            1000,
            // seedOccurrenceMaxMemory
            0,
            // seedOccurrencePercentileCutoff
            0,
            // Expected.
            1000
        },

        {   // Test 4.
            "'Percentile' cutoff is lower than 1000, set to 500. Histogram from a real low-complexity dataset: Target = 'm64004_191221_060015/1441/134152_156483', query = 'm64004_191221_060015/1441/156528_178747'",
            // seedHitHistogram
            {
                {0, 341}, {1, 76}, {2, 48}, {3, 9}, {29, 33}, {30, 33}, {31, 64},
                {32, 36}, {38, 170}, {39, 61}, {40, 76}, {41, 101}, {42, 113},
                {43, 117}, {79, 84}, {80, 127}, {82, 114}, {8899, 6518},
            },
            // seedOccurrenceMin
            0,
            // seedOccurrenceMax
            1000,
            // seedOccurrenceMaxMemory
            0,
            // seedOccurrencePercentileCutoff
            500,
            // Expected.
            500
        },

        {   // Test 5.
            "The maximum allowed memory for the seed hits has the threshold set to 83 (82 + 1, because the next seed occurrence is 8899 which happens for 6518 seeds, and exceeeds the memory limit). Histogram from a real low-complexity dataset: Target = 'm64004_191221_060015/1441/134152_156483', query = 'm64004_191221_060015/1441/156528_178747'",
            // seedHitHistogram
            {
                {0, 341}, {1, 76}, {2, 48}, {3, 9}, {29, 33}, {30, 33}, {31, 64},
                {32, 36}, {38, 170}, {39, 61}, {40, 76}, {41, 101}, {42, 113},
                {43, 117}, {79, 84}, {80, 127}, {82, 114}, {8899, 6518},
            },
            // seedOccurrenceMin
            0,
            // seedOccurrenceMax
            1000,
            // seedOccurrenceMaxMemory
            1'600'000, // bytes
            // seedOccurrencePercentileCutoff
            500,
            // Expected.
            83
        },

        {   // Test 6.
            "Minimum occurrence cutoff threshold is set to prevent filtering seeds with too low occurrence. Histogram from a real low-complexity dataset: Target = 'm64004_191221_060015/1441/134152_156483', query = 'm64004_191221_060015/1441/156528_178747'",
            // seedHitHistogram
            {
                {0, 341}, {1, 76}, {2, 48}, {3, 9}, {29, 33}, {30, 33}, {31, 64},
                {32, 36}, {38, 170}, {39, 61}, {40, 76}, {41, 101}, {42, 113},
                {43, 117}, {79, 84}, {80, 127}, {82, 114}, {8899, 6518},
            },
            // seedOccurrenceMin
            100,
            // seedOccurrenceMax
            1000,
            // seedOccurrenceMaxMemory
            1'600'000, // bytes
            // seedOccurrencePercentileCutoff
            500,
            // Expected.
            100
        },

        {   // Test 7.
            "No minimum threshold, and use a very large maximum allowed memory, so that the cutoff is bound at the seedOccurrenceMax (1000). Histogram from a real low-complexity dataset: Target = 'm64004_191221_060015/1441/134152_156483', query = 'm64004_191221_060015/1441/156528_178747'",
            // seedHitHistogram
            {
                {0, 341}, {1, 76}, {2, 48}, {3, 9}, {29, 33}, {30, 33}, {31, 64},
                {32, 36}, {38, 170}, {39, 61}, {40, 76}, {41, 101}, {42, 113},
                {43, 117}, {79, 84}, {80, 127}, {82, 114}, {8899, 6518},
            },
            // seedOccurrenceMin
            0,
            // seedOccurrenceMax
            1000,
            // seedOccurrenceMaxMemory
            10'000'000'000, // bytes
            // seedOccurrencePercentileCutoff
            0,
            // Expected.
            1000
        },
    };
    // clang-format on

    for (const auto& data : testData) {
        SCOPED_TRACE(data.testName);
        std::cerr << "testName = " << data.testName << "\n";

        // Run unit under test.
        const int64_t results = PacBio::Pancake::ComputeOccurrenceThreshold(
            data.seedHitHistogram, data.seedOccurrenceMin, data.seedOccurrenceMax,
            data.seedOccurrenceMaxMemory, data.seedOccurrencePercentileCutoff, false);

        // Evaluate.
        EXPECT_EQ(data.expected, results);
    }
}
