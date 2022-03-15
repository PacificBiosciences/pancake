#include <gtest/gtest.h>
#include <pancake/Twobit.hpp>
#include <tuple>

void HelperRoundTrip_DecompressCPPStyle(const std::string& inBases, int32_t numBases,
                                        const std::vector<uint8_t>& expectedTwobit,
                                        const std::vector<PacBio::Pancake::Range>& expectedRanges,
                                        int32_t expectedCompressedBases, bool expectedException)
{
    // Compress the input sequence.
    std::vector<uint8_t> twobit;
    std::vector<PacBio::Pancake::Range> ranges;
    int32_t numCompressedBases = PacBio::Pancake::CompressSequence(inBases, twobit, ranges);

    // Check the compression results.
    EXPECT_EQ(expectedTwobit, twobit);
    EXPECT_EQ(expectedRanges, ranges);
    EXPECT_EQ(expectedCompressedBases, numCompressedBases);

    // Decompress the sequence.
    // Run round trip, expected to see the same sequence or an exception.
    bool exceptionCaught = false;
    try {
        std::string roundResult;
        PacBio::Pancake::DecompressSequence(twobit, numBases, ranges, roundResult);
        EXPECT_EQ(inBases, roundResult);
    } catch (const std::runtime_error& e) {
        exceptionCaught = true;
    }
    EXPECT_EQ(expectedException, exceptionCaught);
}

void HelperRoundTrip_DecompressCStyle(const std::string& inBases, int32_t numBases,
                                      int32_t numBasesToAlloc,
                                      const std::vector<uint8_t>& expectedTwobit,
                                      const std::vector<PacBio::Pancake::Range>& expectedRanges,
                                      int32_t expectedCompressedBases, bool expectedException)
{
    /*
     * The numBasesToAlloc will be used to alocate the storage for decompressedData.
     * In realistic use cases, this should be the same as numBases, but here we intentionally
     * modify them differently so we can test for exceptions.
    */

    // Compress the input sequence.
    std::vector<uint8_t> twobit;
    std::vector<PacBio::Pancake::Range> ranges;
    int32_t numCompressedBases = PacBio::Pancake::CompressSequence(inBases, twobit, ranges);

    // Check the compression results.
    EXPECT_EQ(expectedTwobit, twobit);
    EXPECT_EQ(expectedRanges, ranges);
    EXPECT_EQ(expectedCompressedBases, numCompressedBases);

    // Decompress the sequence.
    // Run round trip, expected to see the same sequence or an exception.
    bool exceptionCaught = false;
    try {
        std::vector<uint8_t> decompressedData(numBasesToAlloc);
        PacBio::Pancake::DecompressSequence(twobit, numBases, ranges, decompressedData.data());
        std::string roundResult(reinterpret_cast<char*>(decompressedData.data()), numBases);
        EXPECT_EQ(inBases, roundResult);
    } catch (const std::runtime_error& e) {
        exceptionCaught = true;
    }
    EXPECT_EQ(expectedException, exceptionCaught);
}

struct TestData
{
    std::string testName = "";
    std::string inBases = "";
    int32_t numBases = 0;
    int32_t numBasesToAlloc = 0;
    std::vector<uint8_t> expectedTwobit;
    std::vector<PacBio::Pancake::Range> expectedRanges;
    int32_t expectedComprBases = 0;
    bool expectedThrow = false;
};

TEST(Twobit_RoundTrip, Decompress)
{
    std::vector<TestData> inputs = {
        TestData{"Empty input", "", 0, 0, {}, {}, 0, false},
        TestData{"Simple_1_Base_NoN_1", "A", 1, 1, {0}, {{0, 1}}, 1, false},
        TestData{"Simple_1_Base_NoN_2", "G", 1, 1, {128}, {{0, 1}}, 1, false},
        TestData{"Simple_2_Bases_NoN", "GC", 2, 2, {144}, {{0, 2}}, 2, false},
        TestData{"Simple_3_Bases_NoN", "GCT", 3, 3, {156}, {{0, 3}}, 3, false},
        TestData{"Simple_4_Bases_NoN", "ACGT", 4, 4, {27}, {{0, 4}}, 4, false},
        TestData{"Simple_WithN_1", "ACGTNTTT", 8, 8, {27, 252}, {{0, 4}, {5, 8}}, 7, false},
        TestData{"Simple_WithN_2", "GCTNTAT", 7, 7, {159, 48}, {{0, 3}, {4, 7}}, 6, false},
        TestData{"Simple_WithN_In_Last_Pos", "GCATNNN", 7, 7, {147}, {{0, 4}}, 4, false},
        TestData{"Simple_WithN_In_FirstPos", "NNNGCAT", 7, 7, {147}, {{3, 7}}, 4, false},
        TestData{"DecompressSanityCheck1", "ACGT", 4, 4, {27}, {{0, 4}}, 4, false},
        TestData{"DecompressSanityCheck2", "ACGT", -1, 0, {27}, {{0, 4}}, 4, true},
        TestData{"Wrong number of bases unpacked", "ACGT", 3, 10, {27}, {{0, 4}}, 4, true},
    };

    for (const auto& testData : inputs) {
        SCOPED_TRACE(testData.testName);
        HelperRoundTrip_DecompressCPPStyle(testData.inBases, testData.numBases,
                                           testData.expectedTwobit, testData.expectedRanges,
                                           testData.expectedComprBases, testData.expectedThrow);
    }
}

TEST(Twobit_RoundTrip, DecompressCStyle_MultipleTests)
{
    std::vector<TestData> inputs = {
        TestData{"Empty input", "", 0, 0, {}, {}, 0, false},
        TestData{"Simple_1_Base_NoN_1", "A", 1, 1, {0}, {{0, 1}}, 1, false},
        TestData{"Simple_1_Base_NoN_2", "G", 1, 1, {128}, {{0, 1}}, 1, false},
        TestData{"Simple_2_Bases_NoN", "GC", 2, 2, {144}, {{0, 2}}, 2, false},
        TestData{"Simple_3_Bases_NoN", "GCT", 3, 3, {156}, {{0, 3}}, 3, false},
        TestData{"Simple_4_Bases_NoN", "ACGT", 4, 4, {27}, {{0, 4}}, 4, false},
        TestData{"Simple_WithN_1", "ACGTNTTT", 8, 8, {27, 252}, {{0, 4}, {5, 8}}, 7, false},
        TestData{"Simple_WithN_2", "GCTNTAT", 7, 7, {159, 48}, {{0, 3}, {4, 7}}, 6, false},
        TestData{"Simple_WithN_In_Last_Pos", "GCATNNN", 7, 7, {147}, {{0, 4}}, 4, false},
        TestData{"Simple_WithN_In_FirstPos", "NNNGCAT", 7, 7, {147}, {{3, 7}}, 4, false},
        TestData{"DecompressSanityCheck1", "ACGT", 4, 4, {27}, {{0, 4}}, 4, false},
        TestData{"DecompressSanityCheck2", "ACGT", -1, 0, {27}, {{0, 4}}, 4, true},
        TestData{"Wrong number of bases unpacked", "ACGT", 3, 10, {27}, {{0, 4}}, 4, true},
    };

    for (const auto& testData : inputs) {
        SCOPED_TRACE(testData.testName);
        HelperRoundTrip_DecompressCStyle(
            testData.inBases, testData.numBases, testData.numBasesToAlloc, testData.expectedTwobit,
            testData.expectedRanges, testData.expectedComprBases, testData.expectedThrow);
    }
}
