// Authors: Ivan Sovic

#include <gtest/gtest.h>
#include <cstdint>
#include <pancake/util/RunLengthEncoding.hpp>
#include <string>
#include <tuple>
#include <vector>

namespace RLECompressTests {

TEST(Test_RunLengthEncoding, EmptyInput)
{
    std::string input("");
    std::string expectedStr("");
    std::vector<int32_t> expectedCounts;

    std::string resultStr;
    std::vector<int32_t> resultCounts;
    PacBio::Pancake::RunLengthEncoding(input, resultStr, resultCounts);

    EXPECT_EQ(expectedStr, resultStr);
    EXPECT_EQ(expectedCounts, resultCounts);
}

TEST(Test_RunLengthEncoding, SingleBase)
{
    std::string input("A");
    std::string expectedStr("A");
    std::vector<int32_t> expectedCounts = {1};

    std::string resultStr;
    std::vector<int32_t> resultCounts;
    PacBio::Pancake::RunLengthEncoding(input, resultStr, resultCounts);

    EXPECT_EQ(expectedStr, resultStr);
    EXPECT_EQ(expectedCounts, resultCounts);
}

TEST(Test_RunLengthEncoding, SimpleTestCase)
{
    std::string input("ACCTGTGGGCAA");
    std::string expectedStr("ACTGTGCA");
    std::vector<int32_t> expectedCounts = {1, 2, 1, 1, 1, 3, 1, 2};

    std::string resultStr;
    std::vector<int32_t> resultCounts;
    PacBio::Pancake::RunLengthEncoding(input, resultStr, resultCounts);

    EXPECT_EQ(expectedStr, resultStr);
    EXPECT_EQ(expectedCounts, resultCounts);
}

TEST(Test_RunLengthEncoding, SeveralTestCases_From_in_2_small_fasta)
{
    //clang-format off
    std::vector<std::tuple<std::string, std::string, std::vector<int32_t>>> testSets = {
        {"ATCGGTTCAGGCCGCG",
         "ATCGTCAGCGCG",
         {
             1,
             1,
             1,
             2,
             2,
             1,
             1,
             2,
             2,
             1,
             1,
             1,
         }},
        {"AGCCGCCAGCAGAACAATACCGAATAATGCCA",
         "AGCGCAGCAGACATACGATATGCA",
         {1, 1, 2, 1, 2, 1, 1, 1, 1, 1, 2, 1, 2, 1, 1, 2, 1, 2, 1, 2, 1, 1, 2, 1}},
        {"CC", "C", {2}},
        {"A", "A", {1}},
        {"ACTTCGCGGCAAGTGGACATCGCCAGAGGTATAAAAGCAGGCGTTAACCCAGCGTAACGGTGG",
         "ACTCGCGCAGTGACATCGCAGAGTATAGCAGCGTACAGCGTACGTG",
         {1, 1, 2, 1, 1, 1, 2, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 2,
          1, 1, 1, 4, 1, 1, 1, 2, 1, 1, 2, 2, 3, 1, 1, 1, 1, 1, 2, 1, 2, 1, 2}},
    };
    //clang-format onf

    for (const auto& singleTest : testSets) {
        const auto& input = std::get<0>(singleTest);
        const auto& expectedStr = std::get<1>(singleTest);
        const auto& expectedRLE = std::get<2>(singleTest);
        SCOPED_TRACE(input);

        std::string resultStr;
        std::vector<int32_t> resultRLE;
        PacBio::Pancake::RunLengthEncoding(input, resultStr, resultRLE);

        EXPECT_EQ(expectedStr, resultStr);
        EXPECT_EQ(expectedRLE, resultRLE);
    }
}

TEST(Test_RunLengthEncoding, SimpleTestCase_WithCoords)
{
    std::string input("ACCTGTGGGCAA");
    std::string expectedStr("ACTGTGCA");
    std::vector<int32_t> expectedSeqToHPC = {0, 1, 1, 2, 3, 4, 5, 5, 5, 6, 7, 7};
    std::vector<int32_t> expectedHPCToSeq = {0, 2, 3, 4, 5, 8, 9, 11};

    std::string resultStr;
    std::vector<int32_t> resultSeqToHPCCoords;
    std::vector<int32_t> resultHPCToSeqCoords;
    PacBio::Pancake::RunLengthEncoding(input, resultStr, resultSeqToHPCCoords,
                                       resultHPCToSeqCoords);

    EXPECT_EQ(expectedStr, resultStr);
    EXPECT_EQ(expectedSeqToHPC, resultSeqToHPCCoords);
    EXPECT_EQ(expectedHPCToSeq, resultHPCToSeqCoords);
}
}  // namespace RLECompressTests
