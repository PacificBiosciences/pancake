// Authors: Ivan Sovic

#include <PancakeTestData.h>
#include <gtest/gtest.h>
#include <pancake/SesAlignBanded.hpp>
#include <pancake/SesDistanceBanded.hpp>
#include <sstream>
#include <tuple>

TEST(SesDistanceBanded, SesDistasnceBanded_EmptyQueryEmptyTarget)
{
    // Inputs.
    std::string query = "";
    std::string target = "";
    int32_t maxDiffs = 100;
    int32_t bandwidth = 30;

    // Order of elements in SesResult: lastQueryPos, lastTargetPos, diffs, minK, maxK, valid
    PacBio::Pancake::Alignment::SesResults expected = {0, 0, 0, true};

    // Run.
    PacBio::Pancake::Alignment::SesResults result = PacBio::Pancake::Alignment::SESDistanceBanded(
        query.c_str(), query.size(), target.c_str(), target.size(), maxDiffs, bandwidth);

    // Evaluate.
    EXPECT_EQ(expected, result);
}

TEST(SesDistanceBanded, SesDistasnceBanded_EmptyQueryNonemptyTarget)
{
    // Inputs.
    std::string query = "";
    std::string target = "ACTG";
    int32_t maxDiffs = 100;
    int32_t bandwidth = 30;

    // Order of elements in SesResult: lastQueryPos, lastTargetPos, diffs, minK, maxK, valid
    PacBio::Pancake::Alignment::SesResults expected = {0, 0, 0, true};

    // Run.
    PacBio::Pancake::Alignment::SesResults result = PacBio::Pancake::Alignment::SESDistanceBanded(
        query.c_str(), query.size(), target.c_str(), target.size(), maxDiffs, bandwidth);

    // Evaluate.
    EXPECT_EQ(expected, result);
}

TEST(SesDistanceBanded, SesDistasnceBanded_NonEmptyQueryEmptyTarget)
{
    // Inputs.
    std::string query = "ACTG";
    std::string target = "";
    int32_t maxDiffs = 100;
    int32_t bandwidth = 30;

    // Order of elements in SesResult: lastQueryPos, lastTargetPos, diffs, minK, maxK, valid
    PacBio::Pancake::Alignment::SesResults expected = {0, 0, 0, true};

    // Run.
    PacBio::Pancake::Alignment::SesResults result = PacBio::Pancake::Alignment::SESDistanceBanded(
        query.c_str(), query.size(), target.c_str(), target.size(), maxDiffs, bandwidth);

    // Evaluate.
    EXPECT_EQ(expected, result);
}

TEST(SesDistanceBanded, SesDistasnceBanded_SimpleExactMatch)
{
    // Inputs.
    std::string query = "ACTG";
    std::string target = "ACTG";
    int32_t maxDiffs = 100;
    int32_t bandwidth = 30;

    // Order of elements in SesResult: lastQueryPos, lastTargetPos, diffs, minK, maxK, valid
    PacBio::Pancake::Alignment::SesResults expected = {4, 4, 0, true};

    // Run.
    PacBio::Pancake::Alignment::SesResults result = PacBio::Pancake::Alignment::SESDistanceBanded(
        query.c_str(), query.size(), target.c_str(), target.size(), maxDiffs, bandwidth);

    // Evaluate.
    EXPECT_EQ(expected, result);
}

TEST(SesDistanceBanded, SesDistasnceBanded_SimpleSingleDiff)
{
    // Inputs.
    std::string query = "ACG";
    std::string target = "ACTG";
    int32_t maxDiffs = 100;
    int32_t bandwidth = 30;

    // Order of elements in SesResult: lastQueryPos, lastTargetPos, diffs, minK, maxK, valid
    PacBio::Pancake::Alignment::SesResults expected = {3, 4, 1, true};

    // Run.
    PacBio::Pancake::Alignment::SesResults result = PacBio::Pancake::Alignment::SESDistanceBanded(
        query.c_str(), query.size(), target.c_str(), target.size(), maxDiffs, bandwidth);

    // Evaluate.
    EXPECT_EQ(expected, result);
}

TEST(SesDistanceBanded, SesDistasnceBanded_SimpleTwoDiffs)
{
    // Inputs.
    std::string query = "AAAAAAAAAAAAAAAAAAAA";
    std::string target = "AAAAAAGGGGGAAAAAAAAAAAAAA";
    int32_t maxDiffs = 100;
    int32_t bandwidth = 30;

    // Order of elements in SesResult: lastQueryPos, lastTargetPos, diffs, minK, maxK, valid
    PacBio::Pancake::Alignment::SesResults expected = {20, 25, 5, true};

    // Run.
    PacBio::Pancake::Alignment::SesResults result = PacBio::Pancake::Alignment::SESDistanceBanded(
        query.c_str(), query.size(), target.c_str(), target.size(), maxDiffs, bandwidth);

    // Evaluate.
    EXPECT_EQ(expected, result);
}

TEST(SesDistanceBanded, SesDistasnceBanded_NormalSmallCase)
{
    /*
     * Note: Mismatches count as 2 diffs.
    */
    // Inputs.
    std::string query = "GGATCAGTT";   // GGAT-CAGTT
                                       // |X|| | |||
    std::string target = "GAATTCGTT";  // GAATTC-GTT
    int32_t maxDiffs = 100;
    int32_t bandwidth = 30;

    // Order of elements in SesResult: lastQueryPos, lastTargetPos, diffs, minK, maxK, valid
    PacBio::Pancake::Alignment::SesResults expected = {9, 9, 4, true};

    // Run.
    PacBio::Pancake::Alignment::SesResults result = PacBio::Pancake::Alignment::SESDistanceBanded(
        query.c_str(), query.size(), target.c_str(), target.size(), maxDiffs, bandwidth);

    // Evaluate.
    EXPECT_EQ(expected, result);
}

TEST(SesDistanceBanded, SesDistasnceBanded_OutOfBandwidth)
{
    /*
     * Note: Mismatches count as 2 diffs.
    */
    // Inputs.
    std::string query = "GGATCAGTT";   // GGAT-CAGTT
                                       // |X|| | |||
    std::string target = "GAATTCGTT";  // GAATTC-GTT
    int32_t maxDiffs = 100;
    int32_t bandwidth = 3;

    // Order of elements in SesResult: lastQueryPos, lastTargetPos, diffs, minK, maxK, valid
    PacBio::Pancake::Alignment::SesResults expected = {3, 2, 2, false};

    // Run.
    PacBio::Pancake::Alignment::SesResults result = PacBio::Pancake::Alignment::SESDistanceBanded(
        query.c_str(), query.size(), target.c_str(), target.size(), maxDiffs, bandwidth);

    // Evaluate.
    EXPECT_EQ(expected, result);
}

TEST(SesDistanceBanded, SesDistasnceBanded_AboveMaxDiffs)
{
    /*
     * Note: Mismatches count as 2 diffs.
    */
    // Inputs.
    std::string query = "GGATCAGTT";   // GGAT-CAGTT
                                       // |X|| | |||
    std::string target = "GAATTCGTT";  // GAATTC-GTT
    int32_t maxDiffs = 3;
    int32_t bandwidth = 30;

    // Order of elements in SesResult: lastQueryPos, lastTargetPos, diffs, minK, maxK, valid
    PacBio::Pancake::Alignment::SesResults expected = {3, 2, 2, false};

    // Run.
    PacBio::Pancake::Alignment::SesResults result = PacBio::Pancake::Alignment::SESDistanceBanded(
        query.c_str(), query.size(), target.c_str(), target.size(), maxDiffs, bandwidth);

    // Evaluate.
    EXPECT_EQ(expected, result);
}
