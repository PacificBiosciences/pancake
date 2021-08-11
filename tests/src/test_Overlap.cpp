#include <gtest/gtest.h>
#include <pacbio/pancake/Overlap.h>
#include <pacbio/pancake/OverlapWriterBase.h>
#include <vector>

std::vector<std::string> WrapRunHeuristicExtendOverlapFlanks(
    const std::vector<std::string>& inOverlaps, int32_t allowedDist)
{
    std::vector<std::string> results;
    for (const auto& ovlStr : inOverlaps) {
        auto ovl = PacBio::Pancake::ParseM4OverlapFromString(ovlStr);
        // ovl->Type = PacBio::Pancake::DetermineOverlapType(ovl, allowedDist);
        PacBio::Pancake::HeuristicExtendOverlapFlanks(ovl, allowedDist);
        std::string resultStr =
            PacBio::Pancake::OverlapWriterBase::PrintOverlapAsM4(*ovl, "", "", true, false);
        results.emplace_back(resultStr);
    }
    return results;
}

TEST(ExtendOverlapFlanks, EmptyInput)
{
    // Inputs.
    std::vector<PacBio::Pancake::OverlapPtr> inOverlaps = {};
    int32_t allowedDist = 50;

    // Expected results.
    std::vector<PacBio::Pancake::Overlap> expected = {};

    std::vector<PacBio::Pancake::Overlap> results;
    for (const auto& ovl : inOverlaps) {
        results.emplace_back(*PacBio::Pancake::HeuristicExtendOverlapFlanks(ovl, allowedDist));
    }

    // Evaluate.
    EXPECT_EQ(expected, results);
}

TEST(ExtendOverlapFlanks, NormalWithAllowedDist50)
{
    // Inputs.
    int32_t allowedDist = 50;
    std::vector<std::string> inOverlaps = {
        "000000001 000000002 -1000 99.00 0 0 1000 5000 0 7000 8000 8000",
        "000000001 000000002 -1000 99.00 0 0 1000 5000 1 0 1000 8000",
        "000000001 000000002 -1000 99.00 0 4000 5000 5000 1 7000 8000 8000",
        "000000001 000000002 -1000 99.00 0 4000 5000 5000 0 0 1000 8000",
        "000000001 000000002 -1000 99.00 0 4000 4998 5000 0 3 1002 5000",
        ///
        "000000001 000000002 -1000 99.00 0 4000 4998 5000 0 3 1002 5000",
        "000000001 000000002 -1000 99.00 0 5 1000 5000 0 3999 4999 5000",
        "000000001 000000002 -1000 99.00 0 3997 5000 5000 1 3998 4999 5000",
        "000000001 000000002 -1000 99.00 0 2000 3000 5000 0 10 1010 1020",
        "000000001 000000002 -1000 99.00 0 10 1010 1020 0 2000 3000 5000",
    };

    // Expected results.
    // All have '*' types because they weren't specified in the input.
    std::vector<std::string> expected = {
        "000000001 000000002 -1000 99.00 0 0 1000 5000 0 7000 8000 8000 *",
        "000000001 000000002 -1000 99.00 0 0 1000 5000 1 0 1000 8000 *",
        "000000001 000000002 -1000 99.00 0 4000 5000 5000 1 7000 8000 8000 *",
        "000000001 000000002 -1000 99.00 0 4000 5000 5000 0 0 1000 8000 *",
        "000000001 000000002 -1000 99.00 0 3997 5000 5000 0 0 1004 5000 *",
        ///
        "000000001 000000002 -1000 99.00 0 3997 5000 5000 0 0 1004 5000 *",
        "000000001 000000002 -1000 99.00 0 0 1001 5000 0 3994 5000 5000 *",
        "000000001 000000002 -1000 99.00 0 3996 5000 5000 1 3998 5000 5000 *",
        "000000001 000000002 -1000 99.00 0 1990 3010 5000 0 0 1020 1020 *",
        "000000001 000000002 -1000 99.00 0 0 1020 1020 0 1990 3010 5000 *",
    };

    // Run and collect results.
    std::vector<std::string> results = WrapRunHeuristicExtendOverlapFlanks(inOverlaps, allowedDist);

    // Evaluate.
    EXPECT_EQ(expected, results);
}

TEST(ExtendOverlapFlanks, NormalWithAllowedDist0)
{
    /*
     * The output should be the same as input because the allowedDist is 0.
    */
    // Inputs.
    int32_t allowedDist = 0;
    std::vector<std::string> inOverlaps = {
        "000000001 000000002 -1000 99.00 0 4000 4998 5000 0 3 1002 5000",
        "000000001 000000002 -1000 99.00 0 5 1000 5000 0 3999 4999 5000",
        "000000001 000000002 -1000 99.00 0 3997 5000 5000 1 3998 4999 5000",
        "000000001 000000002 -1000 99.00 0 2000 3000 5000 0 10 1010 1020",
    };

    // Expected results.
    // All have '*' types because they weren't specified in the input.
    std::vector<std::string> expected = {
        "000000001 000000002 -1000 99.00 0 4000 4998 5000 0 3 1002 5000 *",
        "000000001 000000002 -1000 99.00 0 5 1000 5000 0 3999 4999 5000 *",
        "000000001 000000002 -1000 99.00 0 3997 5000 5000 1 3998 4999 5000 *",
        "000000001 000000002 -1000 99.00 0 2000 3000 5000 0 10 1010 1020 *",
    };

    // Run and collect results.
    std::vector<std::string> results = WrapRunHeuristicExtendOverlapFlanks(inOverlaps, allowedDist);

    // Evaluate.
    EXPECT_EQ(expected, results);
}

TEST(RoundTrip, ParsingOverlapType)
{
    /*
     * If parsing of the overlap types went fine, then the output should be the same
    */
    // Inputs.
    std::vector<std::string> inOverlaps = {
        "000000000 000000001 -1807 100.00 0 181 1988 1988 0 0 1807 1989 3",
        "000000001 000000002 -642 99.84 0 0 642 1989 0 1347 1989 1989 5",
        "000000001 000000002 -1000 100.00 0 0 1000 1000 0 1000 2000 3000 c",
        "000000001 000000002 -1000 100.00 0 0 1000 1000 0 1000 2000 3000 contained",
        "000000001 000000002 -1000 100.00 0 1000 2000 3000 0 0 1000 1000 contains",
        "000000001 000000002 -1000 100.00 0 1000 2000 3000 0 0 1000 1000 C",
        "000000001 000000002 -1000 100.00 0 1000 2000 3000 0 1000 2000 3000 u",
        "000000001 000000002 -1000 100.00 0 1000 2000 3000 0 1000 2000 3000 something_unknown",
    };

    // Expected results.
    // Since we're just printing as M4, the 'type' column uses the legacy long type name.
    std::vector<std::string> expected = {
        "000000000 000000001 -1807 100.00 0 181 1988 1988 0 0 1807 1989 3",
        "000000001 000000002 -642 99.84 0 0 642 1989 0 1347 1989 1989 5",
        "000000001 000000002 -1000 100.00 0 0 1000 1000 0 1000 2000 3000 contained",
        "000000001 000000002 -1000 100.00 0 0 1000 1000 0 1000 2000 3000 contained",
        "000000001 000000002 -1000 100.00 0 1000 2000 3000 0 0 1000 1000 contains",
        "000000001 000000002 -1000 100.00 0 1000 2000 3000 0 0 1000 1000 contains",
        "000000001 000000002 -1000 100.00 0 1000 2000 3000 0 1000 2000 3000 u",
        "000000001 000000002 -1000 100.00 0 1000 2000 3000 0 1000 2000 3000 *",
    };

    // for (size_t i = 0; i < inOverlaps
    std::vector<std::string> results;
    for (const auto& inLine : inOverlaps) {
        auto ovl = PacBio::Pancake::ParseM4OverlapFromString(inLine);
        std::string resultStr =
            PacBio::Pancake::OverlapWriterBase::PrintOverlapAsM4(*ovl, "", "", true, false);
        results.emplace_back(resultStr);
    }

    // Evaluate.
    EXPECT_EQ(expected, results);
}