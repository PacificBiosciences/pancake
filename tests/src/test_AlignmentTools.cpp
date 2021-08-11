#include <gtest/gtest.h>

#include <fstream>
#include <sstream>
#include <tuple>
#include <vector>

#include <pacbio/alignment/AlignmentTools.h>
#include <pbcopper/third-party/edlib.h>

#include "PancakeTestData.h"

using namespace PacBio::Pancake;

namespace AlignmentToolsTests {

TEST(Test_AlignmentTools_EdlibAlignmentToCigar, EmptyInput)
{
    std::vector<unsigned char> input = {};
    PacBio::BAM::Cigar expected;
    PacBio::Pancake::Alignment::DiffCounts expectedDiffs;

    PacBio::Pancake::Alignment::DiffCounts diffs;
    PacBio::BAM::Cigar result = EdlibAlignmentToCigar(input.data(), input.size(), diffs);
    EXPECT_EQ(expected, result);
    EXPECT_EQ(expectedDiffs, diffs);
}

TEST(Test_AlignmentTools_EdlibAlignmentToCigar, SingleOp)
{
    // Edlib move codes: 0: '=', 1: 'I', 2: 'D', 3: 'X'
    std::vector<unsigned char> input = {EDLIB_EDOP_MISMATCH};
    PacBio::BAM::Cigar expected("1X");
    PacBio::Pancake::Alignment::DiffCounts expectedDiffs(0, 1, 0, 0);

    PacBio::Pancake::Alignment::DiffCounts diffs;
    PacBio::BAM::Cigar result = EdlibAlignmentToCigar(input.data(), input.size(), diffs);
    EXPECT_EQ(expected, result);
    EXPECT_EQ(expectedDiffs, diffs);
}

TEST(Test_AlignmentTools_EdlibAlignmentToCigar, SimpleTest)
{
    // Edlib move codes: 0: '=', 1: 'I', 2: 'D', 3: 'X'
    std::vector<unsigned char> input = {EDLIB_EDOP_MATCH,    EDLIB_EDOP_MATCH,  EDLIB_EDOP_MATCH,
                                        EDLIB_EDOP_MISMATCH, EDLIB_EDOP_INSERT, EDLIB_EDOP_DELETE,
                                        EDLIB_EDOP_DELETE,   EDLIB_EDOP_INSERT};
    PacBio::BAM::Cigar expected("3=1X1I2D1I");
    PacBio::Pancake::Alignment::DiffCounts expectedDiffs(3, 1, 2, 2);

    PacBio::Pancake::Alignment::DiffCounts diffs;
    PacBio::BAM::Cigar result = EdlibAlignmentToCigar(input.data(), input.size(), diffs);
    EXPECT_EQ(expected, result);
    EXPECT_EQ(expectedDiffs, diffs);
}

TEST(Test_AlignmentTools_ValidateCigar, ArrayOfTests)
{
    // clang-format off
    std::vector<std::tuple<std::string, std::string, std::string, std::string, bool>> testData = {
        // Good cases.
        {"Empty input", "", "", "", false},
        {"Exact match, good", "ACTG", "ACTG", "4=", false},
        {"Cigar with indels", "ACTT", "ACTG", "3=1I1D", false},
        {"Soft clipping", "AAAACTTAAA", "ACTG", "3S3=1X3S", false},
        {"Hard clipping", "ACTT", "ACTG", "3H3=1X3H", false},
        {"Simple test case", "GGATCAGTTTTATATACAC", "GAGTTCGTTCTATATATACAC", "1=1I1=1D1=1D1=1I3=1X3=2D4=2I2D", false},

        // Bad cases.
        {"One base is incorrect", "ACTT", "ACTG", "4=", true},
        {"Length mismatch", "ACT", "ACTG", "4=", true},
        {"Non-zero CIGAR, but empty seqs", "", "", "4=", true},
        {"Length mismatch with indels", "ACTT", "ACTG", "3=1I1D1I", true},
        {"Length mismatch with mismatches", "ACTT", "ACTG", "3=2X", true},
        {"Simple test case", "GGATCAGTTTTATATACAC", "GAGTTCGTTCTATATATACAC", "1=1I1=1D1=1D1=1I3=1X3=2D4=2I", true},

    };
    // clang-format on

    for (const auto& data : testData) {
        // Inputs.
        const std::string testName = std::get<0>(data);
        const std::string& query = std::get<1>(data);
        const std::string& target = std::get<2>(data);
        PacBio::BAM::Cigar cigar(std::get<3>(data));
        bool shouldThrow = std::get<4>(data);

        // Name the test.
        SCOPED_TRACE("ValidateCigar-" + testName);

        // Run.
        if (shouldThrow) {
            EXPECT_THROW(
                {
                    PacBio::Pancake::ValidateCigar(query.c_str(), query.size(), target.c_str(),
                                                   target.size(), cigar, "");
                },
                std::runtime_error);
        } else {
            PacBio::Pancake::ValidateCigar(query.c_str(), query.size(), target.c_str(),
                                           target.size(), cigar, "");
        }
    }
}

class ExtractVariantStringTestCase
{
public:
    std::string testName;
    std::string query;
    std::string target;
    std::string cigar;
    bool maskHomopolymers;
    bool maskSimpleRepeats;
    bool maskHomopolymerSNPs;
    bool maskHomopolymersArbitrary;
    std::string expectedQueryVariants;
    std::string expectedTargetVariants;
    PacBio::Pancake::Alignment::DiffCounts expectedDiffsPerBase;
    PacBio::Pancake::Alignment::DiffCounts expectedDiffsPerEvent;
};

TEST(Test_AlignmentTools_ExtractVariantString, ArrayOfTests)
{
    // clang-format off
    std::vector<ExtractVariantStringTestCase> testData = {
        ExtractVariantStringTestCase{"Empty input", "", "", "", false, false, false, false, "", "", {}, {}},

        /*
        * Q: GGA-T-CAGTTTT AT--ATACAC
        *    | | | | |||X| ||  ||||
        * T: G-AGTTC-GTTCT ATATATAC--
        * Query: indel G is in a homopolymer, A is not in a HP.
        * Target: indel G is not a HP, and T is.
        * There is one insertion event of 2bp and one deletion event of 2bp. Both are in simple repeats.
        */
        ExtractVariantStringTestCase{"Simple test case",
                                     "GGATCAGTTTTATATACAC", "GAGTTCGTTCTATATATACAC", "1=1I1=1D1=1D1=1I3=1X3=2D4=2I", false, false, false, false, "GATAC", "GTCAT", {14, 1, 4, 4}, {14, 1, 3, 3}},

        ExtractVariantStringTestCase{"Test case with HP masking",
                                     "GGATCAGTTTTATATACAC", "GAGTTCGTTCTATATATACAC", "1=1I1=1D1=1D1=1I3=1X3=2D4=2I", true, false, false, false, "gATAC", "GtCAT", {14, 1, 3, 3}, {14, 1, 2, 2}},

        ExtractVariantStringTestCase{"Test case with simple repeat masking",
                                     "GGATCAGTTTTATATACAC", "GAGTTCGTTCTATATATACAC", "1=1I1=1D1=1D1=1I3=1X3=2D4=2I", false, true, false, false, "GATac", "GTCat", {14, 1, 2, 2}, {14, 1, 2, 2}},

        ExtractVariantStringTestCase{"Test case with HP and simple repeat masking",
                                     "GGATCAGTTTTATATACAC", "GAGTTCGTTCTATATATACAC", "1=1I1=1D1=1D1=1I3=1X3=2D4=2I", true, true, false, false, "gATac", "GtCat", {14, 1, 1, 1}, {14, 1, 1, 1}},

        /*
         * This test case ("GCAC', "GAC") detected a false masking bug in the code. The HP masking used
         * to check forward in target for insertions and backward in target for deletions. This is wrong
         * because the coordinate in the strand with the gap points to the base _after_ the gap (and not before
         * like I originally thought). So it should be vice versa.
        */
        ExtractVariantStringTestCase{"Small HP masking test 1 - insertion, should not masking",
                                     "GCAC", "GAC", "1=1I2=", true, true, false, false, "C", "", {3, 0, 1, 0}, {3, 0, 1, 0}},
        ExtractVariantStringTestCase{"Small HP masking test 1 - deletion, should not masking",
                                     "GAC", "GCAC", "1=1D2=", true, true, false, false, "", "C", {3, 0, 0, 1}, {3, 0, 0, 1}},
        /*
         * This is extracted around the same variant (from the dataset), but in the other direction.
         * There were no issues in this direction.
        */
        ExtractVariantStringTestCase{"Small HP masking test 2 - insertion, should not masking",
                                     "TGCC", "TCC", "1=1I2=", true, true, false, false, "G", "", {3, 0, 1, 0}, {3, 0, 1, 0}},
        ExtractVariantStringTestCase{"Small HP masking test 2 - deletion, should not masking",
                                     "TCC", "TGCC", "1=1D2=", true, true, false, false, "", "G", {3, 0, 0, 1}, {3, 0, 0, 1}},

        /*
         * Test HP insertion masking systematically.
        */
        ExtractVariantStringTestCase{"Simple HP masking test 3a - insertion, no masking",
                                     "AAACAAA", "AAAAAA", "3=1I3=", true, true, false, false, "C", "", {6, 0, 1, 0}, {6, 0, 1, 0}},
        ExtractVariantStringTestCase{"Simple HP masking test 3b - insertion, should mask",
                                     "AAACCAA", "AAACAA", "3=1I3=", true, true, false, false, "c", "", {6, 0, 0, 0}, {6, 0, 0, 0}},
        ExtractVariantStringTestCase{"Simple HP masking test 3c - insertion, no masking",
                                     "AAACACA", "AAAACA", "3=1I3=", true, true, false, false, "C", "", {6, 0, 1, 0}, {6, 0, 1, 0}},
        ExtractVariantStringTestCase{"Simple HP masking test 3d - insertion, no masking",
                                     "AAACAAC", "AAAAAC", "3=1I3=", true, true, false, false, "C", "", {6, 0, 1, 0}, {6, 0, 1, 0}},
        ExtractVariantStringTestCase{"Simple HP masking test 3e - insertion, no masking",
                                     "AACCAAA", "AACAAA", "3=1I3=", true, true, false, false, "c", "", {6, 0, 0, 0}, {6, 0, 0, 0}},
        ExtractVariantStringTestCase{"Simple HP masking test 3f - insertion, no masking",
                                     "ACACAAA", "ACAAAA", "3=1I3=", true, true, false, false, "C", "", {6, 0, 1, 0}, {6, 0, 1, 0}},
        ExtractVariantStringTestCase{"Simple HP masking test 3g - insertion, no masking",
                                     "CAACAAA", "CAAAAA", "3=1I3=", true, true, false, false, "C", "", {6, 0, 1, 0}, {6, 0, 1, 0}},
        ExtractVariantStringTestCase{"Simple HP masking test 3h - no indels, perfect match",
                                     "AAACAAA", "AAACAAA", "7=", true, true, false, false, "", "", {7, 0, 0, 0}, {7, 0, 0, 0}},
        // Test masking of arbitrary indels in HPs.
        ExtractVariantStringTestCase{"Simple HP masking test 3i - insertion, should mask because we now allow arbitrary indels in HP events.",
                                     "AAACAAA", "AAAAAA", "3=1I3=", true, true, false, true, "c", "", {6, 0, 0, 0}, {6, 0, 0, 0}},
        ExtractVariantStringTestCase{"Simple HP masking test 3j - insertion, should mask",
                                     "AAACCAA", "AAACAA", "3=1I3=", true, true, false, true, "c", "", {6, 0, 0, 0}, {6, 0, 0, 0}},
        ExtractVariantStringTestCase{"Simple HP masking test 3k - insertion, should mask because we now allow arbitrary indels in HP events.",
                                     "AAACACA", "AAAACA", "3=1I3=", true, true, false, true, "c", "", {6, 0, 0, 0}, {6, 0, 0, 0}},
        ExtractVariantStringTestCase{"Simple HP masking test 3l - insertion, should mask because we now allow arbitrary indels in HP events.",
                                     "AAACAAC", "AAAAAC", "3=1I3=", true, true, false, true, "c", "", {6, 0, 0, 0}, {6, 0, 0, 0}},
        ExtractVariantStringTestCase{"Simple HP masking test 3m - insertion, should mask",
                                     "AACCAAA", "AACAAA", "3=1I3=", true, true, false, true, "c", "", {6, 0, 0, 0}, {6, 0, 0, 0}},
        ExtractVariantStringTestCase{"Simple HP masking test 3n - insertion, should mask because we now allow arbitrary indels in HP events.",
                                     "ACACAAA", "ACAAAA", "3=1I3=", true, true, false, true, "c", "", {6, 0, 0, 0}, {6, 0, 0, 0}},
        ExtractVariantStringTestCase{"Simple HP masking test 3o - insertion, should mask because we now allow arbitrary indels in HP events.",
                                     "CAACAAA", "CAAAAA", "3=1I3=", true, true, false, true, "c", "", {6, 0, 0, 0}, {6, 0, 0, 0}},
        ExtractVariantStringTestCase{"Simple HP masking test 3p - no indels, perfect match",
                                     "AAACAAA", "AAACAAA", "7=", true, true, false, true, "", "", {7, 0, 0, 0}, {7, 0, 0, 0}},

        /*
         * Test HP deletion masking systematically.
        */
        ExtractVariantStringTestCase{"Simple HP masking test 4a - deletion, no masking",
                                     "AAAAAA", "AAACAAA", "3=1D3=", true, true, false, false, "", "C", {6, 0, 0, 1}, {6, 0, 0, 1}},
        ExtractVariantStringTestCase{"Simple HP masking test 4b - deletion, should mask",
                                     "AAACAA", "AAACCAA", "3=1D3=", true, true, false, false, "", "c", {6, 0, 0, 0}, {6, 0, 0, 0}},
        ExtractVariantStringTestCase{"Simple HP masking test 4c - deletion, no masking",
                                     "AAAACA", "AAACACA", "3=1D3=", true, true, false, false, "", "C", {6, 0, 0, 1}, {6, 0, 0, 1}},
        ExtractVariantStringTestCase{"Simple HP masking test 4d - deletion, no masking",
                                     "AAAAAC", "AAACAAC", "3=1D3=", true, true, false, false, "", "C", {6, 0, 0, 1}, {6, 0, 0, 1}},
        ExtractVariantStringTestCase{"Simple HP masking test 4e - deletion, no masking",
                                     "AACAAA", "AACCAAA", "3=1D3=", true, true, false, false, "", "c", {6, 0, 0, 0}, {6, 0, 0, 0}},
        ExtractVariantStringTestCase{"Simple HP masking test 4f - deletion, no masking",
                                     "ACAAAA", "ACACAAA", "3=1D3=", true, true, false, false, "", "C", {6, 0, 0, 1}, {6, 0, 0, 1}},
        ExtractVariantStringTestCase{"Simple HP masking test 4g - deletion, no masking",
                                     "CAAAAA", "CAACAAA", "3=1D3=", true, true, false, false, "", "C", {6, 0, 0, 1}, {6, 0, 0, 1}},
        ExtractVariantStringTestCase{"Simple HP masking test 4h - no indels, perfect match",
                                     "AAACAAA", "AAACAAA", "7=", true, true, false, false, "", "", {7, 0, 0, 0}, {7, 0, 0, 0}},
        ExtractVariantStringTestCase{"Simple HP masking test 4i - deletion, should mask because we now allow arbitrary indels in HP events.",
                                     "AAAAAA", "AAACAAA", "3=1D3=", true, true, false, true, "", "c", {6, 0, 0, 0}, {6, 0, 0, 0}},
        ExtractVariantStringTestCase{"Simple HP masking test 4j - deletion, should mask",
                                     "AAACAA", "AAACCAA", "3=1D3=", true, true, false, true, "", "c", {6, 0, 0, 0}, {6, 0, 0, 0}},
        ExtractVariantStringTestCase{"Simple HP masking test 4k - deletion, should mask because we now allow arbitrary indels in HP events.",
                                     "AAAACA", "AAACACA", "3=1D3=", true, true, false, true, "", "c", {6, 0, 0, 0}, {6, 0, 0, 0}},
        ExtractVariantStringTestCase{"Simple HP masking test 4l - deletion, should mask because we now allow arbitrary indels in HP events.",
                                     "AAAAAC", "AAACAAC", "3=1D3=", true, true, false, true, "", "c", {6, 0, 0, 0}, {6, 0, 0, 0}},
        ExtractVariantStringTestCase{"Simple HP masking test 4m - deletion, no masking",
                                     "AACAAA", "AACCAAA", "3=1D3=", true, true, false, true, "", "c", {6, 0, 0, 0}, {6, 0, 0, 0}},
        ExtractVariantStringTestCase{"Simple HP masking test 4n - deletion, should mask because we now allow arbitrary indels in HP events.",
                                     "ACAAAA", "ACACAAA", "3=1D3=", true, true, false, true, "", "c", {6, 0, 0, 0}, {6, 0, 0, 0}},
        ExtractVariantStringTestCase{"Simple HP masking test 4o - deletion, should mask because we now allow arbitrary indels in HP events.",
                                     "CAAAAA", "CAACAAA", "3=1D3=", true, true, false, true, "", "c", {6, 0, 0, 0}, {6, 0, 0, 0}},
        ExtractVariantStringTestCase{"Simple HP masking test 4p - no indels, perfect match",
                                     "AAACAAA", "AAACAAA", "7=", true, true, false, true, "", "", {7, 0, 0, 0}, {7, 0, 0, 0}},

        /*
         * Test simple releat insertion masking systematically.
        */
        ExtractVariantStringTestCase{"Simple repeat masking test 5a - no indel, perfect match",
                                     "AAAAAAATCAAAAAAA", "AAAAAAATCAAAAAAA", "16=", true, true, false, false, "", "", {16, 0, 0, 0}, {16, 0, 0, 0}},
        ExtractVariantStringTestCase{"Simple repeat masking test 5b - insertion, no masking",
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "AAAAAAAAAAAAAA",      // AAAAAAA--AAAAAAA
                                     "7=2I7=", true, true, false, false, "TC", "", {14, 0, 2, 0}, {14, 0, 1, 0}},
        /*
         * Sliding window for insertion masking - slide a TC event accross the target, while the query has a TC insertion.
         * Only in some cases should masking pick up.
        */
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 5a - insertion, should not mask",
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "ATCAAAAAAAAAAA",      // ATCAAAA--AAAAAAA
                                     "1=2X4=2I7=", true, true, false, false, "AATC", "TC", {12, 2, 2, 0}, {12, 2, 1, 0}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 5b - insertion, should not mask",
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "AATCAAAAAAAAAA",      // AATCAAA--AAAAAAA
                                     "2=2X3=2I7=", true, true, false, false, "AATC", "TC", {12, 2, 2, 0}, {12, 2, 1, 0}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 5c - insertion, should not mask",
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "AAATCAAAAAAAAA",      // AAATCAA--AAAAAAA
                                     "3=2X2=2I7=", true, true, false, false, "AATC", "TC", {12, 2, 2, 0}, {12, 2, 1, 0}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 5d - insertion, should not mask",
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "AAAATCAAAAAAAA",      // AAAATCA--AAAAAAA
                                     "4=2X1=2I7=", true, true, false, false, "AATC", "TC", {12, 2, 2, 0}, {12, 2, 1, 0}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 5e - insertion, should mask",
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "AAAAATCAAAAAAA",      // AAAAATC--AAAAAAA
                                     "5=2X2I7=", true, true, false, false, "AAtc", "TC", {12, 2, 0, 0}, {12, 2, 0, 0}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 5f - insertion, should not mask",
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "AAAAAATCAAAAAA",      // AAAAAAT--CAAAAAA
                                     "6=1X2I1X6=", true, true, false, false, "ATCA", "TC", {12, 2, 2, 0}, {12, 2, 1, 0}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 5g - insertion, should mask",
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "AAAAAAATCAAAAA",      // AAAAAAA--TCAAAAA
                                     "7=2I2X5=", true, true, false, false, "tcAA", "TC", {12, 2, 0, 0}, {12, 2, 0, 0}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 5h - insertion, should not mask",
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "AAAAAAAATCAAAA",      // AAAAAAA--ATCAAAA
                                     "7=2I1=2X4=", true, true, false, false, "TCAA", "TC", {12, 2, 2, 0}, {12, 2, 1, 0}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 5i - insertion, should not mask",
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "AAAAAAAAATCAAA",      // AAAAAAA--AATCAAA
                                     "7=2I2=2X3=", true, true, false, false, "TCAA", "TC", {12, 2, 2, 0}, {12, 2, 1, 0}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 5j - insertion, should not mask",
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "AAAAAAAAAATCAA",      // AAAAAAA--AAATCAA
                                     "7=2I3=2X2=", true, true, false, false, "TCAA", "TC", {12, 2, 2, 0}, {12, 2, 1, 0}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 5k - insertion, should not mask",
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "AAAAAAAAAAATCA",      // AAAAAAA--AAAATCA
                                     "7=2I4=2X1=", true, true, false, false, "TCAA", "TC", {12, 2, 2, 0}, {12, 2, 1, 0}},
        /*
         * Sliding window for deletion masking - slide a TC event accross the query, while the target has a TC insertion.
         * Only in some cases should masking pick up.
        */
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 6a - insertion, should not mask",
                                     "ATCAAAAAAAAAAA",      // ATCAAAA--AAAAAAA
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "1=2X4=2D7=", true, true, false, false, "TC", "AATC", {12, 2, 0, 2}, {12, 2, 0, 1}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 6b - deletion, should not mask",
                                     "AATCAAAAAAAAAA",      // AATCAAA--AAAAAAA
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "2=2X3=2D7=", true, true, false, false, "TC", "AATC", {12, 2, 0, 2}, {12, 2, 0, 1}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 6c - deletion, should not mask",
                                     "AAATCAAAAAAAAA",      // AAATCAA--AAAAAAA
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "3=2X2=2D7=", true, true, false, false, "TC", "AATC", {12, 2, 0, 2}, {12, 2, 0, 1}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 6d - deletion, should not mask",
                                     "AAAATCAAAAAAAA",      // AAAATCA--AAAAAAA
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "4=2X1=2D7=", true, true, false, false, "TC", "AATC", {12, 2, 0, 2}, {12, 2, 0, 1}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 6e - deletion, should mask",
                                     "AAAAATCAAAAAAA",      // AAAAATC--AAAAAAA
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "5=2X2D7=", true, true, false, false, "TC", "AAtc", {12, 2, 0, 0}, {12, 2, 0, 0}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 6f - deletion, should not mask",
                                     "AAAAAATCAAAAAA",      // AAAAAAT--CAAAAAA
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "6=1X2D1X6=", true, true, false, false, "TC", "ATCA", {12, 2, 0, 2}, {12, 2, 0, 1}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 6g - deletion, should mask",
                                     "AAAAAAATCAAAAA",      // AAAAAAA--TCAAAAA
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "7=2D2X5=", true, true, false, false, "TC", "tcAA", {12, 2, 0, 0}, {12, 2, 0, 0}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 6h - deletion, should not mask",
                                     "AAAAAAAATCAAAA",      // AAAAAAA--ATCAAAA
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "7=2D1=2X4=", true, true, false, false, "TC", "TCAA", {12, 2, 0, 2}, {12, 2, 0, 1}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 6i - deletion, should not mask",
                                     "AAAAAAAAATCAAA",      // AAAAAAA--AATCAAA
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "7=2D2=2X3=", true, true, false, false, "TC", "TCAA", {12, 2, 0, 2}, {12, 2, 0, 1}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 6j - deletion, should not mask",
                                     "AAAAAAAAAATCAA",      // AAAAAAA--AAATCAA
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "7=2D3=2X2=", true, true, false, false, "TC", "TCAA", {12, 2, 0, 2}, {12, 2, 0, 1}},
        ExtractVariantStringTestCase{"Simple repeat masking sliding window test 6k - deletion, should not mask",
                                     "AAAAAAAAAAATCA",      // AAAAAAA--AAAATCA
                                     "AAAAAAATCAAAAAAA",    // AAAAAAATCAAAAAAA
                                     "7=2D4=2X1=", true, true, false, false, "TC", "TCAA", {12, 2, 0, 2}, {12, 2, 0, 1}},
        /*
         * Test simple repeat masking in the same sequence (query).
        */
        ExtractVariantStringTestCase{"Simple repeat masking in same sequence window test 7a - should not mask",
                                     "AAAAAAATCTAAAAAA",    // AAAAAAATCTAAAAAA
                                     "AAAAAAAAAAAAAA",      // AAAAAAA--AAAAAAA
                                     "7=2I1X6=", true, true, false, false, "TCT", "A", {13, 1, 2, 0}, {13, 1, 1, 0}},
        ExtractVariantStringTestCase{"Simple repeat masking in same sequence window test 7b - should mask",
                                     "AAAAAAATCTCAAAAA",    // AAAAAAATCTCAAAAA
                                     "AAAAAAAAAAAAAA",      // AAAAAAA--AAAAAAA
                                     "7=2I2X5=", true, true, false, false, "tcTC", "AA", {12, 2, 0, 0}, {12, 2, 0, 0}},
        ExtractVariantStringTestCase{"Simple repeat masking in same sequence window test 7c - should not mask",
                                     "AAAAAATTCAAAAAAA",    // AAAAAATTCAAAAAA
                                     "AAAAAAAAAAAAAA",      // AAAAAAA--AAAAAAA
                                     "6=1X2I7=", true, true, false, false, "TTC", "A", {13, 1, 2, 0}, {13, 1, 1, 0}},
        ExtractVariantStringTestCase{"Simple repeat masking in same sequence window test 7d - should mask",
                                     "AAAAATCTCAAAAAAA",    // AAAAAAATCTCAAAAA
                                     "AAAAAAAAAAAAAA",      // AAAAAAA--AAAAAAA
                                     "5=2X2I7=", true, true, false, false, "TCtc", "AA", {12, 2, 0, 0}, {12, 2, 0, 0}},
        /*
         * Test simple repeat masking in the same sequence (target).
        */
        ExtractVariantStringTestCase{"Simple repeat masking in same sequence window test 8a - should not mask",
                                     "AAAAAAAAAAAAAA",      // AAAAAAA--AAAAAAA
                                     "AAAAAAATCTAAAAAA",    // AAAAAAATCTAAAAAA
                                     "7=2D1X6=", true, true, false, false, "A", "TCT", {13, 1, 0, 2}, {13, 1, 0, 1}},
        ExtractVariantStringTestCase{"Simple repeat masking in same sequence window test 8b - should mask",
                                     "AAAAAAAAAAAAAA",      // AAAAAAA--AAAAAAA
                                     "AAAAAAATCTCAAAAA",    // AAAAAAATCTCAAAAA
                                     "7=2D2X5=", true, true, false, false, "AA", "tcTC", {12, 2, 0, 0}, {12, 2, 0, 0}},
        ExtractVariantStringTestCase{"Simple repeat masking in same sequence window test 8c - should not mask",
                                     "AAAAAAAAAAAAAA",      // AAAAAAA--AAAAAAA
                                     "AAAAAATTCAAAAAAA",    // AAAAAATTCAAAAAA
                                     "6=1X2D7=", true, true, false, false, "A", "TTC", {13, 1, 0, 2}, {13, 1, 0, 1}},
        ExtractVariantStringTestCase{"Simple repeat masking in same sequence window test 8d - should mask",
                                     "AAAAAAAAAAAAAA",      // AAAAAAA--AAAAAAA
                                     "AAAAATCTCAAAAAAA",    // AAAAAAATCTCAAAAA
                                     "5=2X2D7=", true, true, false, false, "AA", "TCtc", {12, 2, 0, 0}, {12, 2, 0, 0}},

        /*
         * Test masking of SNPs in HP regions. The SNP is in the target.
        */
        ExtractVariantStringTestCase{"SNP internal to the homopolymer. Should be masked. (target)",
                                     "AAAAAAAAAAAAAAA",
                                     "AAAAAAATAAAAAAA",
                                     "7=1X7=", false, false, true, false, "a", "t", {14, 0, 0, 0}, {14, 0, 0, 0}},
        ExtractVariantStringTestCase{"SNP at the right boundary of the homopolymer. Should be masked. (target)",
                                     "AAAAAAAACCCCCCC",
                                     "AAAAAAATCCCCCCC",
                                     "7=1X7=", false, false, true, false, "a", "t", {14, 0, 0, 0}, {14, 0, 0, 0}},
        ExtractVariantStringTestCase{"SNP at the left boundary of the homopolymer. Should be masked. (target)",
                                     "CCCCCCCAAAAAAAA",
                                     "CCCCCCCTAAAAAAA",
                                     "7=1X7=", false, false, true, false, "a", "t", {14, 0, 0, 0}, {14, 0, 0, 0}},
        ExtractVariantStringTestCase{"SNP adjacent to the homopolymer, but not a homopolymer. (target)",
                                     "AAAAAAACGCGCGCG",
                                     "AAAAAAATGCGCGCG",
                                     "7=1X7=", false, false, true, false, "C", "T", {14, 1, 0, 0}, {14, 1, 0, 0}},
        ExtractVariantStringTestCase{"Cross-HP SNP. Should be masked. (target)",
                                     "AAAAAAAACCCCCCC",
                                     "AAAAAAACCCCCCCC",
                                     "7=1X7=", false, false, true, false, "a", "c", {14, 0, 0, 0}, {14, 0, 0, 0}},
        ExtractVariantStringTestCase{"SNP next to an insertion. (target)",
                                     "AAAAAAACCCCCCC",  // AAAAAA-ACCCCCCC
                                     "AAAAAAATCCCCCCC", // AAAAAAATCCCCCCC
                                     "6=1D1X7=", false, false, true, false, "a", "At", {13, 0, 0, 1}, {13, 0, 0, 1}},
        ExtractVariantStringTestCase{"SNP next to a deletion. (target)",
                                     "AAAAAAAACCCCCCC", // AAAAAAAACCCCCCC
                                     "AAAAAATCCCCCCC",  // AAAAAA-TCCCCCCC
                                     "6=1I1X7=", false, false, true, false, "Aa", "t", {13, 0, 1, 0}, {13, 0, 1, 0}},
        ExtractVariantStringTestCase{"SNP in a very short HP. (target)",
                                     "ACTTG", // AAAAAAAACCCCCCC
                                     "ACCTG",  // AAAAAA-TCCCCCCC
                                     "2=1X2=", false, false, true, false, "t", "c", {4, 0, 0, 0}, {4, 0, 0, 0}},

        /*
         * Test masking of SNPs in HP regions. The SNP is in the query. Same as above, just swapped query and target values.
        */
        ExtractVariantStringTestCase{"SNP internal to the homopolymer. Should be masked. (query)",
                                     "AAAAAAATAAAAAAA",
                                     "AAAAAAAAAAAAAAA",
                                     "7=1X7=", false, false, true, false, "t", "a", {14, 0, 0, 0}, {14, 0, 0, 0}},
        ExtractVariantStringTestCase{"SNP at the right boundary of the homopolymer. Should be masked. (query)",
                                     "AAAAAAATCCCCCCC",
                                     "AAAAAAAACCCCCCC",
                                     "7=1X7=", false, false, true, false, "t", "a", {14, 0, 0, 0}, {14, 0, 0, 0}},
        ExtractVariantStringTestCase{"SNP at the left boundary of the homopolymer. Should be masked. (query)",
                                     "CCCCCCCTAAAAAAA",
                                     "CCCCCCCAAAAAAAA",
                                     "7=1X7=", false, false, true, false, "t", "a", {14, 0, 0, 0}, {14, 0, 0, 0}},
        ExtractVariantStringTestCase{"SNP adjacent to the homopolymer, but not a homopolymer. (query)",
                                     "AAAAAAATGCGCGCG",
                                     "AAAAAAACGCGCGCG",
                                     "7=1X7=", false, false, true, false, "T", "C", {14, 1, 0, 0}, {14, 1, 0, 0}},
        ExtractVariantStringTestCase{"Cross-HP SNP. Should be masked. (query)",
                                     "AAAAAAACCCCCCCC",
                                     "AAAAAAAACCCCCCC",
                                     "7=1X7=", false, false, true, false, "c", "a", {14, 0, 0, 0}, {14, 0, 0, 0}},
        ExtractVariantStringTestCase{"SNP next to an insertion. (query)",
                                     "AAAAAAATCCCCCCC", // AAAAAAATCCCCCCC
                                     "AAAAAAACCCCCCC",  // AAAAAA-ACCCCCCC
                                     "6=1I1X7=", false, false, true, false, "At", "a", {13, 0, 1, 0}, {13, 0, 1, 0}},
        ExtractVariantStringTestCase{"SNP next to a deletion. (query)",
                                     "AAAAAATCCCCCCC",  // AAAAAA-TCCCCCCC
                                     "AAAAAAAACCCCCCC", // AAAAAAAACCCCCCC
                                     "6=1D1X7=", false, false, true, false, "t", "Aa", {13, 0, 0, 1}, {13, 0, 0, 1}},
        ExtractVariantStringTestCase{"SNP in a very short HP. (query)",
                                     "ACCTG",  // AAAAAA-TCCCCCCC
                                     "ACTTG", // AAAAAAAACCCCCCC
                                     "2=1X2=", false, false, true, false, "c", "t", {4, 0, 0, 0}, {4, 0, 0, 0}},
    };
    // clang-format on

    for (const auto& data : testData) {
        // Name the test.
        SCOPED_TRACE("ExtractVariantString-" + data.testName);
        PacBio::BAM::Cigar cigar(data.cigar);
        std::string resultQueryVariants;
        std::string resultTargetVariants;
        PacBio::Pancake::Alignment::DiffCounts resultDiffsPerBase;
        PacBio::Pancake::Alignment::DiffCounts resultDiffsPerEvent;

        // Run.
        PacBio::Pancake::ExtractVariantString(
            data.query.c_str(), data.query.size(), data.target.c_str(), data.target.size(), cigar,
            data.maskHomopolymers, data.maskSimpleRepeats, data.maskHomopolymerSNPs,
            data.maskHomopolymersArbitrary, resultQueryVariants, resultTargetVariants,
            resultDiffsPerBase, resultDiffsPerEvent);

        // Evaluate.
        EXPECT_EQ(data.expectedQueryVariants, resultQueryVariants);
        EXPECT_EQ(data.expectedTargetVariants, resultTargetVariants);
        EXPECT_EQ(data.expectedDiffsPerBase, resultDiffsPerBase);
        EXPECT_EQ(data.expectedDiffsPerEvent, resultDiffsPerEvent);
    }
}

TEST(Test_AlignmentTools_FindTargetPosFromCigar, ArrayOfTests)
{
    // clang-format off
    std::vector<std::tuple<std::string, std::string, int32_t, int32_t, bool>> testData = {
        {"Empty input", "", 0, -1, true},
        {"Simple matches", "10=", 4, 4, false},
        {"Simple matches", "10=", -1, -1, true},
        {"Simple matches", "10=", 10, -1, true},

        {"Simple CIGAR pos 0", "2=2D1=1X1=2I1=", 0, 0, false},
        {"Simple CIGAR pos 1", "2=2D1=1X1=2I1=", 1, 1, false},
        {"Simple CIGAR pos 2", "2=2D1=1X1=2I1=", 2, 4, false},
        {"Simple CIGAR pos 3", "2=2D1=1X1=2I1=", 3, 5, false},
        {"Simple CIGAR pos 4", "2=2D1=1X1=2I1=", 4, 6, false},
        {"Simple CIGAR pos 5", "2=2D1=1X1=2I1=", 5, 6, false},
        {"Simple CIGAR pos 6", "2=2D1=1X1=2I1=", 6, 6, false},
        {"Simple CIGAR pos 7", "2=2D1=1X1=2I1=", 7, 7, false},

        {"Custom real 1", "59=1I10=1I601=1I48=1D274=1D432=1I84=1I573=1I94=1D1545=1I1=1D1131=1D1042=1I329=1I581=1D452=", 1380, 1379, false},
    };
    // clang-format on

    for (const auto& data : testData) {
        // Get the data.
        const std::string testName = std::get<0>(data);
        PacBio::BAM::Cigar cigar(std::get<1>(data));
        int32_t queryPos = std::get<2>(data);
        int32_t expected = std::get<3>(data);
        bool shouldThrow = std::get<4>(data);

        // Name the test.
        SCOPED_TRACE("FindTargetPosFromCigar-" + testName);

        // Run.
        if (shouldThrow) {
            EXPECT_THROW({ PacBio::Pancake::FindTargetPosFromCigar(cigar, queryPos); },
                         std::runtime_error);
        } else {
            int32_t result = PacBio::Pancake::FindTargetPosFromCigar(cigar, queryPos);
            // Evaluate.
            EXPECT_EQ(expected, result);
        }
    }
}

TEST(Test_AlignmentTools_ComputeDiffCounts, ArrayOfTests)
{
    // clang-format off
    std::vector<std::tuple<std::string, std::string, std::string, std::string, Alignment::DiffCounts, bool, bool>> testData = {
        {"Empty input", "", "", "", Alignment::DiffCounts(), false, false},
        {"Exact match", "4=", "", "", Alignment::DiffCounts(4, 0, 0, 0), false, false},
        {"Variant pos out of bounds (X)", "3=1X1=", "", "", Alignment::DiffCounts(), false, true},
        {"Variant pos out of bounds (I)", "3=1I1=", "", "", Alignment::DiffCounts(), false, true},
        {"Variant pos out of bounds (D)", "3=1D1=", "", "", Alignment::DiffCounts(), false, true},
        {"Valid case, 1X", "3=1X1=", "A", "C", Alignment::DiffCounts(4, 1, 0, 0), false, false},
        {"Valid case, 1X, masked", "3=1X1=", "a", "c", Alignment::DiffCounts(4, 0, 0, 0), false, false},
        {"Valid case, 1X, mixed mask, should throw", "3=1X1=", "a", "C", Alignment::DiffCounts(), false, true},
        {"Valid case with multiple diffs and 1 masked I", "3=1X1=1D2=2I2=1I", "ATAt", "CT", Alignment::DiffCounts(8, 1, 2, 1), false, false},
        {"Valid case with clipping and reference skip", "3S3=1X1=1D2=1I2N2=1I3H", "ATt", "CT", Alignment::DiffCounts(8, 1, 1, 1), false, false},
        {"The 2-base I has mixed masking, but the parameter says not to throw", "3=1X1=1D2=2I2=1I", "ATat", "CT", Alignment::DiffCounts(8, 1, 1, 1), false, false},
        {"The 2-base D has mixed masking, but the parameter says not to throw", "3=1X1=2D2=2I2=1I", "ATat", "CTc", Alignment::DiffCounts(8, 1, 1, 1), false, false},
        {"Bad case, the 2-base I has mixed masking", "3=1X1=1D2=2I2=1I", "ATat", "CT", Alignment::DiffCounts(), true, true},
        {"Bad case, the 2-base B has mixed masking", "3=1X1=2D2=2I2=1I", "ATat", "CTc", Alignment::DiffCounts(), true, true},
    };
    // clang-format on

    for (const auto& data : testData) {
        // Inputs.
        const std::string testName = std::get<0>(data);
        PacBio::BAM::Cigar cigar(std::get<1>(data));
        const std::string& queryVariants = std::get<2>(data);
        const std::string& targetVariants = std::get<3>(data);
        const PacBio::Pancake::Alignment::DiffCounts& expected = std::get<4>(data);
        bool throwOnPartiallyMaskedIndels = std::get<5>(data);
        bool shouldThrow = std::get<6>(data);

        // Name the test.
        SCOPED_TRACE("ComputeDiffCounts-" + testName);

        // Run.
        if (shouldThrow) {
            EXPECT_THROW(
                {
                    PacBio::Pancake::ComputeDiffCounts(cigar, queryVariants, targetVariants,
                                                       throwOnPartiallyMaskedIndels);
                },
                std::runtime_error);
        } else {
            auto result = PacBio::Pancake::ComputeDiffCounts(cigar, queryVariants, targetVariants,
                                                             throwOnPartiallyMaskedIndels);
            // Evaluate.
            EXPECT_EQ(expected, result);
        }
    }
}

TEST(Test_AlignmentTools_NormalizeAlignmentInPlace, ArrayOfTests)
{
    // clang-format off
    std::vector<std::tuple<std::string, std::string, std::string, std::string, std::string, bool>> testData = {
        {"Empty input", "", "", "", "", false},

        {"Exact match",
                        // Input alignment.
                        "ACTG",
                        "ACTG",
                        // Expected alignment.
                        "ACTG",
                        "ACTG",
                        false},

        {"Just a mismatch",
                        // Input alignment.
                        "CAC",
                        "CGC",
                        // Expected alignment.
                        "CAC",
                        "CGC",
                        false},

        {"Simple normalization with 1 indel and 1 mismatch",
                        // Input alignment.
                        "TTGACACT",
                        "TTG-TACT",
                        // Expected alignment.
                        "TTGACACT",
                        "TTGT-ACT",
                        false},
        {"Simple with two deletions in the query",
                        // Input alignment.
                        "AC--TAAC",
                        "ACTATAAC",
                        // Expected alignment.
                        "ACTA-A-C",
                        "ACTATAAC",
                        false},
        {"Test reverse complement alignment of the previous one. Shows that left alignment is not symmetric.",
                        // Input alignment.
                        "GTTATAGT",
                        "GTTA--GT",
                        // Expected alignment.
                        "GTTATAGT",
                        "GTTA--GT",
                        false},

        {"Test shifting of gaps on the query",
                        // Input alignment.
                        "-C--CGT",
                        "CCGAC-T",
                        // Expected alignment.
                        "CCG---T",
                        "CCGAC-T",                  // TODO: Clear meaningles "-":"-" alignments.
                        false},

        {"Test shifting of gaps on the target",
                        // Input alignment.
                        "ATAT-AGCCGGC",
                        "ATATTA---GGC",
                        // Expected alignment.
                        "ATAT-AGCCGGC",
                        "ATATTAG--G-C",
                        false},


    };
    // clang-format on

    for (const auto& data : testData) {
        // Inputs.
        const std::string testName = std::get<0>(data);
        const std::string& queryAln = std::get<1>(data);
        const std::string& targetAln = std::get<2>(data);
        const std::string& expectedQueryAln = std::get<3>(data);
        const std::string& expectedTargetAln = std::get<4>(data);
        const bool shouldThrow = std::get<5>(data);

        // Name the test.
        SCOPED_TRACE("NormalizeAlignmentInPlace-" + testName);

        std::string resultQueryAln = queryAln;
        std::string resultTargetAln = targetAln;

        // Run.
        if (shouldThrow) {
            EXPECT_THROW(
                { PacBio::Pancake::NormalizeAlignmentInPlace(resultQueryAln, resultTargetAln); },
                std::runtime_error);
        } else {
            PacBio::Pancake::NormalizeAlignmentInPlace(resultQueryAln, resultTargetAln);
            // Evaluate.
            EXPECT_EQ(expectedQueryAln, resultQueryAln);
            EXPECT_EQ(expectedTargetAln, resultTargetAln);
        }
    }
}

TEST(Test_AlignmentTools_ConvertCigarToM5, ArrayOfTests)
{
    // clang-format off
    std::vector<std::tuple<std::string, std::string, std::string, std::string, std::string, std::string, bool>> testData = {
        {"Empty input", "", "", "", "", "", false},
        {"Simple match", "ACTG", "ACTG", "4=", "ACTG", "ACTG", false},
        {"Insertion in prefix", "AAAACTG", "ACTG", "3I4=", "AAAACTG", "---ACTG", false},
        {"Insertion in suffix", "ACTGAAA", "ACTG", "4=3I", "ACTGAAA", "ACTG---", false},
        {"Deletion in prefix", "ACTG", "AAAACTG", "3D4=", "---ACTG", "AAAACTG", false},
        {"Deletion in suffix", "ACTG", "ACTGAAA", "4=3D", "ACTG---", "ACTGAAA", false},
        {"Invalid CIGAR string", "ACTG", "ACTGAAA", "4=", "", "", true},
        {"Simple CIGAR with all ops", "ATATAGCCGGC", "ATAGTAGGC", "3=1X1D1=3I3=", "ATAT-AGCCGGC", "ATAGTA---GGC", false},
    };
    // clang-format on

    for (const auto& data : testData) {
        // Inputs.
        const std::string testName = std::get<0>(data);
        const std::string& query = std::get<1>(data);
        const std::string& target = std::get<2>(data);
        const PacBio::BAM::Cigar cigar(std::get<3>(data));
        const std::string& expectedQueryAln = std::get<4>(data);
        const std::string& expectedTargetAln = std::get<5>(data);
        const bool shouldThrow = std::get<6>(data);

        // Name the test.
        SCOPED_TRACE("ConvertCigarToM5-" + testName);

        std::string resultQueryAln;
        std::string resultTargetAln;

        // Run.
        if (shouldThrow) {
            EXPECT_THROW(
                {
                    PacBio::Pancake::ConvertCigarToM5(query.c_str(), query.size(), target.c_str(),
                                                      target.size(), cigar, resultQueryAln,
                                                      resultTargetAln);
                },
                std::runtime_error);
        } else {
            PacBio::Pancake::ConvertCigarToM5(query.c_str(), query.size(), target.c_str(),
                                              target.size(), cigar, resultQueryAln,
                                              resultTargetAln);
            EXPECT_EQ(expectedQueryAln, resultQueryAln);
            EXPECT_EQ(expectedTargetAln, resultTargetAln);
        }
    }
}

TEST(Test_AlignmentTools_ConvertM5ToCigar, ArrayOfTests_RoundTrip)
{
    // clang-format off
    std::vector<std::tuple<std::string, std::string, std::string, std::string, bool>> testData = {
        {"Empty input", "", "", "", false},
        {"Simple match", "ACTG", "ACTG", "4=", false},
        {"Insertion in prefix", "AAAACTG", "ACTG", "3I4=", false},
        {"Insertion in suffix", "ACTGAAA", "ACTG", "4=3I", false},
        {"Deletion in prefix", "ACTG", "AAAACTG", "3D4=", false},
        {"Deletion in suffix", "ACTG", "ACTGAAA", "4=3D", false},
        {"Invalid CIGAR string", "ACTG", "ACTGAAA", "4=", true},
        {"Simple CIGAR with all ops", "ATATAGCCGGC", "ATAGTAGGC", "3=1X1D1=3I3=", false},
        {"Real pair of sequences",
                // >mo_m64001_190914_015449/130352465/ccs:8054-16191
                "TATGAAGAGCATGTGAGCTAGTCCACTCGAGCAACTGAAGACTACACAAAAAGCCTCTTAGTCCGAGTATTTTGGGAATTTTTCACTCAAAATTAATTGCATTTTAATAGTTGCTAATATTTCAAAAAACTTGTGCACTTATTTGCACCGGCATCAGGCTGCTCAGGCTCGCCGCTTTCCCCATCCCCTTCTCCACGAAGCTCCATCAAAAAGTCCTGCTTAAAAATGCAATAACAGGAAATGGCTACGACGGCAACAAAATACAGGAAATTGAGCGCAAAGAAATATGCTATATAGGAGATGAAAAAATGAAAACGCAAACGAGACGAAGTCCAGGCAAACAAACATGGCCGTCAGACATGAACATATGTATACAGCCGCGGTTAACCAATTAGGAGTTTTGAATTCCAGAAAAAAGAATTAATACAAATATGCACTTTATTTACTAAGATTCATTTACATTTGGATTCGTTGGTTGGTGTCAGCAAGCGTACTTGAGCTTTCAGGCAGAAAAGTGCTCATGCTTTAGCATACTGAGTCTTAGAATAGCAAAGTAGCTACAAAAAATTTAAAATCAAAGTTAAGTTAAAAGTTAAAATCAAAGTTAAAAAAAAAAAACTCGAATGAACGTTTTGAAAAAGGGGAAAAGGGAAACCCCATGGAAAGCAACACAACTATTACAATTACTTTCTGTTATATGAGATTAGTATTTTTTTGAAGACACTAAACATGTGGCCCCTGCTGTACACCTGAACCATGCAAAGCGAGCATGAGGCGAAGACAGGTGGAAAGGTGGCGTTACGGGTAAGAGCATTTTGCCCGCTGGAACTTGGCCGTATGTGCTCCGGCAATTTATGTCCAAGGTGCGTGCCACCGAAATGAACCGGTGAACCTGGACGTTACGCTCGTTTCTTGGGGGAAAGTCAGGCAGTTGACCAACCTGTGTCCTTTTTCCTGGTTCTGCGCTTAACGTGTTGCCATTAAATTTTATGGTCAGGTCTTCCTCTCTTTGTCCGTGTTTGCATTTTGCTCAATTTTGTCGCCGCAGAACATAAAAAAGGGAGTCAGGAGGGGCTGAAAAAGAAGCGAAGGAATATTGCACTGTGGACGACATGTGAAAACATATTTTTATGGGACCCCAAAGGGTAATTGGACATTGCACGGGATATCCTAGCTGGCTAATGAGTGTGTCGAGCCTTTTTAGATAATTAATTGAGAGCTTTTCAAATCATCTGTGAAGTATAATGACCACACACACACAATAAGGCAACGAAACGTTGAAGTGTGATTGCACTAATGAGTAGATTAAAGTCAGCTGGAGATAACATCTTAATAAATCCAAGATAATTGAAATATCCTTTCGAGATAGTTGTATAATTTGTAATTCTCCATTCGGGAATATTTGTTTCATTGGCTTCAAGTGAAACCAATTGTGCCGCACTTCATCCAGCGAAACTTTGACGGCTAACTCGCTTAATGCGCGGCTTAAGCCAAACTACGAGTAAGATATTTTGTCTGTGGCAACTTGAAAGTTGCAAGTAATTAGTGGCACTTGAGCAACTACTAGGCTGTCGGGATACGAAAATGGGGCAGGGCATACCCAGATAGCTGGATGCTCCAAAAGTTATGTAGATCACGAAAATCAGTCCTTTGAGAACAAAGAACTTTTAGCAATCTTAAAATATTATTTCTAAGTATTCTAATCGTTATACTTTATTAGTATTTTAGTTCTATGAGTTTCTTAAAAAAAAAAAAACACTGATCTTACAAGCTCGCTTTCCACCTATCCTCTTTCAGACCCATCAGTCTAGTGCTGGAGCCAATCCTGAAGGAGCCCGGCGACATCTCCTTTACGTGGCTCGAGATTCACCGCACCAAGAACGCGCTGCTGCAGCAGTTGGACCTGCATGTGAACGCTAGCGTGACCACGGGCATGGGCTCCCAGACCTCCTCCCAGTCCGCAGCCGGCTCCAACACGGGCATTGGGGGACTGGGCACCCAGTGGGTCAGGATCTGCTGCCAATGCCGGGGGCTCTGCTGGCGGAGCAACGGCCAACGAGACGCTGAACGAGTTCGGATTCTTTCCCAAGGAGTCCGACTGCGAATACAAATGTCCTGAAATTGATGCATGCATTGCGGCCAGCTTGTGGTGCGATGGTGAGTACAGATTGCCAGTCGAGGCCACAGACCATTTCTCAGGGATTTGCCGAATGATGGGTTTTTGGGATGGAACCACTACATGCAAAGCAGTTCATAGAAGACGACGTGAACAAGAAACTTTCGAGCGTTCATCAATTCGGGTCAATGGGCAACTACTTAAGAAAAGGGCTTCTGAATAAAGAACATACATGTTGAATTATCAAGCAAATTTTGCATTTAAAGTAACTTAAGGTACGAATTCTTACTCAAAATTCAGTAGAAACGTTCTTAGGAAAATTAATACCACTAAAGAGCTAGCAATAGGCTTATCTCCTTACATTGGCATGCTTGGTCTCAATGATTGTACTTTGGGAATTTTGATATCACTATGGAGAGTCCTGGATCCCCATCTTCGCCTGGAGCTTTGCACTATTTGTACGTGATTTGGCGCTCAGTCTTGGGGACTTGGGTTACCCATTTGCAGCGTGCCGCCGGCAGAAGGCAACCAAAAATTTCACTTGTGGCTGGCTTGTCTGTTAATTTCGAAGCGGGCTGTCCACGTCCTTGTAATCAGGCAAAGTTTCCTGACTGCTGCACCAGGCTGCAGTGCACCACTCCTGCTGGTTGTCGGTTGCTAGTTGTTAGTTGCTTCGTTGGTCGGTTGTCGTTGTCGGCGGTTGCTGGTTGCTGGTTGTTGGCGATTGTTGGCCGACTGGCGGGATGGCGGGATGGCGGGATTGGGCCGGTTGTTGGACAAACCCGCAGCGGCAAAACCCCAAAATGTCTTGCAGAAAAAGTTTATCCCCCAGTCCGGCACCAGGCCAACACGAAATGGCCGCAACAGCTGTGGCAGACTTATTTTCAATGCACTTGGGAGATGAAGGGCGCAAAAAAGTAGATGCCAAGTGGGAAAGTTTATAGACTAGGAAATACCGTTTTGTTAGGAATCCGCCGTTCATTTTGTATTTACATGAACCTACCAATACATATTTAAAAGTCATATATTATCCTCCAAGTGAATGTAATTATTTTAATTTTGCTGCTGAAATAGTAAATACTTTAATATGTACATTGAAGATTGTGTGATATACAAGTATAAAATTGCATTGATTTTCCTGCATTTTCCTACTTATACTCCAGAAGCCGCACAGACTCTTTAAATCTTACGCTCGGTGGAAAAAATACAAAGGCTAGCCTGGGCTTTGCCTTTTAGCTTCAATTATTTATGGGCAGCCAGGCGCGTCTGCATCTGGAGCTGGTTAACCGCAAAATCCAGGGCATTAAACAAGCTTTGGACCCCACGGGCCACTCGAGTGACACTGCATTCAATTTGTAACCTTAGAATAGGCCACCCAGCTGCGCTGCTCCTCCCAAACCCCAATCTGCTTTCAATTTTTCAATCAGCGCCTGTACCTCTATCTAGGATACTGAATCTGTTCCTACACTGACAAAAAGATTCTTTTCACATATTACCTTCAAAATATCTATTTGAAACTATATGCAAAATAAAAAAGATACATCCAACTAGGGATTGTATTGTTGAAACTTCTTAAATATATATTGCTTGTCACTAAGACTATTTATGGATAATTTGAAAGAGAGTCTTATAGACCTGCTTAGACCTTAAAAATATTTTTAAGCGGTGTTTTGGGGTAATGAACAAATATCTGAATCAAAAGCTGTAACAGATTAATTCCACCTAAGAAGTGACCGCTTTTTTCTAAAAGTGTCCTTATAATCTCTCCTTTTCAGGCCACCACAACTGCCCCAGCGGTTTCGATGAGTCGGAGGAGGAGTGCGGCACCGCTCGCAAGTTGCTGGAACTGCCAGGCGGAGTTTTCGCCGCCCTTGGATGCATTGCAGCAGCCCTGACCGCCTGCCTGATATTCTGCATGTTCGGCCTGATGAGGAAGCGGAAAAAGTCGGTGGTGCAGAAGGGCGCCGGAGTGGGCGGCATGGGCGGAATGGGAGTGGGTGGGATGGGCGTGGGCATCGGGTTCATGAACGGAGCGAGCAACGGCAACGTGTCCGCCGCCAGCATCGGTACCTTGAAGAAGGACTTCAAGAAGGAGGCGCTGTACATCGATCCGGCCTCCTGACACGGACTACAGCCGGTCGAATATATCCACGTGGACCGGGAGACCACCGTGTGATATGTTGCCGCCGGCTGTGGCATTCCCTATAATCATGTGGAGCTTTCCATGGAGCTGGCCTGACCGATGGATCAGTCCGCTGATCCCGTTCAAAATCTGGCCTGGAGATTTTAAAGATACAGAGTCGACGCCTGGAGGTATGAGAACACTAAATGAACCTACATTTTAAGCCGATAACTTTTTCTAAACGCTTCTGGGAGTAACAAATTGCGCACTAGCGAAGTTGTGTGTATGTAGTCATAAGGAAAATACTGTTATTAAGTTGCGTATGGTTTTTAGTTTGTATATATAATTTAATTTTAGGCATAAATAATTTGAATATCGTTGAAGGAGTCAGAACTGTGCACAGAACTTTGTGCAAGAACTTTTTGAATGTAGACTTTGTTTCCTATCATTAATCTTTCAATATCAGAAAATTATCATCGAGTTCAAACTAAAAATCTTTATAATTTCATGCATTAAAGGTTTCCAAATTAAAATCAAAATACATAAACTTTTCACAGTAAATAAAATTTTAAAATTCACGTCTTAAGGAATGTGTATATTAAGCAATAAAATCACTCACAACTTTTTAAAACTCAAAATGTTCGAGTAATATTAATTTATTGCTATTCACATGGTTGCCCGTGCTGACCCTTGCATAATTGTTTAGTTCGATTATGAGAGTTATCAATTTTTCACTTTAAGACTACATACACAAACGAAAAGAGGACATAAAACTATTGTAAATTAGGTACCTTTAAACCAGTTAGATCAGAGAACAACAATTCGTGGGCGGTTCGTAGATGACTTCGGCATATTACCACCCACATAGCACACAGCACCCATTAAAGGAAGTAAACCAATGACATTTAGCAACTAAATGAAAGATCATAAACTACAATGAAATTGTATAGAAGAGCAACTCAAAGAGCAGATCCCTATGCGCAGGGCCACTCTGGCACTAATTTTAAGCGGTTAGACTAGGTAAAAAACTGTTAGCAGCGATTGGAAAGGAGTTGAGGAAGTTCGAGGGTAATCTCTACAGTGTTGTGATCTAAAGCATGCACATACTGAACTATATTTAAACACTATGTTTAACGTCGAAGCATTTACATGACGAGTACTCAGAACCGGGGGCGTATATTGTATAGGTACTAGCTAATAGCTAAATAGTCAACTAAGAGGCAGCACTTTTCTGTTATTATCGTACAGCACAAACACATAAGCAAATGACATCACCTAGTGGTTTATTATATACAAGTCGATACAGCAAGTTAATTAAATGCAATTGTTATTAATGAATTAAGTAGAGCTCCAAACCAGAGATATGTTCACTGCACTTAATCGAGAGAACTGCCAAAGGATTAGACAACAAAACAACAAAACAACAAAACAATCAATGTAAATATTTCTGTCGAGCAATAACTTGTTTTGCGTTAGCATTTATGTTCGTAATTTGAACACAAACCGATTCCACTTAAGTAGCTGCAGTGAGTAATATTCCTAAAGTCAATAGACTTTGTAAAAACTACGAGCGAAGAGTAATAATAAACTGTATGCAAAACTACCATTTAAGACACGGCTCATAACCCATATAAACAATCAAACTTTTGGCAGTTGCAGTAGGATTAAGCAAACAACACCAAGAAGAAGCGGAAAAAGGACATAAGTACTTTATTTTATAGACACTGAACTAGGGCAAAAATTATGAAACTGAATACTGAAAGCAAAGCAAATACTGGCTGCAGTGGGAATTTTAAAAATTTAAGAACCTCTGCCAACGCCGGCAGAGTCCACATCCATAAAAACCCATTAAAATTTCAACTGCAATGAAACTTTGCCGAAAAGCGTGTGTATAAAAAATATATGATTTAGCCAGCATAAAAAGTAAATGAAAATATATACATAAATTAAATCAAATATTTTCGGTTGGCATTTTTTGATAGCAAATAACCTAAATATTTATAGCTAAAAAAGAAAGAAAGAAAAATTACCAAACGCACTTTGTACATTTTATAAAAAACGAGCAGTTTATATTATTTCCAGACACACAGCATAGCATACACATACATATATGATATTTTACACACATTCGAATAAGAGATAAAATAATACCGAAACATATACGAGAGCGAAAAGAAATTAAATTAATTAAAAAAGGCGCATAGGCCGAAAATGCATTTTTACTGCCGCCAATTGAGGACAAAACTAAAATCGAAAACACGGCCGGCAGGGCATTTGAATAGTTAGCTTTAAGCCCGGCTAAAAATGCTGTAACTACTTTGTTGGTCACCATTGTCACTCACACAAACTCACATCTCAATCATCAGTTCGGTTCAGTTCCTTCTCTCAACCAAACCAAACTATGGAAATTTCTGTGCCAACTGTGTTTAAAGTAGAATTCAAAGATAAAATACGAGATAAAAGAGGGGAGGCCGGCGGATGGAGAAATTGAATTATGACCGCACACGAAACATTGCATATAATTTATTTGAAATGTTCAAAAATAAAGTAGACGAAACGAGAAGAAAACGCTGCGCGATGGTAAATAACAATATAATACAGACAGCTGTGTAAATAGTCAAACTGCGTACAGTCTACGCATATTTAATTAAAACTCGAAAAAAACTGATAAAACTGAAAAAACTGAAAAAAGTGACCCAATGAGATGGAATTTGAGAAAGGAGATCATGTCTTAAGTTCCATTTCGGAATTCACCCGATTTTCATTTCCTTTCACCGTTCGATAGAGTAAGTTTCACTTGTATGATAAATAAAACTAAACCAAATATTTATATAGAGATATTATACGTATATGAAAAGGAAACCTTAAGTTAATGACGCACTGAATCAAATAGATGAAACGAGTATTTAAAACCAAGACAATTAAAACCAAATAGGAACTTTAAGCAAAATAAAACCGATGAAAAAATGAAACTGAACACAGAGCACCTTTCTCTATTAGCCAAAATTGCAAAAGTGGGTGGTTGGGCGTTTGATTTTCGGGTGGAATAGATGGGGTGTTTAGTGGGTGGAATGCCAACTGTTCCGTTTGACGCTCGGGGGAAAACCGTTGGCCAAACTAAGGCGAACTCAAATCGTGTCCAAATAGCTGGCAGTCAAGTAGGCCACAGATTTTGATTGCCCGCATCCGGGCAGCAATCGGTGGACATCCTCGAAAATTCGGAACATCGAACGAGATGCTGTTAAAGCCAAGCCAAGCCGGTAGCTCCAGGCGTGTGTTCTCACCGGTCGGAACTCGGATTGGCCGTGGGATCGGGTTGTGGCTATGGGCTTTCCATGTCCCCCGCCTCGCCCACGCTGAGCCCGCATAAAACATAAGCCACCGGCCTGGTGCTCGTGTGGTGTGTGCTCTATATCCCACAGATTCCGCCTCCTCCCTCGGAATGCCCAGTTCTCTGACTCCGTGTTTATGGCCGGCAGCTCTGTTTCTGGACCTGGCAACTCTGGCCCCTGGACTTGGCCAAAATGGCTTTGATGGGCGTGGGCGACAAAAGGACCGGAAATCGCTAACGAGGCGTGGTAAAGGAGTTTTTGCAGTAAGTTTGAATATCAAGGTAGATACTCAAAAAATTGAATAACAAAGCATAATCAAATTTAATGGAGGGATCTACCTTGAATGTATTTTGTTATTATTGTATTACAGCTTATCAAACTAATCACTAGTATGTTTGTTGCTTTTATTTATTATTAACACATTTAAATATCTGGCTATTAAAGTTAATAATTTTGCAGTTCTTTTTATGCTGTGGTCATTAATATACATATAATACATATATCTACTAAAAGTATTCACTCTAGATTAATGCAAATTGATGAAATTATTTATTTATCGTACGCACTT",
                // >fa_m64001_190914_015449/101320295/ccs:10320-18430 reversed
                "TATGAAGAGCATGTGAGCTAGTCCACTCGAGCAACTGAAGACTACACAAAAAGCCTCTTAGGTCGAGTATTTTGGGAATTTTTCACTCAAAATTAATTGCATTTTAATAGTTGCTAATATTTCAAAAAACTTGTGCACTTATTTGCACCGGCATCAGGCTGCTCAGGCTCGCCGCTTTCCCCATCCCCTTCTCCATGAAGCTCCATCAAAAAGTCCTGCTTAAAAATGCAATAACAGGAAATGGCTACGACGGCAACAAAATACAGGAAATTGAGCGCAAAGAAATATGCTATATAGGAGATGAAAAAATGAAAACGCAAACGAGACGAAGTCCAGGCAAACAAACATGGCCGTCAGACATGAACATATGTATACAGCCGCGGTTAACCAATTAGGAGTTTTGAATTCCAGAAAAAAGAATTAATACAAATATGCACTTTATTTACTAAGATTCATTTACATTTGGATTCGTTGGTTGGTTTCAGCAAGCGTACTTGAGCTTTCAGGCAGAAAAGTGCTCATGCTTTAGCATACTGAGTCTTAGAATAGCAAAGTAGCTACAAAAAATTTAAAATCAAAGTTAAGTTAAAAGTTAAAATCAAAGTTAAAAAAAAAAAACTCGAATGAACGTTTTGAAAAAGGGGAAAAGGGAAACCCCATGGAAAGCAACACAACTATTACAATTACTTTCTGTTATATGAGATTAGTATTTTTTTTGAAGACACTAAACATGTGGCCCCTGCTGTACACCTGAACCATGCAAAGCGAGCATGAGGCGAAGACAGGTGGAAAGGTGGCGTTACGGGTAAGAGCATTTTGCCCGCTGGAACTTGGCCGTATGTGCTCCGGCAATTTATGTCCAAGGTGCGTGCCACCGAAATGAACCGGTGAACCTGGACGTTACGCTCGTTTCTTGGGGGAAAGTCCGGCAGTTGACCAACCTGTGTCCTCTTTCCTGGTTCTGCGCTTAACGTGTTGCCATTAAATTTTATGGTCAGGTCTTCCTCTCTTTGTCCGTGTTTGCATTTTGCTCAATTTTGTCGCCGCAGAACATAAAAAAGGGAGTCAGGATGGGGCTGAAAAAGAAGCGAAGGAATATTGCACTGTGGACGACATGTGAAAACATATTTTTATGGGACCCCAAAGGGTAATTGGACATTGCACGGGATATCCTAGCTGGCTAATGAGTGTGTCGAGCCTTTTTAGATAATTAATTGAGAGCTTTTCAAATCATCTGTGAAGTATAATGACCACACACACACAATAAGGCAACGAAACGTTGAAGTGTGATTGCACTAATGAGTAGATTAAAGTCAGCTGGAGATAACATCTTAATAAATCCAAGATAATTGAAATATCCTTTCGAGATAGTTGTATAATTTGTAATTCTCCATTCGGGAATATTTGTTTCATTGGCTTCAAGTGAAACCAATTGTGCCGCACTTCATCCAGCGAAACTTTGACGGCTAACTCGCTTAATGCGCGGCTTAAGCCAAACTACGAGTAAGATATTTTGTCTGTGGCAACTTGAAAGTTGCAAGTAATTAGTGGCACTTGAGCAACTACTAGGCTGTCGGGATACGAAAATGGGGCAGGGCATACCCAGATAGCTGGATGCTCCAAAAGTTATGTAGATCACGAAAATCAGTCCCTTGAGAACAAATAACTTTTAGCAATCTTAAAATATTATTTCTAAGTATTCTAATCGTTATACTTTATTAGTATTTTAGTTCTATGACTTTCTTTAAAAAAAAAAAAAAAAACACTGATCTTACAAGCTCGCTTTCCACCTATCCTCTTTCAGACCCATCAGCCTAGTGCTGGAGCCAATCCTGAAGGAGCCCGGCGACATCTCCTTTACGTGGCTCGAGATTCACCGCACCAAGAACGCGCTGCTGCAGCAGTTGGACCTGCATGTGAACGCTAGCGTGACCACGGGCATGGGCTCCCAGACCTCCTCCCAGTCCGCAGCCGGCTCCAACACGGGCATTGGGGGACTGGGCACCAGTGGGTCAGGATCTGCTGCCAATGCCGGGGCTCTGCTGGCGGAGCAACGGCCAACGAGACGCTGAACGAGTTCGGATTCTTTCCAAGGAGTCCGACTGCGAATACAAATGTCCTGAAATTGATGCATGCATTGCGGCCAGCTTGTGGTGCGATGGTGAGTACAGATTGCCAGTCGAGGCCACAGACCATTTCTCAGGGATTTGCCGAATGATGGGTTTTTGGGATGGAACCACTACATGCAAAGCAGTGCATAGAAGACGACGTGAACAAGAAACTTTCGAGCGTTCATCAATTCGGGTCAATGGGCAACTACTTAAGAAAAGGGCTTCTGAATAAAGAACATACATGTTGAATTATCAAGCAAATTTTGCATTTAAAGTAACTTAAGGTACGAATTCTTACTCAAAATTCAGTAGAAACGTTCTTAGGAAAATTAATACACTAAAGAGCTAGCAATAGGCTTATCTCCTTACATTGGCATGCTTGGTCTCAATGATTGTACTTTGGGAATTTTGATATCACTATGGAGAGTCCTGGATCCCCATCTTCGCCTGGAGCTTTGCACTATTTGTACGTGATTTGGCGCTCAGTCTTGGGGATTTGGGTTACCCATTTGCAGCGTGCCGCCGGCAGAAGGCAACCAAAAATTTCACTTGTGGCTGGCTTGTCTGTTAATTTCGAAGCGGGCTGTCCACGTCCTTGTAATCAGGCAAAGTTTCCTGACTGCTGCACCAGGCTGCAGTGCACCACTCCTGCTGGTTGTCGGTTGCTAGTTGTTAGTTGCTTCGTTGGTCGGTTGTCGTTGTCGGCGGTTGCTGGTTGCTGGTTGTTGGCGATTGTTGGCCGACTGGCGGGATGGCGGGATGGCGGGATTGGGCCGGTTGTTGGACAAACCCGCAGCGGCAAAACCCCAAAATGTCTTGCAGAAAAAGTTTATCCCCCAGTCCGGCACCAGGCCAACACGAAATGGCCGCAACAGCTGTGGCAGACTTATTTTCAATGCACTTGGGAGATGAAGGGCGCAAAAAAGTAGATGCCAAGTGGGAAAGTTTATAGACTAGGAAATACCGTTTTGTTAGGAATCCGCCGTTCATTTTGTATTTACATGAACCTACCAATACATATTTAAAAGTCATATATTATCCTCCAAGTGAATGTAATTATTTTAAATTTGCTGCTGGAATAGTAAATACTTTAATATGTACATTGAAGATTGTGTGATATACAAGTATAAAATTGCATTGATTTTCCTGCATTTTCCTACTTATACTCCAGAAGCCGCACAGACTCTTTAAATCTTACGCTCGGTGGAAAAAATACAAAGGCTAGCCTGGGCTTTGCCTTTTAGCTTCAATTATTTATGGGCAGCCAGGCGCGTCTGCATCTGGAGCTGGTTAACCGCAAAATCCAGGGCATCAAACAAGCTTTGGACCCCACGGGCCACTCGAGTGACACTGCATTCAATTTGTAACCTTAGAATAGGCCACCAGCTGCGCTGCTCCTCCCAAACCCCAATCTGCTTTCAATTTTTCAATCAGCGCCTGTACCTCTATCTAGGATACTGAATCTGTTCCTACACTGACAAAAAGATTCTTTTCAAATATTACCTTCAAAATATCTATTTGAAACTATATGCAAAATAAAAAAGATACATCCAACTAGGGATTGTATTGTTGAAACTTCTTAAATATATATTGCTTGTCACTAAGACTATTTATGGATAATTTGAAAGAGATTCTTATAGACCTGCTTAGACCTTAAAAATATTTTTAAGCGGTGTTTTGGGGTAATGAACAATATCTGAATCAAAAGCTGTAACAGATAATTCCACCTAAGAAGTGACCGCTTTTTTCTAAAAGTGTCCTTATAATCTCTCCTTTTCAGGCCACCACAACTGCCCCAGTGGTTTTCGATGAGTCGGAGGAGGAGTGCGGCACCGCTCGCAAGTTGCTGGAACTGCCAGGCGGAGTTTTCGCCGCCCTTGGATGCATTGCAGCAGCCCTGACCGCCTGCCTGATATTCTGCATGTTCGGCCTGATGAGGAAGCGGAAAAAGTCGGTGGTGCAGAAGGGCGCCGGAGTGGGCGGCATGGGCGGAATGGGAGTGGGTGGGATGGGCGTGGGCATCGGGTTCATGAACGGAGCGAGCAACGGCAACGTGTCCGCCGCCAGCATCGGTACCTTGAAGAAGGACTTCAAGAAGGAGGCGCTGTACATCGATCCGGCCTCCTGACACGGACTACAGCCGGTCGAATATATCCACGTGGACCGGGAGACCACCGTGTGATATGTTGCCGCCGGCTGTGGCATTCCCTATGATCATGTGGAGCTTTCCATGGAGCTGGCCTGACCGATGGATCAGTCCGCTGATCCCGTTCAAAATCTGGCCTGGAGATTTTAAAGATACAGAGTCGACGCCTGGAGGTATGAGAACACTAAATGAACCTACATTTTAAGCCGATAACTTTTTCTAAACGCTTCTGGGAGTAACAAATTGCGCACTAGCGAAGTTGTGTGTATGTAGTCATAAGGAAAATACTGTTATTAAGTTGCGTATGGTTTTTAGTTTGTATATATAATTTAATTTTAGGCATAAATAATTTGAAAATCGTTGAAGGAGTCAGAACTGTGCACAGAACTTTGTGCAAGAACTTTTTGAATGTAGACTTTGTTTCCTATCATTAATCTTTCAATATCAGAAAATTATCTTGGAGTTCAAACTAAAAATCTTTATAATTTCATGCATTAAAGGTTCCAAATTAAAATCAAAATACATAAACTTTTCACAGTAAATAAAATTTTAAAATTCACGTCTTAAGGAATGTGTATATTAAGCAATAAAATCACACACAACTTTTTAAAACTCAAAATGTTCGAGTAATATTAATTTATTGCTAGTCACATGGTTGCCCGTGCTGACCCTTGCATAATTGTTTTAGTTCGATTATGAGAGTTAACAATTTTTCACTTTAAGACTACATACACAAACGAAAAGAGGACATAAAACTATTGTAAATTAGGTACCTTTAAACCAGTTAGATCAGCTAGAGAACAACAATTCGTGGGCGGATCGTAGATGACTTCGGCATATTACCACATAGCACACAGCACCCATTAAAGAAAGTAAACCAATGACATTTAGCAACTAAATGAAAGATCATAAACTACAATGAAATTGTATAGAAGAGCAACTCAAAGAGCAGATCCCTATGCGCAGGGCCACTCTGGCACTAATTTTAAGCGGTTAGACTAGGTAAAAAACTGTTAGCAGCGATTGGAAAGGAGTTGAGGAAGTTCGAGGGTAATCTCTACAGTGTTGTGATCTAAGGCATGCACATACTGAACTATGTTTAAACACTATGTTTAACGTCGAAGCATTTACATGACGAGTACTCAGAACCGGGGGCGTATATTGTATAGGTACTAGCTAGTAGCTAAATAGTCAACTAAGAGGCAGCACTTTTCTGTTATTATCGTACAGCACAAACAAGCAAATGACATCACCTAGTGGTTTATTATATACAAGTCGATACAGCTAGTTAATTAAATGCAATTGTTATTAATGAATTAAGTAGAGCTCCAAACCAGAGATATGCACTTAATCGAGAGAACTGCCAAAGGATTAGACAACAAAACAACAAAACAATCAATGTAAATATTTCTGTCGAGCAATAACTTGTTTTGCGTTAGCATTTATGTTCGTAATTTGAACACAAACCGATTCCACTTAAGTAGCTGCAGTGAGTAATATTCCTAAAGTCAATAGACTTTGTAAAAACTACGAGCGAAGAATAATAATAAACTGTATGCAAAACTACCATTTAAGACACGGCTCATAACCCATATAAACAATCAAACTTTTGGCAGTTGCAGTAGGATTAAGCAAACAACACCAAGAAGAAGCGGAAAAGGACATAAGTACTTTATTTTATAGACACTGAACTAGGGCAAAAATTATGAAACTGAATACTGAAAGCAAAGCAAATACTGGCTGCAGTGGGAATTTTAAAAATTTAAGAACCTCTGCCACGCCGGCAGAGTCCACATCCATAAAAACCCATTAAATTTCAACTGCAATGAAACTTTGCCGAAAAGCGTGTGTATAAAAAATATATGATTTAGCCAGCATAAAAGTAAATGAAAATATATACATAAATTAAATCAAATATTTTCGGTTGGCATTTTTTGATAGCAAATAACCTAAATATTTATAGCTAAAAAAGAAAGAAAGAAAAATTACCAAAACGCACTTTGTACAATTTAATAAAAAACGAGCAGTTTATATTATTTCCAGACACACAGCATAGCATACACATACATATATGATATTTTACACACATTCGAATAAGAGATAAAATAATACCGAAACATATACGAGAGCGAAAAGAAATTAAATTAATTAAAAAGGGCGCATAGGCCGAAAATGCATTTTTACTGCCGCCAATTGAGGACAAAACTAAAAATCGAAACACGGCCGGCAGGGCATTTGAATAGTTAGCTTTAAGCCCGGCTAAAAATGCTGTAACTACTTTGTTGGCCACCATTGTCACTCACACAAACTCACATCTCAATCATCAGTTCGGTTCAGTTCCTTCTCTCAACCAAACCAAACTATGGAAATTTCTGTGCCAACTGTGTTTAAAGTAGAATTCAAAGATAAAATACGAGATAAAAGAGGGGGGGCCGGCGGATGGAGAAATTGAATTATGACCGCACACGAAACATTGCATATAATTTATTTGAAATGTTCAAAAATAAAGTAGACGAAACGAGAAGAAAACGCTGCGCGATGGTAAATAACAATATAATACAGACAGCTGTGTAAATAGTCAAACTGCGTACAGTCTACGCATATTTAATTAAAAACCGAAAAAAAAAACTGATAAAACTGAAAAAAGTGACGCAATGAGATGGAATTTGAGAAAGGAGATCATGTCTTAAGTTCCATTTCGGAATTCACCCGATTTTCATTTCCTTTCACCGTTCGATAGAGTAAGTTTCACTTGTATGATAAATAAAACTAAACCAAATATTTATATAGAGATATTATACGTATATGAAAAGGAAACCTTAAGTTAATGACGCACTGAATCAAATAGATGAAACGAGTATTTAAAACCAAGACAATTAAAACCAAAAGGAACTTTAAGCAAAATAAAACCGATGAAAAAATTAAACTGAACACAGAGCACCTTTCTCTATTAGCCAAAATTGCAAAAGTGGGTGGTTGGGCGTTTGATTTTCGGGTGGAATAGATGGGGTGTTTAGTGGGTGGAATGCCAACTGTTCCGTTTGACGCTCGGGGGAAAACCGTTGGCCAAACTAAGGCGAACTCAAATCGTGTCCAAATAGCTGGCAGTCAAGTAGGCCACAGATTTTGATTGCCCGCATCCGGGCAGCAATCGGTGGACATCCTCGAAAATTCGGAACATCGAACGAGATGCTGTTAAAGCCAAGCCAAGCCGGTAGCTCCAGGCGTGTGTTCTCACCGGTCGGAACTCGGATTGGCCGTGGGATCGGGTTGTGGCTATGGGCTTTCCATGTCCCCCGCCTCGCCCACGCTGAGCCCGCATAAAACATAAGCCACCGGCCTGGTGCTCGTGTGGTGTGTGCTCTATATCCCACAGATTTCGCCTCCTCCCTCGGATTGCCCAGTTCTCTGACTCCGTGTTTATGGCCGGCAGCTCTGTTTCTGGACCTGGCAACTCTGGCCCGTGGACTTGGCCAAAATGGCTTTGATGGGCGTGGGCGACAAAAGGACCGGAAATCGCTAACGAGGCGTGGTAAAGCAGTTTTTGCAGTAAGTTTGAATATCAAGGTAGATACTCAAAAAATTGAATAACAAAGCATAATCAAATTTAATGGAGGGATCTACCTTGAATGTATTTTGTTATTATTGTATTACAGCTTATCAAACTAATCACTAGTATGTTTGTTGCTTTTATTTATTATTAACACATTTAAATATCTGGCTATTAAAGTTAATAATTTTGCAGTTCTTTTTATGCCGTGGTCATTAATATACATATAATACATATATCTACTAAAAGTATTCACTCTAGATTAATGCAAATTGATGAAATTATTTATTTATCGTACGCACTT",
                // Same CIGAR string as produced by the SES2 global aligner (or Edlib).
                "61=1D2=1I131=1X284=1X235=1D209=1X23=1X120=1D579=1X11=1X74=1X6=1D13=4D50=1X192=1I31=1I54=1I164=1X191=1I159=1X569=1X10=1X235=1X71=1I111=1X135=1X61=1I26=1I80=1X5=1D406=1X288=1X101=1X1=1X43=1I93=1X49=1X38=1D19=1X87=2D3=1D1=1D17=1X27=2I1=2I22=1X206=1X20=1X81=1X58=1I1=3I46=1X58=2I4=4I1=1I46=8I145=1X118=1I119=1I36=1I68=1I111=1D13=1X4=1D142=1X55=1D6=1I68=1X142=1X184=1D1=1I9=3I1=1I2=2I4=1X18=1X235=1I35=1X456=1X16=1X66=1X74=1X218=1X96=",
                false},
    };
    // clang-format on

    for (const auto& data : testData) {
        // Inputs.
        const std::string testName = std::get<0>(data);
        const std::string& query = std::get<1>(data);
        const std::string& target = std::get<2>(data);
        const PacBio::BAM::Cigar cigar(std::get<3>(data));
        const bool shouldThrow = std::get<4>(data);

        // Name the test.
        SCOPED_TRACE("ConvertM5ToCigar-" + testName);

        std::string queryAln;
        std::string targetAln;

        // Run.
        if (shouldThrow) {
            EXPECT_THROW(
                {
                    PacBio::Pancake::ConvertCigarToM5(query.c_str(), query.size(), target.c_str(),
                                                      target.size(), cigar, queryAln, targetAln);
                    PacBio::BAM::Cigar resultCigar =
                        PacBio::Pancake::ConvertM5ToCigar(queryAln, targetAln);
                },
                std::runtime_error);
        } else {
            PacBio::Pancake::ConvertCigarToM5(query.c_str(), query.size(), target.c_str(),
                                              target.size(), cigar, queryAln, targetAln);
            PacBio::BAM::Cigar resultCigar = PacBio::Pancake::ConvertM5ToCigar(queryAln, targetAln);
            EXPECT_EQ(cigar, resultCigar);
        }
    }
}

TEST(Test_AlignmentTools_NormalizeCigar, ArrayOfTests)
{
    /*
     * Tests here are identical to Test_AlignmentTools_NormalizeAlignmentInPlace, with the only addition of a test
     * that throws.
    */

    // clang-format off
    std::vector<std::tuple<std::string, std::string, std::string, std::string, std::string, bool>> testData = {
        {"Empty input", "", "", "", "", false},

        {"Exact match",
                        // Query.
                        "ACTG",
                        // Target.
                        "ACTG",
                        // Input CIGAR.
                        "4=",
                        // Normalized CIGAR.
                        "4=",
                        false},

        {"Just a mismatch",
                        "CAC",
                        "CGC",
                        "1=1X1=",
                        "1=1X1=",
                        false},

        {"Simple normalization with 1 indel and 1 mismatch",
                        /*
                         * Input:
                         *  TTGACACT
                         *  ||| X|||
                         *  TTG-TACT
                         * Expected alignment.
                         *  TTGACACT
                         *  |||X |||
                         *  TTGT-ACT
                        */
                        "TTGACACT", "TTGTACT", "3=1I1X3=", "3=1X1I3=", false},
        {"Simple with two deletions in the query",
                        /*
                         * Input:
                         *  AC--TAAC
                         *  ||  ||||
                         *  ACTATAAC
                         * Expected:
                         *  ACTA-A-C
                         *  |||| | |
                         *  ACTATAAC
                        */
                        "ACTAAC", "ACTATAAC", "2=2D4=", "4=1D1=1D1=", false},
        {"Test reverse complement alignment of the previous one. Shows that left alignment is not symmetric.",
                        /*
                         * Input:
                         *  GTTATAGT
                         *  ||||  ||
                         *  GTTA--GT
                         * Expected:
                         *  GTTATAGT
                         *  ||||  ||
                         *  GTTA--GT
                        */
                        "GTTATAGT", "GTTAGT", "4=2I2=", "4=2I2=", false},

        {"Test shifting of gaps on the query",
                        /*
                         * Input:
                         *  -C--CGT
                         *   |  | |
                         *  CCGAC-T
                         * Expected:
                         *  CCG---T
                         *  |||   |
                         *  CCGAC-T
                        */
                        "CCGT", "CCGACT", "1D1=2D1=1I1=", "3=2D1=", false},

        {"Test shifting of gaps on the target",
                        /*
                         * Input:
                         *  ATAT-AGCCGGC
                         *  |||| |   |||
                         *  ATATTA---GGC
                         * Expected:
                         *  ATAT-AGCCGGC
                         *  |||| ||  | |
                         *  ATATTAG--G-C
                        */
                        "ATATAGCCGGC", "ATATTAGGC", "4=1D1=3I3=", "4=1D2=2I1=1I1=", false},

        {"Bad input CIGAR string, should throw",
                        "ACTG", "ACTG", "10=", "", true},
    };

    // clang-format on

    for (const auto& data : testData) {
        // Inputs.
        const std::string testName = std::get<0>(data);
        const std::string& query = std::get<1>(data);
        const std::string& target = std::get<2>(data);
        const PacBio::BAM::Cigar inputCigar(std::get<3>(data));
        const PacBio::BAM::Cigar expectedCigar(std::get<4>(data));
        const bool shouldThrow = std::get<5>(data);

        // Name the test.
        SCOPED_TRACE("NormalizeAlignmentInPlace-" + testName);

        // Run.
        if (shouldThrow) {
            EXPECT_THROW(
                {
                    const PacBio::BAM::Cigar result = PacBio::Pancake::NormalizeCigar(
                        query.c_str(), query.size(), target.c_str(), target.size(), inputCigar);
                },
                std::runtime_error);
        } else {
            const PacBio::BAM::Cigar result = PacBio::Pancake::NormalizeCigar(
                query.c_str(), query.size(), target.c_str(), target.size(), inputCigar);
            // Evaluate.
            EXPECT_EQ(expectedCigar, result);
        }
    }
}

TEST(Test_AlignmentTools_TrimCigar, ArrayOfTests)
{
    // clang-format off
    // <cigar, windowSize, minMatches, expectedTrimmedCigar, expectedClippedFrontQuery, expectedClippedFrontTarget, expectedClippedBackQuery, expectedClippedBackTarget>
    std::vector<std::tuple<std::string, std::string, int32_t, int32_t, bool, std::string, PacBio::Pancake::TrimmingInfo>> testData {
        {"Empty input", "", 0, 0, true, "", {0, 0, 0, 0}},

        {"Completely bad alignment.", "1000X", 30, 15, true, "", {0, 0, 0, 0}},
        {"Remove leading errors.", "1I1000=", 30, 15, true, "1000=", {1, 0, 0, 0}},
        {"Remove trailing errors, trim to match.", "1000=1X", 30, 15, true, "1000=", {0, 0, 1, 1}},
        {"Remove trailing errors, trim to any op.", "1000=1X", 30, 15, false, "1000=1X", {0, 0, 0, 0}},
        {"Remove leading errors, trim to match.", "1X1000=", 30, 15, true, "1000=", {1, 1, 0, 0}},
        {"Remove leading errors, trim to any op.", "1X1000=", 30, 15, false, "1X1000=", {0, 0, 0, 0}},
        {"Remove leading and trailing errors.", "1I1000=1X", 30, 15, true, "1000=", {1, 0, 1, 1}},

        {"Remove trailing errors with mixed matches, trim to match.", "1000=1X1=1X", 30, 15, true, "1000=1X1=", {0, 0, 1, 1}},
        {"Remove trailing errors with mixed matches, trim to any op.", "1000=1X1=1X", 30, 15, false, "1000=1X1=1X", {0, 0, 0, 0}},

        {"Large match portion followed with a trailing 1D. Trim to first match.", "5026=1D", 30, 15, true, "5026=", {0, 0, 0, 1}},
        {"Large match portion followed with a trailing 1D. Trim to non-match ops allowed.", "5026=1D", 30, 15, false, "5026=1D", {0, 0, 0, 0}},
        {"Large match portion prefixed with a leading 1D. Trim to first match.", "1D5026=", 30, 15, true, "5026=", {0, 1, 0, 0}},
        {"Large match portion prefixed with a leading 1D. Trim to non-match ops allowed.", "1D5026=", 30, 15, false, "1D5026=", {0, 0, 0, 0}},

        {"Large match portion followed with a trailing 1I. Trim to first match.", "5026=1I", 30, 15, true, "5026=", {0, 0, 1, 0}},
        {"Large match portion followed with a trailing 1I. Trim to non-match ops allowed.", "5026=1I", 30, 15, false, "5026=1I", {0, 0, 0, 0}},
        {"Large match portion prefixed with a leading 1I. Trim to first match.", "1I5026=", 30, 15, true, "5026=", {1, 0, 0, 0}},
        {"Large match portion prefixed with a leading 1I. Trim to non-match ops allowed.", "1I5026=", 30, 15, false, "1I5026=", {0, 0, 0, 0}},

        /*
            These test clipping on the right side.
        */
        {"Left clipping in a small window, 5bp with required 5bp matches, trim to first match op.",
                    "5X100=1I100=", 5, 5, true, "100=1I100=", {5, 5, 0, 0}},
        {"Left clipping in a small window, 5bp with required 5bp matches, trim to non-match ops allowed.",
                    "5X100=1I100=", 5, 5, false, "100=1I100=", {5, 5, 0, 0}},
        ///
        {"Left clipping in a small window, 5bp with required 3bp matches, trim to first match op.",
                    "5X100=1I100=", 5, 3, true, "100=1I100=", {5, 5, 0, 0}},
        {"Left clipping in a small window, 5bp with required 3bp matches, trim to non-match ops allowed.",
                    "5X100=1I100=", 5, 3, false, "2X100=1I100=", {3, 3, 0, 0}},
        ///
        {"Left clipping in a larger window, 30bp with required 15bp matches, trim to first match op.",
                    "5X100=1I100=", 30, 15, true, "100=1I100=", {5, 5, 0, 0}},
        {"Left clipping in a larger window, 30bp with required 15bp matches, trim to non-match ops allowed. Window is large enough so that the leading 5X is not clipped.",
                    "5X100=1I100=", 30, 15, false, "5X100=1I100=", {0, 0, 0, 0}},
        ///
        {"Left clipping. Leading diff is larger than window size, trim to first match op.",
                    "1000X100=1I100=", 30, 15, true, "100=1I100=", {1000, 1000, 0, 0}},
        {"Left clipping, simple. Clip larger than window size, trim to non-match ops allowed.",
                    "1000X100=1I100=", 30, 15, false, "15X100=1I100=", {985, 985, 0, 0}},
        ///
        {"Left clipping of insertions, trim to first match op.",
                    "100I100=1I100=", 30, 15, true, "100=1I100=", {100, 0, 0, 0}},
        {"Left clipping of insertions, trim to non-match ops allowed.",
                    "100I100=1I100=", 30, 15, false, "15I100=1I100=", {85, 0, 0, 0}},
        ///
        {"Left clipping of deletions, trim to first match op.",
                    "100D100=1I100=", 30, 15, true, "100=1I100=", {0, 100, 0, 0}},
        {"Left clipping of deletions, trim to non-match ops allowed.",
                    "100D100=1I100=", 30, 15, false, "15D100=1I100=", {0, 85, 0, 0}},

        /*
            These are completely symmetric to the left-clipping test case.
        */
        {"Right clipping in a small window, 5bp with required 5bp matches, trim to first match op.",
                    "100=1I100=5X", 5, 5, true, "100=1I100=", {0, 0, 5, 5}},
        {"Right clipping in a small window, 5bp with required 5bp matches, trim to non-match ops allowed.",
                    "100=1I100=5X", 5, 5, false, "100=1I100=", {0, 0, 5, 5}},
        ///
        {"Right clipping in a small window, 5bp with required 3bp matches, trim to first match op.",
                    "100=1I100=5X", 5, 3, true, "100=1I100=", {0, 0, 5, 5}},
        {"Right clipping in a small window, 5bp with required 3bp matches, trim to non-match ops allowed.",
                    "100=1I100=5X", 5, 3, false, "100=1I100=2X", {0, 0, 3, 3}},
        ///
        {"Right clipping in a larger window, 30bp with required 15bp matches, trim to first match op.",
                    "100=1I100=5X", 30, 15, true, "100=1I100=", {0, 0, 5, 5}},
        {"Right clipping in a larger window, 30bp with required 15bp matches, trim to non-match ops allowed. Window is large enough so that the leading 5X is not clipped.",
                    "100=1I100=5X", 30, 15, false, "100=1I100=5X", {0, 0, 0, 0}},
        ///
        {"Right clipping. Leading diff is larger than window size, trim to first match op.",
                    "100=1I100=1000X", 30, 15, true, "100=1I100=", {0, 0, 1000, 1000}},
        {"Right clipping, simple. Clip larger than window size, trim to non-match ops allowed.",
                    "100=1I100=1000X", 30, 15, false, "100=1I100=15X", {0, 0, 985, 985}},
        ///
        {"Right clipping of insertions, trim to first match op.",
                    "100=1I100=100I", 30, 15, true, "100=1I100=", {0, 0, 100, 0}},
        {"Right clipping of insertions, trim to non-match ops allowed.",
                    "100=1I100=100I", 30, 15, false, "100=1I100=15I", {0, 0, 85, 0}},
        ///
        {"Right clipping of deletions, trim to first match op.",
                    "100=1I100=100D", 30, 15, true, "100=1I100=", {0, 0, 0, 100}},
        {"Right clipping of deletions, trim to non-match ops allowed.",
                    "100=1I100=100D", 30, 15, false, "100=1I100=15D", {0, 0, 0, 85}},

        /* Following test runs a more complex example.
           Maches:     1     3     4        5     7     12    14    15    16
           Diffs:   1     2     4     5  6     8     9     11    13    15    18
                    1X 1= 1D 2= 2D 1= 1X 1D 1= 2I 2= 1X 5= 2D 1= 2I 1= 2X 1= 3I

                    0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2
                    0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 3 3 3
                    X = D = = D D = X D = I I = = X = = = = = D D = I I = X X = I I I
                    |___________________________| | | | |
                    Window: 7= and 8!=            | | | |
                      |___________________________| | | |
                      Window: 7= and 8!=            | | |
                        |___________________________| | |
                        Window: 7= and 8!=            | |
                          |___________________________| |
                          Window: 8= and 7!=            |
                            |___________________________|
                            Window: 8= and 7!=

        */
        {"Left clipping, simple, mixed ops.",
                    "1X1=1D2=2D1=1X1D1=2I2=1X5=2D1=2I1=2X1=3I100=1I100=", 15, 8, true, "2=2D1=1X1D1=2I2=1X5=2D1=2I1=2X1=3I100=1I100=", {2, 3, 0, 0}},
        // Same as the left clipping case, but flipped
        {"Right clipping, simple, mixed ops.",
                    "100=1I100=3I1=2X1=2I1=2D5=1X2=2I1=1D1X1=2D2=1D1=1X", 15, 8, true, "100=1I100=3I1=2X1=2I1=2D5=1X2=2I1=1D1X1=2D2=", {0, 0, 2, 3}},

        /*
            This test case contains a bunch of diffs at the front and a bunch of diffs at the back.
            The diffs at front and back are the same and symmetric.
            The actual diffs are the same as in the above test cases.
        */
        {"Clipping on front and back, mixed ops. Smaller window - 15bp with 8bp matches required.",
                    "1X1=1D2=2D1=1X1D1=2I2=1X5=2D1=2I1=2X1=3I100=1I100=1I100=3I1=2X1=2I1=2D5=1X2=2I1=1D1X1=2D2=1D1=1X", 15, 8, true, "2=2D1=1X1D1=2I2=1X5=2D1=2I1=2X1=3I100=1I100=1I100=3I1=2X1=2I1=2D5=1X2=2I1=1D1X1=2D2=", {2, 3, 2, 3}},
        {"Clipping on front and back, mixed ops. Larger window - 80bp with 80bp matches required.",
                    "1X1=1D2=2D1=1X1D1=2I2=1X5=2D1=2I1=2X1=3I100=1I100=1I100=3I1=2X1=2I1=2D5=1X2=2I1=1D1X1=2D2=1D1=1X", 80, 80, true, "100=1I100=1I100=", {27, 26, 27, 26}},
    };
    // clang-format on

    for (const auto& data : testData) {
        // Inputs.
        const std::string testName = std::get<0>(data);
        const PacBio::BAM::Cigar cigar(std::get<1>(data));
        const int32_t windowSize = std::get<2>(data);
        const int32_t minMatches = std::get<3>(data);
        const bool clipOnFirstMatch = std::get<4>(data);
        // Keep the expected Cigar as string so that we can easily see the diff.
        const std::string expectedTrimmedCigar(std::get<5>(data));
        const auto& expectedTrimmingInfo = std::get<6>(data);

        // Name the test.
        SCOPED_TRACE("TrimCigar-" + testName);

        PacBio::BAM::Cigar resultsCigar;
        TrimmingInfo resultsTrimming;

        PacBio::Pancake::TrimCigar(cigar, windowSize, minMatches, clipOnFirstMatch, resultsCigar,
                                   resultsTrimming);

        EXPECT_EQ(expectedTrimmedCigar, resultsCigar.ToStdString());
        EXPECT_EQ(expectedTrimmingInfo, resultsTrimming);
    }
}

TEST(Test_AlignmentTools_TrimCigar, ArrayOfTests_ShouldThrow)
{
    // clang-format off
    std::vector<std::tuple<std::string, std::string, int32_t, int32_t, bool, bool, std::string, PacBio::Pancake::TrimmingInfo>> testData {
        {"Window size is too large, 1000bp and max is 512bp.", "1I1000=", 1000, 15, true, true, "1000=", {1, 0, 0, 0}},
    };
    // clang-format on

    for (const auto& data : testData) {
        // Inputs.
        const std::string testName = std::get<0>(data);
        const PacBio::BAM::Cigar cigar(std::get<1>(data));
        const int32_t windowSize = std::get<2>(data);
        const int32_t minMatches = std::get<3>(data);
        const bool clipOnFirstMatch = std::get<4>(data);
        // Keep the expected Cigar as string so that we can easily see the diff.
        const bool exptectedThrow = std::get<5>(data);
        const std::string expectedTrimmedCigar(std::get<6>(data));
        const auto& expectedTrimmingInfo = std::get<7>(data);

        // Name the test.
        SCOPED_TRACE("TrimCigar-" + testName);

        PacBio::BAM::Cigar resultsCigar;
        TrimmingInfo resultsTrimming;

        if (exptectedThrow) {
            EXPECT_THROW(
                {
                    PacBio::Pancake::TrimCigar(cigar, windowSize, minMatches, clipOnFirstMatch,
                                               resultsCigar, resultsTrimming);
                },
                std::runtime_error);

        } else {
            PacBio::Pancake::TrimCigar(cigar, windowSize, minMatches, clipOnFirstMatch,
                                       resultsCigar, resultsTrimming);
            EXPECT_EQ(expectedTrimmedCigar, resultsCigar.ToStdString());
            EXPECT_EQ(expectedTrimmingInfo, resultsTrimming);
        }
    }
}

TEST(Test_AlignmentTools_ScoreCigarAlignment, ArrayOfTests)
{
    // clang-format off
    struct TestDataStruct {
        std::string name;
        std::string cigar;
        int32_t match = 1;
        int32_t mismatch = 1;
        int32_t gapOpen = 1;
        int32_t gapExtend = 1;
        int32_t expectedScore = 0;
    };
    std::vector<TestDataStruct> testData = {
        {"Empty input", "", 8, 4, 4, 2, 0},
        {"Single match event", "1000=", 8, 4, 4, 2, 8000},
        {"Single mismatch event", "1000X", 8, 4, 4, 2, -4000},
        {"Single insertion event", "1000I", 8, 4, 4, 2, -2002},
        {"Single deletion event", "1000D", 8, 4, 4, 2, -2002},
        {"Simple mixed CIGAR", "8=2X3=1D2=4I8=", 8, 4, 4, 2, 8*8 - 2*4 + 3*8 - 1*4 + 2*8 - (4 + 3*2) + 8*8},
    };
    // clang-format on

    for (const auto& data : testData) {
        // Inputs.
        const PacBio::BAM::Cigar cigar(data.cigar);

        // Name the test.
        SCOPED_TRACE("ScoreCigarAlignment-" + data.name);

        const int32_t result = PacBio::Pancake::ScoreCigarAlignment(
            cigar, data.match, data.mismatch, data.gapOpen, data.gapExtend);
        EXPECT_EQ(data.expectedScore, result);
    }
}

TEST(Test_AlignmentTools_MergeCigars, ArrayOfTests)
{
    // clang-format off
    struct TestDataStruct {
        std::string name;
        std::string destCigar;
        std::string sourceCigar;
        std::string expected;
    };
    std::vector<TestDataStruct> testData = {
        {"Empty input", "", "", ""},
        {"Empty destination", "", "5=1X3I", "5=1X3I"},
        {"Empty source", "5=1X3I", "", "5=1X3I"},
        {"Nonempty source and destination, no events to merge", "1X", "5=1X3I", "1X5=1X3I"},
        {"Nonempty source and destination with merging of same events", "3=1X3=", "5=1X3I", "3=1X8=1X3I"},
    };
    // clang-format on

    for (const auto& data : testData) {
        // Inputs.
        PacBio::BAM::Cigar destCigar(data.destCigar);
        const PacBio::BAM::Cigar sourceCigar(data.sourceCigar);

        // Name the test.
        SCOPED_TRACE("MergeCigars-" + data.name);

        // Run the unit under test.
        PacBio::Pancake::MergeCigars(destCigar, sourceCigar);

        // Convert the result for comparison.
        const std::string result = destCigar.ToStdString();

        EXPECT_EQ(data.expected, result);
    }
}

TEST(Test_AlignmentTools_ComputeSimpleRepeatMask, ArrayOfTests)
{
    // clang-format off
    struct TestDataStruct {
        std::string name;
        std::string seq;
        int32_t maxSpan = 0;
        std::vector<uint8_t> expected;
        bool expectedThrow = false;
    };
    std::vector<TestDataStruct> testData = {
        {"Empty input", "", 0, {}, false},
        {"Simple sequence, no masking", "ACTGAAACACTG", 0, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, false},
        {"Simple sequence, HP masking", "ACTGAAACACTG", 1, {0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0}, false},
        {"Simple sequence, dinuc masking", "ACTGAAACACTG", 2, {0, 0, 0, 0, 1, 1, 3, 2, 2, 2, 0, 0}, false},
        {"Simple sequence, trinuc masking",
                    "ACTG" "AA" "ACAC" "TG" "ATGATGATG" "ACTGACTG", 3,
                    /*
                        Seq:        ACTGAAACACTGATGATGATGACTGACTG
                        No mask:    0000  |  |         | |      |
                        HP:             111  |         | |      |
                        Dinuc:            2222         | |      |
                        Trinuc:               444444444444      |
                        Quadnuc:                       8888888888
                    */
                    {
                        0, 0, 0, 0,                     // ACTG
                        1, 1, 3, 2, 2, 2, 4, 4,         // AAACACTG
                        4, 4, 4, 4, 4, 4, 4, 4, 4,      // ATGATGATG
                        4, 0, 0, 0, 0, 0, 0, 0,         // ACTGACTG
                    },
                    false
        },
        {"Simple sequence, quadnuc masking",
                    "ACTG" "AA" "ACAC" "TG" "ATGATGATG" "ACTGACTG", 4,
                    /*
                        Seq:        ACTGAAACACTGATGATGATGACTGACTG
                        No mask:    0000  |  |         | |      |
                        HP:             111  |         | |      |
                        Dinuc:            2222         | |      |
                        Trinuc:               444444444444      |
                        Quadnuc:                       8888888888
                    */
                    {
                        0, 0, 0, 0,                     // ACTG
                        1, 1, 3, 2, 2, 2, 4, 4,         // AAACACTG
                        4, 4, 4, 4, 4, 4, 4, 12, 12,    // ATGATGATG
                        12, 8, 8, 8, 8, 8, 8, 8,        // ACTGACTG
                    },
                    false
        },
        {"Simple sequence, level 7 masking",
                    "ACTG" "AA" "ACAC" "TG" "ATGATGATG" "ACTGACTG", 7,
                    {
                        0, 0, 0, 0,                             // ACTG
                        1, 1, 3, 2, 2, 2, 36, 36,               // AAACACTG
                        36, 36, 36, 36, 36, 36, 36, 44, 44,     // ATGATGATG
                        44, 8, 8, 8, 8, 8, 8, 8,                // ACTGACTG
                    },
                    false
        },
        {"Too large max span, return everything masked at level 7",
                    "ACTG" "AA" "ACAC" "TG" "ATGATGATG" "ACTGACTG", 8,
                    {
                        0, 0, 0, 0,                             // ACTG
                        1, 1, 3, 2, 2, 2, 36, 36,               // AAACACTG
                        36, 36, 36, 36, 36, 36, 36, 44, 44,     // ATGATGATG
                        44, 8, 8, 8, 8, 8, 8, 8,                // ACTGACTG
                    },
                    false
        },
        {"Two homopolymers", "AAAAAAAACCCCCCC", 1, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, false},
        {"Long HPs are also dinucs and trinucs", "AAAAAAAACCCCCCC", 3, {7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7}, false},
        {"Regular sequence with HPs", "GGATCAGTTTTATATACAC", 3, {1, 1, 0, 0, 0, 0, 0, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2}, false},
        {"Non-ACTG base, simple sequence, HP masking", "ACTGANAACACTG", 1, {0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0}, false},
        {"Multiple non-ACTG bases, simple sequence, HP masking", "ACTGANNNANACACTG", 1, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, false},

    };
    // clang-format on

    for (const auto& data : testData) {

        // Name the test.
        SCOPED_TRACE("ComputeSimpleRepeatMask-" + data.name);

        // std::cerr << "Test: " << data.name << "\n";

        // Run.
        if (data.expectedThrow) {
            EXPECT_THROW(
                {
                    const std::vector<uint8_t> results =
                        ComputeSimpleRepeatMask(data.seq.c_str(), data.seq.size(), data.maxSpan);
                },
                std::runtime_error);
        } else {

            const std::vector<uint8_t> results =
                ComputeSimpleRepeatMask(data.seq.c_str(), data.seq.size(), data.maxSpan);

            // std::cerr << "Len: seq.size() = " << data.seq.size()
            //           << ", results.size() = " << results.size()
            //           << ", expected.size() = " << data.expected.size() << "\n";
            // std::cerr << "Seq: ";
            // for (const auto& v : data.seq) {
            //     std::cerr << v << " ";
            // }
            // std::cerr << "\n";
            // std::cerr << "Res: ";
            // for (const auto& v : results) {
            //     std::cerr << static_cast<uint32_t>(v) << " ";
            // }
            // std::cerr << "\n";
            // std::cerr << "Exp: ";
            // for (const auto& v : data.expected) {
            //     std::cerr << static_cast<uint32_t>(v) << " ";
            // }
            // std::cerr << "\n";

            EXPECT_EQ(data.expected, results);
        }
        // std::cerr << "\n";
    }
}

TEST(Test_AlignmentTools_CheckAlignmentOutOfBand, ArrayOfTests)
{
    // clang-format off
    struct TestDataStruct {
        const std::string name;
        const std::string cigar;
        const int32_t bandwidth;
        const bool expected;
    };
    std::vector<TestDataStruct> testData = {
        // Within band.
        {"Empty input", "", 0, false},
        {"No indels, within band", "100=", 10, false},
        {"Insertions within band", "50=5I50=", 10, false},
        {"Deletions within band", "50=5D50=", 10, false},
        {"Indels within band", "50=5I5D50=", 10, false},
        {"No indels, marginal to bandwidth 0", "100=", 0, false},
        {"No indels, marginal to bandwidth 1", "100=", 1, false},
        {"No indels, within band", "100=", -1, false},

        // Out of band.
        {"Single insertion", "100=1I", 0, true},
        {"Single deletion", "100=1D", 0, true},
        {"Single insertion, bw 1", "100=1I", 1, true},
        {"Single insertion, bw -1", "100=1I", -1, true},
        {"Single deletion, bw 1", "100=1D", 1, true},
        {"Single deletion, bw -1", "100=1D", -1, true},
        {"Insertions, touches the boundary", "50=10I50=", 10, true},
        {"Insertions, touches the diagonal -1 from the boundary (fuzz buffer)", "50=9I50=", 10, true},
        {"Deletions, touches the boundary", "50=10D50=", 10, true},
        {"Deletions, touches the diagonal -1 from the boundary (fuzz buffer)", "50=9D50=", 10, true},
        {"Indels, out of band", "50=10I10D50=", 10, true},
        {"Indels, out of band", "50=5I2D10=3I1D1=5I50=", 10, true},
        {"Indels, out of band", "50=10D10I50=", 10, true},
        {"Indels, out of band", "50=5D2I10=3D1I1=5D50=", 10, true},
    };
    // clang-format on

    for (const auto& data : testData) {
        // Inputs.
        const PacBio::BAM::Cigar cigar(data.cigar);

        // Name the test.
        SCOPED_TRACE("CheckAlignmentOutOfBand-" + data.name);

        // Run the unit under test.
        const bool result = PacBio::Pancake::CheckAlignmentOutOfBand(cigar, data.bandwidth);

        // Evaluate.
        EXPECT_EQ(data.expected, result);
    }
}
}
