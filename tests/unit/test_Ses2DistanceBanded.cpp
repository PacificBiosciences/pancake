// Authors: Ivan Sovic

#include <PancakeTestData.h>
#include <gtest/gtest.h>
#include <pancake/Ses2DistanceBanded.hpp>
#include <pancake/SesAlignBanded.hpp>
#include <sstream>
#include <tuple>

namespace PacBio {
namespace Pancake {
namespace Alignment {
namespace Test {

struct TestData
{
    std::string testName;
    std::string query;
    std::string target;
    int32_t maxDiffs = 0;
    int32_t bandwidth = 0;
    SesResults expected = {0, 0, 0, false};
};

std::vector<TestData> testDataSemiglobal = {
    TestData{"EmptyQueryEmptyTarget", "", "", 100, 30, SesResults(0, 0, 0, true)},
    TestData{"EmptyQueryNonemptyTarget", "", "ACTG", 100, 30, SesResults(0, 0, 0, true)},
    TestData{"NonEmptyQueryEmptyTarget", "ACTG", "", 100, 30, SesResults(0, 0, 0, true)},
    TestData{"SimpleExactMatch", "ACTG", "ACTG", 100, 30, SesResults(4, 4, 0, true)},

    TestData{"SimpleSingleIndelDiff", "ACG", "ACTG", 100, 30, SesResults(3, 4, 1, true)},
    // TestData{"NormalSmallCase_MultipleMismatches", "CCCC", "GGGG", 100, 30, SesResults(0, 4, 4, true)},

    // It could also be (5, 5, 1, true), depends on the preferences of the aligner.
    TestData{"SimpleSingleMismatchDiff", "AAAAA", "AAATA", 100, 30, SesResults(4, 5, 1, true)},

    // Front matches.
    TestData{"TrailingIndelInTarget", "ACTG", "ACTGGGG", 100, 30, SesResults(4, 4, 0, true)},

    // In this case, the opening gap (k = -3) is processed before the diagonal (k = 0) and both ave the same
    // number of diffs. Since smaller k is processed first, we reach a diagonal we can traverse, and zip down to the end.
    TestData{"FrontIndelInTarget", "ACTG", "GGGACTG", 100, 30, SesResults(4, 7, 3, true)},

    // Front matches.
    TestData{"TrailingIndelInQuery", "ACTGGGG", "ACTG", 100, 30, SesResults(4, 4, 0, true)},

    /*
     *  ---GGGACTG
     *     |
     *  ACTG------
    */
    TestData{"FrontIndelInQuery", "GGGACTG", "ACTG", 100, 30, SesResults(1, 4, 3, true)},

    TestData{"SimpleFiveBaseDeletion", "AAAAAAAAAAAAAAAAAAAA", "AAAAAAGGGGGAAAAAAAAAAAAAA", 100, 30,
             SesResults(20, 25, 5, true)},

    /*
     * Note: Even though a global aligner should correctly demarcate the
     * insertion of 5 bases ("GGGGG"), a semiglobal Ses2 aligner here would
     * prefer to make a mismatch call for the last G base which would allow it
     * to stat going down the diagonal faster than it would have if it had
     * waited for a deletion.
     * That's why the X coordinate ends at 24 instead of 25.
     * Should the global aligner be used, X would be 25.
     *
     * Also note that the original SES aligner from 1986 would not result in the
     * same way because it's unaware of mismatches.
    */
    TestData{"SimpleFiveBaseInsertion", "AAAAAAGGGGGAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAA", 100,
             30, SesResults(20, 20, 5, true)},

    /*
     * Note: Semiglobal aligner would produce the alignment on the right in the
     * diagram below. Alignment stops because the target sequence hits the end,
     * and this happens before the last insertion.
     * That's why only 2 diffs are produced.

     * A global aligner would also include the trailing indels.
     *
     * GA-T-GTTT      GA-T-GTTT
     * || | | ||      || | |||
     * GAATTG-TT      GAATTGTT-
    */
    TestData{"SimpleIndelsOnly", "GATGTTT", "GAATTGTT", 100, 30, SesResults(6, 8, 2, true)},

    /*
     * Note: Mismatches count as 1 diff, unlike the 1986 algorithm.
     * The original O(nd) algorithm would find an alignment like this:
     *      G-GAT-CAGTT
     *      |  || | |||
     *      GA-ATTC-GTT
     * In reality, minimum edit distance would have only 3 mismatches:
     *      GGATCAGTT
     *      |X||XX|||
     *      GAATTCGTT
    */
    TestData{"NormalSmallCase", "GGATCAGTT", "GAATTCGTT", 100, 30, SesResults(9, 9, 3, true)},

    /*
     * The maximum number of diffs is within the allowed limit, but
     * the maximum bandwidth for alignment is below the threshold (3).
     * The output should be marked as "invalid" but also report the
     * furthest reaching (f.r.) wave.
     *
     * Note: Mismatches count as 1 diff, unlike the 1986 algorithm.
     * The original O(nd) algorithm would find an alignment like this:
     *      G-GAT-CAGTT
     *      |  || | |||
     *      GA-ATTC-GTT
     * In reality, minimum edit distance would have only 3 mismatches:
     *      GGATCAGTT
     *      |X||XX|||
     *      GAATTCGTT
    */
    TestData{"NarrowBandwidthButOnlyMismatches", "GGATCAGTT", "GAATTCGTT", 100, 2,
             SesResults(5, 6, 3, false)},

    /*
     * The maximum number of diffs is within the allowed limit, but
     * the maximum bandwidth for alignment is below the threshold (3).
     * The output should be marked as "invalid" but also report the
     * furthest reaching (f.r.) wave.
    */
    TestData{"OutOfBandwidth", "GATGTTT", "GAATTGTT", 100, 2, SesResults(3, 4, 2, false)},

    /*
     * The maximum number of diffs is too small, so the alignment
     * should be marked as not valid.
     * the maximum bandwidth for alignment is GOOD here (30 > 3).
     * The output should be marked as "invalid" but also report the
     * furthest reaching (f.r.) wave.
     *
     * The furthest reaching point is (5, 6) which includes
     * 2 mismatches and an deletion.
     *
     * Mismatches count as 1 diff, unlike the 1986 algorithm.
     * The original O(nd) algorithm would find an alignment like this:
     *      G-GAT-CAGTT
     *      |  || | |||
     *      GA-ATTC-GTT
     * In reality, minimum edit distance would have only 3 mismatches:
     *      GGATCAGTT
     *      |X||XX|||
     *      GAATTCGTT
    */
    TestData{"AboveMaxDiffs", "GGATCAGTT", "GAATTCGTT", 3, 30, SesResults(5, 6, 2, false)},
};

std::vector<TestData> testDataGlobal = {
    TestData{"EmptyQueryEmptyTarget", "", "", 100, 30, SesResults(0, 0, 0, true)},
    TestData{"EmptyQueryNonemptyTarget", "", "ACTG", 100, 30, SesResults(0, 0, 0, true)},
    TestData{"NonEmptyQueryEmptyTarget", "ACTG", "", 100, 30, SesResults(0, 0, 0, true)},
    TestData{"SimpleExactMatch", "ACTG", "ACTG", 100, 30, SesResults(4, 4, 0, true)},
    TestData{"SimpleSingleIndelDiff", "ACG", "ACTG", 100, 30, SesResults(3, 4, 1, true)},
    TestData{"SimpleSingleMismatchDiff", "AAAAA", "AAATA", 100, 30, SesResults(5, 5, 1, true)},
    TestData{"NormalSmallCase_MultipleMismatches", "CCCC", "GGGG", 100, 30,
             SesResults(4, 4, 4, true)},
    TestData{"TrailingIndelInTarget", "ACTG", "ACTGGGG", 100, 30, SesResults(4, 7, 3, true)},
    TestData{"FrontIndelInTarget", "ACTG", "GGGACTG", 100, 30, SesResults(4, 7, 3, true)},
    TestData{"TrailingIndelInQuery", "ACTGGGG", "ACTG", 100, 30, SesResults(7, 4, 3, true)},
    TestData{"FrontIndelInQuery", "GGGACTG", "ACTG", 100, 30, SesResults(7, 4, 3, true)},
    TestData{"SimpleFiveBaseDeletion", "AAAAAAAAAAAAAAAAAAAA", "AAAAAAGGGGGAAAAAAAAAAAAAA", 100, 30,
             SesResults(20, 25, 5, true)},

    /*
     * Global aligner will correctly align the gap and align all the way to the end
     * of both sequences.
    */
    TestData{"SimpleFiveBaseInsertion", "AAAAAAGGGGGAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAA", 100,
             30, SesResults(25, 20, 5, true)},

    /*
     * Subtle difference between this case and the semiglobal case for same input data is
     * that this produces 3 diffs instead of 2, because trailing indels are also computed.
     * Also, the endpoint of the alignment is at the end of both sequences.
     * (Alternatively, if another alignment with the same edit distance was actually
     * internally computed, then the indels would be internal. Check out the diagram below.)
     *
     * GA-T-GTTT      GA-T-GTTT
     * || | | ||      || | |||
     * GAATTG-TT      GAATTGTT-
    */
    TestData{"SimpleIndelsOnly", "GATGTTT", "GAATTGTT", 100, 30, SesResults(7, 8, 3, true)},

    /*
     * Note: Mismatches count as 1 diff, unlike the 1986 algorithm.
     * The original O(nd) algorithm would find an alignment like this:
     *      G-GAT-CAGTT
     *      |  || | |||
     *      GA-ATTC-GTT
     * In reality, minimum edit distance would have only 3 mismatches:
     *      GGATCAGTT
     *      |X||XX|||
     *      GAATTCGTT
    */
    TestData{"NormalSmallCase", "GGATCAGTT", "GAATTCGTT", 100, 30, SesResults(9, 9, 3, true)},

    /*
     * The maximum number of diffs is within the allowed limit, but
     * the maximum bandwidth for alignment is below the threshold (3).
     * The output should be marked as "invalid" but also report the
     * furthest reaching (f.r.) wave.
     *
     * Note: Mismatches count as 1 diff, unlike the 1986 algorithm.
     * The original O(nd) algorithm would find an alignment like this:
     *      G-GAT-CAGTT
     *      |  || | |||
     *      GA-ATTC-GTT
     * In reality, minimum edit distance would have only 3 mismatches:
     *      GGATCAGTT
     *      |X||XX|||
     *      GAATTCGTT
    */
    TestData{"NarrowBandwidthButOnlyMismatches", "GGATCAGTT", "GAATTCGTT", 100, 2,
             SesResults(5, 6, 3, false)},

    /*
     * The maximum number of diffs is within the allowed limit, but
     * the maximum bandwidth for alignment is below the threshold (3).
     * The output should be marked as "invalid" but also report the
     * furthest reaching (f.r.) wave.
    */
    TestData{"OutOfBandwidth", "GATGTTT", "GAATTGTT", 100, 2, SesResults(3, 4, 2, false)},

    /*
     * The maximum number of diffs is too small, so the alignment
     * should be marked as not valid.
     * the maximum bandwidth for alignment is GOOD here (30 > 3).
     * The output should be marked as "invalid" but also report the
     * furthest reaching (f.r.) wave.
     *
     * The furthest reaching point is (5, 6) which includes
     * 2 mismatches and an deletion.
     *
     * Mismatches count as 1 diff, unlike the 1986 algorithm.
     * The original O(nd) algorithm would find an alignment like this:
     *      G-GAT-CAGTT
     *      |  || | |||
     *      GA-ATTC-GTT
     * In reality, minimum edit distance would have only 3 mismatches:
     *      GGATCAGTT
     *      |X||XX|||
     *      GAATTCGTT
    */
    TestData{"AboveMaxDiffs", "GGATCAGTT", "GAATTCGTT", 3, 30, SesResults(5, 6, 2, false)},

    /*
     * This is an actual real set of sequences. Only a subportion which aligns end-to-end
     * was extracted for the unit test here, so that the global aligner can be tested without mapping.
    */
    // clang-format off
    TestData{"RealHifiReads_Subsequences",
                // >mo_m64001_190914_015449/130352465/ccs:8054-16191
                "TATGAAGAGCATGTGAGCTAGTCCACTCGAGCAACTGAAGACTACACAAAAAGCCTCTTAGTCCGAGTATTTTGGGAATTTTTCACTCAAAATTAATTGCATTTTAATAGTTGCTAATATTTCAAAAAACTTGTGCACTTATTTGCACCGGCATCAGGCTGCTCAGGCTCGCCGCTTTCCCCATCCCCTTCTCCACGAAGCTCCATCAAAAAGTCCTGCTTAAAAATGCAATAACAGGAAATGGCTACGACGGCAACAAAATACAGGAAATTGAGCGCAAAGAAATATGCTATATAGGAGATGAAAAAATGAAAACGCAAACGAGACGAAGTCCAGGCAAACAAACATGGCCGTCAGACATGAACATATGTATACAGCCGCGGTTAACCAATTAGGAGTTTTGAATTCCAGAAAAAAGAATTAATACAAATATGCACTTTATTTACTAAGATTCATTTACATTTGGATTCGTTGGTTGGTGTCAGCAAGCGTACTTGAGCTTTCAGGCAGAAAAGTGCTCATGCTTTAGCATACTGAGTCTTAGAATAGCAAAGTAGCTACAAAAAATTTAAAATCAAAGTTAAGTTAAAAGTTAAAATCAAAGTTAAAAAAAAAAAACTCGAATGAACGTTTTGAAAAAGGGGAAAAGGGAAACCCCATGGAAAGCAACACAACTATTACAATTACTTTCTGTTATATGAGATTAGTATTTTTTTGAAGACACTAAACATGTGGCCCCTGCTGTACACCTGAACCATGCAAAGCGAGCATGAGGCGAAGACAGGTGGAAAGGTGGCGTTACGGGTAAGAGCATTTTGCCCGCTGGAACTTGGCCGTATGTGCTCCGGCAATTTATGTCCAAGGTGCGTGCCACCGAAATGAACCGGTGAACCTGGACGTTACGCTCGTTTCTTGGGGGAAAGTCAGGCAGTTGACCAACCTGTGTCCTTTTTCCTGGTTCTGCGCTTAACGTGTTGCCATTAAATTTTATGGTCAGGTCTTCCTCTCTTTGTCCGTGTTTGCATTTTGCTCAATTTTGTCGCCGCAGAACATAAAAAAGGGAGTCAGGAGGGGCTGAAAAAGAAGCGAAGGAATATTGCACTGTGGACGACATGTGAAAACATATTTTTATGGGACCCCAAAGGGTAATTGGACATTGCACGGGATATCCTAGCTGGCTAATGAGTGTGTCGAGCCTTTTTAGATAATTAATTGAGAGCTTTTCAAATCATCTGTGAAGTATAATGACCACACACACACAATAAGGCAACGAAACGTTGAAGTGTGATTGCACTAATGAGTAGATTAAAGTCAGCTGGAGATAACATCTTAATAAATCCAAGATAATTGAAATATCCTTTCGAGATAGTTGTATAATTTGTAATTCTCCATTCGGGAATATTTGTTTCATTGGCTTCAAGTGAAACCAATTGTGCCGCACTTCATCCAGCGAAACTTTGACGGCTAACTCGCTTAATGCGCGGCTTAAGCCAAACTACGAGTAAGATATTTTGTCTGTGGCAACTTGAAAGTTGCAAGTAATTAGTGGCACTTGAGCAACTACTAGGCTGTCGGGATACGAAAATGGGGCAGGGCATACCCAGATAGCTGGATGCTCCAAAAGTTATGTAGATCACGAAAATCAGTCCTTTGAGAACAAAGAACTTTTAGCAATCTTAAAATATTATTTCTAAGTATTCTAATCGTTATACTTTATTAGTATTTTAGTTCTATGAGTTTCTTAAAAAAAAAAAAACACTGATCTTACAAGCTCGCTTTCCACCTATCCTCTTTCAGACCCATCAGTCTAGTGCTGGAGCCAATCCTGAAGGAGCCCGGCGACATCTCCTTTACGTGGCTCGAGATTCACCGCACCAAGAACGCGCTGCTGCAGCAGTTGGACCTGCATGTGAACGCTAGCGTGACCACGGGCATGGGCTCCCAGACCTCCTCCCAGTCCGCAGCCGGCTCCAACACGGGCATTGGGGGACTGGGCACCCAGTGGGTCAGGATCTGCTGCCAATGCCGGGGGCTCTGCTGGCGGAGCAACGGCCAACGAGACGCTGAACGAGTTCGGATTCTTTCCCAAGGAGTCCGACTGCGAATACAAATGTCCTGAAATTGATGCATGCATTGCGGCCAGCTTGTGGTGCGATGGTGAGTACAGATTGCCAGTCGAGGCCACAGACCATTTCTCAGGGATTTGCCGAATGATGGGTTTTTGGGATGGAACCACTACATGCAAAGCAGTTCATAGAAGACGACGTGAACAAGAAACTTTCGAGCGTTCATCAATTCGGGTCAATGGGCAACTACTTAAGAAAAGGGCTTCTGAATAAAGAACATACATGTTGAATTATCAAGCAAATTTTGCATTTAAAGTAACTTAAGGTACGAATTCTTACTCAAAATTCAGTAGAAACGTTCTTAGGAAAATTAATACCACTAAAGAGCTAGCAATAGGCTTATCTCCTTACATTGGCATGCTTGGTCTCAATGATTGTACTTTGGGAATTTTGATATCACTATGGAGAGTCCTGGATCCCCATCTTCGCCTGGAGCTTTGCACTATTTGTACGTGATTTGGCGCTCAGTCTTGGGGACTTGGGTTACCCATTTGCAGCGTGCCGCCGGCAGAAGGCAACCAAAAATTTCACTTGTGGCTGGCTTGTCTGTTAATTTCGAAGCGGGCTGTCCACGTCCTTGTAATCAGGCAAAGTTTCCTGACTGCTGCACCAGGCTGCAGTGCACCACTCCTGCTGGTTGTCGGTTGCTAGTTGTTAGTTGCTTCGTTGGTCGGTTGTCGTTGTCGGCGGTTGCTGGTTGCTGGTTGTTGGCGATTGTTGGCCGACTGGCGGGATGGCGGGATGGCGGGATTGGGCCGGTTGTTGGACAAACCCGCAGCGGCAAAACCCCAAAATGTCTTGCAGAAAAAGTTTATCCCCCAGTCCGGCACCAGGCCAACACGAAATGGCCGCAACAGCTGTGGCAGACTTATTTTCAATGCACTTGGGAGATGAAGGGCGCAAAAAAGTAGATGCCAAGTGGGAAAGTTTATAGACTAGGAAATACCGTTTTGTTAGGAATCCGCCGTTCATTTTGTATTTACATGAACCTACCAATACATATTTAAAAGTCATATATTATCCTCCAAGTGAATGTAATTATTTTAATTTTGCTGCTGAAATAGTAAATACTTTAATATGTACATTGAAGATTGTGTGATATACAAGTATAAAATTGCATTGATTTTCCTGCATTTTCCTACTTATACTCCAGAAGCCGCACAGACTCTTTAAATCTTACGCTCGGTGGAAAAAATACAAAGGCTAGCCTGGGCTTTGCCTTTTAGCTTCAATTATTTATGGGCAGCCAGGCGCGTCTGCATCTGGAGCTGGTTAACCGCAAAATCCAGGGCATTAAACAAGCTTTGGACCCCACGGGCCACTCGAGTGACACTGCATTCAATTTGTAACCTTAGAATAGGCCACCCAGCTGCGCTGCTCCTCCCAAACCCCAATCTGCTTTCAATTTTTCAATCAGCGCCTGTACCTCTATCTAGGATACTGAATCTGTTCCTACACTGACAAAAAGATTCTTTTCACATATTACCTTCAAAATATCTATTTGAAACTATATGCAAAATAAAAAAGATACATCCAACTAGGGATTGTATTGTTGAAACTTCTTAAATATATATTGCTTGTCACTAAGACTATTTATGGATAATTTGAAAGAGAGTCTTATAGACCTGCTTAGACCTTAAAAATATTTTTAAGCGGTGTTTTGGGGTAATGAACAAATATCTGAATCAAAAGCTGTAACAGATTAATTCCACCTAAGAAGTGACCGCTTTTTTCTAAAAGTGTCCTTATAATCTCTCCTTTTCAGGCCACCACAACTGCCCCAGCGGTTTCGATGAGTCGGAGGAGGAGTGCGGCACCGCTCGCAAGTTGCTGGAACTGCCAGGCGGAGTTTTCGCCGCCCTTGGATGCATTGCAGCAGCCCTGACCGCCTGCCTGATATTCTGCATGTTCGGCCTGATGAGGAAGCGGAAAAAGTCGGTGGTGCAGAAGGGCGCCGGAGTGGGCGGCATGGGCGGAATGGGAGTGGGTGGGATGGGCGTGGGCATCGGGTTCATGAACGGAGCGAGCAACGGCAACGTGTCCGCCGCCAGCATCGGTACCTTGAAGAAGGACTTCAAGAAGGAGGCGCTGTACATCGATCCGGCCTCCTGACACGGACTACAGCCGGTCGAATATATCCACGTGGACCGGGAGACCACCGTGTGATATGTTGCCGCCGGCTGTGGCATTCCCTATAATCATGTGGAGCTTTCCATGGAGCTGGCCTGACCGATGGATCAGTCCGCTGATCCCGTTCAAAATCTGGCCTGGAGATTTTAAAGATACAGAGTCGACGCCTGGAGGTATGAGAACACTAAATGAACCTACATTTTAAGCCGATAACTTTTTCTAAACGCTTCTGGGAGTAACAAATTGCGCACTAGCGAAGTTGTGTGTATGTAGTCATAAGGAAAATACTGTTATTAAGTTGCGTATGGTTTTTAGTTTGTATATATAATTTAATTTTAGGCATAAATAATTTGAATATCGTTGAAGGAGTCAGAACTGTGCACAGAACTTTGTGCAAGAACTTTTTGAATGTAGACTTTGTTTCCTATCATTAATCTTTCAATATCAGAAAATTATCATCGAGTTCAAACTAAAAATCTTTATAATTTCATGCATTAAAGGTTTCCAAATTAAAATCAAAATACATAAACTTTTCACAGTAAATAAAATTTTAAAATTCACGTCTTAAGGAATGTGTATATTAAGCAATAAAATCACTCACAACTTTTTAAAACTCAAAATGTTCGAGTAATATTAATTTATTGCTATTCACATGGTTGCCCGTGCTGACCCTTGCATAATTGTTTAGTTCGATTATGAGAGTTATCAATTTTTCACTTTAAGACTACATACACAAACGAAAAGAGGACATAAAACTATTGTAAATTAGGTACCTTTAAACCAGTTAGATCAGAGAACAACAATTCGTGGGCGGTTCGTAGATGACTTCGGCATATTACCACCCACATAGCACACAGCACCCATTAAAGGAAGTAAACCAATGACATTTAGCAACTAAATGAAAGATCATAAACTACAATGAAATTGTATAGAAGAGCAACTCAAAGAGCAGATCCCTATGCGCAGGGCCACTCTGGCACTAATTTTAAGCGGTTAGACTAGGTAAAAAACTGTTAGCAGCGATTGGAAAGGAGTTGAGGAAGTTCGAGGGTAATCTCTACAGTGTTGTGATCTAAAGCATGCACATACTGAACTATATTTAAACACTATGTTTAACGTCGAAGCATTTACATGACGAGTACTCAGAACCGGGGGCGTATATTGTATAGGTACTAGCTAATAGCTAAATAGTCAACTAAGAGGCAGCACTTTTCTGTTATTATCGTACAGCACAAACACATAAGCAAATGACATCACCTAGTGGTTTATTATATACAAGTCGATACAGCAAGTTAATTAAATGCAATTGTTATTAATGAATTAAGTAGAGCTCCAAACCAGAGATATGTTCACTGCACTTAATCGAGAGAACTGCCAAAGGATTAGACAACAAAACAACAAAACAACAAAACAATCAATGTAAATATTTCTGTCGAGCAATAACTTGTTTTGCGTTAGCATTTATGTTCGTAATTTGAACACAAACCGATTCCACTTAAGTAGCTGCAGTGAGTAATATTCCTAAAGTCAATAGACTTTGTAAAAACTACGAGCGAAGAGTAATAATAAACTGTATGCAAAACTACCATTTAAGACACGGCTCATAACCCATATAAACAATCAAACTTTTGGCAGTTGCAGTAGGATTAAGCAAACAACACCAAGAAGAAGCGGAAAAAGGACATAAGTACTTTATTTTATAGACACTGAACTAGGGCAAAAATTATGAAACTGAATACTGAAAGCAAAGCAAATACTGGCTGCAGTGGGAATTTTAAAAATTTAAGAACCTCTGCCAACGCCGGCAGAGTCCACATCCATAAAAACCCATTAAAATTTCAACTGCAATGAAACTTTGCCGAAAAGCGTGTGTATAAAAAATATATGATTTAGCCAGCATAAAAAGTAAATGAAAATATATACATAAATTAAATCAAATATTTTCGGTTGGCATTTTTTGATAGCAAATAACCTAAATATTTATAGCTAAAAAAGAAAGAAAGAAAAATTACCAAACGCACTTTGTACATTTTATAAAAAACGAGCAGTTTATATTATTTCCAGACACACAGCATAGCATACACATACATATATGATATTTTACACACATTCGAATAAGAGATAAAATAATACCGAAACATATACGAGAGCGAAAAGAAATTAAATTAATTAAAAAAGGCGCATAGGCCGAAAATGCATTTTTACTGCCGCCAATTGAGGACAAAACTAAAATCGAAAACACGGCCGGCAGGGCATTTGAATAGTTAGCTTTAAGCCCGGCTAAAAATGCTGTAACTACTTTGTTGGTCACCATTGTCACTCACACAAACTCACATCTCAATCATCAGTTCGGTTCAGTTCCTTCTCTCAACCAAACCAAACTATGGAAATTTCTGTGCCAACTGTGTTTAAAGTAGAATTCAAAGATAAAATACGAGATAAAAGAGGGGAGGCCGGCGGATGGAGAAATTGAATTATGACCGCACACGAAACATTGCATATAATTTATTTGAAATGTTCAAAAATAAAGTAGACGAAACGAGAAGAAAACGCTGCGCGATGGTAAATAACAATATAATACAGACAGCTGTGTAAATAGTCAAACTGCGTACAGTCTACGCATATTTAATTAAAACTCGAAAAAAACTGATAAAACTGAAAAAACTGAAAAAAGTGACCCAATGAGATGGAATTTGAGAAAGGAGATCATGTCTTAAGTTCCATTTCGGAATTCACCCGATTTTCATTTCCTTTCACCGTTCGATAGAGTAAGTTTCACTTGTATGATAAATAAAACTAAACCAAATATTTATATAGAGATATTATACGTATATGAAAAGGAAACCTTAAGTTAATGACGCACTGAATCAAATAGATGAAACGAGTATTTAAAACCAAGACAATTAAAACCAAATAGGAACTTTAAGCAAAATAAAACCGATGAAAAAATGAAACTGAACACAGAGCACCTTTCTCTATTAGCCAAAATTGCAAAAGTGGGTGGTTGGGCGTTTGATTTTCGGGTGGAATAGATGGGGTGTTTAGTGGGTGGAATGCCAACTGTTCCGTTTGACGCTCGGGGGAAAACCGTTGGCCAAACTAAGGCGAACTCAAATCGTGTCCAAATAGCTGGCAGTCAAGTAGGCCACAGATTTTGATTGCCCGCATCCGGGCAGCAATCGGTGGACATCCTCGAAAATTCGGAACATCGAACGAGATGCTGTTAAAGCCAAGCCAAGCCGGTAGCTCCAGGCGTGTGTTCTCACCGGTCGGAACTCGGATTGGCCGTGGGATCGGGTTGTGGCTATGGGCTTTCCATGTCCCCCGCCTCGCCCACGCTGAGCCCGCATAAAACATAAGCCACCGGCCTGGTGCTCGTGTGGTGTGTGCTCTATATCCCACAGATTCCGCCTCCTCCCTCGGAATGCCCAGTTCTCTGACTCCGTGTTTATGGCCGGCAGCTCTGTTTCTGGACCTGGCAACTCTGGCCCCTGGACTTGGCCAAAATGGCTTTGATGGGCGTGGGCGACAAAAGGACCGGAAATCGCTAACGAGGCGTGGTAAAGGAGTTTTTGCAGTAAGTTTGAATATCAAGGTAGATACTCAAAAAATTGAATAACAAAGCATAATCAAATTTAATGGAGGGATCTACCTTGAATGTATTTTGTTATTATTGTATTACAGCTTATCAAACTAATCACTAGTATGTTTGTTGCTTTTATTTATTATTAACACATTTAAATATCTGGCTATTAAAGTTAATAATTTTGCAGTTCTTTTTATGCTGTGGTCATTAATATACATATAATACATATATCTACTAAAAGTATTCACTCTAGATTAATGCAAATTGATGAAATTATTTATTTATCGTACGCACTT",
                // >fa_m64001_190914_015449/101320295/ccs:10320-18430 reversed
                "TATGAAGAGCATGTGAGCTAGTCCACTCGAGCAACTGAAGACTACACAAAAAGCCTCTTAGGTCGAGTATTTTGGGAATTTTTCACTCAAAATTAATTGCATTTTAATAGTTGCTAATATTTCAAAAAACTTGTGCACTTATTTGCACCGGCATCAGGCTGCTCAGGCTCGCCGCTTTCCCCATCCCCTTCTCCATGAAGCTCCATCAAAAAGTCCTGCTTAAAAATGCAATAACAGGAAATGGCTACGACGGCAACAAAATACAGGAAATTGAGCGCAAAGAAATATGCTATATAGGAGATGAAAAAATGAAAACGCAAACGAGACGAAGTCCAGGCAAACAAACATGGCCGTCAGACATGAACATATGTATACAGCCGCGGTTAACCAATTAGGAGTTTTGAATTCCAGAAAAAAGAATTAATACAAATATGCACTTTATTTACTAAGATTCATTTACATTTGGATTCGTTGGTTGGTTTCAGCAAGCGTACTTGAGCTTTCAGGCAGAAAAGTGCTCATGCTTTAGCATACTGAGTCTTAGAATAGCAAAGTAGCTACAAAAAATTTAAAATCAAAGTTAAGTTAAAAGTTAAAATCAAAGTTAAAAAAAAAAAACTCGAATGAACGTTTTGAAAAAGGGGAAAAGGGAAACCCCATGGAAAGCAACACAACTATTACAATTACTTTCTGTTATATGAGATTAGTATTTTTTTTGAAGACACTAAACATGTGGCCCCTGCTGTACACCTGAACCATGCAAAGCGAGCATGAGGCGAAGACAGGTGGAAAGGTGGCGTTACGGGTAAGAGCATTTTGCCCGCTGGAACTTGGCCGTATGTGCTCCGGCAATTTATGTCCAAGGTGCGTGCCACCGAAATGAACCGGTGAACCTGGACGTTACGCTCGTTTCTTGGGGGAAAGTCCGGCAGTTGACCAACCTGTGTCCTCTTTCCTGGTTCTGCGCTTAACGTGTTGCCATTAAATTTTATGGTCAGGTCTTCCTCTCTTTGTCCGTGTTTGCATTTTGCTCAATTTTGTCGCCGCAGAACATAAAAAAGGGAGTCAGGATGGGGCTGAAAAAGAAGCGAAGGAATATTGCACTGTGGACGACATGTGAAAACATATTTTTATGGGACCCCAAAGGGTAATTGGACATTGCACGGGATATCCTAGCTGGCTAATGAGTGTGTCGAGCCTTTTTAGATAATTAATTGAGAGCTTTTCAAATCATCTGTGAAGTATAATGACCACACACACACAATAAGGCAACGAAACGTTGAAGTGTGATTGCACTAATGAGTAGATTAAAGTCAGCTGGAGATAACATCTTAATAAATCCAAGATAATTGAAATATCCTTTCGAGATAGTTGTATAATTTGTAATTCTCCATTCGGGAATATTTGTTTCATTGGCTTCAAGTGAAACCAATTGTGCCGCACTTCATCCAGCGAAACTTTGACGGCTAACTCGCTTAATGCGCGGCTTAAGCCAAACTACGAGTAAGATATTTTGTCTGTGGCAACTTGAAAGTTGCAAGTAATTAGTGGCACTTGAGCAACTACTAGGCTGTCGGGATACGAAAATGGGGCAGGGCATACCCAGATAGCTGGATGCTCCAAAAGTTATGTAGATCACGAAAATCAGTCCCTTGAGAACAAATAACTTTTAGCAATCTTAAAATATTATTTCTAAGTATTCTAATCGTTATACTTTATTAGTATTTTAGTTCTATGACTTTCTTTAAAAAAAAAAAAAAAAACACTGATCTTACAAGCTCGCTTTCCACCTATCCTCTTTCAGACCCATCAGCCTAGTGCTGGAGCCAATCCTGAAGGAGCCCGGCGACATCTCCTTTACGTGGCTCGAGATTCACCGCACCAAGAACGCGCTGCTGCAGCAGTTGGACCTGCATGTGAACGCTAGCGTGACCACGGGCATGGGCTCCCAGACCTCCTCCCAGTCCGCAGCCGGCTCCAACACGGGCATTGGGGGACTGGGCACCAGTGGGTCAGGATCTGCTGCCAATGCCGGGGCTCTGCTGGCGGAGCAACGGCCAACGAGACGCTGAACGAGTTCGGATTCTTTCCAAGGAGTCCGACTGCGAATACAAATGTCCTGAAATTGATGCATGCATTGCGGCCAGCTTGTGGTGCGATGGTGAGTACAGATTGCCAGTCGAGGCCACAGACCATTTCTCAGGGATTTGCCGAATGATGGGTTTTTGGGATGGAACCACTACATGCAAAGCAGTGCATAGAAGACGACGTGAACAAGAAACTTTCGAGCGTTCATCAATTCGGGTCAATGGGCAACTACTTAAGAAAAGGGCTTCTGAATAAAGAACATACATGTTGAATTATCAAGCAAATTTTGCATTTAAAGTAACTTAAGGTACGAATTCTTACTCAAAATTCAGTAGAAACGTTCTTAGGAAAATTAATACACTAAAGAGCTAGCAATAGGCTTATCTCCTTACATTGGCATGCTTGGTCTCAATGATTGTACTTTGGGAATTTTGATATCACTATGGAGAGTCCTGGATCCCCATCTTCGCCTGGAGCTTTGCACTATTTGTACGTGATTTGGCGCTCAGTCTTGGGGATTTGGGTTACCCATTTGCAGCGTGCCGCCGGCAGAAGGCAACCAAAAATTTCACTTGTGGCTGGCTTGTCTGTTAATTTCGAAGCGGGCTGTCCACGTCCTTGTAATCAGGCAAAGTTTCCTGACTGCTGCACCAGGCTGCAGTGCACCACTCCTGCTGGTTGTCGGTTGCTAGTTGTTAGTTGCTTCGTTGGTCGGTTGTCGTTGTCGGCGGTTGCTGGTTGCTGGTTGTTGGCGATTGTTGGCCGACTGGCGGGATGGCGGGATGGCGGGATTGGGCCGGTTGTTGGACAAACCCGCAGCGGCAAAACCCCAAAATGTCTTGCAGAAAAAGTTTATCCCCCAGTCCGGCACCAGGCCAACACGAAATGGCCGCAACAGCTGTGGCAGACTTATTTTCAATGCACTTGGGAGATGAAGGGCGCAAAAAAGTAGATGCCAAGTGGGAAAGTTTATAGACTAGGAAATACCGTTTTGTTAGGAATCCGCCGTTCATTTTGTATTTACATGAACCTACCAATACATATTTAAAAGTCATATATTATCCTCCAAGTGAATGTAATTATTTTAAATTTGCTGCTGGAATAGTAAATACTTTAATATGTACATTGAAGATTGTGTGATATACAAGTATAAAATTGCATTGATTTTCCTGCATTTTCCTACTTATACTCCAGAAGCCGCACAGACTCTTTAAATCTTACGCTCGGTGGAAAAAATACAAAGGCTAGCCTGGGCTTTGCCTTTTAGCTTCAATTATTTATGGGCAGCCAGGCGCGTCTGCATCTGGAGCTGGTTAACCGCAAAATCCAGGGCATCAAACAAGCTTTGGACCCCACGGGCCACTCGAGTGACACTGCATTCAATTTGTAACCTTAGAATAGGCCACCAGCTGCGCTGCTCCTCCCAAACCCCAATCTGCTTTCAATTTTTCAATCAGCGCCTGTACCTCTATCTAGGATACTGAATCTGTTCCTACACTGACAAAAAGATTCTTTTCAAATATTACCTTCAAAATATCTATTTGAAACTATATGCAAAATAAAAAAGATACATCCAACTAGGGATTGTATTGTTGAAACTTCTTAAATATATATTGCTTGTCACTAAGACTATTTATGGATAATTTGAAAGAGATTCTTATAGACCTGCTTAGACCTTAAAAATATTTTTAAGCGGTGTTTTGGGGTAATGAACAATATCTGAATCAAAAGCTGTAACAGATAATTCCACCTAAGAAGTGACCGCTTTTTTCTAAAAGTGTCCTTATAATCTCTCCTTTTCAGGCCACCACAACTGCCCCAGTGGTTTTCGATGAGTCGGAGGAGGAGTGCGGCACCGCTCGCAAGTTGCTGGAACTGCCAGGCGGAGTTTTCGCCGCCCTTGGATGCATTGCAGCAGCCCTGACCGCCTGCCTGATATTCTGCATGTTCGGCCTGATGAGGAAGCGGAAAAAGTCGGTGGTGCAGAAGGGCGCCGGAGTGGGCGGCATGGGCGGAATGGGAGTGGGTGGGATGGGCGTGGGCATCGGGTTCATGAACGGAGCGAGCAACGGCAACGTGTCCGCCGCCAGCATCGGTACCTTGAAGAAGGACTTCAAGAAGGAGGCGCTGTACATCGATCCGGCCTCCTGACACGGACTACAGCCGGTCGAATATATCCACGTGGACCGGGAGACCACCGTGTGATATGTTGCCGCCGGCTGTGGCATTCCCTATGATCATGTGGAGCTTTCCATGGAGCTGGCCTGACCGATGGATCAGTCCGCTGATCCCGTTCAAAATCTGGCCTGGAGATTTTAAAGATACAGAGTCGACGCCTGGAGGTATGAGAACACTAAATGAACCTACATTTTAAGCCGATAACTTTTTCTAAACGCTTCTGGGAGTAACAAATTGCGCACTAGCGAAGTTGTGTGTATGTAGTCATAAGGAAAATACTGTTATTAAGTTGCGTATGGTTTTTAGTTTGTATATATAATTTAATTTTAGGCATAAATAATTTGAAAATCGTTGAAGGAGTCAGAACTGTGCACAGAACTTTGTGCAAGAACTTTTTGAATGTAGACTTTGTTTCCTATCATTAATCTTTCAATATCAGAAAATTATCTTGGAGTTCAAACTAAAAATCTTTATAATTTCATGCATTAAAGGTTCCAAATTAAAATCAAAATACATAAACTTTTCACAGTAAATAAAATTTTAAAATTCACGTCTTAAGGAATGTGTATATTAAGCAATAAAATCACACACAACTTTTTAAAACTCAAAATGTTCGAGTAATATTAATTTATTGCTAGTCACATGGTTGCCCGTGCTGACCCTTGCATAATTGTTTTAGTTCGATTATGAGAGTTAACAATTTTTCACTTTAAGACTACATACACAAACGAAAAGAGGACATAAAACTATTGTAAATTAGGTACCTTTAAACCAGTTAGATCAGCTAGAGAACAACAATTCGTGGGCGGATCGTAGATGACTTCGGCATATTACCACATAGCACACAGCACCCATTAAAGAAAGTAAACCAATGACATTTAGCAACTAAATGAAAGATCATAAACTACAATGAAATTGTATAGAAGAGCAACTCAAAGAGCAGATCCCTATGCGCAGGGCCACTCTGGCACTAATTTTAAGCGGTTAGACTAGGTAAAAAACTGTTAGCAGCGATTGGAAAGGAGTTGAGGAAGTTCGAGGGTAATCTCTACAGTGTTGTGATCTAAGGCATGCACATACTGAACTATGTTTAAACACTATGTTTAACGTCGAAGCATTTACATGACGAGTACTCAGAACCGGGGGCGTATATTGTATAGGTACTAGCTAGTAGCTAAATAGTCAACTAAGAGGCAGCACTTTTCTGTTATTATCGTACAGCACAAACAAGCAAATGACATCACCTAGTGGTTTATTATATACAAGTCGATACAGCTAGTTAATTAAATGCAATTGTTATTAATGAATTAAGTAGAGCTCCAAACCAGAGATATGCACTTAATCGAGAGAACTGCCAAAGGATTAGACAACAAAACAACAAAACAATCAATGTAAATATTTCTGTCGAGCAATAACTTGTTTTGCGTTAGCATTTATGTTCGTAATTTGAACACAAACCGATTCCACTTAAGTAGCTGCAGTGAGTAATATTCCTAAAGTCAATAGACTTTGTAAAAACTACGAGCGAAGAATAATAATAAACTGTATGCAAAACTACCATTTAAGACACGGCTCATAACCCATATAAACAATCAAACTTTTGGCAGTTGCAGTAGGATTAAGCAAACAACACCAAGAAGAAGCGGAAAAGGACATAAGTACTTTATTTTATAGACACTGAACTAGGGCAAAAATTATGAAACTGAATACTGAAAGCAAAGCAAATACTGGCTGCAGTGGGAATTTTAAAAATTTAAGAACCTCTGCCACGCCGGCAGAGTCCACATCCATAAAAACCCATTAAATTTCAACTGCAATGAAACTTTGCCGAAAAGCGTGTGTATAAAAAATATATGATTTAGCCAGCATAAAAGTAAATGAAAATATATACATAAATTAAATCAAATATTTTCGGTTGGCATTTTTTGATAGCAAATAACCTAAATATTTATAGCTAAAAAAGAAAGAAAGAAAAATTACCAAAACGCACTTTGTACAATTTAATAAAAAACGAGCAGTTTATATTATTTCCAGACACACAGCATAGCATACACATACATATATGATATTTTACACACATTCGAATAAGAGATAAAATAATACCGAAACATATACGAGAGCGAAAAGAAATTAAATTAATTAAAAAGGGCGCATAGGCCGAAAATGCATTTTTACTGCCGCCAATTGAGGACAAAACTAAAAATCGAAACACGGCCGGCAGGGCATTTGAATAGTTAGCTTTAAGCCCGGCTAAAAATGCTGTAACTACTTTGTTGGCCACCATTGTCACTCACACAAACTCACATCTCAATCATCAGTTCGGTTCAGTTCCTTCTCTCAACCAAACCAAACTATGGAAATTTCTGTGCCAACTGTGTTTAAAGTAGAATTCAAAGATAAAATACGAGATAAAAGAGGGGGGGCCGGCGGATGGAGAAATTGAATTATGACCGCACACGAAACATTGCATATAATTTATTTGAAATGTTCAAAAATAAAGTAGACGAAACGAGAAGAAAACGCTGCGCGATGGTAAATAACAATATAATACAGACAGCTGTGTAAATAGTCAAACTGCGTACAGTCTACGCATATTTAATTAAAAACCGAAAAAAAAAACTGATAAAACTGAAAAAAGTGACGCAATGAGATGGAATTTGAGAAAGGAGATCATGTCTTAAGTTCCATTTCGGAATTCACCCGATTTTCATTTCCTTTCACCGTTCGATAGAGTAAGTTTCACTTGTATGATAAATAAAACTAAACCAAATATTTATATAGAGATATTATACGTATATGAAAAGGAAACCTTAAGTTAATGACGCACTGAATCAAATAGATGAAACGAGTATTTAAAACCAAGACAATTAAAACCAAAAGGAACTTTAAGCAAAATAAAACCGATGAAAAAATTAAACTGAACACAGAGCACCTTTCTCTATTAGCCAAAATTGCAAAAGTGGGTGGTTGGGCGTTTGATTTTCGGGTGGAATAGATGGGGTGTTTAGTGGGTGGAATGCCAACTGTTCCGTTTGACGCTCGGGGGAAAACCGTTGGCCAAACTAAGGCGAACTCAAATCGTGTCCAAATAGCTGGCAGTCAAGTAGGCCACAGATTTTGATTGCCCGCATCCGGGCAGCAATCGGTGGACATCCTCGAAAATTCGGAACATCGAACGAGATGCTGTTAAAGCCAAGCCAAGCCGGTAGCTCCAGGCGTGTGTTCTCACCGGTCGGAACTCGGATTGGCCGTGGGATCGGGTTGTGGCTATGGGCTTTCCATGTCCCCCGCCTCGCCCACGCTGAGCCCGCATAAAACATAAGCCACCGGCCTGGTGCTCGTGTGGTGTGTGCTCTATATCCCACAGATTTCGCCTCCTCCCTCGGATTGCCCAGTTCTCTGACTCCGTGTTTATGGCCGGCAGCTCTGTTTCTGGACCTGGCAACTCTGGCCCGTGGACTTGGCCAAAATGGCTTTGATGGGCGTGGGCGACAAAAGGACCGGAAATCGCTAACGAGGCGTGGTAAAGCAGTTTTTGCAGTAAGTTTGAATATCAAGGTAGATACTCAAAAAATTGAATAACAAAGCATAATCAAATTTAATGGAGGGATCTACCTTGAATGTATTTTGTTATTATTGTATTACAGCTTATCAAACTAATCACTAGTATGTTTGTTGCTTTTATTTATTATTAACACATTTAAATATCTGGCTATTAAAGTTAATAATTTTGCAGTTCTTTTTATGCCGTGGTCATTAATATACATATAATACATATATCTACTAAAAGTATTCACTCTAGATTAATGCAAATTGATGAAATTATTTATTTATCGTACGCACTT",
                // Alen = 8138, Blen = 8111
                static_cast<int32_t>(0.03 * 8138),
                static_cast<int32_t>(0.01 * 8111),
                // Full length alignment.
                // Edit distance was confirmed with both Raptor and Minimap2.
                SesResults(8138, 8111, 105, true)
    },
    // clang-format on
};

TEST(Ses2DistanceBanded, Semiglobal_AllTests)
{
    for (const auto& data : testDataSemiglobal) {
        // Name the test.
        SCOPED_TRACE("Semiglobal-" + data.testName);

        // Run.
        SesResults result = SES2DistanceBanded<SESAlignMode::Semiglobal, SESTrimmingMode::Disabled>(
            data.query, data.target, data.maxDiffs, data.bandwidth);

        // Evaluate.
        EXPECT_EQ(data.expected, result);
    }
}

TEST(Ses2DistanceBanded, Global_AllTests)
{
    for (const auto& data : testDataGlobal) {
        // Name the test.
        SCOPED_TRACE("Global-" + data.testName);

        // Run.
        SesResults result = SES2DistanceBanded<SESAlignMode::Global, SESTrimmingMode::Disabled>(
            data.query, data.target, data.maxDiffs, data.bandwidth);

        // Evaluate.
        EXPECT_EQ(data.expected, result);
    }
}

}  // namespace Test
}  // namespace Alignment
}  // namespace Pancake
}  // namespace PacBio
