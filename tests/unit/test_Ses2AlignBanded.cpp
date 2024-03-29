// Authors: Ivan Sovic

#include <PancakeTestData.h>
#include <gtest/gtest.h>
#include <pancake/Ses2AlignBanded.hpp>
#include <sstream>
#include <tuple>

namespace PacBio {
namespace Pancake {
namespace Alignment {
namespace Tests {

struct TestData
{
    std::string testName;
    std::string query;
    std::string target;
    int32_t maxDiffs = 0;
    int32_t bandwidth = 0;
    SesResults expectedGlobal = {0, 0, 0, false};
    SesResults expectedSemiglobal = {0, 0, 0, false};
};

// clang-format off
std::vector<TestData> testDataGlobal = {
    TestData{"EmptyQueryEmptyTarget", "", "", 100, 30,
                SesResults(0, 0, 0, 0, 0, 0, 0, true, PacBio::Data::Cigar("")),
                SesResults(0, 0, 0, 0, 0, 0, 0, true, PacBio::Data::Cigar("")),
    },
    TestData{"EmptyQueryNonemptyTarget", "", "ACTG", 100, 30,
                SesResults(0, 0, 0, 0, 0, 0, 0, true, PacBio::Data::Cigar("")),
                SesResults(0, 0, 0, 0, 0, 0, 0, true, PacBio::Data::Cigar("")),
    },
    TestData{"NonEmptyQueryEmptyTarget", "ACTG", "", 100, 30,
                SesResults(0, 0, 0, 0, 0, 0, 0, true, PacBio::Data::Cigar("")),
                SesResults(0, 0, 0, 0, 0, 0, 0, true, PacBio::Data::Cigar("")),
    },
    TestData{"SimpleSingleIndelDiff", "ACG", "ACTG", 15, 30,
                SesResults(3, 4, 1, 3, 0, 0, 1, true, PacBio::Data::Cigar("2=1D1=")),
                SesResults(3, 4, 1, 3, 0, 0, 1, true, PacBio::Data::Cigar("2=1D1=")),
    },
    TestData{"SimpleSingleMismatchDiff", "AAAAA", "AAATA", 15, 30,
                SesResults(5, 5, 1, 4, 1, 0, 0, true, PacBio::Data::Cigar("3=1X1=")),
                SesResults(4, 5, 1, 4, 0, 0, 1, true, PacBio::Data::Cigar("3=1D1=")),
    },
    TestData{"SimpleFiveBaseDeletion", "AAAAAAAAAAAAAAAAAAAA", "AAAAAAGGGGGAAAAAAAAAAAAAA", 15, 30,
                SesResults(20, 25, 5, 20, 0, 0, 5, true, PacBio::Data::Cigar("6=5D14=")),
                /*
                 * Deletions are handled before mismatches and insertions, that's why the
                 * semiglobal case is identical to the global one. Otherwise, an alternative would be
                 *      SesResults(20, 20, 5, 15, 5, 0, 5, true, PacBio::Data::Cigar("6=5X9=")),
                 * See the counter example below (SimpleFiveBaseInsertion).
                */
                SesResults(20, 25, 5, 20, 0, 0, 5, true, PacBio::Data::Cigar("6=5D14=")),
    },
    TestData{"SimpleFiveBaseInsertion", "AAAAAAGGGGGAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAA", 15, 30,
                SesResults(25, 20, 5, 20, 0, 5, 0, true, PacBio::Data::Cigar("6=5I14=")),
                // The semiglobal aligner tries to find the shortest alignment to reach an end of a sequence.
                SesResults(20, 20, 5, 15, 5, 0, 0, true, PacBio::Data::Cigar("6=5X9=")),
    },
    TestData{"SimpleIndelsOnly", "GATGTTT", "GAATTGTT", 15, 30,
                SesResults(7, 8, 3, 6, 0, 1, 2, true, PacBio::Data::Cigar("2=1D1=1D3=1I")),
                SesResults(6, 8, 2, 6, 0, 0, 2, true, PacBio::Data::Cigar("2=1D1=1D3=")),
                /*
                 * GA-T-GTTT
                 * || | |||
                 * GAATTGTT-
                */
    },
    TestData{"NarrowBandwidthButOnlyMismatches_NotValid", "GGATCAGTT", "GAATTCGTT", 100, 2,
                /*
                * The maximum number of diffs is within the allowed limit, but
                * the maximum bandwidth for alignment is below the threshold (3).
                * The output should be marked as "invalid" but also report the
                * furthest reaching (f.r.) wave.
                *      GGATCAGTT
                *      |X||XX|||
                *      GAATTCGTT
                */
                SesResults(5, 4, 2, 3, 1, 1, 0, false, PacBio::Data::Cigar("1=1X2=1I")),
                SesResults(5, 4, 2, 3, 1, 1, 0, false, PacBio::Data::Cigar("1=1X2=1I")),
    },
    TestData{"OutOfBandwidth_NotValid", "GATGTTT", "GAATTGTT", 100, 2,
                /*
                * The maximum number of diffs is within the allowed limit, but
                * the maximum bandwidth for alignment is below the threshold (3).
                * The output should be marked as "invalid" but also report the
                * furthest reaching (f.r.) wave.
                * GA-T-GTTT
                * || | |||
                * GAATTGTT-
                */
                SesResults(3, 2, 1, 2, 0, 1, 0, false, PacBio::Data::Cigar("2=1I")),
                SesResults(3, 2, 1, 2, 0, 1, 0, false, PacBio::Data::Cigar("2=1I")),
    },
    TestData{"AboveMaxDiffs_NotValid", "GGATCAGTT", "GAATTCGTT", 3, 30,
                /*
                * The maximum number of diffs is too small, so the alignment
                * should be marked as not valid.
                * The maximum bandwidth for alignment is GOOD here (30 > 3).
                * The output should be marked as "invalid" but also report the
                * furthest reaching (f.r.) wave.
                *      GGATCAGTT
                *      |X||XX|||
                *      GAATTCGTT
                */
                SesResults(5, 4, 2, 3, 1, 1, 0, false, PacBio::Data::Cigar("1=1X2=1I")),
                SesResults(5, 4, 2, 3, 1, 1, 0, false, PacBio::Data::Cigar("1=1X2=1I")),
    },


    TestData{"NormalSmallCase_SingleMatch", "A", "A", 15, 30,
                SesResults(1, 1, 0, 1, 0, 0, 0, true, PacBio::Data::Cigar("1=")),
                SesResults(1, 1, 0, 1, 0, 0, 0, true, PacBio::Data::Cigar("1=")),
    },
    TestData{"NormalSmallCase_MultipleExactMatches", "ACTG", "ACTG", 15, 30,
                SesResults(4, 4, 0, 4, 0, 0, 0, true, PacBio::Data::Cigar("4=")),
                SesResults(4, 4, 0, 4, 0, 0, 0, true, PacBio::Data::Cigar("4=")),
    },
    TestData{"NormalSmallCase_SingleMismatch", "A", "C", 5, 30,
                SesResults(1, 1, 1, 0, 1, 0, 0, true, PacBio::Data::Cigar("1X")),
                SesResults(0, 1, 1, 0, 0, 0, 1, true, PacBio::Data::Cigar("1D")),
    },
    TestData{"NormalSmallCase_MultipleMismatches", "CCCC", "GGGG", 15, 30,
                SesResults(4, 4, 4, 0, 4, 0, 0, true, PacBio::Data::Cigar("4X")),
                SesResults(0, 4, 4, 0, 0, 0, 4, true, PacBio::Data::Cigar("4D")),
    },
    TestData{"NormalSmallCase_CIGAR_5I_in_suffix", "ACTGAAAAA", "ACTG", 15, 30,
                SesResults(9, 4, 5, 4, 0, 5, 0, true, PacBio::Data::Cigar("4=5I")),
                SesResults(4, 4, 0, 4, 0, 0, 0, true, PacBio::Data::Cigar("4=")),
    },
    TestData{"NormalSmallCase_CIGAR_5D_in_suffix", "ACTG", "ACTGAAAAA", 15, 30,
                SesResults(4, 9, 5, 4, 0, 0, 5, true, PacBio::Data::Cigar("4=5D")),
                SesResults(4, 4, 0, 4, 0, 0, 0, true, PacBio::Data::Cigar("4=")),
    },
    TestData{"NormalSmallCase_CIGAR_5D_in_prefix", "ACTG", "CCCCCACTG", 15, 30,
                SesResults(4, 9, 5, 4, 0, 0, 5, true, PacBio::Data::Cigar("5D4=")),
                SesResults(4, 4, 3, 1, 3, 0, 0, true, PacBio::Data::Cigar("1X1=2X")),
    },
    TestData{"NormalSmallCase_CIGAR_5I_in_prefix", "CCCCCACTG", "ACTG", 15, 30,
                SesResults(9, 4, 5, 4, 0, 5, 0, true, PacBio::Data::Cigar("5I4=")),
                SesResults(1, 4, 3, 1, 0, 0, 3, true, PacBio::Data::Cigar("1D1=2D")),
    },
    TestData{"NormalSmallCase_SimpleSeq", "GGATCAGTT", "GAATTCGTT", 5, 30,
                SesResults(9, 9, 3, 7, 1, 1, 1, true, PacBio::Data::Cigar("1=1X2=1D1=1I3=")),
                SesResults(9, 9, 3, 7, 1, 1, 1, true, PacBio::Data::Cigar("1=1X2=1D1=1I3=")),
                /*
                 * GGAT-CAGTT
                 * |X|| | |||
                 * GAATTC-GTT
                */
    },
    TestData{"NormalSmallCase_AnotherSimpleSeq", "GGATTCAGTT", "GAATTCCGTT", 15, 30,
                SesResults(10, 10, 2, 8, 2, 0, 0, true, PacBio::Data::Cigar("1=1X4=1X3=")),
                SesResults(10, 10, 2, 8, 2, 0, 0, true, PacBio::Data::Cigar("1=1X4=1X3=")),
                /*
                * Correct alignment:
                * 1=1X4=1X3=
                * GGATTCAGTT
                * |X||||X|||
                * GAATTCCGTT
                *
                * This is the one that 1986 SES would find:
                * 1=1I1=1D3=1X3=
                * GGA-TTCAGTT
                * | | |||X|||
                * G-AATTCCGTT
                */
    },

    TestData{"RealHifiReads_SmallPieceFrom_sesResultLeft",
                "CATGGTGAGTCACCTCTGACTGAGAGTTTACTCACTTAGCCGCGTGTCCACTATTGCTGGGTAAGATCAGACCGTTATTCTCGACAGCGGAAGTACGACAATGCTTATCGCCGAAGGATTAATGACCGCCAAAAAATATCACGGTGATTACCAACAGTCTCCCGGCAGCGTTTGCCCTTTCCGAAAATAAAGACATTACTCTGGTCGTCTGTGGTGGCACGGTCCGCCATAAAACGCGCTCG",
                "CATGTGAGTCACCTCTGACTGAGAGTTTACTCACTTAGCCGCGTGTCCACTATTGCTGGGTAAGATCAGATTACGGTTGCGCCTGTTACCGCGGCAACGTCCTGTGCACAGAAGCTCTTATGCGTCCCCAGGTAATGAATAATTGCCTCTTTGCCCGTCATACACTTGCTCCTTTCAGTCCGAACTTAGCTTTAATTTCTGCGATCTTCGCCAGAGCCTGTGCACGATTTAGAGGTCTACCGCCCATAACAGGAAGTTGTTTTACTGGTTCAGGTATCGTCTCACCACGGTTAATTCGCGCTGTCATACAGGTCAGTTCATCGGCAGCCTTGCGCCGTAATTCCGCGTCAGCCAGCGCATTGGCCCGCATGTTCTGGTACAAGTTGGTAACCAACCAGTAATGCGCGTTCGATTTCCACGGATAAGACTCTGCATCCGGATACAGGCCACGCTTCCGGCAATACTCGTACCTCCCGGGATTTCATGAAATTCCGGCTCGGTGGTTTCGAGGCAATAAAATCGGCTTACATGGCCCAGGTGCAGTACAGCATGTGGGTGACGCGAAAAGATGCCTGGTACTTTGCCAACTATGACCCGCGCATGAAGCGTGAAGGCCTGCATTATGTCGTGATTGAGCGGAATGAAAAGTACATGGCGAGTTTTGACGAGATGGTGCCGGAGTTCATCGAAAAAATGGACGAGGCACTGGCTGAAATTGGTTTTGTATTTGGGGAGCAATGGCGATGACGCATCCTCACGATAATATCCGGGTACCTCACAACACGGCAAGCCTGCATTGCGGCGCTTCAGTCTCCGCTGCATACTGTCCAGGTGAGCGCGGGTGATGGCATAACAGAGGAAAGAAAATGTCACTCTTCCGCAGAAATGAAATATGGTATGCCTCGTATTCGCTCCCGGGCGGGAAACGAATTAAGGAATCTCTTGGCACAAAGGACAAACGGCAAGCTCAGGAGTTGCACGACAAGCGAAAAGCAGAACTCTGGCGAGTAGAAAAGCTAGGGGATTTACCTGATGTCACTTTTGAAGAGGCCTGCCTAAGATGGCTTGAGGAAAAAGCTGATAAAAAATCTCTCGATTCAGATAAAAGCCGGATTGAGTTCTGGCTTGAACATTTTGAGGGTATAAGGCTTAAAGATATCTCGGAGGCAAAGATTTACTCTGCTGTAAGCAGAATGCATAACAGAAAGACGAAAGAAATATGGAAACAGAAAGTTCAGGCCGCCATCAGGAAAGGTAAAGAACTGCCTGTTTATGAACCAAAGCCAGTATCAACTCAGACAAAGGCAAAGCATCTTGCCATGATAAAGGCCATTCTCCGTGCTGCAGAACGCGACTGGAAGTGGCTGGAAAAAGCGCCTGTCATCAAGATACCAGCGGTCAGAAACAAGCGAGTCAGATGGCTGGAAAAGGAGGAAGCAAAACGCCTTATTGATGAGTGCCCCGAACCACTGAAATCTGTCGTCAAGTTTGCGCTGGCAACTGGTCTGAGAAAGTCGAACATCATAAATCTGGAATGGCAACAAATCGACATGCAGCGACGAGTTGCCTGGGTGAATCCAGAAGAGAGCAAATCAAACCGCGCCATTGGTGTGGCGCTGAACGATACCGCCTGTAAAGTGTTGCGTGATCAAATAGGCAAGCATCACAAATGGGTGTTTGTACATACCAAGGCGGCTAAGCGAGCAGATGGAACATCAACGCCTGCGGTCAGGAAGATGCGCATCGACAGCAAGACATCATGGCTATCAGCTTGTCGTCGTGCAGGAATTGAAGATTTCCGTTTCCATGACCTCAGACACACCTGGGCAAGCTGGCTGATTCAGTCAGGCGTCCCATTATCAGTGCTTCAGGAAATGGGCGGATGGGAGTCCATAGAAATGGTTCGTAGGTATGCTCACCTTGCGCCTAATCATTTGACAGAGCATGCGAGGAAAATAGACGACATTTTTGGTGATAATGTCCCAAATATGTCCCACTCTGAAATTATGGAGGATATAAAGAAGGCGTAACTGATTGAATTGTAATGGCGCGCCCTGCAGGATTCGAACCTGCGGCCCACGACTTAGAAGGTCGTTGCTCTATCCAACTGAGCTAAGGGCGCGTTGATACCGCAATGCGGTGTAATCGCGTGAATTATACGGTCAACCCTTGCTGAGTCAATGGCTTTTGATCTGGTTGCTGAACAAGTGAACGACCGCGTCTGATTTTCTGATTTATTTCGCTATAGCGGCAAACAAACGCACACCGCTGCGCGTCTGAATCAAGAAAACCCGTATTTTCATGTATCAAAGTACAATTTCCCGACCTAACGGAAAATTGTCCGCTCCTATGAGACTGGTAACTATGAAACCAACGTCGGTGATCATTATGGATACTCATCCTATCATCAGAATGTCTATTGAAGTTCTGTTGCAAAAAAACAGTGAATTGCAGATTGTCCTGAAAACGGATGATTATCGCATAACCATCGATTATCTCCGAACCCGTCCTGTTGATTTAATCATTATGGATATAGACTTGCCCGGAACAGACGGTTTTACCTTCCTGAAAAGGATCAAACAAATCCAGAGCACAGTGAAAGTGTTATTTTTATCATCGAAATCAGAATGCTTTTATGCTGGCAGAGCGATACAAGCTGGTGCTAACGGTTTTGTCAGTAAATGCAATGATCAGAATGATATTTTTCATGCCGTTCAGATGATCCTCTCCGGATACACGTTTTTTCCCAGCGAAACGCTTAACTATATAAAAAGCAATAAATGTAGTACGAATAGTTCAACGGTCACTGTGCTATCTAATCGTGAAGTGACCATATTACGTTATCTGGTTAGCGGATTATCTAATAAAGAAATTGCCGATAAGTTATTACTTAGCAATAAAACAGTTAGTGCGCATAAATCTAATATTTATGGCAAGCTAGGTTTGCATTCAATTGTAGAGCTTATCGACTACGCCAAATTATACGAATTAATATAATATTAATTATAATTGATCATAAATATCGCATCCGCTTTCGCCACACCTGGCCGAACACCGCTGGCTAACGCTCGATAATTGGCAAAAAAGTTTAGTGTGACATTGCCATTTGCATCTACTTCCTCAGTCGGGCTCGCCTCCCCCAGTGCGAGCCGGGAGCGATCGCTATTACGTAATTCGATGGCGACGGTTTGTGCCATTGCGGGATCATCCAGAGCCAGCAGGTTGGTATCGGATGCCGGCGTTCCCGTAAATAAAATCGCAACTGAACCCGGAGGACATCCCTCCAGCCGCAGGCTAAAAGGGACGAGTGCCGTGGTATCGCCAGCGTTCAGTAGTTGTGTCGTAGGCCATCTGCCTAAATCTACCGTCTTATCAATATCCGCTGTGTTTACGGTACAGGAGAAATCAACAACGTTACCGTGCAAATTGATATTAATCGTTCCTAAAGGGTCAACTGCCCATCCACTGGAACTCCACAGTAGCCCGCAGAAACAGCTAAAGAGTACTCTTCTCATTCTCCGTACCTCATTGATAATCGACGCGTAAATACCCCAGCGCGCTAAACGGCCCTTCGGTCGGTTTTTGACCGGTAATACTGATAGGCCAGGCGCGAAGTGTGACATTGGCTGCCGCAGCTGCATCCAGACGGAAAGGAATAACGCTATTGAGATCGTTAGGCGTGATCGGCGTATCGTTCTGATCGGCGACAATAAAACCTAAATCCTGATTGTCCGACACCATCGCCTGACCAGAAACGGCACTGGCTTCCAGACGCATTGTTAAATAAGCCTGCGCAGCAACATTCGTACATTTGACCGCAATGCTCTTGGTTTGCGGCATGACACCAGCAGGTCGATTACCCGGCCCTGCCGCACTAAATAACGATGCGCCGATATCACCAAAATCAAATTCAACAATCTGCCCGGCATTTAATTCGCAGTTTTGCGGTACTTCAACCCGGCCACCAAAACTAATGGTATAAACAGTTGT",
                240,    // maxDiffs
                96,     // bandwidth
                SesResults(215, 175, 78, 143, 26, 46, 6, false, PacBio::Data::Cigar("4=1I66=2D1X2=1D2=3X1=1D1=1D2X2=1X4=1X1=2I3=1X1=1X1I2=1I1=1X2I1=1X1=1X4=1D2X3=1I2=1I1=2X2=2X1=2X1I2=2I1=3I3=1X2=2I1=4I1=1I1=1X2=1X1I1=4I7=1X1=1I1=2I1=3I2=4I2=3I2=1I2=1X2=2I2=1I1=2I1=")),
                SesResults(215, 175, 78, 143, 26, 46, 6, false, PacBio::Data::Cigar("4=1I66=2D1X2=1D2=3X1=1D1=1D2X2=1X4=1X1=2I3=1X1=1X1I2=1I1=1X2I1=1X1=1X4=1D2X3=1I2=1I1=2X2=2X1=2X1I2=2I1=3I3=1X2=2I1=4I1=1I1=1X2=1X1I1=4I7=1X1=1I1=2I1=3I2=4I2=3I2=1I2=1X2=2I2=1I1=2I1=")),
                /*
                 * This was failing on a real E. Coli run. For some reason, the match count is wrong.
                 *
                 * EDIT: The issue isn't as simply testable like this, because the problem was in the ScratchSpace which wasn't reset properly.
                 * Upon a rerun, some values weren't initialized properly and resulted in a bad alignment.
                 * Still, it cannot hurt to have this test here too.
                */
    },

    TestData{"RealHifiReads_Subsequences_CheckBandwidthAndMaxDiffsInAllocations",
                // >query len=351
                "GACAGTACAGTTTGAAACACTGTTTTTGTGGAGTCTGCAAGTGGATATTTGGCTAGATTTGAGGATTTCGTTGGAAACGGGATAAGGTATAAAAAGCAGACAGCAGCATTCTCAGCAACTTCTTTGTGATGTTTGCATTCAAGTCACAGAATGGAACATTCCCTTTCACAGAGCAGGTTTGAAACCACCTTTTTGTAGTGTCTGTACTGGACTTTTGGAGCACTTTCCGGCCTAAGGTGAAAAAGGACATATCTTCCCATGAAAACTAGACAGAAGCATTCTCAGAAATTTACTCATGATGTGTGTCCTCAACTGACGGAGTAGAACCTTTCTTTTGATAGAGCAGTTTTG",
                // >target len=702
                "GGAGAGTATAGTTTTGAAACACTCTTTATTGTGGAGTCTGCAAGTGGATATTTGGCTGGATTTGAGGATTTCGTTGGAAACGGGATAAGGTATAAAAAGCAGACAGAAGCATTCTCAGCAATTTCTTTGTGATGTTTGCATTCAAGTCACAGAATTGAACATTCCCTTTCACAGAGCAGGTTTGAAACACTCTTTTTGTAGTGTCTGTAACTGGACTTTTGGAGCGCTTTCCGGCCTAAGGAGAAAAAGGACATATCTTCCATAAAAACTAGACAGAAGCATTGTCAGAAACTTACTCGTGATGTGTGTCTTCAACTGACGGAGTAGAACCTTTCTTTTGATAGAGCAGTTTTGAAACACTCTTTTTGTAGAATCTCCAAGTGGATATTTGGATAGCTTTGAGGATTTCGTTGGAAACGGGAATATCTTCATATAAAACCTAGACAGAAGCATTCTCAGAAACTTCCTTGTGATGGTTGCATTCAAGTCACGGAGTTGAACATTGGCTTTCATAGAGCAGGTTGGAAACACTCTTTTTTCCATTCCTGGAAGTGGACATTTGGAGCGCTTTGAGGCCTATGGTGAAAAAGGAAATATCTTCCCATAAAAACTAGACAGAAGCATTCTCAGAAACTTCTTTGTGATGTGTGTCCTCAACTGACAGAGTTGAACATGTCTTTTGAGAGAGCAGTTTTGAAACAC",
                6, 110,
                SesResults(22, 23, 5, 19, 2, 1, 2, false, PacBio::Data::Cigar("1=1D1=1X4=1X5=1D8=1I")),
                SesResults(22, 23, 5, 19, 2, 1, 2, false, PacBio::Data::Cigar("1=1D1=1X4=1X5=1D8=1I")),
                /*
                 * The max_diffs is much smaller than the given bandwidth.
                 * This tests if the memory allocation will do the right thing, otherwise AddressSanitizer will complain.
                 * This is a real test case from a human dataset.
                */
    },


    TestData{"RealHifiReads_Subsequences",
                // >mo_m64001_190914_015449/130352465/ccs:8054-16191
                "TATGAAGAGCATGTGAGCTAGTCCACTCGAGCAACTGAAGACTACACAAAAAGCCTCTTAGTCCGAGTATTTTGGGAATTTTTCACTCAAAATTAATTGCATTTTAATAGTTGCTAATATTTCAAAAAACTTGTGCACTTATTTGCACCGGCATCAGGCTGCTCAGGCTCGCCGCTTTCCCCATCCCCTTCTCCACGAAGCTCCATCAAAAAGTCCTGCTTAAAAATGCAATAACAGGAAATGGCTACGACGGCAACAAAATACAGGAAATTGAGCGCAAAGAAATATGCTATATAGGAGATGAAAAAATGAAAACGCAAACGAGACGAAGTCCAGGCAAACAAACATGGCCGTCAGACATGAACATATGTATACAGCCGCGGTTAACCAATTAGGAGTTTTGAATTCCAGAAAAAAGAATTAATACAAATATGCACTTTATTTACTAAGATTCATTTACATTTGGATTCGTTGGTTGGTGTCAGCAAGCGTACTTGAGCTTTCAGGCAGAAAAGTGCTCATGCTTTAGCATACTGAGTCTTAGAATAGCAAAGTAGCTACAAAAAATTTAAAATCAAAGTTAAGTTAAAAGTTAAAATCAAAGTTAAAAAAAAAAAACTCGAATGAACGTTTTGAAAAAGGGGAAAAGGGAAACCCCATGGAAAGCAACACAACTATTACAATTACTTTCTGTTATATGAGATTAGTATTTTTTTGAAGACACTAAACATGTGGCCCCTGCTGTACACCTGAACCATGCAAAGCGAGCATGAGGCGAAGACAGGTGGAAAGGTGGCGTTACGGGTAAGAGCATTTTGCCCGCTGGAACTTGGCCGTATGTGCTCCGGCAATTTATGTCCAAGGTGCGTGCCACCGAAATGAACCGGTGAACCTGGACGTTACGCTCGTTTCTTGGGGGAAAGTCAGGCAGTTGACCAACCTGTGTCCTTTTTCCTGGTTCTGCGCTTAACGTGTTGCCATTAAATTTTATGGTCAGGTCTTCCTCTCTTTGTCCGTGTTTGCATTTTGCTCAATTTTGTCGCCGCAGAACATAAAAAAGGGAGTCAGGAGGGGCTGAAAAAGAAGCGAAGGAATATTGCACTGTGGACGACATGTGAAAACATATTTTTATGGGACCCCAAAGGGTAATTGGACATTGCACGGGATATCCTAGCTGGCTAATGAGTGTGTCGAGCCTTTTTAGATAATTAATTGAGAGCTTTTCAAATCATCTGTGAAGTATAATGACCACACACACACAATAAGGCAACGAAACGTTGAAGTGTGATTGCACTAATGAGTAGATTAAAGTCAGCTGGAGATAACATCTTAATAAATCCAAGATAATTGAAATATCCTTTCGAGATAGTTGTATAATTTGTAATTCTCCATTCGGGAATATTTGTTTCATTGGCTTCAAGTGAAACCAATTGTGCCGCACTTCATCCAGCGAAACTTTGACGGCTAACTCGCTTAATGCGCGGCTTAAGCCAAACTACGAGTAAGATATTTTGTCTGTGGCAACTTGAAAGTTGCAAGTAATTAGTGGCACTTGAGCAACTACTAGGCTGTCGGGATACGAAAATGGGGCAGGGCATACCCAGATAGCTGGATGCTCCAAAAGTTATGTAGATCACGAAAATCAGTCCTTTGAGAACAAAGAACTTTTAGCAATCTTAAAATATTATTTCTAAGTATTCTAATCGTTATACTTTATTAGTATTTTAGTTCTATGAGTTTCTTAAAAAAAAAAAAACACTGATCTTACAAGCTCGCTTTCCACCTATCCTCTTTCAGACCCATCAGTCTAGTGCTGGAGCCAATCCTGAAGGAGCCCGGCGACATCTCCTTTACGTGGCTCGAGATTCACCGCACCAAGAACGCGCTGCTGCAGCAGTTGGACCTGCATGTGAACGCTAGCGTGACCACGGGCATGGGCTCCCAGACCTCCTCCCAGTCCGCAGCCGGCTCCAACACGGGCATTGGGGGACTGGGCACCCAGTGGGTCAGGATCTGCTGCCAATGCCGGGGGCTCTGCTGGCGGAGCAACGGCCAACGAGACGCTGAACGAGTTCGGATTCTTTCCCAAGGAGTCCGACTGCGAATACAAATGTCCTGAAATTGATGCATGCATTGCGGCCAGCTTGTGGTGCGATGGTGAGTACAGATTGCCAGTCGAGGCCACAGACCATTTCTCAGGGATTTGCCGAATGATGGGTTTTTGGGATGGAACCACTACATGCAAAGCAGTTCATAGAAGACGACGTGAACAAGAAACTTTCGAGCGTTCATCAATTCGGGTCAATGGGCAACTACTTAAGAAAAGGGCTTCTGAATAAAGAACATACATGTTGAATTATCAAGCAAATTTTGCATTTAAAGTAACTTAAGGTACGAATTCTTACTCAAAATTCAGTAGAAACGTTCTTAGGAAAATTAATACCACTAAAGAGCTAGCAATAGGCTTATCTCCTTACATTGGCATGCTTGGTCTCAATGATTGTACTTTGGGAATTTTGATATCACTATGGAGAGTCCTGGATCCCCATCTTCGCCTGGAGCTTTGCACTATTTGTACGTGATTTGGCGCTCAGTCTTGGGGACTTGGGTTACCCATTTGCAGCGTGCCGCCGGCAGAAGGCAACCAAAAATTTCACTTGTGGCTGGCTTGTCTGTTAATTTCGAAGCGGGCTGTCCACGTCCTTGTAATCAGGCAAAGTTTCCTGACTGCTGCACCAGGCTGCAGTGCACCACTCCTGCTGGTTGTCGGTTGCTAGTTGTTAGTTGCTTCGTTGGTCGGTTGTCGTTGTCGGCGGTTGCTGGTTGCTGGTTGTTGGCGATTGTTGGCCGACTGGCGGGATGGCGGGATGGCGGGATTGGGCCGGTTGTTGGACAAACCCGCAGCGGCAAAACCCCAAAATGTCTTGCAGAAAAAGTTTATCCCCCAGTCCGGCACCAGGCCAACACGAAATGGCCGCAACAGCTGTGGCAGACTTATTTTCAATGCACTTGGGAGATGAAGGGCGCAAAAAAGTAGATGCCAAGTGGGAAAGTTTATAGACTAGGAAATACCGTTTTGTTAGGAATCCGCCGTTCATTTTGTATTTACATGAACCTACCAATACATATTTAAAAGTCATATATTATCCTCCAAGTGAATGTAATTATTTTAATTTTGCTGCTGAAATAGTAAATACTTTAATATGTACATTGAAGATTGTGTGATATACAAGTATAAAATTGCATTGATTTTCCTGCATTTTCCTACTTATACTCCAGAAGCCGCACAGACTCTTTAAATCTTACGCTCGGTGGAAAAAATACAAAGGCTAGCCTGGGCTTTGCCTTTTAGCTTCAATTATTTATGGGCAGCCAGGCGCGTCTGCATCTGGAGCTGGTTAACCGCAAAATCCAGGGCATTAAACAAGCTTTGGACCCCACGGGCCACTCGAGTGACACTGCATTCAATTTGTAACCTTAGAATAGGCCACCCAGCTGCGCTGCTCCTCCCAAACCCCAATCTGCTTTCAATTTTTCAATCAGCGCCTGTACCTCTATCTAGGATACTGAATCTGTTCCTACACTGACAAAAAGATTCTTTTCACATATTACCTTCAAAATATCTATTTGAAACTATATGCAAAATAAAAAAGATACATCCAACTAGGGATTGTATTGTTGAAACTTCTTAAATATATATTGCTTGTCACTAAGACTATTTATGGATAATTTGAAAGAGAGTCTTATAGACCTGCTTAGACCTTAAAAATATTTTTAAGCGGTGTTTTGGGGTAATGAACAAATATCTGAATCAAAAGCTGTAACAGATTAATTCCACCTAAGAAGTGACCGCTTTTTTCTAAAAGTGTCCTTATAATCTCTCCTTTTCAGGCCACCACAACTGCCCCAGCGGTTTCGATGAGTCGGAGGAGGAGTGCGGCACCGCTCGCAAGTTGCTGGAACTGCCAGGCGGAGTTTTCGCCGCCCTTGGATGCATTGCAGCAGCCCTGACCGCCTGCCTGATATTCTGCATGTTCGGCCTGATGAGGAAGCGGAAAAAGTCGGTGGTGCAGAAGGGCGCCGGAGTGGGCGGCATGGGCGGAATGGGAGTGGGTGGGATGGGCGTGGGCATCGGGTTCATGAACGGAGCGAGCAACGGCAACGTGTCCGCCGCCAGCATCGGTACCTTGAAGAAGGACTTCAAGAAGGAGGCGCTGTACATCGATCCGGCCTCCTGACACGGACTACAGCCGGTCGAATATATCCACGTGGACCGGGAGACCACCGTGTGATATGTTGCCGCCGGCTGTGGCATTCCCTATAATCATGTGGAGCTTTCCATGGAGCTGGCCTGACCGATGGATCAGTCCGCTGATCCCGTTCAAAATCTGGCCTGGAGATTTTAAAGATACAGAGTCGACGCCTGGAGGTATGAGAACACTAAATGAACCTACATTTTAAGCCGATAACTTTTTCTAAACGCTTCTGGGAGTAACAAATTGCGCACTAGCGAAGTTGTGTGTATGTAGTCATAAGGAAAATACTGTTATTAAGTTGCGTATGGTTTTTAGTTTGTATATATAATTTAATTTTAGGCATAAATAATTTGAATATCGTTGAAGGAGTCAGAACTGTGCACAGAACTTTGTGCAAGAACTTTTTGAATGTAGACTTTGTTTCCTATCATTAATCTTTCAATATCAGAAAATTATCATCGAGTTCAAACTAAAAATCTTTATAATTTCATGCATTAAAGGTTTCCAAATTAAAATCAAAATACATAAACTTTTCACAGTAAATAAAATTTTAAAATTCACGTCTTAAGGAATGTGTATATTAAGCAATAAAATCACTCACAACTTTTTAAAACTCAAAATGTTCGAGTAATATTAATTTATTGCTATTCACATGGTTGCCCGTGCTGACCCTTGCATAATTGTTTAGTTCGATTATGAGAGTTATCAATTTTTCACTTTAAGACTACATACACAAACGAAAAGAGGACATAAAACTATTGTAAATTAGGTACCTTTAAACCAGTTAGATCAGAGAACAACAATTCGTGGGCGGTTCGTAGATGACTTCGGCATATTACCACCCACATAGCACACAGCACCCATTAAAGGAAGTAAACCAATGACATTTAGCAACTAAATGAAAGATCATAAACTACAATGAAATTGTATAGAAGAGCAACTCAAAGAGCAGATCCCTATGCGCAGGGCCACTCTGGCACTAATTTTAAGCGGTTAGACTAGGTAAAAAACTGTTAGCAGCGATTGGAAAGGAGTTGAGGAAGTTCGAGGGTAATCTCTACAGTGTTGTGATCTAAAGCATGCACATACTGAACTATATTTAAACACTATGTTTAACGTCGAAGCATTTACATGACGAGTACTCAGAACCGGGGGCGTATATTGTATAGGTACTAGCTAATAGCTAAATAGTCAACTAAGAGGCAGCACTTTTCTGTTATTATCGTACAGCACAAACACATAAGCAAATGACATCACCTAGTGGTTTATTATATACAAGTCGATACAGCAAGTTAATTAAATGCAATTGTTATTAATGAATTAAGTAGAGCTCCAAACCAGAGATATGTTCACTGCACTTAATCGAGAGAACTGCCAAAGGATTAGACAACAAAACAACAAAACAACAAAACAATCAATGTAAATATTTCTGTCGAGCAATAACTTGTTTTGCGTTAGCATTTATGTTCGTAATTTGAACACAAACCGATTCCACTTAAGTAGCTGCAGTGAGTAATATTCCTAAAGTCAATAGACTTTGTAAAAACTACGAGCGAAGAGTAATAATAAACTGTATGCAAAACTACCATTTAAGACACGGCTCATAACCCATATAAACAATCAAACTTTTGGCAGTTGCAGTAGGATTAAGCAAACAACACCAAGAAGAAGCGGAAAAAGGACATAAGTACTTTATTTTATAGACACTGAACTAGGGCAAAAATTATGAAACTGAATACTGAAAGCAAAGCAAATACTGGCTGCAGTGGGAATTTTAAAAATTTAAGAACCTCTGCCAACGCCGGCAGAGTCCACATCCATAAAAACCCATTAAAATTTCAACTGCAATGAAACTTTGCCGAAAAGCGTGTGTATAAAAAATATATGATTTAGCCAGCATAAAAAGTAAATGAAAATATATACATAAATTAAATCAAATATTTTCGGTTGGCATTTTTTGATAGCAAATAACCTAAATATTTATAGCTAAAAAAGAAAGAAAGAAAAATTACCAAACGCACTTTGTACATTTTATAAAAAACGAGCAGTTTATATTATTTCCAGACACACAGCATAGCATACACATACATATATGATATTTTACACACATTCGAATAAGAGATAAAATAATACCGAAACATATACGAGAGCGAAAAGAAATTAAATTAATTAAAAAAGGCGCATAGGCCGAAAATGCATTTTTACTGCCGCCAATTGAGGACAAAACTAAAATCGAAAACACGGCCGGCAGGGCATTTGAATAGTTAGCTTTAAGCCCGGCTAAAAATGCTGTAACTACTTTGTTGGTCACCATTGTCACTCACACAAACTCACATCTCAATCATCAGTTCGGTTCAGTTCCTTCTCTCAACCAAACCAAACTATGGAAATTTCTGTGCCAACTGTGTTTAAAGTAGAATTCAAAGATAAAATACGAGATAAAAGAGGGGAGGCCGGCGGATGGAGAAATTGAATTATGACCGCACACGAAACATTGCATATAATTTATTTGAAATGTTCAAAAATAAAGTAGACGAAACGAGAAGAAAACGCTGCGCGATGGTAAATAACAATATAATACAGACAGCTGTGTAAATAGTCAAACTGCGTACAGTCTACGCATATTTAATTAAAACTCGAAAAAAACTGATAAAACTGAAAAAACTGAAAAAAGTGACCCAATGAGATGGAATTTGAGAAAGGAGATCATGTCTTAAGTTCCATTTCGGAATTCACCCGATTTTCATTTCCTTTCACCGTTCGATAGAGTAAGTTTCACTTGTATGATAAATAAAACTAAACCAAATATTTATATAGAGATATTATACGTATATGAAAAGGAAACCTTAAGTTAATGACGCACTGAATCAAATAGATGAAACGAGTATTTAAAACCAAGACAATTAAAACCAAATAGGAACTTTAAGCAAAATAAAACCGATGAAAAAATGAAACTGAACACAGAGCACCTTTCTCTATTAGCCAAAATTGCAAAAGTGGGTGGTTGGGCGTTTGATTTTCGGGTGGAATAGATGGGGTGTTTAGTGGGTGGAATGCCAACTGTTCCGTTTGACGCTCGGGGGAAAACCGTTGGCCAAACTAAGGCGAACTCAAATCGTGTCCAAATAGCTGGCAGTCAAGTAGGCCACAGATTTTGATTGCCCGCATCCGGGCAGCAATCGGTGGACATCCTCGAAAATTCGGAACATCGAACGAGATGCTGTTAAAGCCAAGCCAAGCCGGTAGCTCCAGGCGTGTGTTCTCACCGGTCGGAACTCGGATTGGCCGTGGGATCGGGTTGTGGCTATGGGCTTTCCATGTCCCCCGCCTCGCCCACGCTGAGCCCGCATAAAACATAAGCCACCGGCCTGGTGCTCGTGTGGTGTGTGCTCTATATCCCACAGATTCCGCCTCCTCCCTCGGAATGCCCAGTTCTCTGACTCCGTGTTTATGGCCGGCAGCTCTGTTTCTGGACCTGGCAACTCTGGCCCCTGGACTTGGCCAAAATGGCTTTGATGGGCGTGGGCGACAAAAGGACCGGAAATCGCTAACGAGGCGTGGTAAAGGAGTTTTTGCAGTAAGTTTGAATATCAAGGTAGATACTCAAAAAATTGAATAACAAAGCATAATCAAATTTAATGGAGGGATCTACCTTGAATGTATTTTGTTATTATTGTATTACAGCTTATCAAACTAATCACTAGTATGTTTGTTGCTTTTATTTATTATTAACACATTTAAATATCTGGCTATTAAAGTTAATAATTTTGCAGTTCTTTTTATGCTGTGGTCATTAATATACATATAATACATATATCTACTAAAAGTATTCACTCTAGATTAATGCAAATTGATGAAATTATTTATTTATCGTACGCACTT",
                // >fa_m64001_190914_015449/101320295/ccs:10320-18430 reversed
                "TATGAAGAGCATGTGAGCTAGTCCACTCGAGCAACTGAAGACTACACAAAAAGCCTCTTAGGTCGAGTATTTTGGGAATTTTTCACTCAAAATTAATTGCATTTTAATAGTTGCTAATATTTCAAAAAACTTGTGCACTTATTTGCACCGGCATCAGGCTGCTCAGGCTCGCCGCTTTCCCCATCCCCTTCTCCATGAAGCTCCATCAAAAAGTCCTGCTTAAAAATGCAATAACAGGAAATGGCTACGACGGCAACAAAATACAGGAAATTGAGCGCAAAGAAATATGCTATATAGGAGATGAAAAAATGAAAACGCAAACGAGACGAAGTCCAGGCAAACAAACATGGCCGTCAGACATGAACATATGTATACAGCCGCGGTTAACCAATTAGGAGTTTTGAATTCCAGAAAAAAGAATTAATACAAATATGCACTTTATTTACTAAGATTCATTTACATTTGGATTCGTTGGTTGGTTTCAGCAAGCGTACTTGAGCTTTCAGGCAGAAAAGTGCTCATGCTTTAGCATACTGAGTCTTAGAATAGCAAAGTAGCTACAAAAAATTTAAAATCAAAGTTAAGTTAAAAGTTAAAATCAAAGTTAAAAAAAAAAAACTCGAATGAACGTTTTGAAAAAGGGGAAAAGGGAAACCCCATGGAAAGCAACACAACTATTACAATTACTTTCTGTTATATGAGATTAGTATTTTTTTTGAAGACACTAAACATGTGGCCCCTGCTGTACACCTGAACCATGCAAAGCGAGCATGAGGCGAAGACAGGTGGAAAGGTGGCGTTACGGGTAAGAGCATTTTGCCCGCTGGAACTTGGCCGTATGTGCTCCGGCAATTTATGTCCAAGGTGCGTGCCACCGAAATGAACCGGTGAACCTGGACGTTACGCTCGTTTCTTGGGGGAAAGTCCGGCAGTTGACCAACCTGTGTCCTCTTTCCTGGTTCTGCGCTTAACGTGTTGCCATTAAATTTTATGGTCAGGTCTTCCTCTCTTTGTCCGTGTTTGCATTTTGCTCAATTTTGTCGCCGCAGAACATAAAAAAGGGAGTCAGGATGGGGCTGAAAAAGAAGCGAAGGAATATTGCACTGTGGACGACATGTGAAAACATATTTTTATGGGACCCCAAAGGGTAATTGGACATTGCACGGGATATCCTAGCTGGCTAATGAGTGTGTCGAGCCTTTTTAGATAATTAATTGAGAGCTTTTCAAATCATCTGTGAAGTATAATGACCACACACACACAATAAGGCAACGAAACGTTGAAGTGTGATTGCACTAATGAGTAGATTAAAGTCAGCTGGAGATAACATCTTAATAAATCCAAGATAATTGAAATATCCTTTCGAGATAGTTGTATAATTTGTAATTCTCCATTCGGGAATATTTGTTTCATTGGCTTCAAGTGAAACCAATTGTGCCGCACTTCATCCAGCGAAACTTTGACGGCTAACTCGCTTAATGCGCGGCTTAAGCCAAACTACGAGTAAGATATTTTGTCTGTGGCAACTTGAAAGTTGCAAGTAATTAGTGGCACTTGAGCAACTACTAGGCTGTCGGGATACGAAAATGGGGCAGGGCATACCCAGATAGCTGGATGCTCCAAAAGTTATGTAGATCACGAAAATCAGTCCCTTGAGAACAAATAACTTTTAGCAATCTTAAAATATTATTTCTAAGTATTCTAATCGTTATACTTTATTAGTATTTTAGTTCTATGACTTTCTTTAAAAAAAAAAAAAAAAACACTGATCTTACAAGCTCGCTTTCCACCTATCCTCTTTCAGACCCATCAGCCTAGTGCTGGAGCCAATCCTGAAGGAGCCCGGCGACATCTCCTTTACGTGGCTCGAGATTCACCGCACCAAGAACGCGCTGCTGCAGCAGTTGGACCTGCATGTGAACGCTAGCGTGACCACGGGCATGGGCTCCCAGACCTCCTCCCAGTCCGCAGCCGGCTCCAACACGGGCATTGGGGGACTGGGCACCAGTGGGTCAGGATCTGCTGCCAATGCCGGGGCTCTGCTGGCGGAGCAACGGCCAACGAGACGCTGAACGAGTTCGGATTCTTTCCAAGGAGTCCGACTGCGAATACAAATGTCCTGAAATTGATGCATGCATTGCGGCCAGCTTGTGGTGCGATGGTGAGTACAGATTGCCAGTCGAGGCCACAGACCATTTCTCAGGGATTTGCCGAATGATGGGTTTTTGGGATGGAACCACTACATGCAAAGCAGTGCATAGAAGACGACGTGAACAAGAAACTTTCGAGCGTTCATCAATTCGGGTCAATGGGCAACTACTTAAGAAAAGGGCTTCTGAATAAAGAACATACATGTTGAATTATCAAGCAAATTTTGCATTTAAAGTAACTTAAGGTACGAATTCTTACTCAAAATTCAGTAGAAACGTTCTTAGGAAAATTAATACACTAAAGAGCTAGCAATAGGCTTATCTCCTTACATTGGCATGCTTGGTCTCAATGATTGTACTTTGGGAATTTTGATATCACTATGGAGAGTCCTGGATCCCCATCTTCGCCTGGAGCTTTGCACTATTTGTACGTGATTTGGCGCTCAGTCTTGGGGATTTGGGTTACCCATTTGCAGCGTGCCGCCGGCAGAAGGCAACCAAAAATTTCACTTGTGGCTGGCTTGTCTGTTAATTTCGAAGCGGGCTGTCCACGTCCTTGTAATCAGGCAAAGTTTCCTGACTGCTGCACCAGGCTGCAGTGCACCACTCCTGCTGGTTGTCGGTTGCTAGTTGTTAGTTGCTTCGTTGGTCGGTTGTCGTTGTCGGCGGTTGCTGGTTGCTGGTTGTTGGCGATTGTTGGCCGACTGGCGGGATGGCGGGATGGCGGGATTGGGCCGGTTGTTGGACAAACCCGCAGCGGCAAAACCCCAAAATGTCTTGCAGAAAAAGTTTATCCCCCAGTCCGGCACCAGGCCAACACGAAATGGCCGCAACAGCTGTGGCAGACTTATTTTCAATGCACTTGGGAGATGAAGGGCGCAAAAAAGTAGATGCCAAGTGGGAAAGTTTATAGACTAGGAAATACCGTTTTGTTAGGAATCCGCCGTTCATTTTGTATTTACATGAACCTACCAATACATATTTAAAAGTCATATATTATCCTCCAAGTGAATGTAATTATTTTAAATTTGCTGCTGGAATAGTAAATACTTTAATATGTACATTGAAGATTGTGTGATATACAAGTATAAAATTGCATTGATTTTCCTGCATTTTCCTACTTATACTCCAGAAGCCGCACAGACTCTTTAAATCTTACGCTCGGTGGAAAAAATACAAAGGCTAGCCTGGGCTTTGCCTTTTAGCTTCAATTATTTATGGGCAGCCAGGCGCGTCTGCATCTGGAGCTGGTTAACCGCAAAATCCAGGGCATCAAACAAGCTTTGGACCCCACGGGCCACTCGAGTGACACTGCATTCAATTTGTAACCTTAGAATAGGCCACCAGCTGCGCTGCTCCTCCCAAACCCCAATCTGCTTTCAATTTTTCAATCAGCGCCTGTACCTCTATCTAGGATACTGAATCTGTTCCTACACTGACAAAAAGATTCTTTTCAAATATTACCTTCAAAATATCTATTTGAAACTATATGCAAAATAAAAAAGATACATCCAACTAGGGATTGTATTGTTGAAACTTCTTAAATATATATTGCTTGTCACTAAGACTATTTATGGATAATTTGAAAGAGATTCTTATAGACCTGCTTAGACCTTAAAAATATTTTTAAGCGGTGTTTTGGGGTAATGAACAATATCTGAATCAAAAGCTGTAACAGATAATTCCACCTAAGAAGTGACCGCTTTTTTCTAAAAGTGTCCTTATAATCTCTCCTTTTCAGGCCACCACAACTGCCCCAGTGGTTTTCGATGAGTCGGAGGAGGAGTGCGGCACCGCTCGCAAGTTGCTGGAACTGCCAGGCGGAGTTTTCGCCGCCCTTGGATGCATTGCAGCAGCCCTGACCGCCTGCCTGATATTCTGCATGTTCGGCCTGATGAGGAAGCGGAAAAAGTCGGTGGTGCAGAAGGGCGCCGGAGTGGGCGGCATGGGCGGAATGGGAGTGGGTGGGATGGGCGTGGGCATCGGGTTCATGAACGGAGCGAGCAACGGCAACGTGTCCGCCGCCAGCATCGGTACCTTGAAGAAGGACTTCAAGAAGGAGGCGCTGTACATCGATCCGGCCTCCTGACACGGACTACAGCCGGTCGAATATATCCACGTGGACCGGGAGACCACCGTGTGATATGTTGCCGCCGGCTGTGGCATTCCCTATGATCATGTGGAGCTTTCCATGGAGCTGGCCTGACCGATGGATCAGTCCGCTGATCCCGTTCAAAATCTGGCCTGGAGATTTTAAAGATACAGAGTCGACGCCTGGAGGTATGAGAACACTAAATGAACCTACATTTTAAGCCGATAACTTTTTCTAAACGCTTCTGGGAGTAACAAATTGCGCACTAGCGAAGTTGTGTGTATGTAGTCATAAGGAAAATACTGTTATTAAGTTGCGTATGGTTTTTAGTTTGTATATATAATTTAATTTTAGGCATAAATAATTTGAAAATCGTTGAAGGAGTCAGAACTGTGCACAGAACTTTGTGCAAGAACTTTTTGAATGTAGACTTTGTTTCCTATCATTAATCTTTCAATATCAGAAAATTATCTTGGAGTTCAAACTAAAAATCTTTATAATTTCATGCATTAAAGGTTCCAAATTAAAATCAAAATACATAAACTTTTCACAGTAAATAAAATTTTAAAATTCACGTCTTAAGGAATGTGTATATTAAGCAATAAAATCACACACAACTTTTTAAAACTCAAAATGTTCGAGTAATATTAATTTATTGCTAGTCACATGGTTGCCCGTGCTGACCCTTGCATAATTGTTTTAGTTCGATTATGAGAGTTAACAATTTTTCACTTTAAGACTACATACACAAACGAAAAGAGGACATAAAACTATTGTAAATTAGGTACCTTTAAACCAGTTAGATCAGCTAGAGAACAACAATTCGTGGGCGGATCGTAGATGACTTCGGCATATTACCACATAGCACACAGCACCCATTAAAGAAAGTAAACCAATGACATTTAGCAACTAAATGAAAGATCATAAACTACAATGAAATTGTATAGAAGAGCAACTCAAAGAGCAGATCCCTATGCGCAGGGCCACTCTGGCACTAATTTTAAGCGGTTAGACTAGGTAAAAAACTGTTAGCAGCGATTGGAAAGGAGTTGAGGAAGTTCGAGGGTAATCTCTACAGTGTTGTGATCTAAGGCATGCACATACTGAACTATGTTTAAACACTATGTTTAACGTCGAAGCATTTACATGACGAGTACTCAGAACCGGGGGCGTATATTGTATAGGTACTAGCTAGTAGCTAAATAGTCAACTAAGAGGCAGCACTTTTCTGTTATTATCGTACAGCACAAACAAGCAAATGACATCACCTAGTGGTTTATTATATACAAGTCGATACAGCTAGTTAATTAAATGCAATTGTTATTAATGAATTAAGTAGAGCTCCAAACCAGAGATATGCACTTAATCGAGAGAACTGCCAAAGGATTAGACAACAAAACAACAAAACAATCAATGTAAATATTTCTGTCGAGCAATAACTTGTTTTGCGTTAGCATTTATGTTCGTAATTTGAACACAAACCGATTCCACTTAAGTAGCTGCAGTGAGTAATATTCCTAAAGTCAATAGACTTTGTAAAAACTACGAGCGAAGAATAATAATAAACTGTATGCAAAACTACCATTTAAGACACGGCTCATAACCCATATAAACAATCAAACTTTTGGCAGTTGCAGTAGGATTAAGCAAACAACACCAAGAAGAAGCGGAAAAGGACATAAGTACTTTATTTTATAGACACTGAACTAGGGCAAAAATTATGAAACTGAATACTGAAAGCAAAGCAAATACTGGCTGCAGTGGGAATTTTAAAAATTTAAGAACCTCTGCCACGCCGGCAGAGTCCACATCCATAAAAACCCATTAAATTTCAACTGCAATGAAACTTTGCCGAAAAGCGTGTGTATAAAAAATATATGATTTAGCCAGCATAAAAGTAAATGAAAATATATACATAAATTAAATCAAATATTTTCGGTTGGCATTTTTTGATAGCAAATAACCTAAATATTTATAGCTAAAAAAGAAAGAAAGAAAAATTACCAAAACGCACTTTGTACAATTTAATAAAAAACGAGCAGTTTATATTATTTCCAGACACACAGCATAGCATACACATACATATATGATATTTTACACACATTCGAATAAGAGATAAAATAATACCGAAACATATACGAGAGCGAAAAGAAATTAAATTAATTAAAAAGGGCGCATAGGCCGAAAATGCATTTTTACTGCCGCCAATTGAGGACAAAACTAAAAATCGAAACACGGCCGGCAGGGCATTTGAATAGTTAGCTTTAAGCCCGGCTAAAAATGCTGTAACTACTTTGTTGGCCACCATTGTCACTCACACAAACTCACATCTCAATCATCAGTTCGGTTCAGTTCCTTCTCTCAACCAAACCAAACTATGGAAATTTCTGTGCCAACTGTGTTTAAAGTAGAATTCAAAGATAAAATACGAGATAAAAGAGGGGGGGCCGGCGGATGGAGAAATTGAATTATGACCGCACACGAAACATTGCATATAATTTATTTGAAATGTTCAAAAATAAAGTAGACGAAACGAGAAGAAAACGCTGCGCGATGGTAAATAACAATATAATACAGACAGCTGTGTAAATAGTCAAACTGCGTACAGTCTACGCATATTTAATTAAAAACCGAAAAAAAAAACTGATAAAACTGAAAAAAGTGACGCAATGAGATGGAATTTGAGAAAGGAGATCATGTCTTAAGTTCCATTTCGGAATTCACCCGATTTTCATTTCCTTTCACCGTTCGATAGAGTAAGTTTCACTTGTATGATAAATAAAACTAAACCAAATATTTATATAGAGATATTATACGTATATGAAAAGGAAACCTTAAGTTAATGACGCACTGAATCAAATAGATGAAACGAGTATTTAAAACCAAGACAATTAAAACCAAAAGGAACTTTAAGCAAAATAAAACCGATGAAAAAATTAAACTGAACACAGAGCACCTTTCTCTATTAGCCAAAATTGCAAAAGTGGGTGGTTGGGCGTTTGATTTTCGGGTGGAATAGATGGGGTGTTTAGTGGGTGGAATGCCAACTGTTCCGTTTGACGCTCGGGGGAAAACCGTTGGCCAAACTAAGGCGAACTCAAATCGTGTCCAAATAGCTGGCAGTCAAGTAGGCCACAGATTTTGATTGCCCGCATCCGGGCAGCAATCGGTGGACATCCTCGAAAATTCGGAACATCGAACGAGATGCTGTTAAAGCCAAGCCAAGCCGGTAGCTCCAGGCGTGTGTTCTCACCGGTCGGAACTCGGATTGGCCGTGGGATCGGGTTGTGGCTATGGGCTTTCCATGTCCCCCGCCTCGCCCACGCTGAGCCCGCATAAAACATAAGCCACCGGCCTGGTGCTCGTGTGGTGTGTGCTCTATATCCCACAGATTTCGCCTCCTCCCTCGGATTGCCCAGTTCTCTGACTCCGTGTTTATGGCCGGCAGCTCTGTTTCTGGACCTGGCAACTCTGGCCCGTGGACTTGGCCAAAATGGCTTTGATGGGCGTGGGCGACAAAAGGACCGGAAATCGCTAACGAGGCGTGGTAAAGCAGTTTTTGCAGTAAGTTTGAATATCAAGGTAGATACTCAAAAAATTGAATAACAAAGCATAATCAAATTTAATGGAGGGATCTACCTTGAATGTATTTTGTTATTATTGTATTACAGCTTATCAAACTAATCACTAGTATGTTTGTTGCTTTTATTTATTATTAACACATTTAAATATCTGGCTATTAAAGTTAATAATTTTGCAGTTCTTTTTATGCCGTGGTCATTAATATACATATAATACATATATCTACTAAAAGTATTCACTCTAGATTAATGCAAATTGATGAAATTATTTATTTATCGTACGCACTT",
                // Alen = 8138, Blen = 8111
                static_cast<int32_t>(0.03 * 8138),
                static_cast<int32_t>(0.01 * 8111),
                SesResults(8138, 8111, 105, 8051, 42, 45, 18, true, PacBio::Data::Cigar("61=1D2=1I131=1X284=1X235=1D209=1X23=1X120=1D579=1X11=1X74=1X6=1D13=4D50=1X192=1I31=1I54=1I164=1X191=1I159=1X569=1X10=1X235=1X71=1I111=1X135=1X61=1I26=1I80=1X5=1D406=1X288=1X101=1X1=1X43=1I93=1X49=1X38=1D19=1X87=2D3=1D1=1D17=1X27=2I1=2I22=1X206=1X20=1X81=1X58=1I1=3I46=1X58=2I4=4I1=1I46=8I145=1X118=1I119=1I36=1I68=1I111=1D13=1X4=1D142=1X55=1D6=1I68=1X142=1X184=1D1=1I9=3I1=1I2=2I4=1X18=1X235=1I35=1X456=1X16=1X66=1X74=1X218=1X96=")),
                SesResults(8138, 8111, 105, 8051, 42, 45, 18, true, PacBio::Data::Cigar("61=1D2=1I131=1X284=1X235=1D209=1X23=1X120=1D579=1X11=1X74=1X6=1D13=4D50=1X192=1I31=1I54=1I164=1X191=1I159=1X569=1X10=1X235=1X71=1I111=1X135=1X61=1I26=1I80=1X5=1D406=1X288=1X101=1X1=1X43=1I93=1X49=1X38=1D19=1X87=2D3=1D1=1D17=1X27=2I1=2I22=1X206=1X20=1X81=1X58=1I1=3I46=1X58=2I4=4I1=1I46=8I145=1X118=1I119=1I36=1I68=1I111=1D13=1X4=1D142=1X55=1D6=1I68=1X142=1X184=1D1=1I9=3I1=1I2=2I4=1X18=1X235=1I35=1X456=1X16=1X66=1X74=1X218=1X96=")),
                /*
                * This is an actual real set of sequences. Only a subportion which aligns end-to-end
                * was extracted for the unit test here, so that the global aligner can be tested without mapping.
                *
                * NOTE: I compared this alignment with Edlib, and Edlib produces the identical CIGAR string!
                */
    },
};
// clang-format on

TEST(SES2AlignBanded_Global, AllTests)
{
    for (const auto& data : testDataGlobal) {
        // Global alignment with traceback.
        {
            // Name the test.
            SCOPED_TRACE("Global-" + data.testName);

            // Run.
            SesResults result =
                PacBio::Pancake::SES2AlignBanded<PacBio::Pancake::SESAlignMode::Global,
                                                 PacBio::Pancake::SESTrimmingMode::Disabled,
                                                 PacBio::Pancake::SESTracebackMode::Enabled>(
                    data.query, data.target, data.maxDiffs, data.bandwidth);

            // Evaluate.
            EXPECT_EQ(data.expectedGlobal, result);
        }
        // Semiglobal alignment with traceback.
        {
            // Name the test.
            SCOPED_TRACE("Semiglobal-" + data.testName);

            // Run.
            SesResults result = SES2AlignBanded<SESAlignMode::Semiglobal, SESTrimmingMode::Disabled,
                                                SESTracebackMode::Enabled>(
                data.query, data.target, data.maxDiffs, data.bandwidth);

            // Evaluate.
            EXPECT_EQ(data.expectedSemiglobal, result);
        }
    }
}
}  // namespace Tests
}  // namespace Alignment
}  // namespace Pancake
}  // namespace PacBio
