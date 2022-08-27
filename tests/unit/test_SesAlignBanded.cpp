// Authors: Ivan Sovic

#include <PancakeTestData.h>
#include <gtest/gtest.h>
#include <pancake/SesAlignBanded.hpp>
#include <sstream>
#include <tuple>

TEST(SESAlignBanded_Global, NormalSmallCase_SimpleSeqWithAllDiffs)
{
    /*
     * Note: Mismatches count as 2 diffs.
    */
    // Inputs.
    std::string query = "GGATCAGTT";   // GGAT-CAGTT
                                       // |X|| | |||
    std::string target = "GAATTCGTT";  // GAATTC-GTT
    int32_t maxDiffs = 6;
    int32_t bandwidth = 30;
    std::string cigar = "1=1I1=1D1=1D1=1I3=";

    /*
     * GGA-T-CAGTT
     * | | | | |||
     * G-AATTC-GTT
    */

    // Order of elements in SesResult: lastQueryPos, lastTargetPos, diffs, minK, maxK, valid
    PacBio::Pancake::SesResults expected(9, 9, 4, 7, 0, 2, 2, true, PacBio::Data::Cigar(cigar));

    // Run.
    PacBio::Pancake::SesResults result =
        PacBio::Pancake::SESAlignBanded<PacBio::Pancake::SESAlignMode::Global,
                                        PacBio::Pancake::SESTracebackMode::Enabled>(
            query, target, maxDiffs, bandwidth);

    std::cerr << "CIGAR: " << result.cigar.ToStdString() << "\n";

    // Evaluate.
    EXPECT_EQ(expected, result);
}

TEST(SESAlignBanded_Global, NormalSmallCase_AnotherSimpleSeqWithAllDiffs)
{
    /*
     * Note: Mismatches count as 2 diffs.
    */
    // Inputs.
    std::string query = "GGATTCAGTT";   // GGA-TTCAGTT
                                        // | | |||X|||
    std::string target = "GAATTCCGTT";  // G-AATTCCGTT
    int32_t maxDiffs = 15;
    int32_t bandwidth = 30;
    std::string cigar = "1=1I1=1D3=1X3=";

    /*
     * GGA-T-CAGTT
     * | | | | |||
     * G-AATTC-GTT
    */

    // Order of elements in SesResult: lastQueryPos, lastTargetPos, diffs, minK, maxK, valid
    PacBio::Pancake::SesResults expected(10, 10, 3, 8, 1, 1, 1, true, PacBio::Data::Cigar(cigar));

    // Run.
    PacBio::Pancake::SesResults result =
        PacBio::Pancake::SESAlignBanded<PacBio::Pancake::SESAlignMode::Global,
                                        PacBio::Pancake::SESTracebackMode::Enabled>(
            query, target, maxDiffs, bandwidth);

    std::cerr << "CIGAR: " << result.cigar.ToStdString() << "\n";

    // Evaluate.
    EXPECT_EQ(expected, result);
}

TEST(SESAlignBanded_Global, NormalSmallCase_CIGAR_5D_in_suffix)
{
    /*
     * Note: Mismatches count as 2 diffs.
    */
    // Inputs.
    std::string query = "ACTG";
    std::string target = "ACTGAAAAA";
    int32_t maxDiffs = 15;
    int32_t bandwidth = 30;
    std::string cigar = "4=5D";

    // Order of elements in SesResult: lastQueryPos, lastTargetPos, diffs, minK, maxK, valid
    PacBio::Pancake::SesResults expected(query.size(), target.size(), 5, 4, 0, 0, 5, true,
                                         PacBio::Data::Cigar(cigar));

    // Run.
    PacBio::Pancake::SesResults result =
        PacBio::Pancake::SESAlignBanded<PacBio::Pancake::SESAlignMode::Global,
                                        PacBio::Pancake::SESTracebackMode::Enabled>(
            query, target, maxDiffs, bandwidth);

    std::cerr << "CIGAR: " << result.cigar.ToStdString() << "\n";

    // Evaluate.
    EXPECT_EQ(expected, result);
}

TEST(SESAlignBanded_Global, NormalSmallCase_CIGAR_5I_in_suffix)
{
    /*
     * Note: Mismatches count as 2 diffs.
    */
    // Inputs.
    std::string query = "ACTGAAAAA";
    std::string target = "ACTG";
    int32_t maxDiffs = 15;
    int32_t bandwidth = 30;
    std::string cigar = "4=5I";

    // Order of elements in SesResult: lastQueryPos, lastTargetPos, diffs, minK, maxK, valid
    PacBio::Pancake::SesResults expected(query.size(), target.size(), 5, 4, 0, 5, 0, true,
                                         PacBio::Data::Cigar(cigar));

    // Run.
    PacBio::Pancake::SesResults result =
        PacBio::Pancake::SESAlignBanded<PacBio::Pancake::SESAlignMode::Global,
                                        PacBio::Pancake::SESTracebackMode::Enabled>(
            query, target, maxDiffs, bandwidth);

    std::cerr << "CIGAR: " << result.cigar.ToStdString() << "\n";

    // Evaluate.
    EXPECT_EQ(expected, result);
}

TEST(SESAlignBanded_Global, NormalSmallCase_CIGAR_5D_in_prefix)
{
    /*
     * Note: Mismatches count as 2 diffs.
    */
    // Inputs.
    std::string query = "ACTG";
    std::string target = "CCCCCACTG";
    int32_t maxDiffs = 15;
    int32_t bandwidth = 30;
    std::string cigar = "5D4=";

    // Order of elements in SesResult: lastQueryPos, lastTargetPos, diffs, minK, maxK, valid
    PacBio::Pancake::SesResults expected(query.size(), target.size(), 5, 4, 0, 0, 5, true,
                                         PacBio::Data::Cigar(cigar));

    // Run.
    PacBio::Pancake::SesResults result =
        PacBio::Pancake::SESAlignBanded<PacBio::Pancake::SESAlignMode::Global,
                                        PacBio::Pancake::SESTracebackMode::Enabled>(
            query, target, maxDiffs, bandwidth);

    std::cerr << "CIGAR: " << result.cigar.ToStdString() << "\n";

    // Evaluate.
    EXPECT_EQ(expected, result);
}

TEST(SESAlignBanded_Global, NormalSmallCase_CIGAR_5I_in_prefix)
{
    /*
     * Note: Mismatches count as 2 diffs.
    */
    // Inputs.
    std::string query = "CCCCCACTG";
    std::string target = "ACTG";
    int32_t maxDiffs = 15;
    int32_t bandwidth = 30;
    std::string cigar = "5I4=";

    // Order of elements in SesResult: lastQueryPos, lastTargetPos, diffs, minK, maxK, valid
    PacBio::Pancake::SesResults expected(query.size(), target.size(), 5, 4, 0, 5, 0, true,
                                         PacBio::Data::Cigar(cigar));

    // Run.
    PacBio::Pancake::SesResults result =
        PacBio::Pancake::SESAlignBanded<PacBio::Pancake::SESAlignMode::Global,
                                        PacBio::Pancake::SESTracebackMode::Enabled>(
            query, target, maxDiffs, bandwidth);

    std::cerr << "CIGAR: " << result.cigar.ToStdString() << "\n";

    // Evaluate.
    EXPECT_EQ(expected, result);
}

TEST(SESAlignBanded_Global, NormalSmallCase_SingleMatch)
{
    /*
     * Note: Mismatches count as 2 diffs.
    */
    // Inputs.
    std::string query = "A";
    std::string target = "A";
    int32_t maxDiffs = 5;
    int32_t bandwidth = 30;
    std::string cigar = "1=";

    // Order of elements in SesResult: lastQueryPos, lastTargetPos, diffs, minK, maxK, valid
    PacBio::Pancake::SesResults expected(query.size(), target.size(), 0, 1, 0, 0, 0, true,
                                         PacBio::Data::Cigar(cigar));

    // Run.
    PacBio::Pancake::SesResults result =
        PacBio::Pancake::SESAlignBanded<PacBio::Pancake::SESAlignMode::Global,
                                        PacBio::Pancake::SESTracebackMode::Enabled>(
            query, target, maxDiffs, bandwidth);

    std::cerr << "CIGAR: " << result.cigar.ToStdString() << "\n";

    // Evaluate.
    EXPECT_EQ(expected, result);
}

TEST(SESAlignBanded_Global, NormalSmallCase_MultipleExactMatches)
{
    /*
     * Note: Mismatches count as 2 diffs.
    */
    // Inputs.
    std::string query = "ACTG";
    std::string target = "ACTG";
    int32_t maxDiffs = 10;
    int32_t bandwidth = 30;
    std::string cigar = "4=";

    // Order of elements in SesResult: lastQueryPos, lastTargetPos, diffs, minK, maxK, valid
    PacBio::Pancake::SesResults expected(query.size(), target.size(), 0, 4, 0, 0, 0, true,
                                         PacBio::Data::Cigar(cigar));

    // Run.
    PacBio::Pancake::SesResults result =
        PacBio::Pancake::SESAlignBanded<PacBio::Pancake::SESAlignMode::Global,
                                        PacBio::Pancake::SESTracebackMode::Enabled>(
            query, target, maxDiffs, bandwidth);

    std::cerr << "CIGAR: " << result.cigar.ToStdString() << "\n";

    // Evaluate.
    EXPECT_EQ(expected, result);
}

TEST(SESAlignBanded_Global, NormalSmallCase_SingleMismatch)
{
    /*
     * Note: Mismatches count as 2 diffs.
    */
    // Inputs.
    std::string query = "A";
    std::string target = "C";
    int32_t maxDiffs = 5;
    int32_t bandwidth = 30;
    std::string cigar = "1X";

    // Order of elements in SesResult: lastQueryPos, lastTargetPos, diffs, minK, maxK, valid
    PacBio::Pancake::SesResults expected(query.size(), target.size(), 1, 0, 1, 0, 0, true,
                                         PacBio::Data::Cigar(cigar));

    // Run.
    PacBio::Pancake::SesResults result =
        PacBio::Pancake::SESAlignBanded<PacBio::Pancake::SESAlignMode::Global,
                                        PacBio::Pancake::SESTracebackMode::Enabled>(
            query, target, maxDiffs, bandwidth);

    std::cerr << "CIGAR: " << result.cigar.ToStdString() << "\n";

    // Evaluate.
    EXPECT_EQ(expected, result);
}

TEST(SESAlignBanded_Global, NormalSmallCase_MultipleMismatches)
{
    /*
     * Note: Mismatches count as 2 diffs.
    */
    // Inputs.
    std::string query = "CCCC";
    std::string target = "GGGG";
    int32_t maxDiffs = 30;
    int32_t bandwidth = 30;
    std::string cigar = "4X";

    // Order of elements in SesResult: lastQueryPos, lastTargetPos, diffs, minK, maxK, valid
    PacBio::Pancake::SesResults expected(query.size(), target.size(), 4, 0, 4, 0, 0, true,
                                         PacBio::Data::Cigar(cigar));

    // Run.
    PacBio::Pancake::SesResults result =
        PacBio::Pancake::SESAlignBanded<PacBio::Pancake::SESAlignMode::Global,
                                        PacBio::Pancake::SESTracebackMode::Enabled>(
            query, target, maxDiffs, bandwidth);

    std::cerr << "CIGAR: " << result.cigar.ToStdString() << "\n";

    // Evaluate.
    EXPECT_EQ(expected, result);
}

/////////////////////////////////////////////
/////////////////////////////////////////////
/////////////////////////////////////////////

TEST(SESAlignBanded_Semiglobal, NormalSmallCase_SimpleSeqWithAllDiffs)
{
    /*
     * Note: Mismatches count as 2 diffs.
    */
    // Inputs.
    std::string query = "GGATCAGTT";   // GGAT-CAGTT
                                       // |X|| | |||
    std::string target = "GAATTCGTT";  // GAATTC-GTT
    int32_t maxDiffs = 6;
    int32_t bandwidth = 30;
    std::string cigar = "1=1I1=1D1=1D1=1I3=";

    /*
     * GGA-T-CAGTT
     * | | | | |||
     * G-AATTC-GTT
    */

    // Order of elements in SesResult: lastQueryPos, lastTargetPos, diffs, minK, maxK, valid
    PacBio::Pancake::SesResults expected(9, 9, 4, 7, 0, 2, 2, true, PacBio::Data::Cigar(cigar));

    // Run.
    PacBio::Pancake::SesResults result =
        PacBio::Pancake::SESAlignBanded<PacBio::Pancake::SESAlignMode::Semiglobal,
                                        PacBio::Pancake::SESTracebackMode::Enabled>(
            query, target, maxDiffs, bandwidth);

    std::cerr << "CIGAR: " << result.cigar.ToStdString() << "\n";

    // Evaluate.
    EXPECT_EQ(expected, result);
}

TEST(SESAlignBanded_Semiglobal, NormalSmallCase_AnotherSimpleSeqWithAllDiffs)
{
    /*
     * Note: Mismatches count as 2 diffs.
    */
    // Inputs.
    std::string query = "GGATTCAGTT";   // GGA-TTCAGTT
                                        // | | |||X|||
    std::string target = "GAATTCCGTT";  // G-AATTCCGTT
    int32_t maxDiffs = 15;
    int32_t bandwidth = 30;
    std::string cigar = "1=1I1=1D3=1X3=";

    /*
     * GGA-T-CAGTT
     * | | | | |||
     * G-AATTC-GTT
    */

    // Order of elements in SesResult: lastQueryPos, lastTargetPos, diffs, minK, maxK, valid
    PacBio::Pancake::SesResults expected(10, 10, 3, 8, 1, 1, 1, true, PacBio::Data::Cigar(cigar));

    // Run.
    PacBio::Pancake::SesResults result =
        PacBio::Pancake::SESAlignBanded<PacBio::Pancake::SESAlignMode::Semiglobal,
                                        PacBio::Pancake::SESTracebackMode::Enabled>(
            query, target, maxDiffs, bandwidth);

    std::cerr << "CIGAR: " << result.cigar.ToStdString() << "\n";

    // Evaluate.
    EXPECT_EQ(expected, result);
}

TEST(SESAlignBanded_Semiglobal, NormalSmallCase_CIGAR_5D_in_suffix)
{
    /*
     * Note: Mismatches count as 2 diffs.
    */
    // Inputs.
    std::string query = "ACTG";
    std::string target = "ACTGAAAAA";
    int32_t maxDiffs = 15;
    int32_t bandwidth = 30;
    std::string cigar = "4=";

    // Order of elements in SesResult: lastQueryPos, lastTargetPos, diffs, numEq, numX, numI, numD, valid, cigar
    PacBio::Pancake::SesResults expected(4, 4, 0, 4, 0, 0, 0, true, PacBio::Data::Cigar(cigar));

    // Run.
    PacBio::Pancake::SesResults result =
        PacBio::Pancake::SESAlignBanded<PacBio::Pancake::SESAlignMode::Semiglobal,
                                        PacBio::Pancake::SESTracebackMode::Enabled>(
            query, target, maxDiffs, bandwidth);

    std::cerr << "CIGAR: " << result.cigar.ToStdString() << "\n";

    // Evaluate.
    EXPECT_EQ(expected, result);
}

TEST(SESAlignBanded_Semiglobal, NormalSmallCase_CIGAR_5I_in_suffix)
{
    /*
     * Note: Mismatches count as 2 diffs.
    */
    // Inputs.
    std::string query = "ACTGAAAAA";
    std::string target = "ACTG";
    int32_t maxDiffs = 15;
    int32_t bandwidth = 30;
    std::string cigar = "4=";

    // Order of elements in SesResult: lastQueryPos, lastTargetPos, diffs, numEq, numX, numI, numD, valid, cigar
    PacBio::Pancake::SesResults expected(4, 4, 0, 4, 0, 0, 0, true, PacBio::Data::Cigar(cigar));

    // Run.
    PacBio::Pancake::SesResults result =
        PacBio::Pancake::SESAlignBanded<PacBio::Pancake::SESAlignMode::Semiglobal,
                                        PacBio::Pancake::SESTracebackMode::Enabled>(
            query, target, maxDiffs, bandwidth);

    std::cerr << "CIGAR: " << result.cigar.ToStdString() << "\n";

    // Evaluate.
    EXPECT_EQ(expected, result);
}

TEST(SESAlignBanded_Semiglobal, NormalSmallCase_CIGAR_5D_in_prefix)
{
    /*
     * Note: Mismatches count as 2 diffs.
     * IMPORTANT: In the semiglobal alignment mode where the beginning is locked
     * and the end can be open ended, the alignment of a small query with a target
     * that has a large prefix indel will likely result in an alignment with
     * a lot of diffs.
     *
     * In this particular case, one option would be "1X1=2X" because one "C" matches
     * the deleted prefix, and it's much less expensive than introducing 5D and then
     * matching.
     *
     * A valid alternative is also: "1I1=2I".
     *
    */
    // Inputs.
    std::string query = "ACTG";
    std::string target = "CCCCCACTG";
    int32_t maxDiffs = 15;
    int32_t bandwidth = 30;
    std::string cigar = "1I1=2I";

    // Order of elements in SesResult: lastQueryPos, lastTargetPos, diffs, numEq, numX, numI, numD, valid, cigar
    PacBio::Pancake::SesResults expected(4, 1, 3, 1, 0, 3, 0, true, PacBio::Data::Cigar(cigar));

    // Run.
    PacBio::Pancake::SesResults result =
        PacBio::Pancake::SESAlignBanded<PacBio::Pancake::SESAlignMode::Semiglobal,
                                        PacBio::Pancake::SESTracebackMode::Enabled>(
            query, target, maxDiffs, bandwidth);

    std::cerr << "CIGAR: " << result.cigar.ToStdString() << "\n";

    // Evaluate.
    EXPECT_EQ(expected, result);
}

TEST(SESAlignBanded_Semiglobal, NormalSmallCase_CIGAR_5I_in_prefix)
{
    /*
     * Note: Mismatches count as 2 diffs.
    */
    // Inputs.
    std::string query = "CCCCCACTG";
    std::string target = "ACTG";
    int32_t maxDiffs = 15;
    int32_t bandwidth = 30;
    std::string cigar = "1D1=2D";

    // Order of elements in SesResult: lastQueryPos, lastTargetPos, diffs, numEq, numX, numI, numD, valid, cigar
    PacBio::Pancake::SesResults expected(1, 4, 3, 1, 0, 0, 3, true, PacBio::Data::Cigar(cigar));

    // Run.
    PacBio::Pancake::SesResults result =
        PacBio::Pancake::SESAlignBanded<PacBio::Pancake::SESAlignMode::Semiglobal,
                                        PacBio::Pancake::SESTracebackMode::Enabled>(
            query, target, maxDiffs, bandwidth);

    std::cerr << "CIGAR: " << result.cigar.ToStdString() << "\n";

    // Evaluate.
    EXPECT_EQ(expected, result);
}

TEST(SESAlignBanded_Semiglobal, NormalSmallCase_SingleMatch)
{
    /*
     * Note: Mismatches count as 2 diffs.
    */
    // Inputs.
    std::string query = "A";
    std::string target = "A";
    int32_t maxDiffs = 5;
    int32_t bandwidth = 30;
    std::string cigar = "1=";

    // Order of elements in SesResult: lastQueryPos, lastTargetPos, diffs, numEq, numX, numI, numD, valid, cigar
    PacBio::Pancake::SesResults expected(query.size(), target.size(), 0, 1, 0, 0, 0, true,
                                         PacBio::Data::Cigar(cigar));

    // Run.
    PacBio::Pancake::SesResults result =
        PacBio::Pancake::SESAlignBanded<PacBio::Pancake::SESAlignMode::Semiglobal,
                                        PacBio::Pancake::SESTracebackMode::Enabled>(
            query, target, maxDiffs, bandwidth);

    std::cerr << "CIGAR: " << result.cigar.ToStdString() << "\n";

    // Evaluate.
    EXPECT_EQ(expected, result);
}

TEST(SESAlignBanded_Semiglobal, NormalSmallCase_MultipleExactMatches)
{
    /*
     * Note: Mismatches count as 2 diffs.
    */
    // Inputs.
    std::string query = "ACTG";
    std::string target = "ACTG";
    int32_t maxDiffs = 10;
    int32_t bandwidth = 30;
    std::string cigar = "4=";

    // Order of elements in SesResult: lastQueryPos, lastTargetPos, diffs, numEq, numX, numI, numD, valid, cigar
    PacBio::Pancake::SesResults expected(query.size(), target.size(), 0, 4, 0, 0, 0, true,
                                         PacBio::Data::Cigar(cigar));

    // Run.
    PacBio::Pancake::SesResults result =
        PacBio::Pancake::SESAlignBanded<PacBio::Pancake::SESAlignMode::Semiglobal,
                                        PacBio::Pancake::SESTracebackMode::Enabled>(
            query, target, maxDiffs, bandwidth);

    std::cerr << "CIGAR: " << result.cigar.ToStdString() << "\n";

    // Evaluate.
    EXPECT_EQ(expected, result);
}

TEST(SESAlignBanded_Semiglobal, NormalSmallCase_SingleMismatch)
{
    /*
     * Note: Mismatches count as 2 diffs.
    */
    // Inputs.
    std::string query = "A";
    std::string target = "C";
    int32_t maxDiffs = 5;
    int32_t bandwidth = 30;
    std::string cigar = "1D";

    // Order of elements in SesResult: lastQueryPos, lastTargetPos, diffs, numEq, numX, numI, numD, valid, cigar
    PacBio::Pancake::SesResults expected(0, 1, 1, 0, 0, 0, 1, true, PacBio::Data::Cigar(cigar));

    // Run.
    PacBio::Pancake::SesResults result =
        PacBio::Pancake::SESAlignBanded<PacBio::Pancake::SESAlignMode::Semiglobal,
                                        PacBio::Pancake::SESTracebackMode::Enabled>(
            query, target, maxDiffs, bandwidth);

    std::cerr << "CIGAR: " << result.cigar.ToStdString() << "\n";

    // Evaluate.
    EXPECT_EQ(expected, result);
}

TEST(SESAlignBanded_Semiglobal, NormalSmallCase_MultipleMismatches)
{
    /*
     * Note: Mismatches count as 2 diffs.
    */
    // Inputs.
    std::string query = "CCCC";
    std::string target = "GGGG";
    int32_t maxDiffs = 15;
    int32_t bandwidth = 30;
    std::string cigar = "4D";

    // Order of elements in SesResult: lastQueryPos, lastTargetPos, diffs, numEq, numX, numI, numD, valid, cigar
    PacBio::Pancake::SesResults expected(0, 4, 4, 0, 0, 0, 4, true, PacBio::Data::Cigar(cigar));

    // Run.
    PacBio::Pancake::SesResults result =
        PacBio::Pancake::SESAlignBanded<PacBio::Pancake::SESAlignMode::Semiglobal,
                                        PacBio::Pancake::SESTracebackMode::Enabled>(
            query, target, maxDiffs, bandwidth);

    std::cerr << "CIGAR: " << result.cigar.ToStdString() << "\n";

    // Evaluate.
    EXPECT_EQ(expected, result);
}
