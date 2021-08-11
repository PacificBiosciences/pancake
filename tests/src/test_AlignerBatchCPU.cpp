#include <PancakeTestData.h>
#include <gtest/gtest.h>
#include <pacbio/pancake/AlignerBatchCPU.h>
#include <pacbio/pancake/AlignerFactory.h>
#include <pacbio/pancake/AlignmentParameters.h>
#include <iostream>
#include "TestHelperUtils.h"

TEST(AlignerBatchCPU, ArrayOfTests_Small)
{
    using namespace PacBio::Pancake;

    struct TestData
    {
        std::string testName;
        // Tuple: <query, target, isGlobal>
        std::vector<std::tuple<std::string, std::string, bool>> batchData;
        int32_t numThreads;
        std::vector<PacBio::Pancake::AlignmentResult> expectedAlns;
    };

    // clang-format off
    std::vector<TestData> testData = {
        {
            "Batch of multiple query/target pairs.",
            {
                // Global alignment.
                {"ACTGACTGAC", "ACTGTCTGAC", true},
                {"ACTG", "ACTG", true},
                {"A", "T", true},
                // Extension alignment.
                {"AAAAAAAAAAAAAAAAAAAACCCCCACCCCCCCCCCCCCCCCCCCCCCCCC", "AAAAAAAAAAAAAAAAAAAAGGGGGAGGGGGGGGGGGGGGGGGGGGGGGGG", false},
            },
            4,
            // Expected results.
            {
                PacBio::Pancake::AlignmentResult{PacBio::BAM::Cigar("4=1X5="), 10, 10, 10, 10, true, 14, 14, false, Alignment::DiffCounts(9, 1, 0, 0)},
                PacBio::Pancake::AlignmentResult{PacBio::BAM::Cigar("4="), 4, 4, 4, 4, true, 8, 8, false, Alignment::DiffCounts(4, 0, 0, 0)},
                PacBio::Pancake::AlignmentResult{PacBio::BAM::Cigar("1X"), 1, 1, 1, 1, true, -4, -4, false, Alignment::DiffCounts(0, 1, 0, 0)},
                PacBio::Pancake::AlignmentResult{PacBio::BAM::Cigar("20="), 20, 20, 19, 19, true, 40, 40, true, Alignment::DiffCounts(20, 0, 0, 0)},
            },
        },
        {
            "Empty batch.",
            {
            },
            4,
            // Expected results.
            {
            },
        },
        {
            "Another batch of edge cases.",
            {
                // Global alignment.
                {"", "", true},
                {"A", "", true},
                {"", "A", true},
                {"A", "T", true},
            },
            4,
            // Expected results.
            {
                PacBio::Pancake::AlignmentResult{PacBio::BAM::Cigar(""), 0, 0, 0, 0, true, 0, 0, false, Alignment::DiffCounts(0, 0, 0, 0)},
                PacBio::Pancake::AlignmentResult{PacBio::BAM::Cigar("1I"), 1, 0, 1, 0, true, -4, -4, false, Alignment::DiffCounts(0, 0, 1, 0)},
                PacBio::Pancake::AlignmentResult{PacBio::BAM::Cigar("1D"), 0, 1, 0, 1, true, -4, -4, false, Alignment::DiffCounts(0, 0, 0, 1)},
                PacBio::Pancake::AlignmentResult{PacBio::BAM::Cigar("1X"), 1, 1, 1, 1, true, -4, -4, false, Alignment::DiffCounts(0, 1, 0, 0)},
            },
        },
    };
    // clang-format on

    // Alignment parameters.
    PacBio::Pancake::AlignerType alignerTypeGlobal = AlignerType::EDLIB;
    PacBio::Pancake::AlignmentParameters alnParamsGlobal;
    PacBio::Pancake::AlignerType alignerTypeExt = AlignerType::KSW2;
    PacBio::Pancake::AlignmentParameters alnParamsExt;
    // These parameters are the same as current defaults, but they are specified here in case
    // the defaults ever change, so that the test doesn't have to be updated.
    alnParamsGlobal.zdrop = 100;
    alnParamsGlobal.zdrop2 = 500;
    alnParamsGlobal.alignBandwidth = 500;
    alnParamsGlobal.endBonus = 50;
    alnParamsGlobal.matchScore = 2;
    alnParamsGlobal.mismatchPenalty = 4;
    alnParamsGlobal.gapOpen1 = 4;
    alnParamsGlobal.gapExtend1 = 2;
    alnParamsGlobal.gapOpen2 = 24;
    alnParamsGlobal.gapExtend2 = 1;
    alnParamsExt = alnParamsGlobal;

    for (const auto& data : testData) {
        // Debug info.
        SCOPED_TRACE(data.testName);
        std::cerr << "testName = " << data.testName << "\n";

        PacBio::Pancake::AlignerBatchCPU aligner(data.numThreads, alignerTypeGlobal,
                                                 alnParamsGlobal, alignerTypeExt, alnParamsExt);

        for (const auto& seqPair : data.batchData) {
            const auto& query = std::get<0>(seqPair);
            const auto& target = std::get<1>(seqPair);
            const bool isGlobal = std::get<2>(seqPair);
            if (isGlobal) {
                aligner.AddSequencePairForGlobalAlignment(query.c_str(), query.size(),
                                                          target.c_str(), target.size());
            } else {
                aligner.AddSequencePairForExtensionAlignment(query.c_str(), query.size(),
                                                             target.c_str(), target.size());
            }
        }

        // Run alignment.
        aligner.AlignAll();

        const std::vector<PacBio::Pancake::AlignmentResult>& results = aligner.GetAlnResults();

        // std::cerr << "results.size() = " << results.size() << "\n";
        // for (size_t i = 0; i < results.size(); ++i) {
        //     const auto& aln = results[i];
        //     std::cerr << "[result " << i << "] " << aln << "\n";
        // }

        // Evaluate.
        ASSERT_EQ(data.expectedAlns, results);
    }
}
