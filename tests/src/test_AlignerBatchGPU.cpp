#include <PancakeTestData.h>
#include <gtest/gtest.h>
#include <pacbio/pancake/AlignerBatchGPU.h>
#include <claraparabricks/genomeworks/cudaaligner/aligner.hpp>
#include <claraparabricks/genomeworks/cudaaligner/alignment.hpp>
#include <claraparabricks/genomeworks/cudaaligner/cudaaligner.hpp>
#include <iostream>
#include <string>
#include <vector>
#include "TestHelperUtils.h"

namespace TestsAlignerBatchGPU {

struct TestData
{
    std::string testName;
    std::vector<std::pair<std::string, std::string>> batchData;
    std::vector<PacBio::Pancake::AlignmentResult> expectedAlns;
};

// clang-format off
std::vector<TestData> testData = {
    {
        "Batch of multiple query/target pairs.",
        {
            {"ACTGACTGAC", "ACTGTCTGAC"},
            {"ACTG", "ACTG"},
            {"A", "T"},
        },
        // Expected results.
        {
            PacBio::Pancake::AlignmentResult{PacBio::BAM::Cigar("4=1X5="), 10, 10, 10, 10, true, 14, 14, false, PacBio::Pancake::Alignment::DiffCounts(9, 1, 0, 0)},
            PacBio::Pancake::AlignmentResult{PacBio::BAM::Cigar("4="), 4, 4, 4, 4, true, 8, 8, false, PacBio::Pancake::Alignment::DiffCounts(4, 0, 0, 0)},
            PacBio::Pancake::AlignmentResult{PacBio::BAM::Cigar("1X"), 1, 1, 1, 1, true, -4, -4, false, PacBio::Pancake::Alignment::DiffCounts(0, 1, 0, 0)},
        },
    },
    {
        "Empty batch.",
        {
        },
        // Expected results.
        {
        },
    },
    {
        "Another batch.",
        {
            {"AAAAA", "AAAAA"},
        },
        // Expected results.
        {
            PacBio::Pancake::AlignmentResult{PacBio::BAM::Cigar("5="), 5, 5, 5, 5, true, 10, 10, false, PacBio::Pancake::Alignment::DiffCounts(5, 0, 0, 0)},
        },
    },
    {   // Cudaaligner can not handle zero-length query or target sequences at the moment. These are handled in the AlignerBatchGPU class.
        "Another batch of edge cases.",
        {
            {"", ""},
            {"A", ""},
            {"", "A"},
            {"A", "T"},
        },
        // Expected results.
        {
            PacBio::Pancake::AlignmentResult{PacBio::BAM::Cigar(""), 0, 0, 0, 0, false, 0, 0, false, PacBio::Pancake::Alignment::DiffCounts(0, 0, 0, 0)},
            PacBio::Pancake::AlignmentResult{PacBio::BAM::Cigar("1I"), 1, 0, 1, 0, true, -4, -4, false, PacBio::Pancake::Alignment::DiffCounts(0, 0, 1, 0)},
            PacBio::Pancake::AlignmentResult{PacBio::BAM::Cigar("1D"), 0, 1, 0, 1, true, -4, -4, false, PacBio::Pancake::Alignment::DiffCounts(0, 0, 0, 1)},
            PacBio::Pancake::AlignmentResult{PacBio::BAM::Cigar("1X"), 1, 1, 1, 1, true, -4, -4, false, PacBio::Pancake::Alignment::DiffCounts(0, 1, 0, 0)},
        },
    },
};
// clang-format on

TEST(AlignerBatchGPU, ArrayOfTests_Small)
{
    using namespace PacBio::Pancake;

    uint32_t maxBandwidth = 2000;
    uint32_t deviceId = 0;
    PacBio::Pancake::AlignmentParameters alnParams;

    const int32_t maxMemoryCap = 100 * 1024 * 1024;
    auto aligner =
        PacBio::Pancake::AlignerBatchGPU(alnParams, maxBandwidth, deviceId, maxMemoryCap);

    for (const auto& data : testData) {
        // Debug info.
        SCOPED_TRACE(data.testName);
        // std::cerr << "testName = " << data.testName << "\n";

        // Reuse the aligner in multiple batches.
        aligner.Clear();

        for (const auto& seqPair : data.batchData) {
            const auto& query = seqPair.first;
            const auto& target = seqPair.second;
            aligner.AddSequencePair(query.c_str(), query.size(), target.c_str(), target.size());
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
        EXPECT_EQ(data.expectedAlns, results);
    }
}
}  // namespace TestsAlignerBatchGPU
