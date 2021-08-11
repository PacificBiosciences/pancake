#include <PancakeTestData.h>
#include <gtest/gtest.h>
#include <pacbio/pancake/MapperHiFi.h>
#include <pacbio/pancake/OverlapWriterBase.h>
#include <iostream>
#include <memory>
#include "TestHelperUtils.h"

TEST(MapperHiFi, CheckMappping_LoadFromFile)
{
    struct TestData
    {
        const std::string testName;
        const std::string targetFile;
        const std::string queryFile;
        const PacBio::Pancake::SeedDB::SeedDBParameters seedParams;
        const std::string expectedOverlapFile;
    };

    // clang-format off
    std::vector<TestData> testData = {
        {
            "Empty input",
            // Target.
            "",
            // Query.
            "",
            // SeedParams.
            PacBio::Pancake::SeedDB::SeedDBParameters{28, 80, 0, false, false, true},
            // Expected results.
            ""
        },
        {
            "Small synthetic test case with 2 refs and one query.",
            // Target.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapping/test-1-no-secondary-aln.ref.fasta",
            // Query.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapping/test-1-no-secondary-aln.reads.fasta",
            // SeedParams.
            PacBio::Pancake::SeedDB::SeedDBParameters{28, 80, 0, false, false, true},
            // Expected results.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapping/test-1-no-secondary-aln.out.ovl"
        },
        {
            "A set of real reads, used as both the query and target sequences",
            // Target.
            PacBio::PancakeTestsConfig::Data_Dir + "/hifi-ovl/reads.pile1-5prime.fasta",
            // Query.
            PacBio::PancakeTestsConfig::Data_Dir + "/hifi-ovl/reads.pile1-5prime.fasta",
            // SeedParams.
            PacBio::Pancake::SeedDB::SeedDBParameters{28, 80, 0, false, false, true},
            // Expected results.
            PacBio::PancakeTestsConfig::Data_Dir + "/hifi-ovl/reads.pile1-5prime.out.all_vs_all.ovl"
        },

    };
    // clang-format on

    PacBio::Pancake::OverlapHifiSettings settings;

    for (const auto& data : testData) {
        // Load the sequences from files, and take only the first one.
        const std::vector<std::string> allTargetSeqs =
            data.targetFile.empty()
                ? std::vector<std::string>()
                : PacBio::PancakeTests::HelperLoadFastaAsStringVector(data.targetFile);
        const std::vector<std::string> allQuerySeqs =
            data.queryFile.empty()
                ? std::vector<std::string>()
                : PacBio::PancakeTests::HelperLoadFastaAsStringVector(data.queryFile);

        std::vector<std::string> expectedOverlaps =
            data.expectedOverlapFile.empty()
                ? std::vector<std::string>()
                : PacBio::PancakeTests::HelperLoadFile(data.expectedOverlapFile);

        SCOPED_TRACE(data.testName);

        // std::cerr << "testName = " << data.testName << "\n";

        std::vector<PacBio::Pancake::OverlapHiFi::MapperResult> result =
            PacBio::Pancake::OverlapHiFi::MapHiFi(allTargetSeqs, allQuerySeqs, data.seedParams,
                                                  settings);

        std::vector<std::string> resultsStr;
        for (const auto& queryMappings : result) {
            for (size_t i = 0; i < queryMappings.overlaps.size(); ++i) {
                if (queryMappings.overlaps[i] == nullptr) {
                    continue;
                }
                const auto& mapping = queryMappings.overlaps[i];

                // std::cerr << PacBio::Pancake::OverlapWriterBase::PrintOverlapAsM4(*mapping, "", "",
                //                                                                   true, false)
                //           << "\n";

                resultsStr.emplace_back(PacBio::Pancake::OverlapWriterBase::PrintOverlapAsM4(
                    *mapping, "", "", true, false));
            }
        }

        ASSERT_EQ(expectedOverlaps, resultsStr);
    }
}

TEST(MapperHiFi, ArbitrarySequenceIDs)
{
    struct TestData
    {
        const std::string testName;
        const std::string targetFile;
        const std::string queryFile;
        const int32_t sequenceIdOffset = 0;
        const PacBio::Pancake::SeedDB::SeedDBParameters seedParams;
        const std::string expectedOverlapFile;
    };

    // clang-format off
    std::vector<TestData> testData = {
        {
            "Overlap the same set of reads with itself. No ID offset (the first query and target have the ID == 0).",
            // Target.
            PacBio::PancakeTestsConfig::Data_Dir + "/hifi-ovl/reads.pile1-5prime.fasta",
            // Query.
            PacBio::PancakeTestsConfig::Data_Dir + "/hifi-ovl/reads.pile1-5prime.fasta",
            // Sequence ID offset.
            0,
            // SeedParams.
            PacBio::Pancake::SeedDB::SeedDBParameters{28, 80, 0, false, false, true},
            // Expected results.
            PacBio::PancakeTestsConfig::Data_Dir + "/hifi-ovl/reads.pile1-5prime.out.all_vs_all.m4",
        },
        {
            "Overlap the same set of reads with itself. Offset the IDs by 10000 to check that IDs can be arbitrary.",
            // Target.
            PacBio::PancakeTestsConfig::Data_Dir + "/hifi-ovl/reads.pile1-5prime.fasta",
            // Query.
            PacBio::PancakeTestsConfig::Data_Dir + "/hifi-ovl/reads.pile1-5prime.fasta",
            // Sequence ID offset.
            10000,
            // SeedParams.
            PacBio::Pancake::SeedDB::SeedDBParameters{28, 80, 0, false, false, true},
            // Expected results.
            PacBio::PancakeTestsConfig::Data_Dir + "/hifi-ovl/reads.pile1-5prime.out.all_vs_all.id_offset_10000.m4",
        },
    };
    // clang-format on

    PacBio::Pancake::OverlapHifiSettings settings;

    for (const auto& data : testData) {
        // Load the sequences from files, and take only the first one.
        const std::vector<PacBio::Pancake::FastaSequenceId> allTargetSeqs =
            data.targetFile.empty() ? std::vector<PacBio::Pancake::FastaSequenceId>()
                                    : PacBio::PancakeTests::HelperLoadFastaWithId(
                                          data.targetFile, data.sequenceIdOffset);
        const std::vector<PacBio::Pancake::FastaSequenceId> allQuerySeqs =
            data.queryFile.empty() ? std::vector<PacBio::Pancake::FastaSequenceId>()
                                   : PacBio::PancakeTests::HelperLoadFastaWithId(
                                         data.queryFile, data.sequenceIdOffset);
        const PacBio::Pancake::FastaSequenceCachedStore targetCacheStore(allTargetSeqs);
        const PacBio::Pancake::FastaSequenceCachedStore queryCacheStore(allQuerySeqs);

        std::vector<std::string> expectedOverlaps =
            data.expectedOverlapFile.empty()
                ? std::vector<std::string>()
                : PacBio::PancakeTests::HelperLoadFile(data.expectedOverlapFile);

        SCOPED_TRACE(data.testName);

        std::vector<PacBio::Pancake::OverlapHiFi::MapperResult> result =
            PacBio::Pancake::OverlapHiFi::MapHiFi(targetCacheStore, queryCacheStore,
                                                  data.seedParams, settings);

        std::vector<std::string> resultsStr;
        for (const auto& queryMappings : result) {
            for (size_t i = 0; i < queryMappings.overlaps.size(); ++i) {
                if (queryMappings.overlaps[i] == nullptr) {
                    continue;
                }
                const auto& mapping = queryMappings.overlaps[i];

                resultsStr.emplace_back(
                    PacBio::Pancake::OverlapWriterBase::PrintOverlapAsM4(*mapping));
            }
        }

        // Sort so that the comparison is good.
        std::sort(expectedOverlaps.begin(), expectedOverlaps.end());
        std::sort(resultsStr.begin(), resultsStr.end());

        // std::cerr << "Expected:\n";
        // for (size_t i = 0; i < expectedOverlaps.size(); ++i) {
        //     const auto& ovlStr = expectedOverlaps[i];
        //     std::cerr << "[i = " << i << "] " << ovlStr << "\n";
        // }

        // std::cerr << "Results:\n";
        // for (size_t i = 0; i < resultsStr.size(); ++i) {
        //     const auto& ovlStr = resultsStr[i];
        //     std::cerr << "[i = " << i << "] " << ovlStr << "\n";
        // }

        // for (size_t i = 0; i < std::min(expectedOverlaps.size(), resultsStr.size()); ++i) {
        //     if (resultsStr[i] != expectedOverlaps[i]) {
        //         std::cerr << "Different [i = " << i << "]:\n\tExpected: " << expectedOverlaps[i] << "\n\tResult:   " << resultsStr[i] << "\n";
        //     }
        // }

        ASSERT_EQ(expectedOverlaps, resultsStr);
    }
}
