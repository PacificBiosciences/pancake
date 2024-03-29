#include "TestHelperUtils.hpp"

#include <PancakeTestData.h>
#include <pancake/MapperBatchCPU.hpp>
#include <pancake/OverlapWriterBase.hpp>

#include <gtest/gtest.h>
#include <iostream>

TEST(MapperBatchCPU, BatchMapping_ArrayOfTests)
{
    using namespace PacBio::Pancake;

    struct TestData
    {
        const std::string testName;
        const std::vector<std::pair<std::string, std::string>> batchData;
        const int32_t sequenceIdOffset = 0;
        const PacBio::Pancake::AlignerType alignerTypeGlobal;
        const PacBio::Pancake::SeedDBParameters seedParamsPrimary;
        const PacBio::Pancake::SeedDBParameters seedParamsFallback;
        const std::vector<std::vector<std::string>> expectedOverlaps;
    };

    // clang-format off
    std::vector<TestData> testData = {
        {
            "Batch 1 of multiple query/target vectors.",
            {
                {
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-2-real-insert-target.fasta",
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-2-real-insert-query.fasta",
                },
                {
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-3-real-subseq-full-global-match-target.fasta",
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-3-real-subseq-full-global-match-query.fasta",
                },
                {
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-4-small-revcmp-perfect-aln-target.fasta",
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-4-small-revcmp-perfect-aln-query.fasta",
                },
                {
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-5-revcmp-aln-target.fasta",
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-5-revcmp-aln-query.fasta",
                },
                {
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-6-small-revcmp-imperfect-aln-target.fasta",
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-6-small-revcmp-imperfect-aln-query.fasta",
                },
                {
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-1-poor-aln-overlapping-seeds.target.fasta",
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-1-poor-aln-overlapping-seeds.query.fasta",
                },
                {
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-8-no-back-flank-extension-target.fasta",
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-8-no-back-flank-extension-query.fasta",
                },

            },
            // Sequence ID offset. Zero means that the first query/target after loading has the ID == 0.
            0,
            // Aligner type for global alignment.
            AlignerType::KSW2,
            // SeedParams - primary.
            PacBio::Pancake::SeedDBParameters{19, 10, 0, false, false, true},
            // SeedParams - fallback.
            PacBio::Pancake::SeedDBParameters{10, 5, 0, false, false, true},
            // Expected results.
            {
                {
                    "000000000 000000000 16402 81.12 0 0 18779 18779 0 0 18864 18865 *",
                },
                {
                    "000000000 000000000 15684 98.71 0 0 8111 8111 0 0 8138 8138 *"
                },
                {
                    "000000000 000000000 300 100.00 0 0 150 150 1 0 150 150 *"
                },
                {
                    "000000000 000000000 12717 74.59 0 0 21382 22015 1 701 21992 22028 *"
                },
                {
                    "000000000 000000000 270 96.77 0 0 155 155 1 0 150 150 *"
                },
                {
                    "000000000 000000000 284 68.58 0 6209 7938 43446 0 7261 8999 46238 *"
                },
                {
                    "000000000 000000000 9288 75.87 0 0 15753 15753 1 2 15953 15953 *"
                },
            },
        },
        {
            "Batch 2 of multiple query/target vectors. Using different seeding parameters.",
            {
                {
                    // This pair will result in a seed hit to be placed at the very end of the target sequence, which means
                    // that there will be no flank sequence to align (extend). This is a useful test to verify that
                    // a stitching end condition works well.
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-8-no-back-flank-extension-target.fasta",
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-8-no-back-flank-extension-query.fasta",
                },
                {
                    // This pair will result in a seed hit to be placed at the very beginning of the query sequence, which means
                    // that there will be no flank sequence to align (extend). This is a useful test to verify that
                    // a stitching end condition works well.
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-9-no-front-flank-extension-target.fasta",
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-9-no-front-flank-extension-query.fasta",
                },
            },
            // Sequence ID offset. Zero means that the first query/target after loading has the ID == 0.
            0,
            // Aligner type for global alignment.
            AlignerType::EDLIB,
            // SeedParams - primary.
            PacBio::Pancake::SeedDBParameters{15, 5, 0, false, true, true},
            // SeedParams - fallback.
            PacBio::Pancake::SeedDBParameters{10, 5, 0, false, false, true},
            // Expected results.
            {
                {
                    "000000000 000000000 12370 76.04 0 0 15753 15753 1 2 15953 15953 *"
                },
                {
                    "000000000 000000000 8312 78.47 0 0 9230 9230 0 8372 17577 17578 *"
                },
            },
        },
        {
            "Batch 3 - same as Batch 1, but the query/target IDs start at an arbitrary value > 0.",
            {
                {
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-2-real-insert-target.fasta",
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-2-real-insert-query.fasta",
                },
                {
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-3-real-subseq-full-global-match-target.fasta",
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-3-real-subseq-full-global-match-query.fasta",
                },
                {
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-4-small-revcmp-perfect-aln-target.fasta",
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-4-small-revcmp-perfect-aln-query.fasta",
                },
                {
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-5-revcmp-aln-target.fasta",
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-5-revcmp-aln-query.fasta",
                },
                {
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-6-small-revcmp-imperfect-aln-target.fasta",
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-6-small-revcmp-imperfect-aln-query.fasta",
                },
                {
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-1-poor-aln-overlapping-seeds.target.fasta",
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-1-poor-aln-overlapping-seeds.query.fasta",
                },
                {
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-8-no-back-flank-extension-target.fasta",
                    PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-8-no-back-flank-extension-query.fasta",
                },

            },
            // Sequence ID offset. Zero means that the first query/target after loading has the ID == 0.
            1234567,
            // Aligner type for global alignment.
            AlignerType::KSW2,
            // SeedParams - primary.
            PacBio::Pancake::SeedDBParameters{19, 10, 0, false, false, true},
            // SeedParams - fallback.
            PacBio::Pancake::SeedDBParameters{10, 5, 0, false, false, true},
            // Expected results.
            {
                {
                    "001234567 001234567 16402 81.12 0 0 18779 18779 0 0 18864 18865 *",
                },
                {
                    "001234567 001234567 15684 98.71 0 0 8111 8111 0 0 8138 8138 *"
                },
                {
                    "001234567 001234567 300 100.00 0 0 150 150 1 0 150 150 *"
                },
                {
                    "001234567 001234567 12717 74.59 0 0 21382 22015 1 701 21992 22028 *"
                },
                {
                    "001234567 001234567 270 96.77 0 0 155 155 1 0 150 150 *"
                },
                {
                    "001234567 001234567 284 68.58 0 6209 7938 43446 0 7261 8999 46238 *"
                },
                {
                    "001234567 001234567 9288 75.87 0 0 15753 15753 1 2 15953 15953 *"
                },
            },
        },
    };
    // clang-format on

    for (const auto& data : testData) {
        // Debug info.
        SCOPED_TRACE(data.testName);
        std::cerr << "testName = " << data.testName << "\n";

        // Load the batch sequence data. The helper function takes
        // a vector of target-query filename pairs.
        std::vector<MapperBatchChunk> batchData;
        std::vector<PacBio::BAM::FastaSequence> allSeqs;
        PacBio::PancakeTests::HelperLoadBatchData(data.batchData, data.sequenceIdOffset, 0.000,
                                                  data.seedParamsPrimary, data.seedParamsFallback,
                                                  batchData, allSeqs);

        // Set the seed parameter settings and create a mapper.
        PacBio::Pancake::MapperCLRAlignSettings alignSettings;
        alignSettings.alignerTypeGlobal = data.alignerTypeGlobal;
        PacBio::Pancake::MapperBatchCPU mapper(alignSettings, 1);

        // Run the unit under test.
        // std::vector<std::vector<MapperBaseResult>> results = mapper.DummyMapAndAlign(batchData);
        std::vector<std::vector<MapperBaseResult>> results = mapper.MapAndAlign(batchData);

        // Format the results for comparison.
        std::vector<std::vector<std::string>> resultsStr =
            PacBio::PancakeTests::HelperFormatBatchMappingResults(results);

        // Evaluate.
        ASSERT_EQ(data.expectedOverlaps, resultsStr);
    }
}

TEST(MapperBatchCPU, CheckSelfHitPolicyAndSkippingSymmetrical)
{
    PacBio::Pancake::MapperCLRSettings settingsDefaultPolicy;
    {
        auto& settings = settingsDefaultPolicy;
        settings.map.bestNSecondary = -1;
        settings.map.secondaryMinScoreFraction = 0.0;
        settings.map.secondaryAllowedOverlapFractionTarget = 0.0;
        settings.map.secondaryAllowedOverlapFractionQuery = 0.0;
        settings.map.skipSymmetricOverlaps = false;
        settings.map.selfHitPolicy = PacBio::Pancake::MapperSelfHitPolicy::DEFAULT;
        settings.align.selfHitPolicy = PacBio::Pancake::MapperSelfHitPolicy::DEFAULT;
        settings.map.freqPercentile = 0.000;
        settings.map.seedParams = PacBio::Pancake::SeedDBParameters{19, 10, 0, false, false, true};
        settings.map.seedParamsFallback =
            PacBio::Pancake::SeedDBParameters{10, 5, 0, false, false, true};
        settings.align.alignerTypeGlobal = PacBio::Pancake::AlignerType::EDLIB;
    }

    PacBio::Pancake::MapperCLRSettings settingsSkipSelfHitsInBothMapAndAlign;
    {
        auto& settings = settingsSkipSelfHitsInBothMapAndAlign;
        settings.map.bestNSecondary = -1;
        settings.map.secondaryMinScoreFraction = 0.0;
        settings.map.secondaryAllowedOverlapFractionTarget = 0.0;
        settings.map.secondaryAllowedOverlapFractionQuery = 0.0;
        settings.map.skipSymmetricOverlaps = false;
        settings.map.selfHitPolicy = PacBio::Pancake::MapperSelfHitPolicy::SKIP;
        settings.align.selfHitPolicy = PacBio::Pancake::MapperSelfHitPolicy::SKIP;
        settings.map.freqPercentile = 0.000;
        settings.map.seedParams = PacBio::Pancake::SeedDBParameters{19, 10, 0, false, false, true};
        settings.map.seedParamsFallback =
            PacBio::Pancake::SeedDBParameters{10, 5, 0, false, false, true};
        settings.align.alignerTypeGlobal = PacBio::Pancake::AlignerType::EDLIB;
    }

    PacBio::Pancake::MapperCLRSettings settingsPerfectAlignSelfHitsInBothMapAndAlign;
    {
        auto& settings = settingsPerfectAlignSelfHitsInBothMapAndAlign;
        settings.map.bestNSecondary = -1;
        settings.map.secondaryMinScoreFraction = 0.0;
        settings.map.secondaryAllowedOverlapFractionTarget = 0.0;
        settings.map.secondaryAllowedOverlapFractionQuery = 0.0;
        settings.map.skipSymmetricOverlaps = false;
        settings.map.selfHitPolicy = PacBio::Pancake::MapperSelfHitPolicy::PERFECT_ALIGNMENT;
        settings.align.selfHitPolicy = PacBio::Pancake::MapperSelfHitPolicy::PERFECT_ALIGNMENT;
        settings.map.freqPercentile = 0.000;
        settings.map.seedParams = PacBio::Pancake::SeedDBParameters{19, 10, 0, false, false, true};
        settings.map.seedParamsFallback =
            PacBio::Pancake::SeedDBParameters{10, 5, 0, false, false, true};
        settings.align.alignerTypeGlobal = PacBio::Pancake::AlignerType::EDLIB;
    }

    PacBio::Pancake::MapperCLRSettings settingsSkipSymmetricOverlaps;
    {
        auto& settings = settingsSkipSymmetricOverlaps;
        settings.map.bestNSecondary = -1;
        settings.map.secondaryMinScoreFraction = 0.0;
        settings.map.secondaryAllowedOverlapFractionTarget = 0.0;
        settings.map.secondaryAllowedOverlapFractionQuery = 0.0;
        settings.map.skipSymmetricOverlaps = true;
        settings.map.selfHitPolicy = PacBio::Pancake::MapperSelfHitPolicy::DEFAULT;
        settings.align.selfHitPolicy = PacBio::Pancake::MapperSelfHitPolicy::DEFAULT;
        settings.map.freqPercentile = 0.000;
        settings.map.seedParams = PacBio::Pancake::SeedDBParameters{19, 10, 0, false, false, true};
        settings.map.seedParamsFallback =
            PacBio::Pancake::SeedDBParameters{10, 5, 0, false, false, true};
        settings.align.alignerTypeGlobal = PacBio::Pancake::AlignerType::EDLIB;
    }

    PacBio::Pancake::MapperCLRSettings settingsSkipSelfAndSymmetricOverlaps;
    {
        auto& settings = settingsSkipSelfAndSymmetricOverlaps;
        settings.map.bestNSecondary = -1;
        settings.map.secondaryMinScoreFraction = 0.0;
        settings.map.secondaryAllowedOverlapFractionTarget = 0.0;
        settings.map.secondaryAllowedOverlapFractionQuery = 0.0;
        settings.map.skipSymmetricOverlaps = true;
        settings.map.selfHitPolicy = PacBio::Pancake::MapperSelfHitPolicy::SKIP;
        settings.align.selfHitPolicy = PacBio::Pancake::MapperSelfHitPolicy::SKIP;
        settings.map.freqPercentile = 0.000;
        settings.map.seedParams = PacBio::Pancake::SeedDBParameters{19, 10, 0, false, false, true};
        settings.map.seedParamsFallback =
            PacBio::Pancake::SeedDBParameters{10, 5, 0, false, false, true};
        settings.align.alignerTypeGlobal = PacBio::Pancake::AlignerType::EDLIB;
    }

    PacBio::Pancake::MapperCLRSettings settingsSkipSelfInMappingButDefaultInAlignment;
    {
        auto& settings = settingsSkipSelfInMappingButDefaultInAlignment;
        settings.map.bestNSecondary = -1;
        settings.map.secondaryMinScoreFraction = 0.0;
        settings.map.secondaryAllowedOverlapFractionTarget = 0.0;
        settings.map.secondaryAllowedOverlapFractionQuery = 0.0;
        settings.map.skipSymmetricOverlaps = false;
        settings.map.selfHitPolicy = PacBio::Pancake::MapperSelfHitPolicy::SKIP;
        settings.align.selfHitPolicy = PacBio::Pancake::MapperSelfHitPolicy::DEFAULT;
        settings.map.freqPercentile = 0.000;
        settings.map.seedParams = PacBio::Pancake::SeedDBParameters{19, 10, 0, false, false, true};
        settings.map.seedParamsFallback =
            PacBio::Pancake::SeedDBParameters{10, 5, 0, false, false, true};
        settings.align.alignerTypeGlobal = PacBio::Pancake::AlignerType::EDLIB;
    }

    PacBio::Pancake::MapperCLRSettings settingsDefaultSelfInMappingButSkipInAlignment;
    {
        auto& settings = settingsDefaultSelfInMappingButSkipInAlignment;
        settings.map.bestNSecondary = 100;
        settings.map.secondaryMinScoreFraction = 0.0;
        settings.map.secondaryAllowedOverlapFractionTarget = 0.0;
        settings.map.secondaryAllowedOverlapFractionQuery = 0.0;
        settings.map.skipSymmetricOverlaps = false;
        settings.map.selfHitPolicy = PacBio::Pancake::MapperSelfHitPolicy::DEFAULT;
        settings.align.selfHitPolicy = PacBio::Pancake::MapperSelfHitPolicy::SKIP;
        settings.map.freqPercentile = 0.000;
        settings.map.seedParams = PacBio::Pancake::SeedDBParameters{19, 10, 0, false, false, true};
        settings.map.seedParamsFallback =
            PacBio::Pancake::SeedDBParameters{10, 5, 0, false, false, true};
        settings.align.alignerTypeGlobal = PacBio::Pancake::AlignerType::EDLIB;
    }

    PacBio::Pancake::MapperCLRSettings settingsMockSelfInMappingButDefaultInAlignment;
    {
        auto& settings = settingsMockSelfInMappingButDefaultInAlignment;
        settings.map.bestNSecondary = -1;
        settings.map.secondaryMinScoreFraction = 0.0;
        settings.map.secondaryAllowedOverlapFractionTarget = 0.0;
        settings.map.secondaryAllowedOverlapFractionQuery = 0.0;
        settings.map.skipSymmetricOverlaps = false;
        settings.map.selfHitPolicy = PacBio::Pancake::MapperSelfHitPolicy::PERFECT_ALIGNMENT;
        settings.align.selfHitPolicy = PacBio::Pancake::MapperSelfHitPolicy::DEFAULT;
        settings.map.freqPercentile = 0.000;
        settings.map.seedParams = PacBio::Pancake::SeedDBParameters{19, 10, 0, false, false, true};
        settings.map.seedParamsFallback =
            PacBio::Pancake::SeedDBParameters{10, 5, 0, false, false, true};
        settings.align.alignerTypeGlobal = PacBio::Pancake::AlignerType::EDLIB;
    }

    PacBio::Pancake::MapperCLRSettings settingsDefaultSelfInMappingButMockInAlignment;
    {
        auto& settings = settingsDefaultSelfInMappingButMockInAlignment;
        settings.map.bestNSecondary = -1;
        settings.map.secondaryMinScoreFraction = 0.0;
        settings.map.secondaryAllowedOverlapFractionTarget = 0.0;
        settings.map.secondaryAllowedOverlapFractionQuery = 0.0;
        settings.map.skipSymmetricOverlaps = false;
        settings.map.selfHitPolicy = PacBio::Pancake::MapperSelfHitPolicy::DEFAULT;
        settings.align.selfHitPolicy = PacBio::Pancake::MapperSelfHitPolicy::PERFECT_ALIGNMENT;
        settings.map.freqPercentile = 0.000;
        settings.map.seedParams = PacBio::Pancake::SeedDBParameters{19, 10, 0, false, false, true};
        settings.map.seedParamsFallback =
            PacBio::Pancake::SeedDBParameters{10, 5, 0, false, false, true};
        settings.align.alignerTypeGlobal = PacBio::Pancake::AlignerType::EDLIB;
    }

    struct TestData
    {
        const std::string testName;
        // Map settings are passed for each batch separately.
        const std::vector<
            std::tuple<std::string, std::string, PacBio::Pancake::MapperCLRMapSettings>>
            batchData;
        const int32_t sequenceIdOffset = 0;
        const PacBio::Pancake::MapperCLRAlignSettings alignSettings;
        const std::vector<std::string> expectedOverlapsPaths;
    };

    // clang-format off
    std::vector<TestData> testData = {
        {
            "Overlap the same set of reads with itself.",
            // Input batch data.
            {
                {
                    PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.fasta",
                    PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.fasta",
                    settingsDefaultPolicy.map,
                },
            },
            // Sequence ID offset. Zero means that the first query/target after loading has the ID == 0.
            0,
            // Input alignment settings.
            settingsDefaultPolicy.align,
            // Expected results.
            {
                PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.out.all_vs_all.edlib.m4",
            },
        },
        {
            "Overlap the same set of reads with itself. Offset the IDs by 10000 to check that IDs can be arbitrary.",
            // Input batch data.
            {
                {
                    PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.fasta",
                    PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.fasta",
                    settingsDefaultPolicy.map,
                },
            },
            // Sequence ID offset. Zero means that the first query/target after loading has the ID == 0.
            10000,
            // Input alignment settings.
            settingsDefaultPolicy.align,
            // Expected results.
            {
                PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.out.all_vs_all.edlib.id_offset_10000.m4",
            },
        },
        {
            "Skip self hits.",
            // Input batch data.
            {
                {
                    PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.fasta",
                    PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.fasta",
                    settingsSkipSelfHitsInBothMapAndAlign.map,
                },
            },
            // Sequence ID offset. Zero means that the first query/target after loading has the ID == 0.
            0,
            // Input alignment settings.
            settingsSkipSelfHitsInBothMapAndAlign.align,
            // Expected results.
            {
                PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.out.all_vs_all_no_self_hits.edlib.m4",
            },
        },
        {
            "Mock perfect overlaps",
            // Input batch data.
            {
                {
                    PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.fasta",
                    PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.fasta",
                    settingsPerfectAlignSelfHitsInBothMapAndAlign.map,
                },
            },
            // Sequence ID offset. Zero means that the first query/target after loading has the ID == 0.
            0,
            // Input alignment settings.
            settingsPerfectAlignSelfHitsInBothMapAndAlign.align,
            // Expected results.
            {
                PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.out.all_vs_all.edlib.m4",
            },
        },
        {
            "Skip symmetric overlaps.",
            // Input batch data.
            {
                {
                    PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.fasta",
                    PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.fasta",
                    settingsSkipSymmetricOverlaps.map,
                },
            },
            // Sequence ID offset. Zero means that the first query/target after loading has the ID == 0.
            0,
            // Input alignment settings.
            settingsSkipSymmetricOverlaps.align,
            // Expected results.
            {
                PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.out.all_vs_all_no_symmetric.edlib.m4",
            },
        },

        {
            "Skip self and symmetric overlaps.",
            // Input batch data.
            {
                {
                    PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.fasta",
                    PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.fasta",
                    settingsSkipSelfAndSymmetricOverlaps.map,
                },
            },
            // Sequence ID offset. Zero means that the first query/target after loading has the ID == 0.
            0,
            // Input alignment settings.
            settingsSkipSelfAndSymmetricOverlaps.align,
            // Expected results.
            {
                PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.out.all_vs_all_no_self_hits_no_symmetric.edlib.m4",
            },
        },

        {
            "Skip self hits in the mapping stage, but use the default policy during alignment. This should skip the self hits completely.",
            // Input batch data.
            {
                {
                    PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.fasta",
                    PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.fasta",
                    settingsSkipSelfInMappingButDefaultInAlignment.map,
                },
            },
            // Sequence ID offset. Zero means that the first query/target after loading has the ID == 0.
            0,
            // Input alignment settings.
            settingsSkipSelfInMappingButDefaultInAlignment.align,
            // Expected results.
            {
                PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.out.all_vs_all_no_self_hits.edlib.m4",
            },
        },

        // This test case discovered a bug - because alignment stage cleared self-hits, the WrapFlagSecondaryAndSupplementary function
        // caused a skew in the IDs of its input overlaps and the internal tmpOverlaps which don't contain nullptr overlaps.
        {
            "Skip self hits in the alignment stage, but use the default policy during mapping. This should skip the self hits completely.",
            // Input batch data.
            {
                {
                    PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.fasta",
                    PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.fasta",
                    settingsDefaultSelfInMappingButSkipInAlignment.map,
                },
            },
            // Sequence ID offset. Zero means that the first query/target after loading has the ID == 0.
            0,
            // Input alignment settings.
            settingsDefaultSelfInMappingButSkipInAlignment.align,
            // Expected results.
            {
                PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.out.all_vs_all_no_self_hits.edlib.m4",
            },
        },
        {
            "Mock perfect overlaps in the mapping stage, but use the default policy during alignment. This should report proper alignments, like everything was default.",
            // Input batch data.
            {
                {
                    PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.fasta",
                    PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.fasta",
                    settingsMockSelfInMappingButDefaultInAlignment.map,
                },
            },
            // Sequence ID offset. Zero means that the first query/target after loading has the ID == 0.
            0,
            // Input alignment settings.
            settingsMockSelfInMappingButDefaultInAlignment.align,
            // Expected results.
            {
                PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.out.all_vs_all.edlib.m4",
            },
        },

        {
            "Mock perfect overlaps in the alignment stage, but use the default policy during mapping. This should report proper alignments, like everything was default.",
            // Input batch data.
            {
                {
                    PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.fasta",
                    PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.fasta",
                    settingsDefaultSelfInMappingButMockInAlignment.map,
                },
            },
            // Sequence ID offset. Zero means that the first query/target after loading has the ID == 0.
            0,
            // Input alignment settings.
            settingsDefaultSelfInMappingButMockInAlignment.align,
            // Expected results.
            {
                PacBio::PancakeTestsConfig::Data_Dir + "/ovl-clr/reads.pile1-5prime.out.all_vs_all.edlib.m4",
            },
        },
    };
    // clang-format on

    for (const auto& data : testData) {
        // Debug info.
        SCOPED_TRACE(data.testName);
        std::cerr << "testName = " << data.testName << "\n";

        // Load the batch sequence data. The helper function takes
        // a vector of target-query filename pairs.
        std::vector<PacBio::Pancake::MapperBatchChunk> batchData;
        std::vector<PacBio::BAM::FastaSequence> allSeqs;
        PacBio::PancakeTests::HelperLoadBatchData(data.batchData, data.sequenceIdOffset, batchData,
                                                  allSeqs);

        // Create a mapper.
        const PacBio::Pancake::MapperCLRAlignSettings& alignSettings = data.alignSettings;
        PacBio::Pancake::MapperBatchCPU mapper(alignSettings, 1);

        // Run the unit under test.
        // std::vector<std::vector<MapperBaseResult>> results = mapper.DummyMapAndAlign(batchData);
        std::vector<std::vector<PacBio::Pancake::MapperBaseResult>> results =
            mapper.MapAndAlign(batchData);

        // Format the results for comparison.
        std::vector<std::vector<std::string>> resultsStr =
            PacBio::PancakeTests::HelperFormatBatchMappingResults(results);

        // Sort the results for comparison.
        for (auto& chunkResults : resultsStr) {
            std::sort(chunkResults.begin(), chunkResults.end());
        }

        // Prepare expected results, and sort them.
        std::vector<std::vector<std::string>> expectedOverlaps;
        for (const auto& singleExpectedPath : data.expectedOverlapsPaths) {
            std::vector<std::string> temp =
                PacBio::PancakeTests::HelperLoadFile(singleExpectedPath);
            std::sort(temp.begin(), temp.end());
            expectedOverlaps.emplace_back(std::move(temp));
        }

        // std::cerr << "Expected:\n";
        // for (size_t i = 0; i < data.expectedOverlapsPaths.size(); ++i) {
        //     const auto& singleExpectedOverlaps = expectedOverlaps[i];
        //     std::cerr << "  - Expected path: " << data.expectedOverlapsPaths[i] << "\n";
        //     for (const auto& ovlStr : singleExpectedOverlaps) {
        //         std::cerr << "    " << ovlStr << "\n";
        //     }
        //     std::cerr << "  - Results:\n";
        //     for (const auto& ovlStr : resultsStr[i]) {
        //         std::cerr << "    " << ovlStr << "\n";
        //     }
        //     if (singleExpectedOverlaps.size() != resultsStr[i].size()) {
        //         std::cerr << "  - Sizes differ!\n";
        //     } else {
        //         for (size_t j = 0; j < resultsStr[i].size(); ++j) {
        //             if (resultsStr[i][j] != singleExpectedOverlaps[j]) {
        //                 std::cerr << "  - [j = " << j << " / " << resultsStr[i].size()
        //                           << "] Different result:\n";
        //                 std::cerr << "      Expected: " << singleExpectedOverlaps[j] << "\n";
        //                 std::cerr << "      Result:   " << resultsStr[i][j] << "\n";
        //             }
        //         }
        //     }
        //     std::cerr << "\n";
        // }

        // Evaluate.
        EXPECT_EQ(expectedOverlaps, resultsStr);

        // Test the mocking flags.
        {
            // Prepare the data.
            const auto [expectedMocked, resultsMocked] =
                PacBio::PancakeTests::HelperFormatBatchMappingResultsForMockingFlags(
                    results, data.batchData, data.alignSettings.selfHitPolicy);

            // Evaluate.
            EXPECT_EQ(expectedMocked, resultsMocked);
        }
    }
}

TEST(MapperBatchCPU, EdelweissEdgeCases_ArrayOfTests)
{
    using namespace PacBio::Pancake;

    struct TestData
    {
        const std::string testName;
        std::vector<std::pair<std::string, std::string>> batchData;
        const int32_t sequenceIdOffset = 0;
        const PacBio::Pancake::AlignerType alignerTypeGlobal;
        const PacBio::Pancake::SeedDBParameters seedParamsPrimary;
        const PacBio::Pancake::SeedDBParameters seedParamsFallback;
        const std::vector<std::vector<std::string>> expectedOverlaps;
    };

    // Create 12 batches of 1 pair.
    std::vector<TestData> testData = {
        {
            "Batch of 12 sequence pairs to test Edelweiss edge cases.",
            {},
            // Sequence ID offset. Zero means that the first query/target after loading has the ID == 0.
            1234567,
            // Aligner type for global alignment.
            AlignerType::KSW2,
            // SeedParams - primary.
            PacBio::Pancake::SeedDBParameters{15, 5, 0, false, false, true},
            // SeedParams - fallback.
            PacBio::Pancake::SeedDBParameters{10, 5, 0, false, false, true},
            // Expected results.
            {
                {"001234567 001234567 1142 69.12 0 0 2102 2102 0 0 1741 1741 *"},
                {"001234567 001234567 8490 75.94 0 0 11821 11821 0 0 11316 11316 *"},
                {"001234567 001234567 11524 79.13 0 0 13333 13333 0 0 13151 13151 *"},
                {"001234567 001234567 11724 82.36 0 0 11152 11152 0 0 10992 10992 *"},
                {"001234567 001234567 4769 76.75 0 0 5799 5799 0 0 5444 5444 *"},
                {"001234567 001234567 10698 78.69 0 0 11800 11800 0 0 11540 11540 *"},
                {"001234567 001234567 7461 76.19 0 0 11050 11050 0 0 10567 10567 *"},
                {"001234567 001234567 5203 70.09 0 0 11879 11879 0 0 10568 10568 *"},
                {"001234567 001234567 9215 75.03 0 0 13551 13551 0 0 13408 13408 *"},
                {"001234567 001234567 10026 76.84 0 0 13787 13787 0 0 13406 13406 *"},
                {"001234567 001234567 11188 80.94 0 0 12167 12167 0 0 12032 12032 *"},
                {"001234567 001234567 8378 76.91 0 0 12616 12616 0 0 12032 12032 *"},
            },
        },
    };
    for (size_t i = 1; i < 13; ++i) {
        const std::pair<std::string, std::string> inPair = {
            PacBio::PancakeTestsConfig::Data_Dir + "/edelweiss/edge-case-1-bw/test." +
                std::to_string(i) + ".target.fasta",
            PacBio::PancakeTestsConfig::Data_Dir + "/edelweiss/edge-case-1-bw/test." +
                std::to_string(i) + ".query.fasta",
        };
        testData[0].batchData.emplace_back(inPair);
    }

    // Initialize settings and parameters.
    const PacBio::Pancake::MapperCLRSettings settings =
        PacBio::PancakeTests::HelperInitPancakeSettingsSubread();
    const int32_t numThreads = 1;

    for (const auto& data : testData) {
        // Debug info.
        SCOPED_TRACE(data.testName);
        std::cerr << "testName = " << data.testName << "\n";

        // Load the batch sequence data. The helper function takes
        // a vector of target-query filename pairs.
        std::vector<MapperBatchChunk> batchData;
        std::vector<PacBio::BAM::FastaSequence> allSeqs;
        PacBio::PancakeTests::HelperLoadBatchData(data.batchData, data.sequenceIdOffset,
                                                  settings.map, batchData, allSeqs);

        PacBio::Pancake::MapperBatchCPU mapper(settings.align, numThreads);

        // Run the unit under test.
        std::vector<std::vector<MapperBaseResult>> results = mapper.MapAndAlign(batchData);

        // Format the results for comparison.
        std::vector<std::vector<std::string>> resultsStr =
            PacBio::PancakeTests::HelperFormatBatchMappingResults(results);

        // Evaluate.
        EXPECT_EQ(data.expectedOverlaps, resultsStr);
    }
}