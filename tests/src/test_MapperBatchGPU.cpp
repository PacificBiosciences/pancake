#include <PancakeTestData.h>
#include <gtest/gtest.h>
#include <pacbio/pancake/MapperBatchGPU.h>
#include <pacbio/pancake/OverlapWriterBase.h>
#include <iostream>
#include "TestHelperUtils.h"

TEST(MapperBatchGPU, BatchMapping_ArrayOfTests)
{
    using namespace PacBio::Pancake;

    struct TestData
    {
        const std::string testName;
        const std::vector<std::pair<std::string, std::string>> batchData;
        const int32_t sequenceIdOffset = 0;
        const PacBio::Pancake::AlignerType alignerTypeGlobal;
        const PacBio::Pancake::SeedDB::SeedDBParameters seedParamsPrimary;
        const PacBio::Pancake::SeedDB::SeedDBParameters seedParamsFallback;
        const std::vector<std::vector<std::string>> expectedOverlaps;
    };

    // clang-format off
    std::vector<TestData> testData = {
        {
            "Batch of multiple query/target vectors.",
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
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, false, true},
            // SeedParams - fallback.
            PacBio::Pancake::SeedDB::SeedDBParameters{10, 5, 0, false, false, true},
            // Expected results.
            {
                {
                    "000000000 000000000 19072 81.17 0 0 18779 18779 0 0 18864 18865 *",
                },
                {
                    "000000000 000000000 15726 98.71 0 0 8111 8111 0 0 8138 8138 *"
                },
                {
                    "000000000 000000000 300 100.00 0 0 150 150 1 0 150 150 *"
                },
                {
                    "000000000 000000000 16112 75.43 0 0 21382 22015 1 701 21992 22028 *"
                },
                {
                    "000000000 000000000 280 96.77 0 0 155 155 1 0 150 150 *"
                },
                {
                    "000000000 000000000 284 68.58 0 6209 7938 43446 0 7261 8999 46238 *"
                },
                {
                    "000000000 000000000 11218 75.98 0 0 15753 15753 1 2 15953 15953 *"
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
            PacBio::Pancake::SeedDB::SeedDBParameters{15, 5, 0, false, true, true},
            // SeedParams - fallback.
            PacBio::Pancake::SeedDB::SeedDBParameters{10, 5, 0, false, false, true},
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
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, false, true},
            // SeedParams - fallback.
            PacBio::Pancake::SeedDB::SeedDBParameters{10, 5, 0, false, false, true},
            // Expected results.
            {
                {
                    "001234567 001234567 19072 81.17 0 0 18779 18779 0 0 18864 18865 *",
                },
                {
                    "001234567 001234567 15726 98.71 0 0 8111 8111 0 0 8138 8138 *"
                },
                {
                    "001234567 001234567 300 100.00 0 0 150 150 1 0 150 150 *"
                },
                {
                    "001234567 001234567 16112 75.43 0 0 21382 22015 1 701 21992 22028 *"
                },
                {
                    "001234567 001234567 280 96.77 0 0 155 155 1 0 150 150 *"
                },
                {
                    "001234567 001234567 284 68.58 0 6209 7938 43446 0 7261 8999 46238 *"
                },
                {
                    "001234567 001234567 11218 75.98 0 0 15753 15753 1 2 15953 15953 *"
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

        // Set the alignment parameters.
        PacBio::Pancake::MapperCLRAlignSettings alignSettings;
        alignSettings.alignerTypeGlobal = data.alignerTypeGlobal;

        const uint32_t gpuDeviceId = 0;
        const int64_t gpuMaxMemoryCap =
            static_cast<int64_t>(100) * static_cast<int64_t>(1024 * 1024);
        const int32_t numThreads = 1;
        const int32_t startBandwidth = 500;
        const int32_t maxBandwidth = 2000;
        bool alignRemainingOnCpu = false;
        PacBio::Pancake::MapperBatchGPU mapper(alignSettings, numThreads, startBandwidth,
                                               maxBandwidth, gpuDeviceId, gpuMaxMemoryCap,
                                               alignRemainingOnCpu);

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

TEST(MapperBatchGPU, CheckSelfHitPolicyAndSkippingSymmetrical)
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
        settings.map.seedParams =
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, false, true};
        settings.map.seedParamsFallback =
            PacBio::Pancake::SeedDB::SeedDBParameters{10, 5, 0, false, false, true};
        // settings.align.alignerTypeGlobal = AlignerType::EDLIB;
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
        settings.map.seedParams =
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, false, true};
        settings.map.seedParamsFallback =
            PacBio::Pancake::SeedDB::SeedDBParameters{10, 5, 0, false, false, true};
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
        settings.map.seedParams =
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, false, true};
        settings.map.seedParamsFallback =
            PacBio::Pancake::SeedDB::SeedDBParameters{10, 5, 0, false, false, true};
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
        settings.map.seedParams =
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, false, true};
        settings.map.seedParamsFallback =
            PacBio::Pancake::SeedDB::SeedDBParameters{10, 5, 0, false, false, true};
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
        settings.map.seedParams =
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, false, true};
        settings.map.seedParamsFallback =
            PacBio::Pancake::SeedDB::SeedDBParameters{10, 5, 0, false, false, true};
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
        settings.map.seedParams =
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, false, true};
        settings.map.seedParamsFallback =
            PacBio::Pancake::SeedDB::SeedDBParameters{10, 5, 0, false, false, true};
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
        settings.map.seedParams =
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, false, true};
        settings.map.seedParamsFallback =
            PacBio::Pancake::SeedDB::SeedDBParameters{10, 5, 0, false, false, true};
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
        settings.map.seedParams =
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, false, true};
        settings.map.seedParamsFallback =
            PacBio::Pancake::SeedDB::SeedDBParameters{10, 5, 0, false, false, true};
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
        settings.map.seedParams =
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, false, true};
        settings.map.seedParamsFallback =
            PacBio::Pancake::SeedDB::SeedDBParameters{10, 5, 0, false, false, true};
    }

    struct TestData
    {
        const std::string testName;
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
            "Overlap the same set of reads with itself. Offset the IDs by 10000.",
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

        const PacBio::Pancake::MapperCLRAlignSettings& alignSettings = data.alignSettings;

        // Create the mapper.
        const uint32_t gpuDeviceId = 0;
        const int64_t gpuMaxMemoryCap =
            static_cast<int64_t>(100) * static_cast<int64_t>(1024 * 1024);
        const int32_t numThreads = 2;
        const int32_t startBandwidth = 500;
        const int32_t maxBandwidth = 2000;
        bool alignRemainingOnCpu = false;
        PacBio::Pancake::MapperBatchGPU mapper(alignSettings, numThreads, startBandwidth,
                                               maxBandwidth, gpuDeviceId, gpuMaxMemoryCap,
                                               alignRemainingOnCpu);

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
    }
}
