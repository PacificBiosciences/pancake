#include <PancakeTestData.h>
#include <gtest/gtest.h>
#include <pacbio/pancake/MapperCLR.h>
#include <pacbio/pancake/OverlapWriterBase.h>
#include <iostream>
#include "TestHelperUtils.h"

TEST(MapperCLR, CheckMappping_LoadFromFile)
{
    struct TestData
    {
        std::string testName;
        std::string targetFile;
        std::string queryFile;
        PacBio::Pancake::SeedDB::SeedDBParameters seedParamsPrimary;
        PacBio::Pancake::SeedDB::SeedDBParameters seedParamsFallback;
        std::vector<std::string> expectedOverlaps;
    };

    // clang-format off
    std::vector<TestData> testData = {
        {
            "Pair of real subreads. CCS.",
            // Target.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-2-real-insert-target.fasta",
            // Query.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-2-real-insert-query.fasta",
            // SeedParams - primary.
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, false, 255, true},
            // SeedParams - fallback.
            PacBio::Pancake::SeedDB::SeedDBParameters{10, 5, 0, false, false, 255, true},
            // Expected results.
            {
                "000000000 000000000 -16391 81.12 0 0 18779 18779 0 0 18864 18865 *"
            },
        },
        {
            "RealHifiReads, subsequences with a full global match.",
            // Target.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-3-real-subseq-full-global-match-target.fasta",
            // Query.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-3-real-subseq-full-global-match-query.fasta",
            // SeedParams - primary.
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, false, 255, true},
            // SeedParams - fallback.
            PacBio::Pancake::SeedDB::SeedDBParameters{10, 5, 0, false, false, 255, true},
            // Expected results.
            {
                "000000000 000000000 -8049 98.71 0 0 8111 8111 0 0 8138 8138 *"
            },
        },
        {
            "Small reverse complement perfect alignment",
            // Target.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-4-small-revcmp-perfect-aln-target.fasta",
            // Query.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-4-small-revcmp-perfect-aln-query.fasta",
            // SeedParams - primary.
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, true, 255, true},
            // SeedParams - fallback.
            PacBio::Pancake::SeedDB::SeedDBParameters{10, 5, 0, false, true, 255, true},
            // Expected results.
            {
                "000000000 000000000 -150 100.00 0 0 150 150 1 0 150 150 *"
            },
        },
        {
            "Reverse complement alignment. This one was challenging for the AlignmentSeeded implementation, because offseting the reverse coordinates had a bug.",
            // Target.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-5-revcmp-aln-target.fasta",
            // Query.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-5-revcmp-aln-query.fasta",
            // SeedParams - primary.
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, false, 255, true},
            // SeedParams - fallback.
            PacBio::Pancake::SeedDB::SeedDBParameters{10, 5, 0, false, false, 255, true},
            // Expected results.
            {
                "000000000 000000000 -17688 74.58 0 0 21382 22015 1 701 22028 22028 *"
            },
        },
        {
            "Small reverse complement imperfect alignment. There are 5 insertions sprinkled around, all 'A' bases.",
            // Target.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-6-small-revcmp-imperfect-aln-target.fasta",
            // Query.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-6-small-revcmp-imperfect-aln-query.fasta",
            // SeedParams - primary.
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, true, 255, true},
            // SeedParams - fallback.
            PacBio::Pancake::SeedDB::SeedDBParameters{10, 5, 0, false, true, 255, true},
            // Expected results.
            {
                "000000000 000000000 -150 96.77 0 0 155 155 1 0 150 150 *"
            },
        },
        {
            "Test seed fallback. Without fallback, using large parameters (k = 19, w = 10, spacing = 0, HPC = true) will not map.",
            // Target.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-7-seed-fallback-target.fasta",
            // Query.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-7-seed-fallback-query.fasta",
            // SeedParams - primary.
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, true, 255, true},
            // SeedParams - fallback.
            PacBio::Pancake::SeedDB::SeedDBParameters{10, 5, 0, false, true, 255, true},
            // Expected results.
            {
                "000000000 000000000 -984 71.07 0 0 1241 1241 0 0 1225 1225 *"
            },
        },
        {
            "The same test but without seed fallback. Without fallback, using large parameters (k = 19, w = 10, spacing = 0, HPC = true) will not map.",
            // Target.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-7-seed-fallback-target.fasta",
            // Query.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-7-seed-fallback-query.fasta",
            // SeedParams - primary.
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, true, 255, true},
            // SeedParams - fallback - same as primary.
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, true, 255, true},
            // Expected results.
            {
            },
        },
        {
            "RealHifiReads, poor matching pair of reads. CCS.",
            // Target.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-1-poor-aln-overlapping-seeds.target.fasta",
            // Query.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-1-poor-aln-overlapping-seeds.query.fasta",
            // SeedParams - primary.
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, false, 255, true},
            // SeedParams - fallback.
            PacBio::Pancake::SeedDB::SeedDBParameters{10, 5, 0, false, false, 255, true},
            // Expected results.
            {
                "000000000 000000000 -1345 68.58 0 6209 7938 43446 0 7261 8999 46238 *"
            },
        },

    };
    // clang-format on

    PacBio::Pancake::MapperCLRSettings settings;
    settings.freqPercentile = 0.000;

    for (const auto& data : testData) {
        // Load the sequences from files, and take only the first one.
        const std::vector<PacBio::BAM::FastaSequence> allTargetSeqs =
            PacBio::PancakeTests::HelperLoadFasta(data.targetFile);
        const std::vector<PacBio::BAM::FastaSequence> allQuerySeqs =
            PacBio::PancakeTests::HelperLoadFasta(data.queryFile);
        const std::string& target = allTargetSeqs[0].Bases();
        const std::string& query = allQuerySeqs[0].Bases();

        SCOPED_TRACE(data.testName);

        std::cerr << "testName = " << data.testName << "\n";

        settings.seedParams = data.seedParamsPrimary;
        settings.seedParamsFallback = data.seedParamsFallback;
        PacBio::Pancake::MapperCLR mapper(settings);

        std::vector<PacBio::Pancake::MapperBaseResult> result =
            mapper.MapAndAlign({target}, {query});
        std::vector<std::string> resultsStr;
        for (const auto& queryMappings : result) {
            for (const auto& mapping : queryMappings.mappings) {
                std::cerr << PacBio::Pancake::OverlapWriterBase::PrintOverlapAsM4(
                                 mapping->mapping, "", "", true, false)
                          << "\n";

                resultsStr.emplace_back(PacBio::Pancake::OverlapWriterBase::PrintOverlapAsM4(
                    mapping->mapping, "", "", true, false));
            }
        }

        ASSERT_EQ(data.expectedOverlaps, resultsStr);
    }
}

void DebugPrintChainedRegion(std::ostream& oss, int32_t regionId,
                             const PacBio::Pancake::ChainedRegion& cr)
{
    oss << "[regionId " << regionId << "] chain.hits = " << cr.chain.hits.size()
        << ", chain.score = " << cr.chain.score << ", chain.covQ = " << cr.chain.coveredBasesQuery
        << ", chain.covT = " << cr.chain.coveredBasesTarget << ", priority = " << cr.priority
        << ", isSuppl = " << (cr.isSupplementary ? "true" : "false") << ", ovl: "
        << PacBio::Pancake::OverlapWriterBase::PrintOverlapAsM4(cr.mapping, "", "", true, false)
        << ", diagStart = " << (cr.mapping->Astart - cr.mapping->Bstart)
        << ", diagEnd = " << (cr.mapping->Aend - cr.mapping->Bend);

    oss << "\n";
    for (size_t j = 0; j < cr.chain.hits.size(); ++j) {
        const auto& hit = cr.chain.hits[j];
        std::cerr << "    [hit " << j << "] " << hit << "\n";
    }
}

TEST(MapperCLR, CheckMappingAndSeedHits)
{
    using namespace PacBio::Pancake;

    struct TestData
    {
        std::string testName;
        std::string targetFile;
        std::string queryFile;
        PacBio::Pancake::SeedDB::SeedDBParameters seedParamsPrimary;
        PacBio::Pancake::SeedDB::SeedDBParameters seedParamsFallback;
        std::vector<std::vector<SeedHit>> expectedSeedHits;
        std::vector<std::string> expectedMappings;
    };

    // clang-format off
    std::vector<TestData> testData = {
        {
            "RealHifiReads_Subsequences",
            // Target
            PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-1-poor-aln-overlapping-seeds.target.fasta",
            // Query.
            PacBio::PancakeTestsConfig::Data_Dir + "/mapper-clr/test-1-poor-aln-overlapping-seeds.query.fasta",
            // SeedParams - primary
            PacBio::Pancake::SeedDB::SeedDBParameters{19, 10, 0, false, true, 100, true},
            // SeedParams - fallback
            PacBio::Pancake::SeedDB::SeedDBParameters{15, 5, 0, false, true, 100, true},
            ///// These were the old 'wrong' mappings, where we had two seed hits which had the
            ///// same X-coordinate.
            // // Expected seed hits.
            // {
            //     {
            //         SeedHit(0, false, 38493, 6585, 30, 25, 0),
            //         SeedHit(0, false, 40327, 6585, 24, 25, 0),
            //     },
            //     {
            //         SeedHit(0, false, 11779, 29565, 28, 29, 0),
            //         SeedHit(0, false, 11791, 29577, 25, 25, 0),
            //         SeedHit(0, false, 11849, 29637, 24, 27, 0),
            //     },
            // },
            // // Expected mappings.
            // {
            //     "000000000 000000000 -862 25.72 0 6468 7589 43446 0 38370 41325 46238 *",
            //     "000000000 000000000 -183 73.89 0 29527 29753 43446 0 11747 11966 46238 *"
            // },
            /////
            ///// This is what the seed hits and mappings look like after the LongMergeChains_ strictly
            ///// disallowed the overlap, even at the same position.
            // Expected seed hits.
            {
                {
                    SeedHit(0, false, 7919, 6832, 28, 28, 0),
                    SeedHit(0, false, 11021, 9922, 29, 30, 0),
                    SeedHit(0, false, 11401, 10312, 25, 24, 0),
                },
                {
                    SeedHit(0, false, 37479, 5606, 26, 24, 0),
                    SeedHit(0, false, 37486, 5612, 26, 25, 0),
                },
            },
            // Expected mappings.
            {
                "000000000 000000000 -4710 67.66 0 6209 12358 43446 0 7261 13474 46238 *",
                "000000000 000000000 -525 71.66 0 5287 5940 43446 0 37146 37820 46238 *",
            },


        },
    };
    // clang-format on

    PacBio::Pancake::MapperCLRSettings settings;
    settings.freqPercentile = 0.000;

    for (const auto& data : testData) {
        // Load the sequences from files, and take only the first one.
        const std::vector<PacBio::BAM::FastaSequence> allTargetSeqs =
            PacBio::PancakeTests::HelperLoadFasta(data.targetFile);
        const std::vector<PacBio::BAM::FastaSequence> allQuerySeqs =
            PacBio::PancakeTests::HelperLoadFasta(data.queryFile);
        const std::string& target = allTargetSeqs[0].Bases();
        const std::string& query = allQuerySeqs[0].Bases();

        SCOPED_TRACE(data.testName);

        std::cerr << "testName = " << data.testName << "\n";

        settings.seedParams = data.seedParamsPrimary;
        settings.seedParamsFallback = data.seedParamsFallback;
        // settings.align = false;
        PacBio::Pancake::MapperCLR mapper(settings);

        std::vector<PacBio::Pancake::MapperBaseResult> result =
            mapper.MapAndAlign({target}, {query});
        std::vector<std::vector<PacBio::Pancake::SeedHit>> resultSeedHits;
        std::vector<std::string> resultsMappings;
        for (const auto& queryMappings : result) {
            for (const auto& mapping : queryMappings.mappings) {
                std::cerr << PacBio::Pancake::OverlapWriterBase::PrintOverlapAsM4(
                                 mapping->mapping, "", "", true, false)
                          << "\n";

                // DebugPrintChainedRegion(std::cerr, 0, *mapping);
                resultSeedHits.emplace_back(mapping->chain.hits);

                resultsMappings.emplace_back(PacBio::Pancake::OverlapWriterBase::PrintOverlapAsM4(
                    mapping->mapping, "", "", true, false));
            }
        }

        EXPECT_EQ(data.expectedSeedHits, resultSeedHits);
        EXPECT_EQ(data.expectedMappings, resultsMappings);
    }

    // exit(1);
}
