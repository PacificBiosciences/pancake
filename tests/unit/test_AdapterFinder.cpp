#include <pancake/AdapterFinder.hpp>

#include <PancakeTestData.h>
#include "TestHelperUtils.hpp"

#include <gtest/gtest.h>
#include <pancake/OverlapWriterBase.hpp>

#include <pbcopper/logging/Logging.h>

#include <iostream>
#include <memory>
#include <string>

TEST(AdapterFinder, ArrayOfTests_LoadFromFile)
{
    struct TestData
    {
        const std::string TestName;
        const std::string InReadsFn;
        const int32_t KmerSize = 15;
        const int32_t WindowSize = 2;
        const int32_t MedianLen = 0;
        const std::vector<PacBio::Pancake::Adapter::AdapterPosition> ExpectedAdapters;
    };

    const PacBio::Pancake::SeedDBParameters seedParams1{15, 2, 0, false, false, true};

    // clang-format off
    std::vector<TestData> testData = {
        {
            "Normal Case 1 - missed adapter.",
            // Input sequences for adapter finding.
            PacBio::PancakeTestsConfig::Data_Dir + "/adapter-finder/true-positives/test.01.missed_adapter.good_antidiag.m64011_190228_190319-4064130-208993_227621.fasta",
            // KmerSize, WindowSize
            seedParams1.KmerSize, seedParams1.MinimizerWindow,
            // Median insert length.
            7885,
            // Expected results.
            {
                {9247, 9247, PacBio::Pancake::Adapter::AdapterDetectionType::PALINDROMIC_FINE},
            },
        },
        {
            "Normal Case 2 - missed adapter.",
            // Input sequences for adapter finding.
            PacBio::PancakeTestsConfig::Data_Dir + "/adapter-finder/true-positives/test.02.missed_adapter.good_antitiag.m64011_190228_190319-4064055-202294_220031.fasta",
            // KmerSize, WindowSize
            seedParams1.KmerSize, seedParams1.MinimizerWindow,
            // Median insert length.
            9133,
            // Expected results.
            {
                {8666, 8666, PacBio::Pancake::Adapter::AdapterDetectionType::PALINDROMIC_FINE},
            },
        },
        {
            "Normal Case 3 - missed adapter.",
            // Input sequences for adapter finding.
            PacBio::PancakeTestsConfig::Data_Dir + "/adapter-finder/true-positives/test.03.missed_adapter.good_antitiag.m64011_190228_190319-4064055-131140_149079.fasta",
            // KmerSize, WindowSize
            seedParams1.KmerSize, seedParams1.MinimizerWindow,
            // Median insert length.
            9133,
            // Expected results.
            {
                {8875, 8875, PacBio::Pancake::Adapter::AdapterDetectionType::PALINDROMIC_FINE},
            },
        },
        {
            "Normal Case 4 - missed adapter.",
            // Input sequences for adapter finding.
            PacBio::PancakeTestsConfig::Data_Dir + "/adapter-finder/true-positives/test.04.missed_adapter.weak_antidiag.m64011_190228_190319-4063950-161516_179997.fasta",
            // KmerSize, WindowSize
            seedParams1.KmerSize, seedParams1.MinimizerWindow,
            // Median insert length.
            9246,
            // Expected results.
            {
                {9199, 9199, PacBio::Pancake::Adapter::AdapterDetectionType::PALINDROMIC_FINE},
            },
        },
        {
            "Normal Case 5 - multiple missed adapters.",
            // Input sequences for adapter finding.
            PacBio::PancakeTestsConfig::Data_Dir + "/adapter-finder/true-positives/test.05.missed_adapters_multiple.m64011_190228_190319-4065675-230077_272103.fasta",
            // KmerSize, WindowSize
            seedParams1.KmerSize, seedParams1.MinimizerWindow,
            // Median insert length.
            10679,
            // Expected results.
            {
                {10471, 10471, PacBio::Pancake::Adapter::AdapterDetectionType::PALINDROMIC_FINE},
                {21119, 21119, PacBio::Pancake::Adapter::AdapterDetectionType::PALINDROMIC_FINE},
                {31565, 31565, PacBio::Pancake::Adapter::AdapterDetectionType::PALINDROMIC_FINE},
            },
        },
        {
            "Normal Case 6 - multiple missed adapters.",
            // Input sequences for adapter finding.
            PacBio::PancakeTestsConfig::Data_Dir + "/adapter-finder/true-positives/test.06.missed_adapters_multiple.weak_antidiag.m64011_190228_190319-4259851-61111_108768.fasta",
            // KmerSize, WindowSize
            seedParams1.KmerSize, seedParams1.MinimizerWindow,
            // Median insert length.
            11704,
            // Expected results.
            {
            },
        },

        // Dificult cases.
        {
            "Difficult Case 1.",
            // Input sequences for adapter finding.
            PacBio::PancakeTestsConfig::Data_Dir + "/adapter-finder/difficult-positives/test.01.difficult_antidiag.m64011_190228_190319-4064355-152070_168593.fasta",
            // KmerSize, WindowSize
            seedParams1.KmerSize, seedParams1.MinimizerWindow,
            // Median insert length.
            9138,
            // Expected results.
            {
                {8885, 8885, PacBio::Pancake::Adapter::AdapterDetectionType::PALINDROMIC_FINE},
            },
        },

        // Weird cases.
        {
            "Weird Case 1. Twice longer than median, but no antidiag nor low-complexity nor repeat.",
            // Input sequences for adapter finding.
            PacBio::PancakeTestsConfig::Data_Dir + "/adapter-finder/weird-cases/test.01.twice_longer_than_median_but_no_antidiag_nor_lowcmplx_burst_nor_repeats.m64011_190228_190319-4064632-228963_254288.fasta",
            // KmerSize, WindowSize
            seedParams1.KmerSize, seedParams1.MinimizerWindow,
            // Median insert length.
            10143,
            // Expected results.
            {
            },
        },

        // Synthetic cases.
        {
            "Synthetic Case 1. Single adapter, perfect palindrome of 5000bp. Forward adapter sequence.",
            // Input sequences for adapter finding.
            PacBio::PancakeTestsConfig::Data_Dir + "/adapter-finder/synthetic/test.01.5000bp-fwd_adapter-5000bp_palindrome.fasta",
            // KmerSize, WindowSize
            seedParams1.KmerSize, seedParams1.MinimizerWindow,
            // Median insert length.
            5000,
            // Expected results.
            {
                {5009, 5036, PacBio::Pancake::Adapter::AdapterDetectionType::ALIGNED},
            },
        },
        {
            "Synthetic Case 2. Single adapter, perfect palindrome of 5000bp. Reverse adapter sequence.",
            // Input sequences for adapter finding.
            PacBio::PancakeTestsConfig::Data_Dir + "/adapter-finder/synthetic/test.02.5000bp-rev_adapter-5000bp_palindrome.fasta",
            // KmerSize, WindowSize
            seedParams1.KmerSize, seedParams1.MinimizerWindow,
            // Median insert length.
            5000,
            // Expected results.
            {
                {5009, 5036, PacBio::Pancake::Adapter::AdapterDetectionType::ALIGNED},
            },
        },
        {
            "Synthetic Case 3. Single adapter, forward adapter sequence. Palindromic copy is a partial length. 5000bp+adapter+3000bp.",
            // Input sequences for adapter finding.
            PacBio::PancakeTestsConfig::Data_Dir + "/adapter-finder/synthetic/test.03.5000bp-fwd_adapter-3000bp_palindrome.fasta",
            // KmerSize, WindowSize
            seedParams1.KmerSize, seedParams1.MinimizerWindow,
            // Median insert length.x
            5000,
            // Expected results.
            {
                {5009, 5036, PacBio::Pancake::Adapter::AdapterDetectionType::ALIGNED},
            },
        },
        {
            "Synthetic Case 4. Two adapters: 5000bp seq + fwd adapter + rev 5000bp seq + rev adapter + fwd 3000bp seq.",
            // Input sequences for adapter finding.
            PacBio::PancakeTestsConfig::Data_Dir + "/adapter-finder/synthetic/test.04.5000bp-fwd_adapter-5000bp-rev_adapter-3000bp.fasta",
            // KmerSize, WindowSize
            seedParams1.KmerSize, seedParams1.MinimizerWindow,
            // Median insert length.
            5000,
            // Expected results.
            {
                {5009, 5036, PacBio::Pancake::Adapter::AdapterDetectionType::ALIGNED},
                {10054, 10081, PacBio::Pancake::Adapter::AdapterDetectionType::ALIGNED},
            },
        },
        {
            "Synthetic Case 5. fwd 5000bp + no adapter + rev 5000bp.",
            // Input sequences for adapter finding.
            PacBio::PancakeTestsConfig::Data_Dir + "/adapter-finder/synthetic/test.05.5000bp-no_adapter-5000bp.fasta",
            // KmerSize, WindowSize
            seedParams1.KmerSize, seedParams1.MinimizerWindow,
            // Median insert length.
            5000,
            // Expected results.
            {
                {5000, 5000, PacBio::Pancake::Adapter::AdapterDetectionType::PALINDROMIC_FINE},
            },
        },
        {
            "Synthetic Case 6. fwd 5000bp + no adapter + rev 3000bp.",
            // Input sequences for adapter finding.
            PacBio::PancakeTestsConfig::Data_Dir + "/adapter-finder/synthetic/test.06.5000bp-no_adapter-3000bp.fasta",
            // KmerSize, WindowSize
            seedParams1.KmerSize, seedParams1.MinimizerWindow,
            // Median insert length.
            5000,
            // Expected results.
            {
                {5000, 5000, PacBio::Pancake::Adapter::AdapterDetectionType::PALINDROMIC_FINE},
            },
        },

        // Celegans dataset (very old, RSII).
        {
            "Celegans true positive case 1.",
            // Input sequences for adapter finding.
            PacBio::PancakeTestsConfig::Data_Dir + "/adapter-finder/true-positives/celegans/01.m140928_184123_42139_c100719602550000001823155305141590_s1_p0-966-0_22057.fasta",
            // KmerSize, WindowSize
            seedParams1.KmerSize, seedParams1.MinimizerWindow,
            // Median inssert length. Value <= 0 means unknown.
            16000,
            // Expected results.
            {
                {16107, 16107, PacBio::Pancake::Adapter::AdapterDetectionType::PALINDROMIC_FINE},
            },
        },
        {
            "Celegans true positive case 2.",
            // Input sequences for adapter finding.
            PacBio::PancakeTestsConfig::Data_Dir + "/adapter-finder/true-positives/celegans/02.m140928_184123_42139_c100719602550000001823155305141590_s1_p0-1132-0_25753.fasta",
            // KmerSize, WindowSize
            seedParams1.KmerSize, seedParams1.MinimizerWindow,
            // Median insert length. Value <= 0 means unknown.
            25000,
            // Expected results.
            {
                {24949, 24949, PacBio::Pancake::Adapter::AdapterDetectionType::PALINDROMIC_FINE},
            },
        },
        {
            "Celegans true positive case 3.",
            // Input sequences for adapter finding.
            PacBio::PancakeTestsConfig::Data_Dir + "/adapter-finder/true-positives/celegans/03.m140928_184123_42139_c100719602550000001823155305141590_s1_p0-1459-0_5452.fasta",
            // KmerSize, WindowSize
            seedParams1.KmerSize, seedParams1.MinimizerWindow,
            // Median insert length. Value <= 0 means unknown.
            3000,
            // Expected results.
            {
                {2969, 2969, PacBio::Pancake::Adapter::AdapterDetectionType::PALINDROMIC_FINE},
            },
        },
        {
            "Celegans true positive case 4.",
            // Input sequences for adapter finding.
            PacBio::PancakeTestsConfig::Data_Dir + "/adapter-finder/true-positives/celegans/04.m140928_184123_42139_c100719602550000001823155305141590_s1_p0-1607-0_22586.fasta",
            // KmerSize, WindowSize
            seedParams1.KmerSize, seedParams1.MinimizerWindow,
            // Median insert length. Value <= 0 means unknown.
            17700,
            // Expected results.
            {
                {17706, 17706, PacBio::Pancake::Adapter::AdapterDetectionType::PALINDROMIC_FINE},
            },
        },
        {
            "Celegans true positive case 5.",
            // Input sequences for adapter finding.
            PacBio::PancakeTestsConfig::Data_Dir + "/adapter-finder/true-positives/celegans/05.m140928_184123_42139_c100719602550000001823155305141590_s1_p0-2194-0_32522.fasta",
            // KmerSize, WindowSize
            seedParams1.KmerSize, seedParams1.MinimizerWindow,
            // Median insert length. Value <= 0 means unknown.
            23000,
            // Expected results.
            {
                {23002, 23002, PacBio::Pancake::Adapter::AdapterDetectionType::PALINDROMIC_FINE},
            },
        },
        /////
        {
            "False positive test case 1.",
            // Input sequences for adapter finding.
            PacBio::PancakeTestsConfig::Data_Dir + "/adapter-finder/false-positives/test.01.false_positive_subread.fasta",
            // KmerSize, WindowSize
            seedParams1.KmerSize, seedParams1.MinimizerWindow,
            // Median insert length. Value <= 0 means unknown.
            23000,
            // Expected results.
            {
                {5832, 5832, PacBio::Pancake::Adapter::AdapterDetectionType::PALINDROMIC_FINE},
            },
        },
    };
    // clang-format on

    // Reverse: TTTTCCTCCTCCTCCGTTGTTGTTGTT , Fwd. loop+stem: ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT , Rev loop+stem: ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT
    const std::string adapterSeq("AACAACAACAACGGAGGAGGAGGAAAA");
    // const std::string adapterSeq("TTTTCCTCCTCCTCCGTTGTTGTTGTT");

    // PacBio::Logging::Logger logger{std::cerr, PacBio::Logging::LogLevel::DEBUG};
    // PacBio::Logging::Logger::Current(&logger);

    const int32_t minNumSeeds = 10;
    const double minIdentity = 0.60;
    const double minDiagonalSpanFraction = 0.60;
    const int32_t maxAllowedFlank = -1;  // 300;
    const int32_t maxAdapterEditDist = 5;

    for (const auto& data : testData) {
        // Load the sequences from files, and take only the first one.
        const std::vector<PacBio::BAM::FastaSequence> allSeqs =
            PacBio::PancakeTests::HelperLoadFasta(data.InReadsFn);
        const PacBio::Pancake::FastaSequenceCached seq(
            allSeqs[0].Name(), allSeqs[0].Bases().c_str(), allSeqs[0].Bases().size(), 0);

        SCOPED_TRACE(data.TestName);

        // std::cerr << "TestName = " << data.TestName << "\n";
        // std::cerr.flush();

        // Run unit under test.
        const std::vector<PacBio::Pancake::Adapter::AdapterPosition> result =
            PacBio::Pancake::Adapter::FindPotentialAdapterLocations(
                seq, adapterSeq, data.KmerSize, data.WindowSize, minNumSeeds, minIdentity,
                minDiagonalSpanFraction, maxAllowedFlank, maxAdapterEditDist);

        std::cerr.flush();

        // std::cerr << "Median length: " << data.MedianLen << "\n";
        // std::cerr << "Test results:\n";
        // for (size_t i = 0; i < result.size(); ++i) {
        //     std::cerr << "[" << i << "] " << result[i] << "\n";
        // }
        // std::cerr << "\n-----------------------\n";

        EXPECT_EQ(data.ExpectedAdapters, result);
    }
}
