#include <gtest/gtest.h>
#include <pacbio/pancake/Minimizers.h>
#include <pacbio/pancake/Seed.h>
#include <pacbio/util/CommonTypes.h>
// #include <iostream>

using namespace PacBio::Pancake;

void HelperTestGenerateMinimizers(const std::string& seq, int32_t seqId, int32_t k, int32_t w,
                                  int32_t space, bool useHPC, bool useRC, int32_t expectedRv,
                                  const std::vector<PacBio::Pancake::Int128t>& expectedSeeds,
                                  bool expectedThrow)
{
    // Run the unit under test.
    const uint8_t* seqData = reinterpret_cast<const uint8_t*>(seq.data());
    int32_t seqLen = seq.size();
    std::vector<PacBio::Pancake::Int128t> seeds;

    if (expectedThrow) {
        EXPECT_THROW(
            {
                PacBio::Pancake::SeedDB::GenerateMinimizers(seeds, seqData, seqLen, 0, seqId, k, w,
                                                            space, useRC, useHPC);
            },
            std::runtime_error);

    } else {
        const int rv = PacBio::Pancake::SeedDB::GenerateMinimizers(seeds, seqData, seqLen, 0, seqId,
                                                                   k, w, space, useRC, useHPC);

        // std::cerr << "Results:\n";
        // for (const auto& val : seeds) {
        //     auto s = PacBio::Pancake::SeedDB::Seed(val);
        //     // std::cerr << s.Verbose() << "\n";
        //     std::cerr << "PacBio::Pancake::SeedDB::Seed::Encode(" << s.key << ", " << s.span << ", "
        //               << s.seqID << ", " << s.pos << ", " << s.seqRev << "),\n";
        // }

        // std::cerr << "Expected:\n";
        // for (const auto& val : expectedSeeds) {
        //     auto s = PacBio::Pancake::SeedDB::Seed(val);
        //     // std::cerr << s.Verbose() << "\n";
        //     std::cerr << "PacBio::Pancake::SeedDB::Seed::Encode(" << s.key << ", " << s.span << ", "
        //               << s.seqID << ", " << s.pos << ", " << s.seqRev << "),\n";
        // }

        // Compare the results.
        EXPECT_EQ(expectedRv, rv);
        EXPECT_EQ(expectedSeeds, seeds);
    }
}

TEST(GenerateMinimizers, SmallTest1)
{
    // Inputs.
    const std::string seq = "AAAAAAAAAA";
    const int32_t seqId = 123;
    const int32_t k = 5;
    const int32_t w = 1;
    const int32_t space = 0;
    const bool useHPC = false;
    const bool useRC = true;

    // Helper function.
    const uint64_t mask = SeedDB::ComputeKmerMask(k);
    auto Hash = [&](uint64_t val) { return SeedDB::InvertibleHash(val, mask); };

    // Expected results.
    const int32_t expectedRv = 0;
    const bool expectedThrow = false;
    const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0), k, seqId, 0, false),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0), k, seqId, 1, false),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0), k, seqId, 2, false),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0), k, seqId, 3, false),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0), k, seqId, 4, false),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0), k, seqId, 5, false),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, useRC, expectedRv, expectedSeeds,
                                 expectedThrow);
}

TEST(GenerateMinimizers, SmallTest2)
{
    // Inputs.
    const std::string seq = "CAAAAAAAAA";
    const int32_t seqId = 123;
    const int32_t k = 5;
    const int32_t w = 1;
    const int32_t space = 0;
    const bool useHPC = false;
    const bool useRC = true;

    // Helper function.
    const uint64_t mask = SeedDB::ComputeKmerMask(k);
    auto Hash = [&](uint64_t val) { return SeedDB::InvertibleHash(val, mask); };

    // Expected results.
    const int32_t expectedRv = 0;
    const bool expectedThrow = false;
    const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(256), k, seqId, 0, false),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0), k, seqId, 1, false),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0), k, seqId, 2, false),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0), k, seqId, 3, false),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0), k, seqId, 4, false),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0), k, seqId, 5, false),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, useRC, expectedRv, expectedSeeds,
                                 expectedThrow);
}

TEST(GenerateMinimizers, SmallTest3)
{
    // Inputs.
    const std::string seq = "AGCTTTTCATTCTGACTGCA";
    const int32_t seqId = 123;
    const int32_t k = 5;
    const int32_t w = 1;
    const int32_t space = 0;
    const bool useHPC = false;
    const bool useRC = true;

    // Helper function.
    const uint64_t mask = SeedDB::ComputeKmerMask(k);
    auto Hash = [&](uint64_t val) { return SeedDB::InvertibleHash(val, mask); };

    // Expected results.
    const int32_t expectedRv = 0;
    const bool expectedThrow = false;
    const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(39), k, seqId, 0, true),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(9), k, seqId, 1, true),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(2), k, seqId, 2, true),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(512), k, seqId, 3, true),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(896), k, seqId, 4, true),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(224), k, seqId, 5, true),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(56), k, seqId, 6, true),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(317), k, seqId, 7, false),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(131), k, seqId, 8, true),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(288), k, seqId, 9, true),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(840), k, seqId, 10, true),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(481), k, seqId, 11, false),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(180), k, seqId, 12, true),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(301), k, seqId, 13, true),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(121), k, seqId, 14, false),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(484), k, seqId, 15, false),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, useRC, expectedRv, expectedSeeds,
                                 expectedThrow);
}

TEST(GenerateMinimizers, SmallTest4)
{
    // Inputs.
    const std::string seq = "TCTGACTGCA";
    const int32_t seqId = 123;
    const int32_t k = 5;
    const int32_t w = 1;
    const int32_t space = 0;
    const bool useHPC = false;
    const bool useRC = true;

    // Helper function.
    const uint64_t mask = SeedDB::ComputeKmerMask(k);
    auto Hash = [&](uint64_t val) { return SeedDB::InvertibleHash(val, mask); };

    // Expected results.
    const int32_t expectedRv = 0;
    const bool expectedThrow = false;
    const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(840), k, seqId, 0, true),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(481), k, seqId, 1, false),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(180), k, seqId, 2, true),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(301), k, seqId, 3, true),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(121), k, seqId, 4, false),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(484), k, seqId, 5, false),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, useRC, expectedRv, expectedSeeds,
                                 expectedThrow);
}

TEST(GenerateMinimizers, SmallTest5)
{
    // Inputs.
    const std::string seq =
        "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC";
    const int32_t seqId = 123;
    const int32_t k = 15;
    const int32_t w = 5;
    const int32_t space = 0;
    const bool useHPC = false;
    const bool useRC = true;

    // Expected results.
    const int32_t expectedRv = 0;
    const bool expectedThrow = false;
    /*
     * Without the invertible hash function, these would be the generated seeds.
    */
    // const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(167726968), k, seqId, 0, false),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(189282306), k, seqId, 2, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(315756032), k, seqId, 3, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(259553542), k, seqId, 8, false),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(364792648), k, seqId, 10, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(126904899), k, seqId, 14, false),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(80713150), k, seqId, 18, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(27856109), k, seqId, 19, false),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(76041465), k, seqId, 24, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(54228923), k, seqId, 26, false),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(55644705), k, seqId, 31, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(502196464), k, seqId, 33, false),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(518950656), k, seqId, 35, false),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(536855377), k, seqId, 39, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(503315509), k, seqId, 41, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(125828877), k, seqId, 42, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(74973168), k, seqId, 44, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(35767), k, seqId, 46, false),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(143070), k, seqId, 47, false),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(572280), k, seqId, 48, false),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(2289123), k, seqId, 49, false),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(9156492), k, seqId, 50, false),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(36625970), k, seqId, 51, false),
    // };

    /*
     * WITH the invertible hash function, these would be the generated seeds.
    */
    const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(132144311, 15, 123, 4, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(33567509, 15, 123, 5, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(18410160, 15, 123, 6, 0),
        PacBio::Pancake::SeedDB::Seed::Encode(457467302, 15, 123, 10, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(60328715, 15, 123, 12, 0),
        PacBio::Pancake::SeedDB::Seed::Encode(35576091, 15, 123, 15, 0),
        PacBio::Pancake::SeedDB::Seed::Encode(65050698, 15, 123, 20, 0),
        PacBio::Pancake::SeedDB::Seed::Encode(107140825, 15, 123, 23, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(330443040, 15, 123, 24, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(668958223, 15, 123, 27, 0),
        PacBio::Pancake::SeedDB::Seed::Encode(152777239, 15, 123, 30, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(235365058, 15, 123, 34, 0),
        PacBio::Pancake::SeedDB::Seed::Encode(39491807, 15, 123, 37, 0),
        PacBio::Pancake::SeedDB::Seed::Encode(156927623, 15, 123, 38, 0),
        PacBio::Pancake::SeedDB::Seed::Encode(511189515, 15, 123, 41, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(394814865, 15, 123, 46, 0),
        PacBio::Pancake::SeedDB::Seed::Encode(222756012, 15, 123, 47, 0),
        PacBio::Pancake::SeedDB::Seed::Encode(423331484, 15, 123, 52, 0),
        PacBio::Pancake::SeedDB::Seed::Encode(70557206, 15, 123, 54, 0),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, useRC, expectedRv, expectedSeeds,
                                 expectedThrow);
}

TEST(GenerateMinimizers, SmallTest6)
{
    // Inputs.
    const std::string seq = "AGCTTTTCATTCTGACTGCA";
    const int32_t seqId = 123;
    const int32_t k = 5;
    const int32_t w = 4;
    const int32_t space = 0;
    const bool useHPC = false;
    const bool useRC = true;

    // Expected results.
    const int32_t expectedRv = 0;
    const bool expectedThrow = false;
    /*
     * Without the invertible hash function, these would be the generated seeds.
    */
    // const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(2), k, seqId, 2, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(56), k, seqId, 6, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(131), k, seqId, 8, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(180), k, seqId, 12, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(Hash(121), k, seqId, 14, false),
    // };
    /*
     * WITH the invertible hash function, these would be the generated seeds.
    */
    const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(67, 5, 123, 3, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(195, 5, 123, 4, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(227, 5, 123, 5, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(235, 5, 123, 6, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(419, 5, 123, 9, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(351, 5, 123, 12, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(239, 5, 123, 15, 0),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, useRC, expectedRv, expectedSeeds,
                                 expectedThrow);
}

TEST(GenerateMinimizers, SmallTest7WithNBases)
{
    // Inputs.
    const std::string seq = "AGCTTTTCATTCTGACTGCANNNACTNNNNNAGCTTTTCATTCTGACTGCA";
    const int32_t seqId = 123;
    const int32_t k = 5;
    const int32_t w = 4;
    const int32_t space = 0;
    const bool useHPC = false;
    const bool useRC = true;

    // Expected results.
    const int32_t expectedRv = 0;
    const bool expectedThrow = false;
    /*
     * Without the invertible hash function, these would be the generated seeds.
    */
    // const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
    //     PacBio::Pancake::SeedDB::Seed::Encode(2, seqId, 2, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(56, seqId, 6, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(131, seqId, 8, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(180, seqId, 12, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(121, seqId, 14, false),
    //     PacBio::Pancake::SeedDB::Seed::Encode(2, seqId, 33, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(56, seqId, 37, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(131, seqId, 39, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(180, seqId, 43, true),
    //     PacBio::Pancake::SeedDB::Seed::Encode(121, seqId, 45, false),
    // };
    /*
     * WITH the invertible hash function, these would be the generated seeds.
    */
    const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(67, 5, 123, 3, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(195, 5, 123, 4, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(227, 5, 123, 5, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(235, 5, 123, 6, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(419, 5, 123, 9, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(351, 5, 123, 12, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(239, 5, 123, 15, 0),
        PacBio::Pancake::SeedDB::Seed::Encode(67, 5, 123, 34, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(195, 5, 123, 35, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(227, 5, 123, 36, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(235, 5, 123, 37, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(419, 5, 123, 40, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(351, 5, 123, 43, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(239, 5, 123, 46, 0),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, useRC, expectedRv, expectedSeeds,
                                 expectedThrow);
}

TEST(GenerateMinimizers, SeedSize32BasePairs_PolyA)
{
    // Inputs.
    const std::string seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    const int32_t seqId = 0;
    const int32_t k = 32;
    const int32_t w = 1;
    const int32_t space = 0;
    const bool useHPC = false;
    const bool useRC = true;

    // Expected results.
    const int32_t expectedRv = 0;
    const bool expectedThrow = true;
    const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {};

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, useRC, expectedRv, expectedSeeds,
                                 expectedThrow);
}

TEST(GenerateMinimizers, SeedSize28BasePairs_PolyA)
{
    // Inputs.
    const std::string seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    const int32_t seqId = 0;
    const int32_t k = 28;
    const int32_t w = 1;
    const int32_t space = 0;
    const bool useHPC = false;
    const bool useRC = true;

    // Helper function.
    const uint64_t mask = SeedDB::ComputeKmerMask(k);
    auto Hash = [&](uint64_t val) { return SeedDB::InvertibleHash(val, mask); };

    // Expected results.
    const int32_t expectedRv = 0;
    const bool expectedThrow = false;
    const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0), k, seqId, 0, false),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, useRC, expectedRv, expectedSeeds,
                                 expectedThrow);
}

TEST(GenerateMinimizers, SeedSize28BasePairs_PolyT_WithRC)
{
    // Inputs.
    const std::string seq = "TTTTTTTTTTTTTTTTTTTTTTTTTTTT";
    const int32_t seqId = 0;
    const int32_t k = 28;
    const int32_t w = 1;
    const int32_t space = 0;
    const bool useHPC = false;
    const bool useRC = true;

    // Helper function.
    const uint64_t mask = SeedDB::ComputeKmerMask(k);
    auto Hash = [&](uint64_t val) { return SeedDB::InvertibleHash(val, mask); };

    // Expected results.
    const int32_t expectedRv = 0;
    const bool expectedThrow = false;
    const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0), k, seqId, 0, true),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, useRC, expectedRv, expectedSeeds,
                                 expectedThrow);
}

TEST(GenerateMinimizers, SeedSize28BasePairs_PolyT_OnlyFWD)
{
    // Inputs.
    const std::string seq = "TTTTTTTTTTTTTTTTTTTTTTTTTTTT";
    const int32_t seqId = 0;
    const int32_t k = 28;
    const int32_t w = 1;
    const int32_t space = 0;
    const bool useHPC = false;
    const bool useRC = false;

    // Helper function.
    const uint64_t mask = SeedDB::ComputeKmerMask(k);
    auto Hash = [&](uint64_t val) { return SeedDB::InvertibleHash(val, mask); };

    // Expected results.
    const int32_t expectedRv = 0;
    const bool expectedThrow = false;
    const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0xFFFFFFFFFFFFFFFF), k, seqId, 0, false),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, useRC, expectedRv, expectedSeeds,
                                 expectedThrow);
}

TEST(GenerateMinimizers, HPC1)
{
    // Inputs.
    const std::string seq = "ACGTTTTG";
    const int32_t seqId = 123;
    const int32_t k = 5;
    const int32_t w = 1;
    const int32_t space = 0;
    const bool useHPC = true;
    const bool useRC = true;

    // Helper function.
    const uint64_t mask = SeedDB::ComputeKmerMask(k);
    auto Hash = [&](uint64_t val) { return SeedDB::InvertibleHash(val, mask); };

    // Expected results.
    const int32_t expectedRv = 0;
    const bool expectedThrow = false;
    const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(110), 8, seqId, 0, false),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, useRC, expectedRv, expectedSeeds,
                                 expectedThrow);
}

TEST(GenerateMinimizers, HPC2)
{
    // This test is a bit more comprehensive. It does not explicitly specify inputs and outputs, but performs
    // two runs with different settings and inputs, and expects the same results.
    // First run - use a sequence with homopolymer runs, and use the homopolymer compression options to generate the seeds.
    // Second run - use a sequence with manually compressed homopolymers, and compute minimizers without the compression option.
    // Both runs should result in identical sets of seeds, except for the positions (because one sequence is shorter than the other.

    // Inputs.
    const std::string seqWithHP =
        "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC";
    // Manually removed homopolymers from the sequence (manual HP compression).
    const std::string seqWithoutHP = "AGCTCATCTGACTGCACGCATATGTCTCTGTGTGATAGAGTGTCTGATAGCAGC";
    const int32_t seqId = 123;
    const int32_t k = 5;
    const int32_t w = 1;
    const int32_t space = 0;
    const bool useRC = true;

    // Run the unit under test WITH homopolymer compression for the sequence with homopolymers.
    std::vector<PacBio::Pancake::Int128t> seedsWithHP;
    int32_t rvWithHP = 0;
    {
        const uint8_t* seqData = reinterpret_cast<const uint8_t*>(seqWithHP.data());
        int32_t seqLen = seqWithHP.size();
        rvWithHP = PacBio::Pancake::SeedDB::GenerateMinimizers(seedsWithHP, seqData, seqLen, 0,
                                                               seqId, k, w, space, useRC, true);

        // Reset the span because this will differ between the no-hpc and with-hpc runs.
        for (auto& val : seedsWithHP) {
            auto s = PacBio::Pancake::SeedDB::Seed(val);
            s.span = 0;
            val = s.To128t();
        }
    }

    // Run the unit under test WITHOUT homopolymer compression on the sequence with no homopolymers.
    std::vector<PacBio::Pancake::Int128t> seedsNoHP;
    int32_t rvNoHP = 0;
    {
        const uint8_t* seqData = reinterpret_cast<const uint8_t*>(seqWithoutHP.data());
        int32_t seqLen = seqWithoutHP.size();
        rvNoHP = PacBio::Pancake::SeedDB::GenerateMinimizers(seedsNoHP, seqData, seqLen, 0, seqId,
                                                             k, w, space, useRC, false);

        // Reset the span because this will differ between the no-hpc and with-hpc runs.
        for (auto& val : seedsNoHP) {
            auto s = PacBio::Pancake::SeedDB::Seed(val);
            s.span = 0;
            val = s.To128t();
        }
    }

    // Check the return values.
    EXPECT_EQ(rvWithHP, 0);
    EXPECT_EQ(rvNoHP, 0);

    // There should be the same amount of seeds in both vectors.
    EXPECT_EQ(seedsWithHP.size(), seedsNoHP.size());

    // Need to check each seed manually, because positions are different.
    // Both seed vectors should have the same seeds, just at different positions.
    for (size_t i = 0; i < seedsWithHP.size(); ++i) {
        auto seedWithHP = PacBio::Pancake::SeedDB::Seed(seedsWithHP[i]);
        auto seedNoHP = PacBio::Pancake::SeedDB::Seed(seedsNoHP[i]);
        EXPECT_EQ(seedWithHP.key, seedNoHP.key);
    }
}

TEST(GenerateMinimizers, SpacedSeed_Space1_31bp_JustOneSeed)
{
    /*
     * Here we test the spaced seed construction, with 1 skipped base
     * in between every 2 inclusive bases (i.e. space = 1).
     * This should skip every 'T' base and leave only 'A' bases in the test.
    */
    // Inputs.
    const std::string seq = "ATATATATATATATATATATATATATATATA";
    const int32_t seqId = 0;
    const int32_t k = 16;
    const int32_t w = 1;
    const int32_t space = 1;
    const bool useHPC = false;
    const bool useRC = false;

    // Helper function.
    const uint64_t mask = SeedDB::ComputeKmerMask(k);
    auto Hash = [&](uint64_t val) { return SeedDB::InvertibleHash(val, mask); };

    // Expected results.
    const int32_t expectedRv = 0;
    const bool expectedThrow = false;
    const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0), 31, seqId, 0, false),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, useRC, expectedRv, expectedSeeds,
                                 expectedThrow);
}

TEST(GenerateMinimizers, SpacedSeed_Space1_32bp_TwoSeeds)
{
    /*
     * Here we test the spaced seed construction, with 1 skipped base
     * in between every 2 inclusive bases (i.e. space = 1).
     * This should skip every 'T' base and leave only 'A' bases in the test.
    */
    // Inputs.
    const std::string seq = "TATATATATATATATATATATATATATATATA";
    const int32_t seqId = 0;
    const int32_t k = 16;
    const int32_t w = 1;
    const int32_t space = 1;
    const bool useHPC = false;
    const bool useRC = false;

    // Helper function.
    const uint64_t mask = SeedDB::ComputeKmerMask(k);
    auto Hash = [&](uint64_t val) { return SeedDB::InvertibleHash(val, mask); };

    // Expected results.
    const int32_t expectedRv = 0;
    const bool expectedThrow = false;
    const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0x0FFFFFFFF), 31, seqId, 0,
                                              false),  // 16 bases of 'T's.
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0), 31, seqId, 1, false),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, useRC, expectedRv, expectedSeeds,
                                 expectedThrow);
}

TEST(GenerateMinimizers, SpacedSeed_Space1_33bp_ThreeSeeds)
{
    /*
     * Here we test the spaced seed construction, with 1 skipped base
     * in between every 2 inclusive bases (i.e. space = 1).
     * This should skip every 'T' base and leave only 'A' bases in the test.
    */
    // Inputs.
    const std::string seq = "TATATATATATATATATATATATATATATATAG";
    const int32_t seqId = 0;
    const int32_t k = 16;
    const int32_t w = 1;
    const int32_t space = 1;
    const bool useHPC = false;
    const bool useRC = false;

    // Helper function.
    const uint64_t mask = SeedDB::ComputeKmerMask(k);
    auto Hash = [&](uint64_t val) { return SeedDB::InvertibleHash(val, mask); };

    // Expected results.
    const int32_t expectedRv = 0;
    const bool expectedThrow = false;
    const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0x0FFFFFFFF), 31, seqId, 0,
                                              false),  // 16 bases of 'T's.
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0), 31, seqId, 1, false),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0x0FFFFFFFE), 31, seqId, 2, false),
    };

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, useRC, expectedRv, expectedSeeds,
                                 expectedThrow);
}

TEST(GenerateMinimizers, FromStrings)
{
    // Inputs.
    const int32_t k = 5;
    const int32_t w = 4;
    const int32_t space = 0;
    const bool useHPC = false;
    const bool useRC = true;
    const std::vector<std::string> seqs = {
        "AAAAAAAAAA", "AGCTTTTCATTCTGACTGCANNNACTNNNNNAGCTTTTCATTCTGACTGCA",
    };

    // Helper function.
    const uint64_t mask = SeedDB::ComputeKmerMask(k);
    auto Hash = [&](uint64_t val) { return SeedDB::InvertibleHash(val, mask); };

    // Expected results.
    const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0), k, 0, 3, false),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0), k, 0, 4, false),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(0), k, 0, 5, false),

        PacBio::Pancake::SeedDB::Seed::Encode(67, 5, 1, 3, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(195, 5, 1, 4, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(227, 5, 1, 5, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(235, 5, 1, 6, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(419, 5, 1, 9, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(351, 5, 1, 12, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(239, 5, 1, 15, 0),
        PacBio::Pancake::SeedDB::Seed::Encode(67, 5, 1, 34, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(195, 5, 1, 35, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(227, 5, 1, 36, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(235, 5, 1, 37, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(419, 5, 1, 40, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(351, 5, 1, 43, 1),
        PacBio::Pancake::SeedDB::Seed::Encode(239, 5, 1, 46, 0),
    };
    std::vector<int32_t> expectedSequenceLengths = {10, 51};

    // Run unit under test.
    std::vector<PacBio::Pancake::Int128t> results;
    std::vector<int32_t> sequenceLengths;
    PacBio::Pancake::SeedDB::GenerateMinimizers(results, sequenceLengths, seqs, k, w, space, useRC,
                                                useHPC);

    // Evaluate.
    EXPECT_EQ(expectedSeeds, results);
    EXPECT_EQ(expectedSequenceLengths, sequenceLengths);
}

TEST(GenerateMinimizers, HPCompression_SmallExampleThatFitsInSeedSpan)
{
    // clang-format off
    // Inputs.
    const std::string seq = "CAAAAACTCTC";
    const int32_t seqId = 123;
    const int32_t k = 5;
    const int32_t w = 1;
    const int32_t space = 0;
    const bool useHPC = true;
    const bool useRC = false;

    // Helper function.
    const uint64_t mask = SeedDB::ComputeKmerMask(k);
    auto Hash = [&](uint64_t val) { return SeedDB::InvertibleHash(val, mask); };

    // Expected results.
    /*
        Pos:        0 1 2 3 4 5 6 7 8 9 10
        Seq:        C A A A A A C T C T C
                   .               . . . .
        Windows:   .               . . . .
                   |C A - - - - C T C| . .      Window 0, pos = 0, span = 9. (C, 5*A, C, T, C, T)
                   |_________________| . .
                     |A - - - - C T C T| .      Window 1, pos = 1, span = 9. (5*A, C, T, C, T, C)
                     |_________________| .
                               |C T C T C|      Window 2, pos = 6, span = 5.
                               |_________|
    */
    const uint64_t CACTC = 0b0001'0001'1101;
    const uint64_t ACTCT = 0b0000'0111'0111;
    const uint64_t CTCTC = 0b0001'1101'1101;

    const int32_t expectedRv = 0;
    const bool expectedThrow = false;
    const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {

        PacBio::Pancake::SeedDB::Seed::Encode(Hash(CACTC), 9, seqId, 0, 0),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(ACTCT), 9, seqId, 1, 0),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(CTCTC), 5, seqId, 6, 0),
    };
    // clang-format on

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, useRC, expectedRv, expectedSeeds,
                                 expectedThrow);
}

TEST(GenerateMinimizers, HPCompression_VeryLongHPC_150bp)
{
    // clang-format off
    // Inputs.
    // There are 150*A bases in the seq. There need to be more bases before/after the homopolymer because we need to be
    // able to fill a 15-mer (with HP compression, all 150 bases are considered as one).
    const std::string seq = "CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACTCTCTCTCTCTCTCTCTCTCT";
    const int32_t seqId = 123;
    const int32_t k = 15;
    const int32_t w = 5;
    const int32_t space = 0;
    const bool useHPC = true;
    const bool useRC = false;

    // Helper function.
    const uint64_t mask = SeedDB::ComputeKmerMask(k);
    auto Hash = [&](uint64_t val) { return SeedDB::InvertibleHash(val, mask); };

    // Expected results.
    /*
                            C  A  C  T  C  T  C  T  C  T  C  T  C  T  C
        CACTCTCTCTCTCTC => 01 00 01 11 01 11 01 11 01 11 01 11 01 11 01  => 0b0001'0001'1101'1101'1101'1101'1101'1101

                            A  C  T  C  T  C  T  C  T  C  T  C  T  C  T
        ACTCTCTCTCTCTCT => 00 01 11 01 11 01 11 01 11 01 11 01 11 01 11  => 0b0000'0111'0111'0111'0111'0111'0111'0111

                            C  T  C  T  C  T  C  T  C  T  C  T  C  T  C
        CTCTCTCTCTCTCTC => 01 11 01 11 01 11 01 11 01 11 01 11 01 11 01  => 0b0001'1101'1101'1101'1101'1101'1101'1101
    */
    // const uint64_t ACTCTCTCTCTCTCT = 0b0000'0111'0111'0111'0111'0111'0111'0111;
    const uint64_t CACTCTCTCTCTCTC = 0b0001'0001'1101'1101'1101'1101'1101'1101;
    const uint64_t CTCTCTCTCTCTCTC = 0b0001'1101'1101'1101'1101'1101'1101'1101;
    const int32_t expectedRv = 0;
    const bool expectedThrow = false;
    const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(CACTCTCTCTCTCTC), 164, seqId, 0, 0),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(CTCTCTCTCTCTCTC), 15, seqId, 153, 0),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(CTCTCTCTCTCTCTC), 15, seqId, 155, 0),
        PacBio::Pancake::SeedDB::Seed::Encode(Hash(CTCTCTCTCTCTCTC), 15, seqId, 157, 0),
    };
    // clang-format on

    HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, useRC, expectedRv, expectedSeeds,
                                 expectedThrow);
}

TEST(GenerateMinimizers, HPCompression_VeryLongHPC_300bp)
{
    /*
        This tests a case where the HP span stretches out of the possible addressable area for a seed span (256 bp).
        The HP length here is 300bp.
        Any seed which spans > 255 bp should simply not be reported in the output.
        This includes any seed that covers the large HP in this example.
    */

    // clang-format off
    // Inputs.
    const std::string seq = "CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACTCTCTCTCTCTCTCTCTCTCT";

    /*
        CA...CTCTCTCTCTCTCTCTCTCTCT
        This is a list of all possible seeds. Any seed that contains the 'A' homopolymer is 299bp longer than shown here.
        These seeds cannot be encoded, so they will be skipped below.
        Win 0:          CACTCTCTCTCTCTC     pos =   0, span = 314
        Win 1:          ACTCTCTCTCTCTCT     pos =   1, span = 314
        Win 2:          CTCTCTCTCTCTCTC     pos = 301, span = 15
        Win 3:          TCTCTCTCTCTCTCT     pos = 302, span = 15
        Win 4:          CTCTCTCTCTCTCTC     pos = 303, span = 15
        Win 5:          TCTCTCTCTCTCTCT     pos = 304, span = 15
        Win 6:          CTCTCTCTCTCTCTC     pos = 305, span = 15
        Win 7:          TCTCTCTCTCTCTCT     pos = 306, span = 15
        Win 8:          CTCTCTCTCTCTCTC     pos = 307, span = 15
        Win 9:          TCTCTCTCTCTCTCT     pos = 308, span = 15
    */

    {   // TEST 1: All valid seeds (minimizer window of 1 bases).
        const int32_t seqId = 123;
        const int32_t k = 15;
        const int32_t w = 1;
        const int32_t space = 0;
        const bool useHPC = true;
        const bool useRC = false;

        // Helper function.
        const uint64_t mask = SeedDB::ComputeKmerMask(k);
        auto Hash = [&](uint64_t val) { return SeedDB::InvertibleHash(val, mask); };

        // Expected results.
        /*
                                C  T  C  T  C  T  C  T  C  T  C  T  C  T  C
            CTCTCTCTCTCTCTC => 01 11 01 11 01 11 01 11 01 11 01 11 01 11 01  => 0b0001'1101'1101'1101'1101'1101'1101'1101
                                T  C  T  C  T  C  T  C  T  C  T  C  T  C  T
            TCTCTCTCTCTCTCC => 11 01 11 01 11 01 11 01 11 01 11 01 11 01 11  => 0b0011'0111'0111'0111'0111'0111'0111'0111
        */
        const uint64_t CTCTCTCTCTCTCTC = 0b0001'1101'1101'1101'1101'1101'1101'1101;
        const uint64_t TCTCTCTCTCTCTCT = 0b0011'0111'0111'0111'0111'0111'0111'0111;
        const int32_t expectedRv = 0;
        const bool expectedThrow = false;
        const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
            // Any seed with a span out of max range (256bp) will simply not be reported.
            PacBio::Pancake::SeedDB::Seed::Encode(Hash(CTCTCTCTCTCTCTC), 15, 123, 301, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(Hash(TCTCTCTCTCTCTCT), 15, 123, 302, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(Hash(CTCTCTCTCTCTCTC), 15, 123, 303, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(Hash(TCTCTCTCTCTCTCT), 15, 123, 304, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(Hash(CTCTCTCTCTCTCTC), 15, 123, 305, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(Hash(TCTCTCTCTCTCTCT), 15, 123, 306, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(Hash(CTCTCTCTCTCTCTC), 15, 123, 307, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(Hash(TCTCTCTCTCTCTCT), 15, 123, 308, 0),
        };
        // clang-format on

        HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, useRC, expectedRv,
                                     expectedSeeds, expectedThrow);
    }

    {  // TEST 2: Minimizer window of 5 bases.
        const int32_t seqId = 123;
        const int32_t k = 15;
        const int32_t w = 5;
        const int32_t space = 0;
        const bool useHPC = true;
        const bool useRC = false;

        // Helper function.
        const uint64_t mask = SeedDB::ComputeKmerMask(k);
        auto Hash = [&](uint64_t val) { return SeedDB::InvertibleHash(val, mask); };

        // Expected results.
        /*
                                C  T  C  T  C  T  C  T  C  T  C  T  C  T  C
            CTCTCTCTCTCTCTC => 01 11 01 11 01 11 01 11 01 11 01 11 01 11 01  => 0b0001'1101'1101'1101'1101'1101'1101'1101
        */
        const uint64_t CTCTCTCTCTCTCTC = 0b0001'1101'1101'1101'1101'1101'1101'1101;
        const int32_t expectedRv = 0;
        const bool expectedThrow = false;
        const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
            // Any seed with a span out of max range (256bp) will simply not be reported.
            PacBio::Pancake::SeedDB::Seed::Encode(Hash(CTCTCTCTCTCTCTC), 15, 123, 303, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(Hash(CTCTCTCTCTCTCTC), 15, 123, 305, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(Hash(CTCTCTCTCTCTCTC), 15, 123, 307, 0),
        };
        // clang-format on

        HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, useRC, expectedRv,
                                     expectedSeeds, expectedThrow);
    }
}

TEST(GenerateMinimizers, HPCompression_TwoVeryLongHPC_300bp_plus_300bp)
{
    /*
     * This tests a case where the HP span stretches out of the possible addressable area for a seed span (256 bp).
     * There are two HPs of 300bp in length. Any seed which spans > 255 should simply not be reported in the output.
     * This includes any seed that covers the large HP.

     * This is a great test case because it exposed an issue where after a span of non-valid seeds
     * we end up with the minimizer window not correctly producing minimizers.
     * Concretely, these seeds were reported for w = 5:
            PacBio::Pancake::SeedDB::Seed::Encode(622765704, 15, 123, 303, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(622765704, 15, 123, 305, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(622765704, 15, 123, 307, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(622765704, 15, 123, 309, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(622765704, 15, 123, 624, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(996900086, 15, 123, 625, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(622765704, 15, 123, 626, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(996900086, 15, 123, 627, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(622765704, 15, 123, 628, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(622765704, 15, 123, 630, 0),
     * Seeds at positions {303, 305, 307, 309} are correct, but starting at 624 (the second CT stretch)
     * it looks like the minimizer generator did not remove the non-minimizer seeds.
    */

    // clang-format off
    // Inputs.
    const std::string seq =
                            // first copy
                            "CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                            "AAAAAAAAAAAAAAAAAAAAAACTCTCTCTCTCTCTCTCTCTCT"
                            // second copy
                            "CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                            "AAAAAAAAAAAAAAAAAAAAAACTCTCTCTCTCTCTCTCTCTCT";

    /*
        CA...CTCTCTCTCTCTCTCTCTCTCT

        --- (first copy begins here
        // These two seeds need to be skipped because they contain the first large HP.
        Win 0:          CACTCTCTCTCTCTC     pos = 0, span = 314
        Win 1:          ACTCTCTCTCTCTCT     pos = 1

        // Valid seeds for minimizer selection.
        Win 2:          CTCTCTCTCTCTCTC     pos = 301
        Win 3:          TCTCTCTCTCTCTCT     pos = 302
        Win 4:          CTCTCTCTCTCTCTC     pos = 303
        Win 5:          TCTCTCTCTCTCTCT     pos = 304
        Win 6:          CTCTCTCTCTCTCTC     pos = 305
        Win 7:          TCTCTCTCTCTCTCT     pos = 306
        Win 8:          CTCTCTCTCTCTCTC     pos = 307
        Win 9:          TCTCTCTCTCTCTCT     pos = 308

        --- (second copy begins here
        Win 10:         CTCTCTCTCTCTCTC     pos = 309

        // The following seeds should be skipped because they contain the second large HP.
        Win 11:         TCTCTCTCTCTCTCA     pos = 310
        Win 12:         CTCTCTCTCTCTCAC     pos = 311
        Win 13:         TCTCTCTCTCTCACT     pos = 312
        Win 14:         CTCTCTCTCTCACTC     pos = 313
        Win 15:         TCTCTCTCTCACTCT     pos = 314
        Win 16:         CTCTCTCTCACTCTC     pos = 315
        Win 17:         TCTCTCTCACTCTCT     pos = 316
        Win 18:         CTCTCTCACTCTCTC     pos = 317
        Win 19:         TCTCTCACTCTCTCT     pos = 318
        Win 20:         CTCTCACTCTCTCTC     pos = 319
        Win 21:         TCTCACTCTCTCTCT     pos = 320
        Win 22:         CTCACTCTCTCTCTC     pos = 321
        Win 23:         TCACTCTCTCTCTCT     pos = 322
        Win 24:         CACTCTCTCTCTCTC     pos = 323
        Win 25:         ACTCTCTCTCTCTCT     pos = 324

        // Valid seeds for minimizer selection.
        Win 26:         CTCTCTCTCTCTCTC     pos = 624
        Win 27:         TCTCTCTCTCTCTCT     pos = 625
        Win 28:         CTCTCTCTCTCTCTC     pos = 626
        Win 29:         TCTCTCTCTCTCTCT     pos = 627
        Win 30:         CTCTCTCTCTCTCTC     pos = 628
        Win 31:         TCTCTCTCTCTCTCT     pos = 629
        Win 32:         CTCTCTCTCTCTCTC     pos = 630
        Win 33:         TCTCTCTCTCTCTCT     pos = 631
    */

    {   // TEST 1:Minimizer window 1, so all valid seeds should be captured.
        const int32_t seqId = 123;
        const int32_t k = 15;
        const int32_t w = 1;
        const int32_t space = 0;
        const bool useHPC = true;
        const bool useRC = false;

        // Helper function.
        const uint64_t mask = SeedDB::ComputeKmerMask(k);
        auto Hash = [&](uint64_t val) { return SeedDB::InvertibleHash(val, mask); };

        // Expected results.
        /*
                                C  T  C  T  C  T  C  T  C  T  C  T  C  T  C
            CTCTCTCTCTCTCTC => 01 11 01 11 01 11 01 11 01 11 01 11 01 11 01  => 0b0001'1101'1101'1101'1101'1101'1101'1101
                                T  C  T  C  T  C  T  C  T  C  T  C  T  C  T
            TCTCTCTCTCTCTCC => 11 01 11 01 11 01 11 01 11 01 11 01 11 01 11  => 0b0011'0111'0111'0111'0111'0111'0111'0111
        */
        const uint64_t CTCTCTCTCTCTCTC = 0b0001'1101'1101'1101'1101'1101'1101'1101;
        const uint64_t TCTCTCTCTCTCTCT = 0b0011'0111'0111'0111'0111'0111'0111'0111;
        const int32_t expectedRv = 0;
        const bool expectedThrow = false;
        const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
            PacBio::Pancake::SeedDB::Seed::Encode(Hash(CTCTCTCTCTCTCTC), 15, seqId, 301, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(Hash(TCTCTCTCTCTCTCT), 15, seqId, 302, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(Hash(CTCTCTCTCTCTCTC), 15, seqId, 303, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(Hash(TCTCTCTCTCTCTCT), 15, seqId, 304, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(Hash(CTCTCTCTCTCTCTC), 15, seqId, 305, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(Hash(TCTCTCTCTCTCTCT), 15, seqId, 306, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(Hash(CTCTCTCTCTCTCTC), 15, seqId, 307, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(Hash(TCTCTCTCTCTCTCT), 15, seqId, 308, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(Hash(CTCTCTCTCTCTCTC), 15, seqId, 309, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(Hash(CTCTCTCTCTCTCTC), 15, seqId, 624, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(Hash(TCTCTCTCTCTCTCT), 15, seqId, 625, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(Hash(CTCTCTCTCTCTCTC), 15, seqId, 626, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(Hash(TCTCTCTCTCTCTCT), 15, seqId, 627, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(Hash(CTCTCTCTCTCTCTC), 15, seqId, 628, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(Hash(TCTCTCTCTCTCTCT), 15, seqId, 629, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(Hash(CTCTCTCTCTCTCTC), 15, seqId, 630, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(Hash(TCTCTCTCTCTCTCT), 15, seqId, 631, 0),
        };

        HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, useRC, expectedRv,
                                    expectedSeeds, expectedThrow);
    }

    {   // TEST 2: Minimizer window 5.
        const int32_t seqId = 123;
        const int32_t k = 15;
        const int32_t w = 5;
        const int32_t space = 0;
        const bool useHPC = true;
        const bool useRC = false;

        // Helper function.
        const uint64_t mask = SeedDB::ComputeKmerMask(k);
        auto Hash = [&](uint64_t val) { return SeedDB::InvertibleHash(val, mask); };

        // Expected results.
        /*
                                C  T  C  T  C  T  C  T  C  T  C  T  C  T  C
            CTCTCTCTCTCTCTC => 01 11 01 11 01 11 01 11 01 11 01 11 01 11 01  => 0b0001'1101'1101'1101'1101'1101'1101'1101
        */
        const uint64_t CTCTCTCTCTCTCTC = 0b0001'1101'1101'1101'1101'1101'1101'1101;
        const int32_t expectedRv = 0;
        const bool expectedThrow = false;
        const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
            PacBio::Pancake::SeedDB::Seed::Encode(Hash(CTCTCTCTCTCTCTC), 15, seqId, 303, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(Hash(CTCTCTCTCTCTCTC), 15, seqId, 305, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(Hash(CTCTCTCTCTCTCTC), 15, seqId, 307, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(Hash(CTCTCTCTCTCTCTC), 15, seqId, 309, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(Hash(CTCTCTCTCTCTCTC), 15, seqId, 624, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(Hash(CTCTCTCTCTCTCTC), 15, seqId, 626, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(Hash(CTCTCTCTCTCTCTC), 15, seqId, 628, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(Hash(CTCTCTCTCTCTCTC), 15, seqId, 630, 0),
        };

        HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, useRC, expectedRv,
                                    expectedSeeds, expectedThrow);
    }
    // clang-format on
}

TEST(GenerateMinimizers, HPCompression_TwoCloseHPsOf150bp)
{
    /*
     * Kmer window should span the first HP fully (there are 12 bases before the HP, HP is counted as 1 base, and after the HP
     * there are 2 more bases; in total that's 15 bases).
     * When it slides one base down, it will pick up on the second HP, and the total span should be longer than MAX_SPAN, and these
     * seeds should be ignored.
     * After that, it will slide down after the HP, and this should produce more valid seeds.
     *
     * Note: each HP is 150bp in length, so each one separately can fit into the Seed's span, but together they go out of range.
    */

    // clang-format off
    // Inputs.
    const std::string seq =
                            // first copy
                            "CTCTCTCTCTCTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                            ""
                            // second copy
                            "CTCTCTCTCTCTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                            "CT";

    {   // TEST 1:Minimizer window 1, so all valid seeds should be captured.
        const int32_t seqId = 123;
        const int32_t k = 15;
        const int32_t w = 1;
        const int32_t space = 0;
        const bool useHPC = true;
        const bool useRC = false;

        // Expected results.
        const int32_t expectedRv = 0;
        const bool expectedThrow = false;
        const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
            // PacBio::Pancake::SeedDB::Seed::Encode(Hash(CTCTCTCTCTCTCTC), 15, seqId, 301, 0),

            PacBio::Pancake::SeedDB::Seed::Encode(493075847, 164, 123, 0, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(757015818, 164, 123, 1, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(460313955, 164, 123, 2, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(630753600, 164, 123, 3, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(1028693452, 164, 123, 4, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(750559860, 164, 123, 5, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(432649153, 164, 123, 6, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(517395662, 164, 123, 7, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(576155398, 164, 123, 8, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(15530431, 164, 123, 9, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(710990505, 164, 123, 10, 0),
            // Position 11 is one base before the first HP, and the HP begins at position 12. Seeds which
            // begin at these positions also cover the second HP, so they are skipped.
            // First base after the first HP is at position 162. There are in total 14 regular bases and 1 HP left from this
            // position on, so only one seed is left.
            PacBio::Pancake::SeedDB::Seed::Encode(493075847, 164, 123, 162, 0),
        };

        HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, useRC, expectedRv,
                                    expectedSeeds, expectedThrow);
    }

    {   // TEST 2: Minimizer window 5.
        const int32_t seqId = 123;
        const int32_t k = 15;
        const int32_t w = 5;
        const int32_t space = 0;
        const bool useHPC = true;
        const bool useRC = false;

        // Expected results.
        const int32_t expectedRv = 0;
        const bool expectedThrow = false;
        const std::vector<PacBio::Pancake::Int128t> expectedSeeds = {
            PacBio::Pancake::SeedDB::Seed::Encode(460313955, 164, 123, 2, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(432649153, 164, 123, 6, 0),
            PacBio::Pancake::SeedDB::Seed::Encode(15530431, 164, 123, 9, 0),
        };

        HelperTestGenerateMinimizers(seq, seqId, k, w, space, useHPC, useRC, expectedRv,
                                    expectedSeeds, expectedThrow);
    }
    // clang-format on
}
