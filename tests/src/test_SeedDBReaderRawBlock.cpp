// Authors: Ivan Sovic

#include <PancakeTestData.h>
#include <gtest/gtest.h>
#include <pacbio/pancake/SeedDBReader.h>
#include <pacbio/pancake/SeedDBReaderRawBlock.h>
#include <pacbio/pancake/SequenceSeeds.h>
#include <sstream>

TEST(SeedDBReaderRawBlock, GetSeedDBContiguousParts_NormalSingleBlock)
{
    /*
     * Fetch the byte span of a single small block of 1 sequence.
     * It doesn't span more than 1 file or sequences.
    */

    const std::string inSeedDB =
        R"(V	0.1.0
P	k=30,w=80,hpc=0,hpc_len=10,rc=1
F	0	test-1a.seeddb.0.seeds	5	35936
S	0	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852	0	0	2976	5852	186
S	1	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983	0	2976	7664	11983	479
S	2	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292	0	10640	12384	24292	774
S	3	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105	0	23024	2608	5105	163
S	4	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001	0	25632	10304	19001	644
B	0	0	1	2976
B	1	1	2	7664
B	2	2	3	12384
B	3	3	4	2608
B	4	4	5	10304
)";
    const int32_t blockId = 0;

    // Expected.
    // Tuple: (file_id, file_offset_start, file_offset_end)
    const std::vector<PacBio::Pancake::ContiguousFilePart> expected = {
        PacBio::Pancake::ContiguousFilePart{0, 0, 2976, {0}}};

    // Load the SeedDB.
    std::istringstream is(inSeedDB);
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(is, "filename.seeddb");

    // Run.
    const auto results = PacBio::Pancake::GetSeedDBContiguousParts(seedDBCache, blockId);

    // Evaluate.
    EXPECT_EQ(results, expected);
}

TEST(SeedDBReaderRawBlock, GetSeedDBContiguousParts_NormalMultipleFiles)
{
    /*
     * Fetch the byte span of a block of 4 sequences, where each sequence is
     * stored in a different file.
    */

    const std::string inSeedDB =
        R"(V	0.1.0
P	k=30,w=80,hpc=0,hpc_len=10,rc=1
F	0	test-1b.seeddb.0.seeds	1	2976
F	1	test-1b.seeddb.1.seeds	1	7664
F	2	test-1b.seeddb.2.seeds	1	12384
F	3	test-1b.seeddb.3.seeds	1	2608
F	4	test-1b.seeddb.4.seeds	1	10304
S	0	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852	0	0	2976	5852	186
S	1	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983	1	0	7664	11983	479
S	2	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292	2	0	12384	24292	774
S	3	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105	3	0	2608	5105	163
S	4	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001	4	0	10304	19001	644
B	0	0	1	2976
B	1	1	5	32960
)";
    const int32_t blockId = 1;

    // Expected.
    // Tuple: (file_id, file_offset_start, file_offset_end)
    const std::vector<PacBio::Pancake::ContiguousFilePart> expected = {
        {1, 0, 7664, {1}},
        {2, 0, 12384, {2}},
        {3, 0, 2608, {3}},
        {4, 0, 10304, {4}},
    };

    // Load the SeedDB.
    std::istringstream is(inSeedDB);
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(is, "filename.seeddb");

    // Run.
    const auto results = PacBio::Pancake::GetSeedDBContiguousParts(seedDBCache, blockId);

    // Evaluate.
    EXPECT_EQ(results, expected);
}

TEST(SeedDBReaderRawBlock, GetSeedDBContiguousParts_NormalTwoBlocksWithGap)
{
    /*
     * Fetch the byte span of a block of 4 sequences, where the first two and last
     * two are separated by a filtered sequence. There should be two parts because
     * of that.
    */

    const std::string inSeedDB =
        R"(V	0.1.0
P	k=30,w=80,hpc=0,hpc_len=10,rc=1
F	0	test-1a.seeddb.0.seeds	4	23552
S	0	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852	0	0	2976	5852	186
S	1	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983	0	2976	7664	11983	479
S	2	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105	0	23024	2608	5105	163
S	3	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001	0	25632	10304	19001	644
B	0	0	4	23552
)";
    const int32_t blockId = 0;

    // Expected.
    // Tuple: (file_id, file_offset_start, file_offset_end)
    const std::vector<PacBio::Pancake::ContiguousFilePart> expected = {{0, 0, 10640, {0, 1}},
                                                                       {0, 23024, 35936, {2, 3}}};

    // Load the SeedDB.
    std::istringstream is(inSeedDB);
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(is, "filename.seeddb");

    // Run.
    const auto results = PacBio::Pancake::GetSeedDBContiguousParts(seedDBCache, blockId);

    // Evaluate.
    EXPECT_EQ(results, expected);
}

TEST(SeedDBReaderRawBlock, GetSeedDBContiguousParts_BlockOutOfBoundsThrows)
{
    /*
     * Block is out of bounds, it should throw.
    */

    const std::string inSeedDB =
        R"(V	0.1.0
P	k=30,w=80,hpc=0,hpc_len=10,rc=1
F	0	test-1a.seeddb.0.seeds	1	2976
S	0	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852	0	0	2976	5852	186
B	0	0	1	2976
)";
    const int32_t blockId = 123;

    // Load the SeedDB.
    std::istringstream is(inSeedDB);
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(is, "filename.seeddb");

    // Run and evaluate.
    EXPECT_THROW(
        { auto results = PacBio::Pancake::GetSeedDBContiguousParts(seedDBCache, blockId); },
        std::runtime_error);
}

TEST(SeedDBReaderRawBlock, GetSeedDBContiguousParts_MalformedBlockThrows)
{
    /*
     * Block is malformed, referencing sequences which do not exist in the index.
     * It should throw.
    */

    const std::string inSeedDB =
        R"(V	0.1.0
P	k=30,w=80,hpc=0,hpc_len=10,rc=1
F	0	test-1a.seeddb.0.seeds	1	2976
S	0	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852	0	0	2976	5852	186
B	0	1	5	2976
)";
    const int32_t blockId = 0;

    // Load the SeedDB.
    std::istringstream is(inSeedDB);
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(is, "filename.seeddb");

    // Run and evaluate.
    EXPECT_THROW(
        { auto results = PacBio::Pancake::GetSeedDBContiguousParts(seedDBCache, blockId); },
        std::runtime_error);
}

TEST(SeedDBReaderRawBlock, GetSeedDBContiguousParts_OverlappingBytesThrows)
{
    /*
     * One 'S' line is malformed (S3), it begins and overlaps S2. This should throw
     * in the GetSeedDBContiguousParts becauss seeds should be distinct byte blocks for each sequence.
    */

    const std::string inSeedDB =
        R"(V	0.1.0
P	k=30,w=80,hpc=0,hpc_len=10,rc=1
F	0	test-1a.seeddb.0.seeds	4	23552
S	0	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852	0	0	2976	5852	186
S	1	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983	0	2976	7664	11983	479
S	2	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105	0	23024	2608	5105	163
S	3	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001	0	20000	10304	19001	644
B	0	0	4	23552
)";
    const int32_t blockId = 0;

    // Load the SeedDB.
    std::istringstream is(inSeedDB);
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(is, "filename.seeddb");

    // Run and evaluate.
    EXPECT_THROW(
        { auto results = PacBio::Pancake::GetSeedDBContiguousParts(seedDBCache, blockId); },
        std::runtime_error);
}

TEST(SeedDBReaderRawBlock, GetSeedDBContiguousParts_OutOfOrder)
{
    /*
     * This is a valid case, where the order of sequences permuted in the SeedDB. This
     * should not throw.
    */

    const std::string inSeedDB =
        R"(V	0.1.0
P	k=30,w=80,hpc=0,hpc_len=10,rc=1
F	0	test-1a.seeddb.0.seeds	4	23552
S	0	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852	0	0	2976	5852	186
S	1	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983	0	2976	7664	11983	479
S	2	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001	0	25632	10304	19001	644
S	3	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105	0	23024	2608	5105	163
B	0	0	4	23552
)";
    const int32_t blockId = 0;

    // Expected.
    // Tuple: (file_id, file_offset_start, file_offset_end)
    const std::vector<PacBio::Pancake::ContiguousFilePart> expected = {
        {0, 0, 10640, {0, 1}}, {0, 25632, 35936, {2}}, {0, 23024, 25632, {3}}};

    // Load the SeedDB.
    std::istringstream is(inSeedDB);
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(is, "filename.seeddb");

    // Run.
    const auto results = PacBio::Pancake::GetSeedDBContiguousParts(seedDBCache, blockId);

    // Evaluate.
    EXPECT_EQ(results, expected);
}

TEST(SeedDBReaderRawBlock, GetBlock_NormalWithMultipleFiles)
{
    /*
     * This will load a block with the new function and return a set of seeds.
     * To verify, each of the sequences in the block will be loaded with the
     * previously tested reader and accumulated, then compared.
     *
     * The difference between this SeedDB specified here in the string and the one
     * stored in "test-data/seeddb-writer/test-1b.seeddb" is only in the block groupings.
    */
    const std::string inSeedDB =
        R"(V	0.1.0
P	k=20,w=80,s=0,hpc=0,hpc_len=10,rc=1
F	0	test-1b.seeddb.0.seeds	1	2464
F	1	test-1b.seeddb.1.seeds	1	4752
F	2	test-1b.seeddb.2.seeds	1	9488
F	3	test-1b.seeddb.3.seeds	1	1904
F	4	test-1b.seeddb.4.seeds	1	7456
S	0	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852	0	0	2464	5852	154
S	1	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983	1	0	4752	11983	297
S	2	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292	2	0	9488	24292	593
S	3	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105	3	0	1904	5105	119
S	4	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001	4	0	7456	19001	466
B	0	0	1	2976
B	1	1	5	23600
)";
    const int32_t blockId = 1;

    // Load the SeedDB index cache.
    std::istringstream is(inSeedDB);
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(
            is, PacBio::PancakeTestsConfig::Data_Dir + "/seeddb-writer/test-1b.seeddb");

    // Expected results. Construct these by loading and flattening the block using
    // the SeedDBReader.
    std::vector<PacBio::Pancake::SeedDB::SeedRaw> expected;
    {
        // Create a SeedDB reader.
        PacBio::Pancake::SeedDBReader readerExp(seedDBCache);
        std::vector<PacBio::Pancake::SequenceSeeds> records;
        readerExp.GetBlock(records, blockId);
        for (const auto& record : records) {
            expected.insert(expected.end(), record.Seeds().begin(), record.Seeds().end());
        }
    }

    // Run.
    PacBio::Pancake::SeedDBReaderRawBlock reader(seedDBCache);
    const auto results = reader.GetBlock(blockId);

    // Evaluate.
    EXPECT_EQ(results, expected);
}

TEST(SeedDBReaderRawBlock, GetBlock_NormalWithSingleFileAndAGap)
{
    /*
     * This will load a block with the new function and return a set of seeds.
     * To verify, each of the sequences in the block will be loaded with the
     * previously tested reader and accumulated, then compared.
     *
     * Here, the block has a contiguity gap, and will be a concatenation
     * of two parts (sequences 0 and 1 as one part, and 2 and 3 as the other part).
    */

    const std::string inSeedDB =
        R"(V	0.1.0
P	k=20,w=80,s=0,hpc=0,hpc_len=10,rc=1
F	0	test-1a.seeddb.0.seeds	5	26064
S	0	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852	0	0	2464	5852	154
S	1	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983	0	2464	4752	11983	297
S	2	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105	0	16704	1904	5105	119
S	3	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001	0	18608	7456	19001	466
B	0	0	4	16576
)";
    const int32_t blockId = 0;

    // Load the SeedDB index cache.
    std::istringstream is(inSeedDB);
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(
            is, PacBio::PancakeTestsConfig::Data_Dir + "/seeddb-writer/test-1.seeddb");

    // Expected results. Construct these by loading and flattening the block using
    // the SeedDBReader.
    std::vector<PacBio::Pancake::SeedDB::SeedRaw> expected;
    {
        // Create a SeedDB reader.
        PacBio::Pancake::SeedDBReader readerExp(seedDBCache);
        std::vector<PacBio::Pancake::SequenceSeeds> records;
        readerExp.GetBlock(records, blockId);
        for (const auto& record : records) {
            expected.insert(expected.end(), record.Seeds().begin(), record.Seeds().end());
        }
    }

    // Run.
    PacBio::Pancake::SeedDBReaderRawBlock reader(seedDBCache);
    const auto results = reader.GetBlock(blockId);

    // Evaluate.
    EXPECT_EQ(results, expected);
}
