// Authors: Ivan Sovic

#include <PancakeTestData.h>
#include <gtest/gtest.h>
#include <pancake/SeedDBReader.hpp>
#include <pancake/SeedDBReaderCachedBlock.hpp>
#include <pancake/SequenceSeeds.hpp>
#include <pancake/util/CommonTypes.hpp>
#include <sstream>

TEST(SeedDBReaderCachedBlock, IterateThroughFirstBlock)
{
    /*
     * Normal test case, test loading the first block into cache.
     * This tests iteration over the loaded records.
    */

    // Input values.
    const std::string inSeedDB =
        PacBio::PancakeTestsConfig::Data_Dir + "/seeddb-writer/test-1c.seeddb";
    const int32_t blockId = 0;

    // Expected values.
    std::vector<std::pair<int64_t, std::string>> expected = {
        {0, "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852"},
        {1, "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983"},
        {2, "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292"},
        {3, "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105"}};

    // Load the SeedDB.
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(inSeedDB);

    // Create a SeedDB reader.
    PacBio::Pancake::SeedDBReaderCachedBlock reader(seedDBCache, {blockId});

    // Iterate through all records in the cache and accumulate just the ID and the header.
    // We want to check whether all sequences are fetched properly.
    // Validity of the seeds is tested elsewhere.
    std::vector<std::pair<int64_t, std::string>> results;
    for (const auto& record : reader.records()) {
        results.emplace_back(std::make_pair(record.Id(), record.Name()));
    }

    // Evaluate.
    EXPECT_EQ(results, expected);
}

TEST(SeedDBReaderCachedBlock, IterateThroughSecondBlock)
{
    /*
     * Normal test case, test loading the second block into cache.
     * This tests iteration over the loaded records.
    */

    // Input values.
    const std::string inSeedDB =
        PacBio::PancakeTestsConfig::Data_Dir + "/seeddb-writer/test-1c.seeddb";
    const int32_t blockId = 1;

    // Expected values.
    std::vector<std::pair<int64_t, std::string>> expected = {
        {4, "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001"}};

    // Load the SeedDB.
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(inSeedDB);

    // Create a SeedDB reader.
    PacBio::Pancake::SeedDBReaderCachedBlock reader(seedDBCache, {blockId});

    // Iterate through all records in the cache and accumulate just the ID and the header.
    // We want to check whether all sequences are fetched properly.
    // Validity of the seeds is tested elsewhere.
    std::vector<std::pair<int64_t, std::string>> results;
    for (const auto& record : reader.records()) {
        results.emplace_back(std::make_pair(record.Id(), record.Name()));
    }

    // Evaluate.
    EXPECT_EQ(results, expected);
}

TEST(SeedDBReaderCachedBlock, GetSeedsForSequenceRandomAccess)
{
    /*
     * Test random access.
    */

    // Input values.
    const std::string inSeedDB =
        PacBio::PancakeTestsConfig::Data_Dir + "/seeddb-writer/test-1c.seeddb";
    const int32_t blockId = 0;

    // Load the SeedDB.
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(inSeedDB);

    // Create a SeedDB reader.
    PacBio::Pancake::SeedDBReaderCachedBlock reader(seedDBCache, {blockId});

    // Access each one of the loaded records by it's ID. The records are iterated one
    // by one in reverse order.
    // The fetched ID and Name should match the current
    // values in the "expected" vector.
    for (auto it = reader.records().rbegin(); it != reader.records().rend(); ++it) {
        const auto& record = reader.GetSeedsForSequence(it->Id());
        EXPECT_EQ(record.Id(), it->Id());
        EXPECT_EQ(record.Name(), it->Name());
    }

    // The same as above, but test fetching by Name instead of ID.
    for (auto it = reader.records().rbegin(); it != reader.records().rend(); ++it) {
        const auto& record = reader.GetSeedsForSequence(it->Name());
        EXPECT_EQ(record.Id(), it->Id());
        EXPECT_EQ(record.Name(), it->Name());
    }
}

TEST(SeedDBReaderCachedBlock, GetSeedsForSequenceWhichDontExist)
{
    /*
     * Try to access sequences which do not exist in the block.
    */

    // Input values.
    const std::string inSeedDB =
        PacBio::PancakeTestsConfig::Data_Dir + "/seeddb-writer/test-1c.seeddb";
    const int32_t blockId = 0;

    // Load the SeedDB.
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(inSeedDB);

    // Create a SeedDB reader.
    PacBio::Pancake::SeedDBReaderCachedBlock reader(seedDBCache, {blockId});

    // Evaluate.
    EXPECT_THROW({ reader.GetSeedsForSequence(12345); }, std::runtime_error);
    EXPECT_THROW({ reader.GetSeedsForSequence("some_nonexistent_name"); }, std::runtime_error);
}

TEST(SeedDBReaderCachedBlock, ConstructFromNonexistentBlock)
{
    /*
     * Try to construct the SeedDBReaderCachedBlock from a block which does not
     * exist in the index.
    */

    // Input values.
    const std::string inSeedDB =
        PacBio::PancakeTestsConfig::Data_Dir + "/seeddb-writer/test-1c.seeddb";

    // Load the SeedDB.
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(inSeedDB);

    // Evaluate.
    EXPECT_THROW({ PacBio::Pancake::SeedDBReaderCachedBlock(seedDBCache, {-1}); },
                 std::runtime_error);
    EXPECT_THROW({ PacBio::Pancake::SeedDBReaderCachedBlock(seedDBCache, {123}); },
                 std::runtime_error);
}

TEST(SeedDBReaderCachedBlock, MultipleInputBlocks)
{
    /*
     * Test loading of multiple blocks at once.
     * Evaluate using a previously tested SeedDBReader as the truth.
    */

    const std::string inSeqDB =
        PacBio::PancakeTestsConfig::Data_Dir + "/seeddb-writer/test-1a.seeddb";
    const std::vector<int32_t> inBlocks = {1, 2, 3};

    // Load the SeedDB.
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> indexCache =
        PacBio::Pancake::LoadSeedDBIndexCache(inSeqDB);

    // Collect all expected seeds lines for the specified input blocks
    // using an orthogonal reader. These are treated as the truth data.
    std::vector<PacBio::Pancake::SequenceSeeds> expected;
    for (const auto& blockId : inBlocks) {
        PacBio::Pancake::SeedDBReader reader(indexCache);
        std::vector<PacBio::Pancake::SequenceSeeds> blockRecords;
        reader.GetBlock(blockRecords, blockId);
        expected.insert(expected.end(), blockRecords.begin(), blockRecords.end());
    }
    std::sort(expected.begin(), expected.end(),
              [](const auto& a, const auto& b) { return a.Id() < b.Id(); });

    // Create a unit under test.
    // Read the sequences for the specified blocks, and convert all
    // FastaSequenceCached to FastaSequenceId for easier comparison.
    PacBio::Pancake::SeedDBReaderCachedBlock readerTest(indexCache, {inBlocks});
    std::vector<PacBio::Pancake::SequenceSeeds> results;
    for (const auto& record : readerTest.records()) {
        std::vector<PacBio::Pancake::Int128t> seeds(record.Seeds(), record.Seeds() + record.Size());
        results.emplace_back(PacBio::Pancake::SequenceSeeds(record.Name(), seeds, record.Id()));
    }
    std::sort(results.begin(), results.end(),
              [](const auto& a, const auto& b) { return a.Id() < b.Id(); });

    // Evaluate the current block.
    EXPECT_EQ(expected, results);
}

TEST(SeedDBReaderCachedBlock, BatchCompareWithSeedDBReader)
{
    /*
     * Test loading of multiple blocks at once.
     * Evaluate using a previously tested SeedDBReader as the truth.
    */

    const std::vector<std::string> inDBs = {
        PacBio::PancakeTestsConfig::Data_Dir + "/seeddb-writer/test-1a.seeddb",
        PacBio::PancakeTestsConfig::Data_Dir + "/seeddb-writer/test-1b.seeddb",
        PacBio::PancakeTestsConfig::Data_Dir + "/seeddb-writer/test-1c.seeddb",
        PacBio::PancakeTestsConfig::Data_Dir + "/seeddb-writer/test-1e-permuted-file-offset.seeddb",
    };

    for (const auto& inSeqDB : inDBs) {
        // Load the SeedDB.
        std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> indexCache =
            PacBio::Pancake::LoadSeedDBIndexCache(inSeqDB);

        const int32_t numBlocks = indexCache->blockLines.size();

        for (int32_t blockId = 0; blockId < numBlocks; ++blockId) {
            SCOPED_TRACE(inSeqDB + ", blockId = " + std::to_string(blockId));

            // Collect all expected seeds lines for the specified input blocks
            // using an orthogonal reader. These are treated as the truth data.
            std::vector<PacBio::Pancake::SequenceSeeds> expected;
            PacBio::Pancake::SeedDBReader readerTruth(indexCache);
            std::vector<PacBio::Pancake::SequenceSeeds> blockRecords;
            readerTruth.GetBlock(blockRecords, blockId);
            expected.insert(expected.end(), blockRecords.begin(), blockRecords.end());
            std::sort(expected.begin(), expected.end(),
                      [](const auto& a, const auto& b) { return a.Id() < b.Id(); });

            // Create a unit under test.
            // Read the sequences for the specified blocks, and convert all
            // vvectors to SequenceSeeds objects for easier comparison.
            PacBio::Pancake::SeedDBReaderCachedBlock readerTest(indexCache, {blockId});
            std::vector<PacBio::Pancake::SequenceSeeds> results;
            for (const auto& record : readerTest.records()) {
                std::vector<PacBio::Pancake::Int128t> seeds(record.Seeds(),
                                                            record.Seeds() + record.Size());
                results.emplace_back(
                    PacBio::Pancake::SequenceSeeds(record.Name(), seeds, record.Id()));
            }
            std::sort(results.begin(), results.end(),
                      [](const auto& a, const auto& b) { return a.Id() < b.Id(); });

            // Evaluate the current block.
            EXPECT_EQ(expected, results);
        }
    }
}
