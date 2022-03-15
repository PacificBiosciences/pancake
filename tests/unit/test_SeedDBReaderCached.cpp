// Authors: Ivan Sovic

#include <PancakeTestData.h>
#include <gtest/gtest.h>
#include <pancake/SeedDBReaderCached.hpp>
#include <pancake/SequenceSeeds.hpp>
#include <sstream>

TEST(SeedDBReaderCached, IterateThroughFirstBlock)
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
    PacBio::Pancake::SeedDBReaderCached reader(seedDBCache, blockId);

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

TEST(SeedDBReaderCached, IterateThroughSecondBlock)
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
    PacBio::Pancake::SeedDBReaderCached reader(seedDBCache, blockId);

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

TEST(SeedDBReaderCached, GetSeedsForSequenceRandomAccess)
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
    PacBio::Pancake::SeedDBReaderCached reader(seedDBCache, blockId);

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

TEST(SeedDBReaderCached, GetSeedsForSequenceWhichDontExist)
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
    PacBio::Pancake::SeedDBReaderCached reader(seedDBCache, blockId);

    // Evaluate.
    EXPECT_THROW({ reader.GetSeedsForSequence(12345); }, std::runtime_error);
    EXPECT_THROW({ reader.GetSeedsForSequence("some_nonexistent_name"); }, std::runtime_error);
}

TEST(SeedDBReaderCached, ConstructFromNonexistentBlock)
{
    /*
     * Try to construct the SeedDBReaderCached from a block which does not
     * exist in the index.
    */

    // Input values.
    const std::string inSeedDB =
        PacBio::PancakeTestsConfig::Data_Dir + "/seeddb-writer/test-1c.seeddb";

    // Load the SeedDB.
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(inSeedDB);

    // Evaluate.
    EXPECT_THROW({ PacBio::Pancake::SeedDBReaderCached(seedDBCache, -1); }, std::runtime_error);
    EXPECT_THROW({ PacBio::Pancake::SeedDBReaderCached(seedDBCache, 123); }, std::runtime_error);
}
