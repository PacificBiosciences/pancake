// Authors: Ivan Sovic

#include <PancakeTestData.h>
#include <gtest/gtest.h>
#include <pancake/SeqDBReaderCached.hpp>
#include <sstream>

TEST(SeqDBReaderCached, IterateThroughSecondBlock)
{
    /*
     * Normal test case, test loading a block into cache.
     * This tests iteration over the loaded records.
    */

    // Input values.
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir +
                                "/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb";
    const int32_t blockId = 2;

    // Expected values.
    std::vector<std::pair<int64_t, std::string>> expected = {
        {2, "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292"}};

    // Load the SeedDB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

    // Create a SeedDB reader.
    PacBio::Pancake::SeqDBReaderCached reader(seqDBCache, blockId);

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

TEST(SeqDBReaderCached, IterateThroughALargerBlock)
{
    /*
     * Normal test case, test loading a larger block into cache.
     * This tests iteration over the loaded records.
    */

    // Input values.
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-6.seqdb";
    const int32_t blockId = 0;

    // Expected values.
    std::vector<std::pair<int64_t, std::string>> expected = {
        {0, "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852"},
        {1, "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983"},
        {2, "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292"},
        {3, "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105"},
        {4, "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001"},
    };

    // Load the SeedDB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

    // Create a SeedDB reader.
    PacBio::Pancake::SeqDBReaderCached reader(seqDBCache, blockId);

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

TEST(SeqDBReaderCached, GetSequenceRandomAccess)
{
    /*
     * Test random access.
    */

    // Input values.
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-6.seqdb";
    const int32_t blockId = 0;

    // Load the SeedDB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

    // Create a SeedDB reader.
    PacBio::Pancake::SeqDBReaderCached reader(seqDBCache, blockId);

    // Access each one of the loaded records by it's ID. The records are iterated one
    // by one in reverse order.
    // The fetched ID and Name should match the current
    // values in the "expected" vector.
    for (auto it = reader.records().rbegin(); it != reader.records().rend(); ++it) {
        const auto& record = reader.GetSequence(it->Id());
        EXPECT_EQ(record.Id(), it->Id());
        EXPECT_EQ(record.Name(), it->Name());
    }

    // The same as above, but test fetching by Name instead of ID.
    for (auto it = reader.records().rbegin(); it != reader.records().rend(); ++it) {
        const auto& record = reader.GetSequence(it->Name());
        EXPECT_EQ(record.Id(), it->Id());
        EXPECT_EQ(record.Name(), it->Name());
    }
}

TEST(SeqDBReaderCached, GetSequenceWhichDoesntExist)
{
    /*
     * Try to access sequences which do not exist in the block.
    */

    // Input values.
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir +
                                "/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb";
    const int32_t blockId = 0;

    // Load the SeedDB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

    // Create a SeedDB reader.
    PacBio::Pancake::SeqDBReaderCached reader(seqDBCache, blockId);

    // Evaluate.
    EXPECT_THROW({ reader.GetSequence(12345); }, std::runtime_error);
    EXPECT_THROW({ reader.GetSequence("some_nonexistent_name"); }, std::runtime_error);
}

TEST(SeqDBReaderCached, ConstructFromNonexistentBlock)
{
    /*
     * Try to construct the SeqDBReaderCached from a block which does not
     * exist in the index.
    */

    // Input values.
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir +
                                "/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb";

    // Load the SeedDB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

    // Evaluate.
    EXPECT_THROW({ PacBio::Pancake::SeqDBReaderCached(seqDBCache, -1); }, std::runtime_error);
    EXPECT_THROW({ PacBio::Pancake::SeqDBReaderCached(seqDBCache, 123); }, std::runtime_error);
}
