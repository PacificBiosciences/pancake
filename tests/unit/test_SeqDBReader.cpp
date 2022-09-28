// Authors: Ivan Sovic

#include <PancakeTestData.h>
#include <gtest/gtest.h>
#include <pbbam/FastaReader.h>
#include <pbbam/FastaSequence.h>
#include <pancake/FastaSequenceId.hpp>
#include <pancake/SeqDBReader.hpp>
#include <sstream>
#include "TestHelperUtils.hpp"

TEST(SeqDBReaderCompressed, GetNext)
{
    // Input values;
    const std::string inSeqFasta = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/in.fasta";
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir +
                                "/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb";

    // Load the expected sequences.
    // const std::string expected = PacBio::PancakeTests::HelperLoadFastaAsString(inSeqFasta);
    const auto expected = PacBio::PancakeTests::HelperLoadFasta(inSeqFasta);

    // Load the DB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

    // Create a reader.
    PacBio::Pancake::SeqDBReader reader(seqDBCache);

    // Collect all sequences as a single string for comparison.
    std::vector<PacBio::BAM::FastaSequence> result;
    PacBio::Pancake::FastaSequenceId record;
    while (reader.GetNext(record)) {
        result.emplace_back(record);
    }

    EXPECT_EQ(expected, result);
}

TEST(SeqDBReaderCompressed, GetSequenceByID1)
{
    /*
     * This tests accessing the sequences by their sequence ID.
    */

    // Input values;
    const std::string inSeqFasta = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/in.fasta";
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir +
                                "/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb";

    // Load the expected sequences.
    const auto expected = PacBio::PancakeTests::HelperLoadFasta(inSeqFasta);

    // Load the DB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

    // Create a reader.
    PacBio::Pancake::SeqDBReader reader(seqDBCache);

    // Reverse access all sequences by ID and check if they're the same.
    for (int32_t i = static_cast<int32_t>(expected.size()) - 1; i >= 0; --i) {
        const auto& expRecord = expected[i];
        PacBio::Pancake::FastaSequenceId record;
        reader.GetSequence(record, i);
        EXPECT_EQ(expRecord, record);
    }
}

TEST(SeqDBReaderCompressed, GetSequenceByID2)
{
    /*
     * Unlike the previous test (GetSequenceByID1) this loads a DB which as
     * permuted seqId values compared to their ordinal IDs. That means,
     * sequence with seqId = 3 appears as the first sequence in the file,
     * record with seqId = 0 is the second one, record with seqId = 4 is the third
     * one, and so on.
     *
     * We enforce that the sequences have IDs in the order of appearance in the file.
    */

    // Input values;
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-2.seqdb";

    EXPECT_THROW(
        {
            // Load the DB.
            std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
                PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);
        },
        std::runtime_error);
}

TEST(SeqDBReaderCompressed, GetSequenceByID3)
{
    /*
     * This tests accessing the same sequence multiple times in a row.
     * There was a bug in SeqDBReader - when a sequence was fread, the fileHandler.pos was
     * not being updated. This failed when the same sequence was being read multiple times in a row,
     * because fseek wouldn't run due to an if statement.
    */

    // Input values;
    const std::string inSeqFasta = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/in.fasta";
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir +
                                "/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb";

    // Load the expected sequences.
    const auto expected = PacBio::PancakeTests::HelperLoadFasta(inSeqFasta);

    // Load the DB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

    // Create a reader.
    PacBio::Pancake::SeqDBReader reader(seqDBCache);

    // Load the last sequence multiple times.
    std::vector<int32_t> loadOrder = {0,
                                      0,
                                      0,
                                      static_cast<int32_t>(expected.size()) - 1,
                                      static_cast<int32_t>(expected.size()) - 1,
                                      static_cast<int32_t>(expected.size()) - 1};

    // Arbitrary set of sequences.
    for (const auto& seqId : loadOrder) {
        const auto& expRecord = expected[seqId];
        PacBio::Pancake::FastaSequenceId record;
        reader.GetSequence(record, seqId);
        EXPECT_EQ(expRecord, record);
    }
}

TEST(SeqDBReaderCompressed, GetSequenceByName)
{
    // Input values;
    const std::string inSeqFasta = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/in.fasta";
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir +
                                "/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb";

    // Load the expected sequences.
    const auto expected = PacBio::PancakeTests::HelperLoadFasta(inSeqFasta);

    // Load the DB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

    // Create a reader.
    PacBio::Pancake::SeqDBReader reader(seqDBCache);

    // Reverse access all sequences by name and check if they're the same.
    for (int32_t i = static_cast<int32_t>(expected.size()) - 1; i >= 0; --i) {
        const auto& expRecord = expected[i];
        PacBio::Pancake::FastaSequenceId record;
        reader.GetSequence(record, expRecord.Name());
        EXPECT_EQ(expRecord, record);
    }
}

TEST(SeqDBReaderCompressed, GetSequenceByIdNotExisting)
{
    // Input values;
    const std::string inSeqFasta = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/in.fasta";
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir +
                                "/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb";

    // Load the DB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

    // Create a reader.
    PacBio::Pancake::SeqDBReader reader(seqDBCache);

    // Fetchin non-existent sequence names throws.
    {
        PacBio::Pancake::FastaSequenceId record;
        try {
            reader.GetSequence(record, 10);
            ASSERT_TRUE(false);
        } catch (const std::runtime_error& e) {
        }
    }
    {
        PacBio::Pancake::FastaSequenceId record;
        try {
            reader.GetSequence(record, -1);
            ASSERT_TRUE(false);
        } catch (const std::runtime_error& e) {
        }
    }
}

TEST(SeqDBReaderCompressed, GetSequenceByNameNotExisting)
{
    // Input values;
    const std::string inSeqFasta = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/in.fasta";
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir +
                                "/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb";

    // Load the DB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

    // Create a reader.
    PacBio::Pancake::SeqDBReader reader(seqDBCache);

    // Fetching non-existent sequence names does no throw but returns false.
    // This is because it's hard for a user to know if the sequence is existing
    // or not, so hard to prevent throws.
    {
        PacBio::Pancake::FastaSequenceId record;
        bool rv = reader.GetSequence(record, "some_nonexistent_name");
        ASSERT_FALSE(rv);
    }
}

TEST(SeqDBReaderCompressed, GetNextBatchAllAtOnce)
{
    // Input values;
    const int32_t batchSize = -1;
    const std::string inSeqFasta = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/in.fasta";
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir +
                                "/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb";

    // Load the expected sequences.
    const auto fastaSeqs = PacBio::PancakeTests::HelperLoadFasta(inSeqFasta);
    // Expect all sequences in the same block.
    const std::vector<std::vector<PacBio::Pancake::FastaSequenceId>> expected = {
        {PacBio::Pancake::FastaSequenceId{fastaSeqs[0], 0},
         PacBio::Pancake::FastaSequenceId{fastaSeqs[1], 1},
         PacBio::Pancake::FastaSequenceId{fastaSeqs[2], 2},
         PacBio::Pancake::FastaSequenceId{fastaSeqs[3], 3},
         PacBio::Pancake::FastaSequenceId{fastaSeqs[4], 4}}};

    // Load the DB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

    // Create a reader.
    PacBio::Pancake::SeqDBReader reader(seqDBCache);

    // Collect all sequences as a single string for comparison.
    std::vector<PacBio::Pancake::FastaSequenceId> records;
    std::vector<std::vector<PacBio::Pancake::FastaSequenceId>> results;
    while (reader.GetNextBatch(records, batchSize)) {
        results.emplace_back(records);
    }
    EXPECT_EQ(expected, results);
}

TEST(SeqDBReaderCompressed, GetNextBatchEverySeqSeparately)
{
    // Input values;
    const int32_t batchSize = 0;
    const std::string inSeqFasta = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/in.fasta";
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir +
                                "/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb";

    // Load the expected sequences.
    const auto fastaSeqs = PacBio::PancakeTests::HelperLoadFasta(inSeqFasta);
    // Expect all sequences in the same block.
    const std::vector<std::vector<PacBio::Pancake::FastaSequenceId>> expected = {
        {PacBio::Pancake::FastaSequenceId{fastaSeqs[0], 0}},
        {PacBio::Pancake::FastaSequenceId{fastaSeqs[1], 1}},
        {PacBio::Pancake::FastaSequenceId{fastaSeqs[2], 2}},
        {PacBio::Pancake::FastaSequenceId{fastaSeqs[3], 3}},
        {PacBio::Pancake::FastaSequenceId{fastaSeqs[4], 4}}};

    // Load the DB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

    // Create a reader.
    PacBio::Pancake::SeqDBReader reader(seqDBCache);

    // Collect all sequences as a single string for comparison.
    std::vector<PacBio::Pancake::FastaSequenceId> records;
    std::vector<std::vector<PacBio::Pancake::FastaSequenceId>> results;
    while (reader.GetNextBatch(records, batchSize)) {
        results.emplace_back(records);
    }
    EXPECT_EQ(expected, results);
}

TEST(SeqDBReaderCompressed, GetNextBatch10kbp)
{
    // Input values;
    const int32_t batchSize = 10000;
    const std::string inSeqFasta = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/in.fasta";
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir +
                                "/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb";

    // Load the expected sequences.
    const auto fastaSeqs = PacBio::PancakeTests::HelperLoadFasta(inSeqFasta);
    // Expect 3 blocks of sequences.
    const std::vector<std::vector<PacBio::Pancake::FastaSequenceId>> expected = {
        {PacBio::Pancake::FastaSequenceId{fastaSeqs[0], 0},
         PacBio::Pancake::FastaSequenceId{fastaSeqs[1], 1}},
        {PacBio::Pancake::FastaSequenceId{fastaSeqs[2], 2}},
        {PacBio::Pancake::FastaSequenceId{fastaSeqs[3], 3},
         PacBio::Pancake::FastaSequenceId{fastaSeqs[4], 4}}};

    // Load the DB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

    // Create a reader.
    PacBio::Pancake::SeqDBReader reader(seqDBCache);

    // Load sequences in batches of ~10kbp.
    std::vector<PacBio::Pancake::FastaSequenceId> records;
    std::vector<std::vector<PacBio::Pancake::FastaSequenceId>> results;
    while (reader.GetNextBatch(records, batchSize)) {
        results.emplace_back(records);
    }

    EXPECT_EQ(expected, results);
}

TEST(SeqDBReaderCompressed, GetBlockOutOfBounds)
{
    // Input values;
    const std::string inSeqFasta = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/in.fasta";
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir +
                                "/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb";
    const std::vector<int32_t> testBlocks = {-1, 10};

    // Load the DB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

    // Create a reader.
    PacBio::Pancake::SeqDBReader reader(seqDBCache);

    // Load sequences in batches of ~10kbp.
    std::vector<PacBio::Pancake::FastaSequenceId> records;
    for (const auto& blockId : testBlocks) {
        try {
            reader.GetBlock(records, blockId);
            ASSERT_TRUE(false);
        } catch (const std::runtime_error& e) {
        }
    }
}

TEST(SeqDBReaderCompressed, JumpToByName)
{
    // Input values;
    const std::string inSeqFasta = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/in.fasta";
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir +
                                "/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb";
    const std::vector<std::string> testNames = {
        "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105",   // 3
        "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292",  // 2
        "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001",  // 4
        "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852",   // 0
        "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983",  // 1
    };

    // Load the expected sequences.
    const auto fastaSeqs = PacBio::PancakeTests::HelperLoadFasta(inSeqFasta);
    const std::vector<PacBio::Pancake::FastaSequenceId> expected = {
        PacBio::Pancake::FastaSequenceId{fastaSeqs[3], 3},
        PacBio::Pancake::FastaSequenceId{fastaSeqs[2], 2},
        PacBio::Pancake::FastaSequenceId{fastaSeqs[4], 4},
        PacBio::Pancake::FastaSequenceId{fastaSeqs[0], 0},
        PacBio::Pancake::FastaSequenceId{fastaSeqs[1], 1}};

    // Load the DB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

    // Create a reader.
    PacBio::Pancake::SeqDBReader reader(seqDBCache);

    // Set the positon to the next sequence in the order specified by testNames.
    // The GetNext should then load that particular sequence.
    std::vector<PacBio::Pancake::FastaSequenceId> results;
    PacBio::Pancake::FastaSequenceId record;
    for (const auto& name : testNames) {
        reader.JumpTo(name);
        reader.GetNext(record);
        results.emplace_back(record);
    }

    EXPECT_EQ(expected, results);
}

TEST(SeqDBReaderCompressed, JumpToById)
{
    // Input values;
    const std::string inSeqFasta = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/in.fasta";
    const std::string inSeqDB = PacBio::PancakeTestsConfig::Data_Dir +
                                "/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb";
    const std::vector<int32_t> testIds = {3, 2, 4, 0, 1};

    // Load the expected sequences.
    const auto fastaSeqs = PacBio::PancakeTests::HelperLoadFasta(inSeqFasta);
    const std::vector<PacBio::Pancake::FastaSequenceId> expected = {
        PacBio::Pancake::FastaSequenceId{fastaSeqs[3], 3},
        PacBio::Pancake::FastaSequenceId{fastaSeqs[2], 2},
        PacBio::Pancake::FastaSequenceId{fastaSeqs[4], 4},
        PacBio::Pancake::FastaSequenceId{fastaSeqs[0], 0},
        PacBio::Pancake::FastaSequenceId{fastaSeqs[1], 1}};

    // Load the DB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

    // Create a reader.
    PacBio::Pancake::SeqDBReader reader(seqDBCache);

    // Set the positon to the next sequence in the order specified by testNames.
    // The GetNext should then load that particular sequence.
    std::vector<PacBio::Pancake::FastaSequenceId> results;
    PacBio::Pancake::FastaSequenceId record;
    for (const auto& id : testIds) {
        reader.JumpTo(id);
        reader.GetNext(record);
        results.emplace_back(record);
    }

    EXPECT_EQ(expected, results);
}

TEST(SeqDBReaderUncompressed, GetNext)
{
    /*
     * This is an equivalent test to TEST(SeqDBReaderCompressed, GetNext), but the
     * test-3.seqdb contains uncompressed sequences.
     * This is enough to prove that uncompressed loading works, because the rest
     * of the functionality is the same in both cases.
    */

    // Input values;
    const std::string inSeqFasta = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/in.fasta";
    const std::string inSeqDB =
        PacBio::PancakeTestsConfig::Data_Dir +
        "/seqdb-writer/test-3-uncompressed-each-seq-one-block-and-file.seqdb";

    // Load the expected sequences.
    // const std::string expected = PacBio::PancakeTests::HelperLoadFastaAsString(inSeqFasta);
    const auto expected = PacBio::PancakeTests::HelperLoadFasta(inSeqFasta);

    // Load the DB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

    // Create a reader.
    PacBio::Pancake::SeqDBReader reader(seqDBCache);

    // Collect all sequences as a single string for comparison.
    std::vector<PacBio::BAM::FastaSequence> result;
    PacBio::Pancake::FastaSequenceId record;
    while (reader.GetNext(record)) {
        result.emplace_back(record);
    }

    EXPECT_EQ(expected, result);
}

TEST(SeqDBReader, MalformedBytes)
{
    /*
     * The sequence length in these DBs is wrong - the DBs specify 10000 bp long
     * sequences, but in reality, the sequence is of length 5852 bases.
     * The test-4.seqdb runs the parser with the decompression step, and
     * test-5.seqdb without it.
     * Both DBs reuse existing data from previous tests (test-1 for compressed, and
     * test-5 for uncompressed).
    */

    // Input values;
    const std::string inSeqFasta = PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/in.fasta";
    const std::vector<std::string> inSeqDBList = {
        PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-4.seqdb",
        PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-5.seqdb"};

    // Run tests.
    for (const auto& inSeqDB : inSeqDBList) {
        // Load the DB.
        std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
            PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

        // Create a reader.
        PacBio::Pancake::SeqDBReader reader(seqDBCache);

        // Collect all sequences as a single string for comparison.
        PacBio::Pancake::FastaSequenceId record;
        while (true) {
            bool rv = false;
            try {
                rv = reader.GetNext(record);
                ASSERT_TRUE(false);
            } catch (const std::runtime_error& e) {
            }
            if (!rv) {
                break;
            }
        }
    }
}
