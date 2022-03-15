// Authors: Ivan Sovic

#include <PancakeTestData.h>
#include <gtest/gtest.h>
#include <iostream>
#include <pancake/SeqDBReader.hpp>
#include <pancake/SeqDBReaderCachedBlock.hpp>

TEST(SeqDBReaderCachedBlock, BatchCompareWithSeqDBReader_UncompressedInput)
{
    /*
     * Batch test loading several SeqDBs in a row.
     * This test compares the same DB loaded with the Unit-Under-Test and with
     * a previously tested SeqDBReader. The loaded sequences should be the same.
     * The only difference is the format in which they return the sequences,
     * so they need to be collected in a comparable format first.
    */

    const std::vector<std::string> inDBs = {
        PacBio::PancakeTestsConfig::Data_Dir +
            "/seqdb-writer/test-3-uncompressed-each-seq-one-block-and-file.seqdb",
        PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-7-uncompressed-2blocks.seqdb",
        PacBio::PancakeTestsConfig::Data_Dir +
            "/seqdb-writer/test-9b-uncompressed-reversed-offsets.seqdb",
    };

    for (const auto& inSeqDB : inDBs) {
        // Load the SeedDB.
        std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
            PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

        const int32_t numBlocks = seqDBCache->blockLines.size();

        for (int32_t blockId = 0; blockId < numBlocks; ++blockId) {
            SCOPED_TRACE(inSeqDB + ", blockId = " + std::to_string(blockId));

            // "Reference" (or "truth") reader. This was tested earlier, thoroughly.
            // Collect the expected results for this block using a trusty reader.
            PacBio::Pancake::SeqDBReader readerTruth(seqDBCache);
            std::vector<PacBio::Pancake::FastaSequenceId> expected;
            readerTruth.GetBlock(expected, blockId);
            std::sort(expected.begin(), expected.end(),
                      [](const auto& a, const auto& b) { return a.Id() < b.Id(); });

            // Create a SeedDB reader under test.
            // Convert all FastaSequenceCached to FastaSequenceId for easier comparison.
            PacBio::Pancake::SeqDBReaderCachedBlock readerTest(seqDBCache, false);
            readerTest.LoadBlocks({blockId});
            std::vector<PacBio::Pancake::FastaSequenceId> results;
            for (const auto& record : readerTest.records()) {
                results.emplace_back(PacBio::Pancake::FastaSequenceId(
                    record.Name(), std::string(record.Bases(), record.Size()), record.Id()));
            }
            std::sort(results.begin(), results.end(),
                      [](const auto& a, const auto& b) { return a.Id() < b.Id(); });

            // Evaluate the current block.
            EXPECT_EQ(expected, results);
        }
    }
}

TEST(SeqDBReaderCachedBlock, BatchCompareWithSeqDBReader_CompressedInput)
{
    /*
     * Same as before, but the input DBs are compressed.
    */

    const std::vector<std::string> inDBs = {
        PacBio::PancakeTestsConfig::Data_Dir +
            "/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb",
        PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-6.seqdb",
        PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-8-compressed-2blocks.seqdb",
        PacBio::PancakeTestsConfig::Data_Dir +
            "/seqdb-writer/test-9a-compressed-reversed-offsets.seqdb",
    };

    for (const auto& inSeqDB : inDBs) {
        // Load the SeqDB.
        std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
            PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

        const int32_t numBlocks = seqDBCache->blockLines.size();

        for (int32_t blockId = 0; blockId < numBlocks; ++blockId) {
            SCOPED_TRACE(inSeqDB + ", blockId = " + std::to_string(blockId));

            // "Reference" (or "truth") reader. This was tested earlier, thoroughly.
            // Collect the expected results for this block using a trusty reader.
            PacBio::Pancake::SeqDBReader readerTruth(seqDBCache);
            std::vector<PacBio::Pancake::FastaSequenceId> expected;
            readerTruth.GetBlock(expected, blockId);
            std::sort(expected.begin(), expected.end(),
                      [](const auto& a, const auto& b) { return a.Id() < b.Id(); });

            // Create a SeedDB reader under test.
            // Convert all FastaSequenceCached to FastaSequenceId for easier comparison.
            PacBio::Pancake::SeqDBReaderCachedBlock readerTest(seqDBCache, false);
            readerTest.LoadBlocks({blockId});
            std::vector<PacBio::Pancake::FastaSequenceId> results;
            for (const auto& record : readerTest.records()) {
                results.emplace_back(PacBio::Pancake::FastaSequenceId(
                    record.Name(), std::string(record.Bases(), record.Size()), record.Id()));
            }
            std::sort(results.begin(), results.end(),
                      [](const auto& a, const auto& b) { return a.Id() < b.Id(); });

            // Evaluate the current block.
            EXPECT_EQ(expected, results);
        }
    }
}

TEST(SeqDBReaderCachedBlock, MultipleInputBlocks)
{
    /*
     * Load multiple blocks at once.
     * Tests both the compressed and uncompressed inputs. The SeqDBs are the same
     * except for the compression.
     * We're loading the same three blocks in all cases.
    */

    const std::vector<std::string> inDBs = {
        PacBio::PancakeTestsConfig::Data_Dir +
            "/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb",
        PacBio::PancakeTestsConfig::Data_Dir +
            "/seqdb-writer/test-3-uncompressed-each-seq-one-block-and-file.seqdb",
    };

    for (const auto& inSeqDB : inDBs) {
        SCOPED_TRACE(inSeqDB);

        const std::vector<int32_t> inBlocks = {1, 2, 3};

        // Load the SeqDB.
        std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
            PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

        // Collect all expected sequences for the specified input blocks
        // using an orthogonal reader. These are treated as the truth sequences.
        std::vector<PacBio::Pancake::FastaSequenceId> expected;
        for (const auto& blockId : inBlocks) {
            // "Reference" (or "truth") reader. This was tested earlier, thoroughly.
            // Collect the expected results for this block using a trusty reader.
            std::vector<PacBio::Pancake::FastaSequenceId> currExpected;
            PacBio::Pancake::SeqDBReader readerTruth(seqDBCache);
            readerTruth.GetBlock(currExpected, blockId);

            for (const auto& val : currExpected) {
                expected.emplace_back(val);
            }
        }
        std::sort(expected.begin(), expected.end(),
                  [](const auto& a, const auto& b) { return a.Id() < b.Id(); });

        // Create a unit under test.
        // Read the sequences for the specified blocks, and convert all
        // FastaSequenceCached to FastaSequenceId for easier comparison.
        PacBio::Pancake::SeqDBReaderCachedBlock readerTest(seqDBCache, false);
        readerTest.LoadBlocks(inBlocks);
        std::vector<PacBio::Pancake::FastaSequenceId> results;
        for (const auto& record : readerTest.records()) {
            results.emplace_back(PacBio::Pancake::FastaSequenceId(
                record.Name(), std::string(record.Bases(), record.Size()), record.Id()));
        }
        std::sort(results.begin(), results.end(),
                  [](const auto& a, const auto& b) { return a.Id() < b.Id(); });

        // Evaluate the current block.
        EXPECT_EQ(expected, results);
    }
}

TEST(SeqDBReaderCachedBlock, LoadSequences_SeqId)
{
    /*
     * Load sequences by their sequence ID.
     * Compare with the previously tested SeqDBReader.
    */

    const std::vector<std::string> inDBs = {
        PacBio::PancakeTestsConfig::Data_Dir +
            "/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb",
        PacBio::PancakeTestsConfig::Data_Dir +
            "/seqdb-writer/test-3-uncompressed-each-seq-one-block-and-file.seqdb",
    };

    for (const auto& inSeqDB : inDBs) {
        SCOPED_TRACE(inSeqDB);

        const std::vector<int32_t> seqIds = {1, 2, 3};  // Load the same sequence multiple times.

        // Load the SeqDB.
        std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
            PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

        // Collect all expected sequences for the specified input blocks
        // using an orthogonal reader. These are treated as the truth sequences.
        PacBio::Pancake::SeqDBReader readerTruth(seqDBCache);
        std::vector<PacBio::Pancake::FastaSequenceId> expected;
        for (const auto& seqId : seqIds) {
            PacBio::Pancake::FastaSequenceId record;
            readerTruth.GetSequence(record, seqId);
            expected.emplace_back(record);
        }
        std::sort(expected.begin(), expected.end(),
                  [](const auto& a, const auto& b) { return a.Id() < b.Id(); });

        // Create a unit under test.
        // Read the sequences for the specified blocks, and convert all
        // FastaSequenceCached to FastaSequenceId for easier comparison.
        PacBio::Pancake::SeqDBReaderCachedBlock readerTest(seqDBCache, false);

        // Run the same exact test multiple times, to see that there aren't any sideffects to
        // reusing the same reader to load another set of sequences.
        for (int32_t iteration = 0; iteration < 3; ++iteration) {
            readerTest.LoadSequences(seqIds);
            std::vector<PacBio::Pancake::FastaSequenceId> results;
            for (const auto& record : readerTest.records()) {
                results.emplace_back(PacBio::Pancake::FastaSequenceId(
                    record.Name(), std::string(record.Bases(), record.Size()), record.Id()));
            }
            std::sort(results.begin(), results.end(),
                      [](const auto& a, const auto& b) { return a.Id() < b.Id(); });

            // Evaluate the current block.
            EXPECT_EQ(expected, results);
        }
    }
}

TEST(SeqDBReaderCachedBlock, LoadSequences_SeqId_DuplicateIds_ShouldFail)
{
    /*
     * Load sequences by their sequence ID.
     * Compare with the previously tested SeqDBReader.
    */

    const std::vector<std::string> inDBs = {
        PacBio::PancakeTestsConfig::Data_Dir +
            "/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb",
        PacBio::PancakeTestsConfig::Data_Dir +
            "/seqdb-writer/test-3-uncompressed-each-seq-one-block-and-file.seqdb",
    };

    for (const auto& inSeqDB : inDBs) {
        SCOPED_TRACE(inSeqDB);

        const std::vector<int32_t> seqIds = {1, 2, 3, 2,
                                             2};  // Load the same sequence multiple times.

        // Load the SeqDB.
        std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
            PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

        // Create a unit under test.
        PacBio::Pancake::SeqDBReaderCachedBlock readerTest(seqDBCache, false);
        EXPECT_THROW({ readerTest.LoadSequences(seqIds); }, std::runtime_error);
    }
}

TEST(SeqDBReaderCachedBlock, LoadSequences_SeqName)
{
    /*
     * Load sequences by their sequence name.
     * Compare with the previously tested SeqDBReader.
    */

    const std::vector<std::string> inDBs = {
        PacBio::PancakeTestsConfig::Data_Dir +
            "/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb",
        PacBio::PancakeTestsConfig::Data_Dir +
            "/seqdb-writer/test-3-uncompressed-each-seq-one-block-and-file.seqdb",
    };

    for (const auto& inSeqDB : inDBs) {
        SCOPED_TRACE(inSeqDB);

        const std::vector<std::string> seqNames = {
            "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983",
            "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292",
            "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105",
        };

        // Load the SeqDB.
        std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
            PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

        // Collect all expected sequences for the specified input blocks
        // using an orthogonal reader. These are treated as the truth sequences.

        PacBio::Pancake::SeqDBReader readerTruth(seqDBCache);
        std::vector<PacBio::Pancake::FastaSequenceId> expected;
        for (const auto& seqName : seqNames) {
            PacBio::Pancake::FastaSequenceId record;
            readerTruth.GetSequence(record, seqName);
            expected.emplace_back(record);
        }
        std::sort(expected.begin(), expected.end(),
                  [](const auto& a, const auto& b) { return a.Id() < b.Id(); });

        // Create a unit under test.
        // Read the sequences for the specified blocks, and convert all
        // FastaSequenceCached to FastaSequenceId for easier comparison.
        PacBio::Pancake::SeqDBReaderCachedBlock readerTest(seqDBCache, false);

        // Run the same exact test multiple times, to see that there aren't any sideffects to
        // reusing the same reader to load another set of sequences.
        for (int32_t iteration = 0; iteration < 3; ++iteration) {
            readerTest.LoadSequences(seqNames);
            std::vector<PacBio::Pancake::FastaSequenceId> results;
            for (const auto& record : readerTest.records()) {
                results.emplace_back(PacBio::Pancake::FastaSequenceId(
                    record.Name(), std::string(record.Bases(), record.Size()), record.Id()));
            }
            std::sort(results.begin(), results.end(),
                      [](const auto& a, const auto& b) { return a.Id() < b.Id(); });

            // Evaluate the current block.
            EXPECT_EQ(expected, results);
        }
    }
}

TEST(SeqDBReaderCachedBlock, LoadSequences_SeqName_DuplicateNames_ShouldFail)
{
    /*
     * Load sequences by their sequence name.
     * Compare with the previously tested SeqDBReader.
    */

    const std::vector<std::string> inDBs = {
        PacBio::PancakeTestsConfig::Data_Dir +
            "/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb",
        PacBio::PancakeTestsConfig::Data_Dir +
            "/seqdb-writer/test-3-uncompressed-each-seq-one-block-and-file.seqdb",
    };

    for (const auto& inSeqDB : inDBs) {
        SCOPED_TRACE(inSeqDB);

        const std::vector<std::string> seqNames = {
            "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/"
            "0_11983",
            "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/"
            "0_24292",
            "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/"
            "0_5105",
            "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/"
            "0_5105",  // Same sequence multiple times.
            "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/"
            "0_5105",
        };

        // Load the SeqDB.
        std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
            PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

        // Create a unit under test.
        PacBio::Pancake::SeqDBReaderCachedBlock readerTest(seqDBCache, false);
        EXPECT_THROW({ readerTest.LoadSequences(seqNames); }, std::runtime_error);
    }
}

TEST(SeqDBReaderCachedBlock, LoadAndHPCompress_From_in_2_small)
{
    /*
     * Load sequences and compress homopolymers in the cache.
    */

    const std::vector<std::string> inDBs = {
        PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/test-10a-compressed-in-2-small.seqdb",
        PacBio::PancakeTestsConfig::Data_Dir +
            "/seqdb-writer/test-10b-uncompressed-in-2-small.seqdb",
    };

    std::vector<PacBio::Pancake::FastaSequenceId> expected = {
        PacBio::Pancake::FastaSequenceId("Seq1", "ATCGTCAGCGCG", 0),
        PacBio::Pancake::FastaSequenceId("Seq2", "AGCGCAGCAGACATACGATATGCA", 1),
        PacBio::Pancake::FastaSequenceId("Seq3", "C", 2),
        PacBio::Pancake::FastaSequenceId("Seq4", "A", 3),
        PacBio::Pancake::FastaSequenceId("Seq5", "ACTCGCGCAGTGACATCGCAGAGTATAGCAGCGTACAGCGTACGTG",
                                         4),
    };

    int32_t blockId = 0;

    for (const auto& inSeqDB : inDBs) {
        SCOPED_TRACE(inSeqDB);

        // Load the SeqDB.
        std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
            PacBio::Pancake::LoadSeqDBIndexCache(inSeqDB);

        // Create a unit under test.
        // Read the sequences for the specified blocks and compress homopolymers.
        // Convert all FastaSequenceCached to FastaSequenceId for easier comparison.
        PacBio::Pancake::SeqDBReaderCachedBlock readerTest(seqDBCache, true);
        readerTest.LoadBlocks({blockId});
        std::vector<PacBio::Pancake::FastaSequenceId> results;
        for (const auto& record : readerTest.records()) {
            results.emplace_back(PacBio::Pancake::FastaSequenceId(
                record.Name(), std::string(record.Bases(), record.Size()), record.Id()));
        }
        std::sort(results.begin(), results.end(),
                  [](const auto& a, const auto& b) { return a.Id() < b.Id(); });

        for (const auto& record : results) {
            std::cerr << ">" << record.Name() << "-" << record.Id() << "\n"
                      << record.Bases() << "\n";
        }

        // Evaluate the current block.
        EXPECT_EQ(expected, results);
    }
}

TEST(SeqDBReaderCachedBlock, GetSeqDBContiguousParts_NormalSingleBlock)
{
    /*
     * Fetch the byte span of a single small block of 1 sequence.
     * It doesn't span more than 1 file or sequences.
    */

    const std::string inSeqDB =
        R"(V	0.1.0
C	1
F	0	test-1.seqdb.0.seq	1	1463	5852
F	1	test-1.seqdb.1.seq	1	2996	11983
F	2	test-1.seqdb.2.seq	1	6073	24292
F	3	test-1.seqdb.3.seq	1	1277	5105
F	4	test-1.seqdb.4.seq	1	4751	19001
S	0	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852	0	0	1463	5852	1	0	5852
S	1	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983	1	0	2996	11983	1	0	11983
S	2	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292	2	0	6073	24292	1	0	24292
S	3	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105	3	0	1277	5105	1	0	5105
S	4	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001	4	0	4751	19001	1	0	19001
B	0	0	1	1463	5852
B	1	1	2	2996	11983
B	2	2	3	6073	24292
B	3	3	4	1277	5105
B	4	4	5	4751	19001
    )";
    const int32_t blockId = 0;

    // Expected.
    // Tuple: (file_id, startOffset, endOffset, startId, endId, file_offset_end)
    const std::vector<PacBio::Pancake::ContiguousFilePart> expected = {
        PacBio::Pancake::ContiguousFilePart{0, 0, 1463, {0}}};

    // Load the SeedDB.
    std::istringstream is(inSeqDB);
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(is, "filename.seqdb");

    // Run.
    const auto results = PacBio::Pancake::GetSeqDBContiguousParts(seqDBCache, blockId);

    // Evaluate.
    EXPECT_EQ(expected, results);
}

TEST(SeqDBReaderCachedBlock, GetSeqDBContiguousParts_NormalMultipleFilesBlock0)
{
    /*
     * Fetch the byte span of a block of 4 sequences, where each sequence is
     * stored in a different file.
    */

    const std::string inSeqDB =
        R"(V	0.1.0
C	1
F	0	test-1.seqdb.0.seq	1	1463	5852
F	1	test-1.seqdb.1.seq	1	2996	11983
F	2	test-1.seqdb.2.seq	1	6073	24292
F	3	test-1.seqdb.3.seq	1	1277	5105
F	4	test-1.seqdb.4.seq	1	4751	19001
S	0	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852	0	0	1463	5852	1	0	5852
S	1	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983	1	0	2996	11983	1	0	11983
S	2	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292	2	0	6073	24292	1	0	24292
S	3	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105	3	0	1277	5105	1	0	5105
S	4	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001	4	0	4751	19001	1	0	19001
B	0	0	4	11809	47232
B	1	4	5	4751	19001
    )";
    const int32_t blockId = 0;

    // Expected.
    // Tuple: (file_id, startOffset, endOffset, startId, endId, file_offset_end)
    const std::vector<PacBio::Pancake::ContiguousFilePart> expected = {
        {0, 0, 1463, {0}},
        {1, 0, 2996, {1}},
        {2, 0, 6073, {2}},
        {3, 0, 1277, {3}},
    };

    // Load the SeedDB.
    std::istringstream is(inSeqDB);
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(is, "filename.seqdb");

    // Run.
    const auto results = PacBio::Pancake::GetSeqDBContiguousParts(seqDBCache, blockId);

    // Evaluate.
    EXPECT_EQ(expected, results);
}
TEST(SeqDBReaderCachedBlock, GetSeqDBContiguousParts_NormalMultipleFilesBlock1)
{
    /*
     * Fetch the byte span of a block of 4 sequences, where each sequence is
     * stored in a different file.
    */

    const std::string inSeqDB =
        R"(V	0.1.0
C	1
F	0	test-1.seqdb.0.seq	1	1463	5852
F	1	test-1.seqdb.1.seq	1	2996	11983
F	2	test-1.seqdb.2.seq	1	6073	24292
F	3	test-1.seqdb.3.seq	1	1277	5105
F	4	test-1.seqdb.4.seq	1	4751	19001
S	0	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852	0	0	1463	5852	1	0	5852
S	1	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983	1	0	2996	11983	1	0	11983
S	2	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292	2	0	6073	24292	1	0	24292
S	3	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105	3	0	1277	5105	1	0	5105
S	4	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001	4	0	4751	19001	1	0	19001
B	0	0	4	11809	47232
B	1	4	5	4751	19001
    )";
    const int32_t blockId = 1;

    // Expected.
    // Tuple: (file_id, startOffset, endOffset, startId, endId, file_offset_end)
    const std::vector<PacBio::Pancake::ContiguousFilePart> expected = {
        {4, 0, 4751, {4}},
    };

    // Load the SeedDB.
    std::istringstream is(inSeqDB);
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(is, "filename.seqdb");

    // Run.
    const auto results = PacBio::Pancake::GetSeqDBContiguousParts(seqDBCache, blockId);

    // Evaluate.
    EXPECT_EQ(expected, results);
}

TEST(SeqDBReaderCachedBlock, GetSeqDBContiguousParts_NormalTwoBlocksWithGap)
{
    /*
     * Fetch the byte span of a block of 4 sequences, where the first two and last
     * two are separated by a filtered sequence. There should be two parts because
     * of that.
    */

    const std::string inSeqDB =
        R"(V	0.1.0
C	1
F	0	test-6.seqdb.0.seq	4	10487	41941
S	0	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852	0	0	1463	5852	1	0	5852
S	1	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983	0	1463	2996	11983	1	0	11983
S	2	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105	0	10532	1277	5105	1	0	5105
S	3	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001	0	11809	4751	19001	1	0	19001
B	0	0	4	10487	41941
    )";
    const int32_t blockId = 0;

    // Expected.
    // Tuple: (file_id, file_offset_start, file_offset_end)
    const std::vector<PacBio::Pancake::ContiguousFilePart> expected = {{0, 0, 4459, {0, 1}},
                                                                       {0, 10532, 16560, {2, 3}}};

    // Load the SeedDB.
    std::istringstream is(inSeqDB);
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(is, "filename.seqdb");

    // Run.
    const auto results = PacBio::Pancake::GetSeqDBContiguousParts(seqDBCache, blockId);

    // Evaluate.
    EXPECT_EQ(expected, results);
}

TEST(SeqDBReaderCachedBlock, GetSeqDBContiguousParts_BlockOutOfBoundsThrows)
{
    /*
     * Block is out of bounds, it should throw.
    */

    const std::string inSeqDB =
        R"(V	0.1.0
C	1
F	0	test-6.seqdb.0.seq	4	10487	41941
S	0	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852	0	0	1463	5852	1	0	5852
S	1	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983	0	1463	2996	11983	1	0	11983
S	2	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105	0	10532	1277	5105	1	0	5105
S	3	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001	0	11809	4751	19001	1	0	19001
B	0	0	4	10487	41941
    )";
    const int32_t blockId = 123;

    // Load the SeedDB.
    std::istringstream is(inSeqDB);
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(is, "filename.seqdb");

    // Run and evaluate.
    EXPECT_THROW({ auto results = PacBio::Pancake::GetSeqDBContiguousParts(seqDBCache, blockId); },
                 std::runtime_error);
}

TEST(SeqDBReaderCachedBlock, GetSeqDBContiguousParts_MalformedBlockThrows)
{
    /*
     * Block is malformed, referencing sequences which do not exist in the index.
     * It should throw.
    */

    const std::string inSeqDB =
        R"(V	0.1.0
C	1
F	0	test-6.seqdb.0.seq	4	10487	41941
S	0	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852	0	0	1463	5852	1	0	5852
S	1	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983	0	1463	2996	11983	1	0	11983
S	2	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105	0	10532	1277	5105	1	0	5105
S	3	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001	0	11809	4751	19001	1	0	19001
B	0	0	10	10487	41941
    )";
    const int32_t blockId = 0;

    // Load the SeedDB.
    std::istringstream is(inSeqDB);
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(is, "filename.seqdb");

    // Run and evaluate.
    EXPECT_THROW({ auto results = PacBio::Pancake::GetSeqDBContiguousParts(seqDBCache, blockId); },
                 std::runtime_error);
}

TEST(SeqDBReaderCachedBlock, GetSeqDBContiguousParts_OverlappingBytesThrows)
{
    /*
     * One 'S' line is malformed (S4), it begins and overlaps S3. This should throw
     * in the GetSeqDBContiguousParts because sequenes should be distinct byte blocks for each sequence.
    */

    const std::string inSeqDB =
        R"(V	0.1.0
C	1
F	0	test-6.seqdb.0.seq	4	10487	41941
S	0	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852	0	0	1463	5852	1	0	5852
S	1	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983	0	1463	2996	11983	1	0	11983
S	2	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105	0	10532	1277	5105	1	0	5105
S	3	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001	0	11000	4751	19001	1	0	19001
B	0	0	4	10487	41941
    )";
    const int32_t blockId = 0;

    // Load the SeedDB.
    std::istringstream is(inSeqDB);
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(is, "filename.seqdb");

    // Run and evaluate.
    EXPECT_THROW({ auto results = PacBio::Pancake::GetSeqDBContiguousParts(seqDBCache, blockId); },
                 std::runtime_error);
}

TEST(SeqDBReaderCachedBlock, GetSeqDBContiguousParts_OutOfOrder)
{
    /*
     * This is a valid case, where the order of sequences permuted in the SeedDB. This
     * should not throw.
    */

    const std::string inSeqDB =
        R"(V	0.1.0
C	1
F	0	test-6.seqdb.0.seq	5	16560	66233
S	0	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001	0	11809	4751	19001	1	0	19001
S	1	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105	0	10532	1277	5105	1	0	5105
S	2	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292	0	4459	6073	24292	1	0	24292
S	3	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983	0	1463	2996	11983	1	0	11983
S	4	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852	0	0	1463	5852	1	0	5852
B	0	0	5	16560	66233
    )";
    const int32_t blockId = 0;

    // Expected.
    // Tuple: (file_id, file_offset_start, file_offset_end)
    const std::vector<PacBio::Pancake::ContiguousFilePart> expected = {
        {0, 0, 16560, {4, 3, 2, 1, 0}}};

    // Load the SeedDB.
    std::istringstream is(inSeqDB);
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(is, "filename.seqdb");

    // Run.
    const auto results = PacBio::Pancake::GetSeqDBContiguousParts(seqDBCache, blockId);

    // Evaluate.
    EXPECT_EQ(expected, results);
}

TEST(SeqDBReaderCachedBlock, GetSeqDBContiguousParts_InputIsAsListOfSeqIDs)
{
    /*
     * Test construction of contiguous file parts but from a list of sequence IDs
     * instead of a block ID.
    */

    const std::string inSeqDB =
        R"(V	0.1.0
C	1
F	0	test-6.seqdb.0.seq	5	16560	66233
S	0	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852	0	0	1463	5852	1	0	5852
S	1	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983	0	1463	2996	11983	1	0	11983
S	2	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292	0	4459	6073	24292	1	0	24292
S	3	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105	0	10532	1277	5105	1	0	5105
S	4	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001	0	11809	4751	19001	1	0	19001
B	0	0	5	16560	66233
    )";
    const std::vector<int32_t> seqIds = {1, 2, 4};

    // Expected.
    // Tuple: (file_id, file_offset_start, file_offset_end)
    const std::vector<PacBio::Pancake::ContiguousFilePart> expected = {{0, 1463, 10532, {1, 2}},
                                                                       {0, 11809, 16560, {4}}};

    // Load the SeedDB.
    std::istringstream is(inSeqDB);
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(is, "filename.seqdb");

    // Run.
    const auto results = PacBio::Pancake::GetSeqDBContiguousParts(seqDBCache, seqIds);

    // Evaluate.
    EXPECT_EQ(expected, results);
}

TEST(SeqDBReaderCachedBlock, GetSeqDBContiguousParts_InputIsAsListOfSeqNames)
{
    /*
     * Test construction of contiguous file parts but from a list of sequence IDs
     * instead of a block ID.
    */

    const std::string inSeqDB =
        R"(V	0.1.0
C	1
F	0	test-6.seqdb.0.seq	5	16560	66233
S	0	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852	0	0	1463	5852	1	0	5852
S	1	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983	0	1463	2996	11983	1	0	11983
S	2	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292	0	4459	6073	24292	1	0	24292
S	3	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105	0	10532	1277	5105	1	0	5105
S	4	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001	0	11809	4751	19001	1	0	19001
B	0	0	5	16560	66233
    )";
    const std::vector<std::string> seqNames = {
        "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983",
        "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292",
        "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001"};

    // Expected.
    // Tuple: (file_id, file_offset_start, file_offset_end)
    const std::vector<PacBio::Pancake::ContiguousFilePart> expected = {{0, 1463, 10532, {1, 2}},
                                                                       {0, 11809, 16560, {4}}};

    // Load the SeedDB.
    std::istringstream is(inSeqDB);
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(is, "filename.seqdb");

    // Run.
    const auto results = PacBio::Pancake::GetSeqDBContiguousParts(seqDBCache, seqNames);

    // Evaluate.
    EXPECT_EQ(expected, results);
}
