#include <PancakeTestData.h>
#include <gtest/gtest.h>
#include <pacbio/pancake/FastaSequenceCachedStore.h>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

TEST(FastaSequenceCachedStore, ArrayOfTests)
{
    struct TestData
    {
        std::string testName;
        // Tuple: <header, id, seq>
        std::vector<std::tuple<std::string, int32_t, std::string>> seqs;
    };

    // clang-format off
    std::vector<TestData> testData = {
        {
            "Empty input",
            {},
        },
        {
            "Single seq",
            {
                {"name1", 15, "ACTG"},
            },
        },
        {
            "Multiple seqs",
            {
                {"name1", 15, "ACTG"},
                {"name8", 24, "ACTG"},
                {"name123", 123, "ACTG"},
                {"name2", 2, "ACTG"},
            },
        },
    };
    // clang-format on

    for (const auto& data : testData) {
        // Create the input FastaSequenceCached objects.
        std::vector<PacBio::Pancake::FastaSequenceCached> records;
        for (const auto& seqData : data.seqs) {
            const auto& header = std::get<0>(seqData);
            const auto& id = std::get<1>(seqData);
            const auto& seq = std::get<2>(seqData);
            records.emplace_back(
                PacBio::Pancake::FastaSequenceCached(header, seq.c_str(), seq.size(), id));
        }

        // Initialize the record store from the constructor.
        PacBio::Pancake::FastaSequenceCachedStore recordStoreFromConstructor(records);

        // Initialize the record store by using the API to add all records at once.
        PacBio::Pancake::FastaSequenceCachedStore recordStoreAllAtOnce;
        recordStoreAllAtOnce.AddRecords(records);

        // Initialize the record store by adding one record at a time.
        PacBio::Pancake::FastaSequenceCachedStore recordStoreOneByOne;
        for (const auto& record : records) {
            recordStoreOneByOne.AddRecord(record);
        }

        // Test that all initialization methods produce the same result.
        EXPECT_EQ(recordStoreFromConstructor, recordStoreAllAtOnce);
        EXPECT_EQ(recordStoreFromConstructor, recordStoreOneByOne);

        // Compare the actual record data with the input data.
        EXPECT_EQ(recordStoreFromConstructor.Size(), records.size());
        for (size_t i = 0; i < recordStoreFromConstructor.Size(); ++i) {
            const auto& resultRecord = recordStoreFromConstructor.records()[i];
            const auto& expectedRecord = records[i];
            EXPECT_EQ(expectedRecord, resultRecord);
        }

        // Test fetching of sequences.
        for (const auto& expectedRecord : records) {
            // Fetch a record by the header.
            {
                const auto& resultRecord =
                    recordStoreFromConstructor.GetSequence(expectedRecord.Name());
                EXPECT_EQ(expectedRecord, resultRecord);
            }
            // Fetch a record by the ID.
            {
                const auto& resultRecord =
                    recordStoreFromConstructor.GetSequence(expectedRecord.Id());
                EXPECT_EQ(expectedRecord, resultRecord);
            }
            // Fetch a copy of the record by the header.
            {
                PacBio::Pancake::FastaSequenceCached resultRecord;
                bool rv =
                    recordStoreFromConstructor.GetSequence(resultRecord, expectedRecord.Name());
                EXPECT_EQ(true, rv);
                EXPECT_EQ(expectedRecord, resultRecord);
            }
            // Fetch a copy of the record by the ID.
            {
                PacBio::Pancake::FastaSequenceCached resultRecord;
                bool rv = recordStoreFromConstructor.GetSequence(resultRecord, expectedRecord.Id());
                EXPECT_EQ(true, rv);
                EXPECT_EQ(expectedRecord, resultRecord);
            }
        }

        // clang-format off
        // Test fetching of a sequence that does not exist.
        {
            // Fetch a record by the header.
            {
                EXPECT_THROW(
                    {
                        recordStoreFromConstructor.GetSequence("nonexistent-seq");
                    },
                    std::runtime_error
                );
            }

            // Fetch a record by the ID.
            {
                EXPECT_THROW(
                    {
                        recordStoreFromConstructor.GetSequence(1234567);
                    },
                    std::runtime_error
                );
            }

            // Fetch a record by the header.
            {
                PacBio::Pancake::FastaSequenceCached resultRecord;
                bool rv = recordStoreFromConstructor.GetSequence(resultRecord, "nonexistent-seq");
                EXPECT_EQ(false, rv);
            }

            // Fetch a record by the ID.
            {
                PacBio::Pancake::FastaSequenceCached resultRecord;
                bool rv = recordStoreFromConstructor.GetSequence(resultRecord, 1234567);
                EXPECT_EQ(false, rv);
            }
        }
        // clang-format on
    }
}
