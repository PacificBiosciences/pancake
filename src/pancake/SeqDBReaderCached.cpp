// Authors: Ivan Sovic

#include <pancake/SeqDBReaderCached.hpp>

#include <pancake/SeqDBReader.hpp>

#include <sstream>

namespace PacBio {
namespace Pancake {

SeqDBReaderCached::SeqDBReaderCached(std::shared_ptr<PacBio::Pancake::SeqDBIndexCache>& seqDBCache,
                                     int32_t blockId)
    : seqDBIndexCache_(seqDBCache)
{
    ValidateSeqDBIndexCache(seqDBCache);
    LoadData_(blockId);
}

SeqDBReaderCached::SeqDBReaderCached(std::shared_ptr<PacBio::Pancake::SeqDBIndexCache>& seqDBCache,
                                     const std::vector<int32_t>& seqIds)
    : seqDBIndexCache_(seqDBCache)
{
    ValidateSeqDBIndexCache(seqDBCache);
    LoadData_(seqIds);
}

SeqDBReaderCached::SeqDBReaderCached(std::shared_ptr<PacBio::Pancake::SeqDBIndexCache>& seqDBCache,
                                     const std::vector<std::string>& seqNames)
    : seqDBIndexCache_(seqDBCache)
{
    ValidateSeqDBIndexCache(seqDBCache);
    LoadData_(seqNames);
}

SeqDBReaderCached::~SeqDBReaderCached() = default;

void SeqDBReaderCached::LoadData_(int32_t blockId)
{
    // Create a SeedDB reader.
    PacBio::Pancake::SeqDBReader reader(seqDBIndexCache_);

    // Fetch the block of data.
    reader.GetBlock(records_, blockId);

    // Initialize the record store with lookups.
    recordStore_.Clear();
    recordStore_.AddRecords(records_);
}

void SeqDBReaderCached::LoadData_(const std::vector<int32_t>& seqIds)
{
    records_.clear();

    // Create a SeedDB reader.
    PacBio::Pancake::SeqDBReader reader(seqDBIndexCache_);

    // Load sequences one by one.
    for (const auto& seqId : seqIds) {
        // Get the records.
        PacBio::Pancake::FastaSequenceId record;
        reader.GetSequence(record, seqId);
        // Store the record.
        records_.emplace_back(std::move(record));
    }
    recordStore_.Clear();
    recordStore_.AddRecords(records_);
}

void SeqDBReaderCached::LoadData_(const std::vector<std::string>& seqNames)
{
    records_.clear();

    // Create a SeedDB reader.
    PacBio::Pancake::SeqDBReader reader(seqDBIndexCache_);

    // Load sequences one by one.
    for (const auto& seqName : seqNames) {
        // Get the records.
        PacBio::Pancake::FastaSequenceId record;
        reader.GetSequence(record, seqName);
        // Store the record.
        records_.emplace_back(std::move(record));
    }
    recordStore_.Clear();
    recordStore_.AddRecords(records_);
}

const FastaSequenceCached& SeqDBReaderCached::GetSequence(int32_t seqId) const
{
    // The following line throws if the identified is not present in the store.
    return recordStore_.GetSequence(seqId);
}

const FastaSequenceCached& SeqDBReaderCached::GetSequence(const std::string& seqName) const
{
    // The following line throws if the identified is not present in the store.
    return recordStore_.GetSequence(seqName);
}

}  // namespace Pancake
}  // namespace PacBio
