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

    for (int32_t i = 0; i < static_cast<int32_t>(records_.size()); ++i) {
        headerToOrdinalId_[records_[i].Name()] = i;
        seqIdToOrdinalId_[records_[i].Id()] = i;
    }
}

void SeqDBReaderCached::LoadData_(const std::vector<int32_t>& seqIds)
{
    records_.clear();

    // Create a SeedDB reader.
    PacBio::Pancake::SeqDBReader reader(seqDBIndexCache_);

    // Load sequences one by one.
    for (const auto& seqId : seqIds) {
        int32_t numRecords = static_cast<int32_t>(records_.size());
        // Get the records.
        PacBio::Pancake::FastaSequenceId record;
        reader.GetSequence(record, seqId);
        // Add it to the lookups.
        headerToOrdinalId_[record.Name()] = numRecords;
        seqIdToOrdinalId_[record.Id()] = numRecords;
        // Store the record.
        records_.emplace_back(std::move(record));
    }
}

void SeqDBReaderCached::LoadData_(const std::vector<std::string>& seqNames)
{
    records_.clear();

    // Create a SeedDB reader.
    PacBio::Pancake::SeqDBReader reader(seqDBIndexCache_);

    // Load sequences one by one.
    for (const auto& seqName : seqNames) {
        int32_t numRecords = static_cast<int32_t>(records_.size());
        // Get the records.
        PacBio::Pancake::FastaSequenceId record;
        reader.GetSequence(record, seqName);
        // Add it to the lookups.
        headerToOrdinalId_[record.Name()] = numRecords;
        seqIdToOrdinalId_[record.Id()] = numRecords;
        // Store the record.
        records_.emplace_back(std::move(record));
    }
}

const FastaSequenceId& SeqDBReaderCached::GetSequence(int32_t seqId) const
{
    auto it = seqIdToOrdinalId_.find(seqId);
    if (it == seqIdToOrdinalId_.end()) {
        std::ostringstream oss;
        oss << "(SeqDBReaderCached) Invalid seqId = " << seqId << ".";
        throw std::runtime_error(oss.str());
    }
    int32_t ordinalId = it->second;
    return records_[ordinalId];
}

void SeqDBReaderCached::GetSequence(FastaSequenceId& record, int32_t seqId)
{
    auto it = seqIdToOrdinalId_.find(seqId);
    if (it == seqIdToOrdinalId_.end()) {
        std::ostringstream oss;
        oss << "(SeqDBReaderCached) Invalid seqId = " << seqId << ".";
        throw std::runtime_error(oss.str());
    }
    int32_t ordinalId = it->second;
    record = records_[ordinalId];
}

const FastaSequenceId& SeqDBReaderCached::GetSequence(const std::string& seqName) const
{
    auto it = headerToOrdinalId_.find(seqName);
    if (it == headerToOrdinalId_.end()) {
        std::ostringstream oss;
        oss << "(SeqDBReaderCached) Invalid seqName: '" << seqName << "'";
        throw std::runtime_error(oss.str());
    }
    int32_t ordinalId = it->second;
    return records_[ordinalId];
}

void SeqDBReaderCached::GetSequence(FastaSequenceId& record, const std::string& seqName)
{
    auto it = headerToOrdinalId_.find(seqName);
    if (it == headerToOrdinalId_.end()) {
        std::ostringstream oss;
        oss << "(SeqDBReaderCached) Invalid seqName: '" << seqName << "'";
        throw std::runtime_error(oss.str());
    }
    int32_t ordinalId = it->second;
    record = records_[ordinalId];
}

}  // namespace Pancake
}  // namespace PacBio
