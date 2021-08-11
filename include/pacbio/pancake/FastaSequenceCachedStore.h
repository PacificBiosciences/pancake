/*
 * FastaSequenceCachedStore.h
 *
 *  Created on: Jan 26, 2021
 *      Author: Ivan Sovic
 */

#ifndef PANCAKE_FASTA_SEQUENCE_CACHED_STORE_H
#define PANCAKE_FASTA_SEQUENCE_CACHED_STORE_H

#include <pacbio/pancake/FastaSequenceCached.h>
#include <pacbio/pancake/FastaSequenceId.h>
#include <cstdint>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace PacBio {
namespace Pancake {

/*
 * This is a container which stores a vector of FastaSequenceCached objects.
 * Each record is associated with a numeric ID and a header. Since the collection
 * of records in a vector might be just one chunk from a larger dataset, the record
 * IDs do not have to correspond to the ordinal IDs of the records in the vector.
 * This container keeps track of the ID relation, and allows a random access
 * based on the sequence ID or sequence header.
*/
class FastaSequenceCachedStore
{
public:
    friend bool operator==(const FastaSequenceCachedStore& lhs,
                           const FastaSequenceCachedStore& rhs);

    FastaSequenceCachedStore() = default;

    FastaSequenceCachedStore(const std::vector<FastaSequenceCached>& records)
    {
        AddRecords(records);
    }

    FastaSequenceCachedStore(const std::vector<FastaSequenceId>& records) { AddRecords(records); }

    void Clear()
    {
        records_.clear();
        headerToOrdinalId_.clear();
        seqIdToOrdinalId_.clear();
    }

    size_t Size() const { return records_.size(); }
    const std::vector<FastaSequenceCached>& records() const { return records_; }
    std::vector<FastaSequenceCached>& records() { return records_; }

    /*
     * Intentionally take a copy of the record, so we can move it.
     * This should allow the clients to move their copy from the outside too.
    */
    void AddRecord(FastaSequenceCached record)
    {
        seqIdToOrdinalId_[record.Id()] = records_.size();
        headerToOrdinalId_[record.Name()] = records_.size();
        records_.emplace_back(std::move(record));
    }

    void AddRecord(const FastaSequenceId& record)
    {
        auto seqCache = PacBio::Pancake::FastaSequenceCached(record.Name(), record.Bases().c_str(),
                                                             record.Bases().size(), record.Id());
        seqIdToOrdinalId_[record.Id()] = records_.size();
        headerToOrdinalId_[record.Name()] = records_.size();
        records_.emplace_back(std::move(seqCache));
    }

    void AddRecords(const std::vector<FastaSequenceCached>& records)
    {
        for (const auto& record : records) {
            AddRecord(record);
        }
    }

    void AddRecords(const std::vector<FastaSequenceId>& records)
    {
        for (const auto& record : records) {
            AddRecord(record);
        }
    }

    const FastaSequenceCached& GetSequence(int32_t seqId) const
    {
        const auto it = seqIdToOrdinalId_.find(seqId);
        if (it == seqIdToOrdinalId_.end()) {
            std::ostringstream oss;
            oss << "(FastaSequenceCachedStore::GetSequence) Invalid seqId, not found in the "
                   "records. seqId = "
                << seqId << ", records_.size() = " << records_.size();
            throw std::runtime_error(oss.str());
        }
        const int32_t ordinalId = it->second;
        return records_[ordinalId];
    }

    const FastaSequenceCached& GetSequence(const std::string& seqName) const
    {
        const auto it = headerToOrdinalId_.find(seqName);
        if (it == headerToOrdinalId_.end()) {
            std::ostringstream oss;
            oss << "(FastaSequenceCachedStore::GetSequence) Invalid seqName, not found in the "
                   "records. seqName = "
                << seqName;
            throw std::runtime_error(oss.str());
        }
        const int32_t ordinalId = it->second;
        return records_[ordinalId];
    }

    bool GetSequence(FastaSequenceCached& record, int32_t seqId) const
    {
        const auto it = seqIdToOrdinalId_.find(seqId);
        if (it == seqIdToOrdinalId_.end()) {
            return false;
        }
        const int32_t ordinalId = it->second;
        record = records_[ordinalId];
        return true;
    }

    bool GetSequence(FastaSequenceCached& record, const std::string& seqName) const
    {
        const auto it = headerToOrdinalId_.find(seqName);
        if (it == headerToOrdinalId_.end()) {
            return false;
        }
        const int32_t ordinalId = it->second;
        record = records_[ordinalId];
        return true;
    }

private:
    std::vector<FastaSequenceCached> records_;
    std::unordered_map<std::string, int32_t> headerToOrdinalId_;
    std::unordered_map<int32_t, int32_t> seqIdToOrdinalId_;
};

inline bool operator==(const FastaSequenceCachedStore& lhs, const FastaSequenceCachedStore& rhs)
{
    return lhs.records_ == rhs.records_ && lhs.headerToOrdinalId_ == rhs.headerToOrdinalId_ &&
           lhs.seqIdToOrdinalId_ == rhs.seqIdToOrdinalId_;
}

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_FASTA_SEQUENCE_CACHED_STORE_H