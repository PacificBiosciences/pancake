// Author: Ivan Sovic

#ifndef PANCAKE_SEQDB_READER_CACHED_H
#define PANCAKE_SEQDB_READER_CACHED_H

#include <pacbio/pancake/FastaSequenceId.h>
#include <pacbio/pancake/SeqDBIndexCache.h>
#include <pacbio/pancake/SeqDBReaderCached.h>
#include <memory>
#include <string>

namespace PacBio {
namespace Pancake {

class SeqDBReaderCached
{
public:
    SeqDBReaderCached(std::shared_ptr<PacBio::Pancake::SeqDBIndexCache>& seedDBCache,
                      int32_t blockId);
    SeqDBReaderCached(std::shared_ptr<PacBio::Pancake::SeqDBIndexCache>& seedDBCache,
                      const std::vector<int32_t>& seqIds);
    SeqDBReaderCached(std::shared_ptr<PacBio::Pancake::SeqDBIndexCache>& seedDBCache,
                      const std::vector<std::string>& seqNames);
    ~SeqDBReaderCached();

    const FastaSequenceId& GetSequence(int32_t seqId) const;
    const FastaSequenceId& GetSequence(const std::string& seqName) const;
    void GetSequence(FastaSequenceId& record, int32_t seqId);
    void GetSequence(FastaSequenceId& record, const std::string& seqName);
    const std::vector<PacBio::Pancake::FastaSequenceId>& records() const { return records_; }

private:
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBIndexCache_;
    std::vector<PacBio::Pancake::FastaSequenceId> records_;

    // Info to allow random access.
    std::unordered_map<std::string, int32_t> headerToOrdinalId_;
    std::unordered_map<int32_t, int32_t> seqIdToOrdinalId_;

    void LoadData_(int32_t blockId_);
    void LoadData_(const std::vector<int32_t>& seqIds);
    void LoadData_(const std::vector<std::string>& seqNames);
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEQDB_READER_CACHED_H