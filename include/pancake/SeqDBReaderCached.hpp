// Author: Ivan Sovic

#ifndef PANCAKE_SEQ_DB_READER_CACHED_HPP
#define PANCAKE_SEQ_DB_READER_CACHED_HPP

#include <pancake/FastaSequenceCachedStore.hpp>
#include <pancake/FastaSequenceId.hpp>
#include <pancake/SeqDBIndexCache.hpp>
#include <pancake/SeqDBReaderCached.hpp>

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

    const FastaSequenceCached& GetSequence(int32_t seqId) const;
    const FastaSequenceCached& GetSequence(const std::string& seqName) const;
    const std::vector<FastaSequenceCached>& records() const { return recordStore_.records(); }
    const FastaSequenceCachedStore& recordStore() const { return recordStore_; }

private:
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBIndexCache_;
    std::vector<PacBio::Pancake::FastaSequenceId> records_;
    FastaSequenceCachedStore recordStore_;

    void LoadData_(int32_t blockId_);
    void LoadData_(const std::vector<int32_t>& seqIds);
    void LoadData_(const std::vector<std::string>& seqNames);
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEQ_DB_READER_CACHED_HPP
