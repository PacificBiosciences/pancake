// Author: Ivan Sovic

/*
 * Similar to SeqDBReaderCached, but it doesn't parse the sequences one by one.
 * Instead, it loads all the data as a block, and provides sequences via C-style
 * pointers to the local data.
 * This works only for uncompressed sequences.
*/

#ifndef PANCAKE_SEQDB_READER_CACHED_BLOCK_H
#define PANCAKE_SEQDB_READER_CACHED_BLOCK_H

#include <pacbio/pancake/FastaSequenceCached.h>
#include <pacbio/pancake/FastaSequenceCachedStore.h>
#include <pacbio/pancake/FastaSequenceId.h>
#include <pacbio/pancake/SeqDBIndexCache.h>
#include <pacbio/pancake/SeqDBReaderCachedBlock.h>
#include <memory>
#include <string>
#include <vector>

namespace PacBio {
namespace Pancake {

class SeqDBReaderCachedBlock
{
public:
    SeqDBReaderCachedBlock(std::shared_ptr<PacBio::Pancake::SeqDBIndexCache>& seedDBCache,
                           bool useHomopolymerCompression);
    ~SeqDBReaderCachedBlock();

    void LoadBlocks(const std::vector<int32_t>& blockIds);
    void LoadSequences(const std::vector<int32_t>& seqIds);
    void LoadSequences(const std::vector<std::string>& seqNames);

    const std::vector<FastaSequenceCached>& records() const { return recordStore_.records(); }
    const FastaSequenceCachedStore& recordStore() const { return recordStore_; }

    const FastaSequenceCached& GetSequence(int32_t seqId) const
    {
        return recordStore_.GetSequence(seqId);
    }

    void GetSequence(FastaSequenceCached& record, int32_t seqId)
    {
        recordStore_.GetSequence(record, seqId);
    }

    const FastaSequenceCached& GetSequence(const std::string& seqName) const
    {
        return recordStore_.GetSequence(seqName);
    }

    void GetSequence(FastaSequenceCached& record, const std::string& seqName)
    {
        recordStore_.GetSequence(record, seqName);
    }

private:
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBIndexCache_;
    bool useHomopolymerCompression_;
    std::vector<uint8_t> data_;
    FastaSequenceCachedStore recordStore_;

    void LoadBlockUncompressed_(const std::vector<ContiguousFilePart>& parts);
    void LoadBlockCompressed_(const std::vector<ContiguousFilePart>& parts);

    void CompressHomopolymers_();
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEQDB_READER_CACHED_BLOCK_H