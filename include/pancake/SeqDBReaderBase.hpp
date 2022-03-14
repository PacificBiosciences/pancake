// Author: Ivan Sovic

#ifndef PANCAKE_SEQDB_READER_BASE_H
#define PANCAKE_SEQDB_READER_BASE_H

#include <pancake/FastaSequenceId.hpp>
#include <pancake/SeqDBIndexCache.hpp>

#include <string>
#include <vector>

namespace PacBio {
namespace Pancake {

/// \brief      Base class that can be specialized for reading a compressed or an
///             uncompressed set of sequences.
///
class SeqDBReaderBase
{
public:
    virtual ~SeqDBReaderBase() {}
    virtual bool GetSequence(Pancake::FastaSequenceId& record, int64_t seqId) = 0;
    virtual bool GetSequence(Pancake::FastaSequenceId& record, const std::string& seqName) = 0;
    virtual bool GetNext(Pancake::FastaSequenceId& record) = 0;
    virtual bool GetNextBatch(std::vector<Pancake::FastaSequenceId>& records,
                              int64_t batchSize) = 0;
    virtual bool GetBlock(std::vector<Pancake::FastaSequenceId>& records, int32_t blockId) = 0;
    virtual bool JumpTo(int64_t seqId) = 0;
    virtual bool JumpTo(const std::string& seqName) = 0;
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEQDB_READER_COMPRESSED_H
