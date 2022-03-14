// Author: Ivan Sovic

#ifndef PANCAKE_SEQDB_WRITER_BASE_H
#define PANCAKE_SEQDB_WRITER_BASE_H

#include <pancake/CompressedSequence.hpp>

#include <ostream>
#include <string>

namespace PacBio {
namespace Pancake {

/// \brief      Base class that can be specialized for writing a compressed or an
///             uncompressed set of sequences.
///
class SeqDBWriterBase
{
public:
    virtual ~SeqDBWriterBase() {}
    virtual void AddSequence(const std::string& header, const std::string& seq) = 0;
    virtual bool WriteSequences() = 0;
    virtual void WriteIndex() = 0;
    virtual void ClearSequenceBuffer() = 0;
    virtual void FlushSequenceBuffer() = 0;
    virtual void CloseFiles() = 0;
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEQDB_WRITER_BASE_H