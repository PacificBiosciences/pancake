// Author: Ivan Sovic

#ifndef PANCAKE_OVERLAP_WRITER_M4_HPP
#define PANCAKE_OVERLAP_WRITER_M4_HPP

#include <pancake/FastaSequenceId.hpp>
#include <pancake/Overlap.hpp>
#include <pancake/OverlapWriterBase.hpp>
#include <pancake/SeqDBReaderCached.hpp>
#include <pancake/SeqDBReaderCachedBlock.hpp>
#include <pancake/util/Util.hpp>

#include <cstdint>
#include <memory>
#include <unordered_map>
#include <vector>

namespace PacBio {
namespace Pancake {

class OverlapWriterM4 : public OverlapWriterBase
{
public:
    OverlapWriterM4(FILE* fpOut, bool writeIds, bool writeCigar);
    ~OverlapWriterM4();

    void WriteHeader(const PacBio::Pancake::SeqDBReaderCached& targetSeqs) override;

    void WriteHeader(const PacBio::Pancake::SeqDBReaderCachedBlock& targetSeqs) override;

    void Write(const Overlap& ovl, const PacBio::Pancake::SeqDBReaderCached& targetSeqs,
               const PacBio::Pancake::FastaSequenceId& querySeq) override;

    void Write(const Overlap& ovl, const PacBio::Pancake::SeqDBReaderCachedBlock& targetSeqs,
               const PacBio::Pancake::FastaSequenceCached& querySeq) override;

private:
    std::string outFile_;
    FILE* fpOut_ = NULL;
    bool shouldClose_ = false;
    bool writeIds_ = false;
    bool writeCigar_ = false;
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_OVERLAP_WRITER_M4_HPP
