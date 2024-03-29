// Author: Ivan Sovic

#ifndef PANCAKE_OVERLAP_WRITER_PAF_HPP
#define PANCAKE_OVERLAP_WRITER_PAF_HPP

#include <pancake/FastaSequenceId.hpp>
#include <pancake/Overlap.hpp>
#include <pancake/OverlapWriterBase.hpp>
#include <pancake/SeqDBReaderCached.hpp>
#include <pancake/SeqDBReaderCachedBlock.hpp>
#include <pancake/util/Util.hpp>

#include <cstdint>
#include <cstdio>
#include <memory>
#include <unordered_map>
#include <vector>

namespace PacBio {
namespace Pancake {

class OverlapWriterPAF : public OverlapWriterBase
{
public:
    OverlapWriterPAF(std::FILE* fpOut, bool writeIds, bool writeCigar);

    ~OverlapWriterPAF();

    void Write(const Overlap& ovl, const PacBio::Pancake::FastaSequenceCachedStore& targetSeqs,
               const PacBio::Pancake::FastaSequenceCached& querySeq);

    void WriteHeader(const PacBio::Pancake::FastaSequenceCachedStore& targetSeqs);

private:
    std::string outFile_;
    std::FILE* fpOut_ = NULL;
    bool shouldClose_ = false;
    bool writeIds_ = false;
    bool writeCigar_ = false;
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_OVERLAP_WRITER_PAF_HPP
