// Authors: Ivan Sovic

#include <pancake/OverlapWriterM4.hpp>

namespace PacBio {
namespace Pancake {

OverlapWriterM4::OverlapWriterM4(std::FILE* fpOut, bool writeIds, bool writeCigar)
    : outFile_(""), fpOut_(fpOut), shouldClose_(false), writeIds_(writeIds), writeCigar_(writeCigar)
{}

OverlapWriterM4::~OverlapWriterM4()
{
    if (shouldClose_) {
        std::fclose(fpOut_);
    }
}

void OverlapWriterM4::WriteHeader(const PacBio::Pancake::FastaSequenceCachedStore& /*targetSeqs*/)
{
    // This format doesn't have a header.
}

void OverlapWriterM4::Write(const Overlap& ovl,
                            const PacBio::Pancake::FastaSequenceCachedStore& targetSeqs,
                            const PacBio::Pancake::FastaSequenceCached& querySeq)
{
    // It's important to know if the overlap was flipped because then the query and target
    // names should be swapped.
    if (ovl.IsFlipped) {
        // Don't look for the actual headers unless required. Saves the cost of a search.
        const auto& qName = writeIds_ ? "" : targetSeqs.GetSequence(ovl.Aid).Name();
        const auto& tName = writeIds_ ? "" : querySeq.Name();
        PrintOverlapAsM4(fpOut_, ovl, qName, tName, writeIds_, writeCigar_);

    } else {
        // Don't look for the actual headers unless required. Saves the cost of a search.
        const auto& qName = writeIds_ ? "" : querySeq.Name();
        const auto& tName = writeIds_ ? "" : targetSeqs.GetSequence(ovl.Bid).Name();
        PrintOverlapAsM4(fpOut_, ovl, qName, tName, writeIds_, writeCigar_);
    }
}

}  // namespace Pancake
}  // namespace PacBio
