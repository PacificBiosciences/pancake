// Authors: Ivan Sovic

#include <pancake/OverlapWriterSAM.hpp>

namespace PacBio {
namespace Pancake {

OverlapWriterSAM::OverlapWriterSAM(std::FILE* fpOut, bool writeIds, bool writeCigar)
    : outFile_(""), fpOut_(fpOut), shouldClose_(false), writeIds_(writeIds), writeCigar_(writeCigar)
{}

OverlapWriterSAM::~OverlapWriterSAM()
{
    if (shouldClose_) {
        std::fclose(fpOut_);
    }
}

void OverlapWriterSAM::WriteHeader(const PacBio::Pancake::FastaSequenceCachedStore& targetSeqs)
{
    std::fprintf(fpOut_, "@HD\tVN:1.5\n");
    if (writeIds_) {
        char buff[100];
        for (const auto& targetSeq : targetSeqs.records()) {
            std::sprintf(buff, "%09d", targetSeq.Id());
            std::fprintf(fpOut_, "@SQ\tSN:%s\tLN:%d\n", buff,
                         static_cast<int32_t>(targetSeq.Size()));
        }
    } else {
        for (const auto& targetSeq : targetSeqs.records()) {
            std::fprintf(fpOut_, "@SQ\tSN:%s\tLN:%d\n", targetSeq.Name().c_str(),
                         static_cast<int32_t>(targetSeq.Size()));
        }
    }
}

void OverlapWriterSAM::Write(const Overlap& ovl,
                             const PacBio::Pancake::FastaSequenceCachedStore& targetSeqs,
                             const PacBio::Pancake::FastaSequenceCached& querySeq)
{
    // It's important to know if the overlap was flipped because then the query and target
    // names should be swapped.
    if (ovl.IsFlipped) {
        const auto& targetSeq = targetSeqs.GetSequence(ovl.Aid);
        // Don't look for the actual headers unless required. Saves the cost of a search.
        const auto& qName = writeIds_ ? "" : targetSeq.Name();
        const auto& tName = writeIds_ ? "" : querySeq.Name();
        PrintOverlapAsSAM(fpOut_, ovl, targetSeq.Bases(), targetSeq.Size(), qName, tName, writeIds_,
                          writeCigar_);

    } else {
        // Don't look for the actual headers unless required. Saves the cost of a search.
        const auto& qName = writeIds_ ? "" : querySeq.Name();
        const auto& tName = writeIds_ ? "" : targetSeqs.GetSequence(ovl.Bid).Name();
        PrintOverlapAsSAM(fpOut_, ovl, querySeq.Bases(), querySeq.Size(), qName, tName, writeIds_,
                          writeCigar_);
    }
}

}  // namespace Pancake
}  // namespace PacBio
