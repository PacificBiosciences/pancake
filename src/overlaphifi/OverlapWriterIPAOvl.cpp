// Authors: Ivan Sovic

#include <pacbio/overlaphifi/OverlapWriterIPAOvl.h>

namespace PacBio {
namespace Pancake {

OverlapWriterIPAOvl::OverlapWriterIPAOvl(FILE* fpOut, bool writeIds, bool writeCigar)
    : outFile_(""), fpOut_(fpOut), shouldClose_(false), writeIds_(writeIds), writeCigar_(writeCigar)
{
}

OverlapWriterIPAOvl::~OverlapWriterIPAOvl()
{
    if (shouldClose_) {
        fclose(fpOut_);
    }
}

void OverlapWriterIPAOvl::Write(const OverlapPtr& ovl,
                                const PacBio::Pancake::SeqDBReaderCached& targetSeqs,
                                const PacBio::Pancake::FastaSequenceId& querySeq, bool isFlipped)
{
    // It's important to know if the overlap was flipped because then the query and target
    // names should be swapped.
    if (isFlipped) {
        // Don't look for the actual headers unless required. Saves the cost of a search.
        const auto& qName = writeIds_ ? "" : targetSeqs.GetSequence(ovl->Aid).Name();
        const auto& tName = writeIds_ ? "" : querySeq.Name();
        PrintOverlapAsIPAOvl(fpOut_, ovl, qName, tName, writeIds_, writeCigar_);
    } else {
        // Don't look for the actual headers unless required. Saves the cost of a search.
        const auto& qName = writeIds_ ? "" : querySeq.Name();
        const auto& tName = writeIds_ ? "" : targetSeqs.GetSequence(ovl->Bid).Name();
        PrintOverlapAsIPAOvl(fpOut_, ovl, qName, tName, writeIds_, writeCigar_);
    }
}

void OverlapWriterIPAOvl::Write(const OverlapPtr& ovl,
                                const PacBio::Pancake::SeqDBReaderCachedBlock& targetSeqs,
                                const PacBio::Pancake::FastaSequenceCached& querySeq,
                                bool isFlipped)
{
    // It's important to know if the overlap was flipped because then the query and target
    // names should be swapped.
    if (isFlipped) {
        // Don't look for the actual headers unless required. Saves the cost of a search.
        const auto& qName = writeIds_ ? "" : targetSeqs.GetSequence(ovl->Aid).Name();
        const auto& tName = writeIds_ ? "" : querySeq.Name();
        PrintOverlapAsIPAOvl(fpOut_, ovl, qName, tName, writeIds_, writeCigar_);

    } else {
        // Don't look for the actual headers unless required. Saves the cost of a search.
        const auto& qName = writeIds_ ? "" : querySeq.Name();
        const auto& tName = writeIds_ ? "" : targetSeqs.GetSequence(ovl->Bid).Name();
        PrintOverlapAsIPAOvl(fpOut_, ovl, qName, tName, writeIds_, writeCigar_);
    }
}

}  // namespace Pancake
}  // namespace PacBio
