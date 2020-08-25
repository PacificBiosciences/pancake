// Authors: Ivan Sovic

#include <pacbio/overlaphifi/OverlapWriterFactory.h>
#include <pacbio/overlaphifi/OverlapWriterM4.h>
#include <pacbio/overlaphifi/OverlapWriterPAF.h>
#include <pacbio/overlaphifi/OverlapWriterSAM.h>

namespace PacBio {
namespace Pancake {

std::unique_ptr<OverlapWriterBase> OverlapWriterFactory(OverlapWriterFormat writerType, FILE* fpOut,
                                                        bool writeIds, bool writeCigar)
{
    if (writerType == OverlapWriterFormat::IPAOvl) {
        return std::unique_ptr<OverlapWriterBase>(
            new OverlapWriterIPAOvl(fpOut, writeIds, writeCigar));
    } else if (writerType == OverlapWriterFormat::M4) {
        return std::unique_ptr<OverlapWriterBase>(new OverlapWriterM4(fpOut, writeIds, writeCigar));
    } else if (writerType == OverlapWriterFormat::PAF) {
        return std::unique_ptr<OverlapWriterBase>(
            new OverlapWriterPAF(fpOut, writeIds, writeCigar));
    } else if (writerType == OverlapWriterFormat::SAM) {
        return std::unique_ptr<OverlapWriterBase>(
            new OverlapWriterSAM(fpOut, writeIds, writeCigar));
    }
    throw std::runtime_error("Unsupported output format!");
    return nullptr;
}

}  // namespace Pancake
}  // namespace PacBio
