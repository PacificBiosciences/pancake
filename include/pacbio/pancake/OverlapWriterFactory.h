// Author: Ivan Sovic

#ifndef PANCAKE_OVERLAPHIFI_OVERLAP_WRITER_FACTORY_H
#define PANCAKE_OVERLAPHIFI_OVERLAP_WRITER_FACTORY_H

#include <pacbio/pancake/OverlapWriterBase.h>
#include <pacbio/pancake/OverlapWriterFormat.h>
#include <pacbio/pancake/OverlapWriterIPAOvl.h>
#include <memory>

namespace PacBio {
namespace Pancake {

std::unique_ptr<OverlapWriterBase> OverlapWriterFactory(OverlapWriterFormat writerType, FILE* fpOut,
                                                        bool writeIds, bool writeCigar);
}
}

#endif  // PANCAKE_OVERLAPHIFI_OVERLAP_WRITER_FACTORY_H
