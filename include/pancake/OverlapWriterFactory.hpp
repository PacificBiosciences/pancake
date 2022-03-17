// Author: Ivan Sovic

#ifndef PANCAKE_OVERLAP_WRITER_FACTORY_HPP
#define PANCAKE_OVERLAP_WRITER_FACTORY_HPP

#include <pancake/OverlapWriterBase.hpp>
#include <pancake/OverlapWriterFormat.hpp>
#include <pancake/OverlapWriterIPAOvl.hpp>

#include <memory>

namespace PacBio {
namespace Pancake {

std::unique_ptr<OverlapWriterBase> OverlapWriterFactory(OverlapWriterFormat writerType, FILE* fpOut,
                                                        bool writeIds, bool writeCigar);
}
}  // namespace PacBio

#endif  // PANCAKE_OVERLAP_WRITER_FACTORY_HPP
