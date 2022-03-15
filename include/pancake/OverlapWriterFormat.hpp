// Author: Ivan Sovic

#ifndef PANCAKE_OVERLAPHIFI_OVERLAP_WRITER_FORMAT_H
#define PANCAKE_OVERLAPHIFI_OVERLAP_WRITER_FORMAT_H

namespace PacBio {
namespace Pancake {

enum class OverlapWriterFormat
{
    IPAOvl,
    M4,
    PAF,
    SAM,
    Unknown
};
}
}  // namespace PacBio

#endif  // PANCAKE_OVERLAPHIFI_OVERLAP_WRITER_FORMAT_H
