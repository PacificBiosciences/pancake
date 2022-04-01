// Author: Ivan Sovic

#ifndef PANCAKE_OVERLAP_WRITER_FORMAT_HPP
#define PANCAKE_OVERLAP_WRITER_FORMAT_HPP

#include <string_view>
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

inline OverlapWriterFormat ParseOverlapWriterFormat(const std::string_view val)
{
    if (val == "ipa") {
        return OverlapWriterFormat::IPAOvl;
    } else if (val == "m4") {
        return OverlapWriterFormat::M4;
    } else if (val == "paf") {
        return OverlapWriterFormat::PAF;
    } else if (val == "sam") {
        return OverlapWriterFormat::SAM;
    }
    return OverlapWriterFormat::Unknown;
}

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_OVERLAP_WRITER_FORMAT_HPP
