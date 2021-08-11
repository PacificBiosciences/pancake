// Author: Ivan Sovic

#ifndef PANCAKE_UTIL_GENOMIC_UNIT_H
#define PANCAKE_UTIL_GENOMIC_UNIT_H

#include <pacbio/pancake/Lookups.h>
#include <cstdio>
#include <memory>
#include <sstream>
#include <vector>

namespace PacBio {
namespace Pancake {

enum class GenomicUnit
{
    Unknown,
    bp,
    kbp,
    Mbp,
    Gbp,
};

inline GenomicUnit GenomicUnitFromString(const std::string& val)
{
    if (val == "bp") {
        return GenomicUnit::bp;
    } else if (val == "kbp") {
        return GenomicUnit::kbp;
    } else if (val == "Mbp") {
        return GenomicUnit::Mbp;
    } else if (val == "Gbp") {
        return GenomicUnit::Gbp;
    }
    return GenomicUnit::Unknown;
}

inline std::string GenomicUnitToString(GenomicUnit unit)
{
    switch (unit) {
        case GenomicUnit::bp:
            return "bp";
        case GenomicUnit::kbp:
            return "kbp";
        case GenomicUnit::Mbp:
            return "Mbp";
        case GenomicUnit::Gbp:
            return "Gbp";
        case GenomicUnit::Unknown:
            return "Unknown";
        default:
            return "Asdfqwe";
    }
}

inline double GenomicUnitToDouble(GenomicUnit unit)
{
    switch (unit) {
        case GenomicUnit::bp:
            return 1.0;
        case GenomicUnit::kbp:
            return 1'000.0;
        case GenomicUnit::Mbp:
            return 1'000'000.0;
        case GenomicUnit::Gbp:
            return 1'000'000'000.0;
        case GenomicUnit::Unknown:
        default:
            return 0.0;
    }
}

class GenomicUnitFromTo
{
public:
    GenomicUnitFromTo(GenomicUnit fromUnit, GenomicUnit toUnit)
        : conversionFactor_(GenomicUnitToDouble(fromUnit) / GenomicUnitToDouble(toUnit))
    {
    }
    double Convert(double fromValue) { return fromValue * conversionFactor_; }
    double conversionFactor() const { return conversionFactor_; }

private:  // until we update ScaleLengthByFactor() elsewhere
    double conversionFactor_;
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_UTIL_GENOMIC_UNIT_H
