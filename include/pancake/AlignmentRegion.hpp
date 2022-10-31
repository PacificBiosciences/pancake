// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNMENT_REGION_HPP
#define PANCAKE_ALIGNMENT_REGION_HPP

#include <cstdint>
#include <ostream>
#include <tuple>
#include <vector>

namespace PacBio {
namespace Pancake {

enum class RegionType
{
    FRONT,
    BACK,
    GLOBAL,
};

inline std::string RegionTypeToString(const RegionType type)
{
    if (type == RegionType::FRONT) {
        return "FRONT";
    } else if (type == RegionType::BACK) {
        return "BACK";
    } else if (type == RegionType::GLOBAL) {
        return "GLOBAL";
    }
    return "UNKNOWN";
}

class AlignmentRegion
{
public:
    int32_t qStart = 0;
    int32_t qSpan = 0;
    int32_t tStart = 0;
    int32_t tSpan = 0;
    bool queryRev = false;
    RegionType type = RegionType::GLOBAL;
    int32_t regionId = -1;
    int32_t maxGap = 0;

    bool operator==(const AlignmentRegion& b) const
    {
        return std::tie(qStart, qSpan, tStart, tSpan, queryRev, type, regionId, maxGap) ==
               std::tie(b.qStart, b.qSpan, b.tStart, b.tSpan, b.queryRev, b.type, b.regionId,
                        b.maxGap);
    }
};
inline std::ostream& operator<<(std::ostream& os, const AlignmentRegion& b)
{
    os << "qStart = " << b.qStart << ", qSpan = " << b.qSpan << ", tStart = " << b.tStart
       << ", tSpan = " << b.tSpan << ", queryRev = " << (b.queryRev ? "true" : "false")
       << ", type = " << RegionTypeToString(b.type) << ", regionId = " << b.regionId;
    return os;
}

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ALIGNMENT_REGION_HPP
