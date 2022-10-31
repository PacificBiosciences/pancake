// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNMENT_REGION_HPP
#define PANCAKE_ALIGNMENT_REGION_HPP

#include <pbcopper/logging/Logging.h>
#include <boost/assign.hpp>
#include <boost/bimap.hpp>

#include <cstdint>
#include <ostream>
#include <string>
#include <string_view>
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

// clang-format off
static const boost::bimap<std::string_view, RegionType> RegionTypeToStringBimap =
    boost::assign::list_of<boost::bimap<std::string_view, RegionType>::relation>
        ("FRONT", Pancake::RegionType::FRONT)
        ("BACK", Pancake::RegionType::BACK)
        ("GLOBAL", Pancake::RegionType::GLOBAL);
// clang-format off

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

inline std::string RegionTypeToString(const RegionType type)
{
    const auto it = RegionTypeToStringBimap.right.find(type);
    if (it != RegionTypeToStringBimap.right.end()) {
        return std::string(it->second);
    }
    PBLOG_DEBUG << "Unknown RegionType given to the RegionTypeToString function.";
    return "UNKNOWN";
}

inline RegionType RegionTypeFromString(const std::string_view type)
{
    const auto it = RegionTypeToStringBimap.left.find(type);
    if (it != RegionTypeToStringBimap.left.end()) {
        return it->second;
    }
    throw std::runtime_error("Unknown RegionType: '" + std::string(type) +
                             "' in RegionTypeFromString.");
}

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
