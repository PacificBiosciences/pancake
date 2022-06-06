// Authors: Robert Vaser

#ifndef PANCAKE_BURST_EXCISION_HPP
#define PANCAKE_BURST_EXCISION_HPP

#include <pancake/FastaSequenceCachedStore.hpp>
#include <pancake/MapperCLR.hpp>

#include <cstdint>
#include <utility>
#include <vector>

namespace PacBio {
namespace Pancake {

/// \brief A work in progress burst detection procedure.
///
/// \param zmw                           subreads of a single zmw (first one will be treated as the backbone for mapping)
/// \param backboneIdx                   subread index which will be used as the backbone
/// \param mapperSettings                MapperCLR settings
/// \param minBurstSize                  minimal distance difference of subread and backbone minimizer hits
/// \param maxHitDistanceRatio           maximal distance ratio of -||-
/// \param lowerInterPulseDistanceRatio  lower ratio of IPD median and per base values (x < (1 - param) * median)
/// \param upperInterPulseDistanceRatio  upper ratio of IPD median and per base values (x > (1 + param) * median)
/// \param lowerPulseWidthRatio          lower ratio of PW -||-
/// \param upperPulseWidthRatio          upper ratio of PW -||-
///
/// \return vector containing vectors of burst [begin, end] positions per subread
///
std::vector<std::vector<std::pair<int32_t, int32_t>>> BurstExcision(
    const FastaSequenceCachedStore& zmw, int32_t backboneIdx,
    const MapperCLRSettings& mapperSettings, int32_t minBurstSize = 32,
    double maxHitDistanceRatio = 0.9, double lowerInterPulseDistanceRatio = 0.6,
    double upperInterPulseDistanceRatio = 2.5, double lowerPulseWidthRatio = 0.45,
    double upperPulseWidthRatio = 0.85);

std::vector<std::vector<std::pair<int32_t, int32_t>>> BurstExcision(
    const FastaSequenceCachedStore& zmw, int32_t backboneIdx,
    const std::vector<MapperBaseResult>& mapperResult, int32_t minBurstSize = 32,
    double maxHitDistanceRatio = 0.9, double lowerInterPulseDistanceRatio = 0.6,
    double upperInterPulseDistanceRatio = 2.5, double lowerPulseWidthRatio = 0.45,
    double upperPulseWidthRatio = 0.85);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_BURST_EXCISION_HPP
