// Author: Ivan Sovic

#ifndef PANCAKE_SES_DISTANCE_BANDED_HPP
#define PANCAKE_SES_DISTANCE_BANDED_HPP

#include <pancake/SesResults.hpp>

#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>

namespace PacBio {
namespace Pancake {
namespace Alignment {

SesResults SESDistanceBanded(const char* query, size_t queryLen, const char* target,
                             size_t targetLen, int32_t maxDiffs, int32_t bandwidth);
}
}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SES_DISTANCE_BANDED_HPP
