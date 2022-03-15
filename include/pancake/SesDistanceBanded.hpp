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

SesResults SESDistanceBanded(std::string_view query, std::string_view target, int32_t maxDiffs,
                             int32_t bandwidth);
}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SES_DISTANCE_BANDED_HPP
