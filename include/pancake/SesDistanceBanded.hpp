// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNMENT_SES_DISTANCE_BANDED_H
#define PANCAKE_ALIGNMENT_SES_DISTANCE_BANDED_H

#include <cstdint>

#include <iostream>
#include <limits>
#include <memory>
#include <sstream>

#include <pancake/SesResults.hpp>

namespace PacBio {
namespace Pancake {
namespace Alignment {

SesResults SESDistanceBanded(const char* query, size_t queryLen, const char* target,
                             size_t targetLen, int32_t maxDiffs, int32_t bandwidth);
}
}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ALIGNMENT_SES_DISTANCE_BANDED_H
