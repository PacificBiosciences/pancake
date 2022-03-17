// Author: Ivan Sovic

#ifndef PANCAKE_UTIL_CONVERSION_HPP
#define PANCAKE_UTIL_CONVERSION_HPP

#include <cstdint>
#include <cstdlib>
#include <string>

namespace PacBio {
namespace Pancake {

/// \brief Parses an int from a string, and returns true if the entire string was parsed.
///         Otherwise, returns false.
inline bool ConvertStringToInt(const std::string& inVal, int32_t& outVal)
{
    char* ptr = NULL;
    outVal = std::strtol(inVal.c_str(), &ptr, 10);
    if (*ptr) {
        return false;
    }
    return true;
}

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_UTIL_CONVERSION_HPP
