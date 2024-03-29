// Author: Ivan Sovic

#ifndef PANCAKE_UTIL_DEBUG_TOOLS_HPP
#define PANCAKE_UTIL_DEBUG_TOOLS_HPP

#include <pancake/util/TicToc.hpp>

#include <cstdint>
#include <string>
#include <unordered_map>

namespace PacBio {
namespace Pancake {

inline void LogTicToc([[maybe_unused]] const std::string& label, [[maybe_unused]] TicToc& tt,
                      [[maybe_unused]] std::unordered_map<std::string, double>& timings)
{
#ifdef PANCAKE_ENABLE_TIMINGS
    tt.Stop();
    timings[label] = tt.GetMicrosecs();
    tt.Start();
#endif
}
inline void LogTicTocAdd([[maybe_unused]] const std::string& label, [[maybe_unused]] TicToc& tt,
                         [[maybe_unused]] std::unordered_map<std::string, double>& timings)
{
#ifdef PANCAKE_ENABLE_TIMINGS
    tt.Stop();
    timings[label] += tt.GetMicrosecs();
    tt.Start();
#endif
}
inline void LogTicTocAdd([[maybe_unused]] const std::string& label,
                         [[maybe_unused]] const double timeToAdd,
                         [[maybe_unused]] std::unordered_map<std::string, double>& timings)
{
#ifdef PANCAKE_ENABLE_TIMINGS
    timings[label] += timeToAdd;
#endif
}

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_UTIL_DEBUG_TOOLS_HPP
