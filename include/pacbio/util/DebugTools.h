// Author: Ivan Sovic

#ifndef PANCAKE_DEBUG_TOOLS_H
#define PANCAKE_DEBUG_TOOLS_H

#include <pacbio/util/TicToc.h>
#include <cstdint>
#include <string>
#include <unordered_map>

namespace PacBio {
namespace Pancake {

#ifdef DEBUG_LOG_TIMINGS
inline void LogTicToc(const std::string& label, TicToc& tt,
                      std::unordered_map<std::string, double>& timings)
{
    tt.Stop();
    timings[label] = tt.GetMicrosecs();
    tt.Start();
}
inline void LogTicTocAdd(const std::string& label, TicToc& tt,
                         std::unordered_map<std::string, double>& timings)
{
    tt.Stop();
    timings[label] += tt.GetMicrosecs();
    tt.Start();
}
inline void LogTicTocAdd(const std::string& label, const double timeToAdd,
                         std::unordered_map<std::string, double>& timings)
{
    timings[label] += timeToAdd;
}
#else
inline void LogTicToc(const std::string& /*label*/, TicToc& /*tt*/,
                      std::unordered_map<std::string, double>& /*timings*/)
{
}
inline void LogTicTocAdd(const std::string& /*label*/, TicToc& /*tt*/,
                         std::unordered_map<std::string, double>& /*timings*/)
{
}
inline void LogTicTocAdd(const std::string& /*label*/, const double /*timeToAdd*/,
                         std::unordered_map<std::string, double>& /*timings*/)
{
}
#endif

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_COMMON_H
