// Author: Ivan Sovic

#ifndef PANCAKE_UTIL_H
#define PANCAKE_UTIL_H

#include <pacbio/pancake/Lookups.h>
#include <cstdio>
#include <memory>
#include <sstream>
#include <vector>

namespace PacBio {
namespace Pancake {

struct FileDeleter
{
    void operator()(std::FILE* fp) const { fclose(fp); }
};

inline std::unique_ptr<FILE, FileDeleter> OpenFile(const std::string& filename,
                                                   const std::string& mode)
{
    auto ret = std::unique_ptr<FILE, FileDeleter>(fopen(filename.c_str(), mode.c_str()));
    if (ret == nullptr) {
        std::ostringstream errOss;
        errOss << "Could not open file '" << filename << "'!";
        throw std::runtime_error(errOss.str());
    }
    return ret;
}

/// \brief  Extracts the parent path and basename of a given path, on Unix platforms.
/// \note   If C++17 is available, this can be replaced by "std::filesystem::path::parent_path".
inline void SplitPathUnix(const std::string& path, std::string& parent, std::string& basename)
{
    // int32_t pathLen = path.size();
    size_t pos = path.size();
    size_t last = 0;
    while ((last = path.find_last_of('/', pos)) != std::string::npos) {
        if (last > 0 && path[last - 1] == '\\') {
            pos = last - 1;
            continue;
        }
        break;
    }
    int32_t posPrefixEnd = (last == std::string::npos) ? 0 : last;
    while (posPrefixEnd > 0 && path[posPrefixEnd - 1] == '.')
        --posPrefixEnd;
    parent = path.substr(0, posPrefixEnd);
    basename =
        (last == std::string::npos) ? path : (last < path.size()) ? path.substr(last + 1) : "";
    basename = (basename == ".") ? "" : basename;
}

/// \brief  Joins two paths.
inline std::string JoinPathUnix(std::string left, const std::string& right)
{
    if (left.empty() && right.empty()) return "";
    if (left.empty()) left = ".";
    return left + "/" + right;
}

/// \brief  Extracts the parent path and basename of a given path, on Windows.
/// \note   If C++17 is available, this can be replaced by "std::filesystem::path::parent_path".
inline void SplitPathWindows(const std::string& path, std::string& parent, std::string& basename)
{
    // Escaping paths in Windows is a bit more complicated than in Linux.
    // This function will only work on simple paths.
    size_t pos = path.size();
    size_t last = path.find_last_of('\\', pos);
    int32_t posPrefixEnd = (last == std::string::npos) ? 0 : last;
    parent = path.substr(0, posPrefixEnd);
    basename =
        (last == std::string::npos) ? path : (last < path.size()) ? path.substr(last + 1) : "";
}

/// \brief  Joins two paths.
inline std::string JoinPathWindows(std::string left, const std::string& right)
{
    if (left.empty() && right.empty()) return "";
    if (left.empty()) left = ".";
    return left + R"(\)" + right;
}

/// \brief  Filesystem tools implementation, currently a substitute for C++17 features.
///         When C++17 can readily be used, these functions can be replaced
///         with Path functions.
#if (defined(_WIN32) || defined(_WIN64))
#define SplitPath SplitPathWindows
#define JoinPath JoinPathWindows
#elif (defined(LINUX) || defined(__linux__) || defined(__unix__) || defined(__unix) || \
       defined(unix) || (defined(__APPLE__) && defined(__MACH__)))
#define SplitPath SplitPathUnix
#define JoinPath JoinPathUnix
#elif

#endif

inline std::string ReverseComplement(const std::string& seq, int64_t start, int64_t end)
{
    int64_t seqLen = seq.size();
    if (seqLen == 0) {
        return {};
    }
    if (start < 0 || end < 0 || start >= end || start > seqLen || end > seqLen) {
        std::ostringstream oss;
        oss << "Invalid start or end in a call to ReverseComplement. start = " << start
            << ", end = " << end << ", seqLen = " << seqLen << ".";
        throw std::runtime_error(oss.str());
    }
    int64_t span = end - start;
    std::string ret(span, '\0');
    int64_t maxI = span / 2;
    for (int64_t i = start, j = 0, k = span - 1; i < maxI; ++i, ++j, --k) {
        ret[j] = PacBio::Pancake::BaseToBaseComplement[static_cast<int32_t>(seq[k + start])];
        ret[k] = PacBio::Pancake::BaseToBaseComplement[static_cast<int32_t>(seq[i])];
    }
    if ((span & 0x01) != 0) {
        ret[maxI] = PacBio::Pancake::BaseToBaseComplement[static_cast<int32_t>(seq[maxI + start])];
    }
    return ret;
}

inline std::string ReverseComplement(const char* seq, int64_t seqLen, int64_t start, int64_t end)
{
    if (seqLen == 0) {
        return {};
    }
    if (start < 0 || end < 0 || start >= end || start > seqLen || end > seqLen) {
        std::ostringstream oss;
        oss << "Invalid start or end in a call to ReverseComplement. start = " << start
            << ", end = " << end << ", seqLen = " << seqLen << ".";
        throw std::runtime_error(oss.str());
    }
    int64_t span = end - start;
    std::string ret(span, '\0');
    int64_t maxI = span / 2;
    for (int64_t i = start, j = 0, k = span - 1; i < maxI; ++i, ++j, --k) {
        ret[j] = PacBio::Pancake::BaseToBaseComplement[static_cast<int32_t>(seq[k + start])];
        ret[k] = PacBio::Pancake::BaseToBaseComplement[static_cast<int32_t>(seq[i])];
    }
    if ((span & 0x01) != 0) {
        ret[maxI] = PacBio::Pancake::BaseToBaseComplement[static_cast<int32_t>(seq[maxI + start])];
    }
    return ret;
}

template <typename T>
std::vector<std::pair<T, T>> DistributeJobLoad(int32_t numThreads, T numJobs)
{
    if (numThreads == 0) {
        return {};
    }
    if (numJobs == 0) {
        return {};
    }
    const int32_t actualThreadCount = (numThreads < numJobs) ? numThreads : numJobs;
    const T minimumJobsPerThreads = (numJobs / actualThreadCount);
    const T remainingJobs = (numJobs % actualThreadCount);
    std::vector<T> jobsPerThread(actualThreadCount, minimumJobsPerThreads);
    for (T i = 0; i < remainingJobs; ++i) {
        ++jobsPerThread[i];
    }
    std::vector<std::pair<T, T>> ranges;
    T totalJobs = 0;
    for (const auto& span : jobsPerThread) {
        ranges.emplace_back(std::make_pair(totalJobs, totalJobs + span));
        totalJobs += span;
    }
    return ranges;
}

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_UTIL_H