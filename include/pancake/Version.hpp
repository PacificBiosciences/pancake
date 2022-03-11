// Author: Armin Töpfer

#ifndef PANCAKE_VERSION_H
#define PANCAKE_VERSION_H

#include <string>
#include <tuple>

namespace PacBio {
namespace Pancake {

std::string PancakeGitSha1();
std::string PancakeVersion();
std::string PancakeFormattedVersion();
std::tuple<int, int, int> PancakeVersionTriple();

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_VERSION_H
