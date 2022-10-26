#include <pancake/LibraryInfo.hpp>

#include "LibraryGitHash.hpp"
#include "LibraryVersion.hpp"

namespace PacBio {
namespace Pancake {

Library::Bundle LibraryBundle()
{
    Library::Bundle bundle{LibraryInfo(), {}};
    return bundle;
}

Library::Info LibraryInfo()
{
    return {"pancake", std::string(RELEASE_VERSION), std::string(LIBRARY_GIT_SHA1)};
}

std::string PancakeFormattedVersion()
{
    return LibraryInfo().Release + " (commit " + LibraryInfo().GitSha1 + ')';
}

}  // namespace Pancake
}  // namespace PacBio
