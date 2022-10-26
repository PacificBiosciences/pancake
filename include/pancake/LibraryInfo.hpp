#ifndef PANCAKE_LIBRARYINFO_HPP
#define PANCAKE_LIBRARYINFO_HPP

#include <pbcopper/library/Bundle.h>
#include <pbcopper/library/Info.h>

namespace PacBio {
namespace Pancake {

Library::Info LibraryInfo();

Library::Bundle LibraryBundle();

std::string PancakeFormattedVersion();

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_LIBRARYINFO_HPP
