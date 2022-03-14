// Authors: Ivan Sovic

#ifndef PANCAKE_SEQDB_INFO_SETTINGS_H
#define PANCAKE_SEQDB_INFO_SETTINGS_H

#include <pancake/util/CommonTypes.hpp>
#include <pancake/util/GenomicUnit.hpp>

#include <pbcopper/cli2/CLI.h>

#include <cstdint>
#include <string>

namespace PacBio {
namespace Pancake {

// clang-format off
struct SeqDBInfoSettings
{
    struct Defaults
    {
        static const bool HumanReadableOutput = false;
        static const GenomicUnit Unit = GenomicUnit::bp;
    };

    std::string InputSeqDB;
    bool HumanReadableOutput = Defaults::HumanReadableOutput;
    GenomicUnit Unit = Defaults::Unit;

    SeqDBInfoSettings();
    SeqDBInfoSettings(const PacBio::CLI_v2::Results& options);
    static PacBio::CLI_v2::Interface CreateCLI();
};
// clang-format on

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEQDB_INFO_SETTINGS_H
