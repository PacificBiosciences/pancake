// Authors: Ivan Sovic

#ifndef PANCAKE_SEQDB_DUMP_SETTINGS_H
#define PANCAKE_SEQDB_DUMP_SETTINGS_H

#include <pacbio/util/CommonTypes.h>
#include <cstdint>
#include <string>

#include <pbcopper/cli2/CLI.h>

namespace PacBio {
namespace Pancake {

// clang-format off
struct SeqDBDumpSettings
{
    struct Defaults
    {
        static const int32_t BlockId = -1;    // -1 to write all blocks, a valid block ID otherwise.
        static const bool WriteIds = false;
        static const bool UseHPC = false;
    };

    std::string InputSeqDB;
    std::string OutputFile;
    int32_t BlockId = Defaults::BlockId;
    bool WriteIds = Defaults::WriteIds;
    bool UseHPC = Defaults::UseHPC;

    SeqDBDumpSettings();
    SeqDBDumpSettings(const PacBio::CLI_v2::Results& options);
    static PacBio::CLI_v2::Interface CreateCLI();
};
// clang-format on

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEQDB_DUMP_SETTINGS_H
