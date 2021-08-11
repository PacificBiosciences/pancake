// Authors: Ivan Sovic

#ifndef PANCAKE_SEEDDB_SETTINGS_H
#define PANCAKE_SEEDDB_SETTINGS_H

#include <cstdint>
#include <string>

#include <pacbio/pancake/SeedDBParameters.h>
#include <pbcopper/cli2/CLI.h>

namespace PacBio {
namespace Pancake {
namespace SeedDB {

struct SeedDBSettings
{
    struct Defaults
    {
        static const size_t NumThreads = 1;
        static const bool SplitBlocks = false;
        static const int32_t KmerSize = 28;
        static const int32_t MinimizerWindow = 80;
        static const int32_t Spacing = 0;
        static const bool UseHPC = false;
        static const bool UseHPCForSeedsOnly = false;
        static const bool NoRevCmp = false;
    };

    std::string InputFile;
    std::string OutputPrefix;
    size_t NumThreads = Defaults::NumThreads;
    bool SplitBlocks = Defaults::SplitBlocks;
    SeedDBParameters SeedParameters{
        Defaults::KmerSize, Defaults::MinimizerWindow,    Defaults::Spacing,
        Defaults::UseHPC,   Defaults::UseHPCForSeedsOnly, !Defaults::NoRevCmp};

    SeedDBSettings();
    SeedDBSettings(const PacBio::CLI_v2::Results& options);
    static PacBio::CLI_v2::Interface CreateCLI();
};

}  // namespace SeedDB
}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEQDB_SETTINGS_H