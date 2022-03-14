// Authors: Ivan Sovic

#ifndef PANCAKE_SEQDB_SETTINGS_H
#define PANCAKE_SEQDB_SETTINGS_H

#include <pbcopper/cli2/CLI.h>

#include <cstdint>
#include <string>

namespace PacBio {
namespace Pancake {

struct SeqDBSettings
{
    struct Defaults
    {
        static const size_t NumThreads = 1;
        static const int32_t CompressionLevel = 1;
        static constexpr float BufferSize = 1000.0f;
        static constexpr float BlockSize = 1000.0f;
        static const bool SplitBlocks = false;
    };

    std::string OutputPrefix;
    std::vector<std::string> InputFiles;
    size_t NumThreads = Defaults::NumThreads;
    int32_t CompressionLevel = Defaults::CompressionLevel;
    float BufferSize = Defaults::BufferSize;
    float BlockSize = Defaults::BlockSize;
    bool SplitBlocks = Defaults::SplitBlocks;

    SeqDBSettings();
    SeqDBSettings(const PacBio::CLI_v2::Results& options);
    static PacBio::CLI_v2::Interface CreateCLI();
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEQDB_SETTINGS_H