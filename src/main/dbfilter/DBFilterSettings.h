// Authors: Ivan Sovic

#ifndef PANCAKE_DBFILTER_SETTINGS_H
#define PANCAKE_DBFILTER_SETTINGS_H

#include <pancake/util/CommonTypes.hpp>

#include <pbcopper/cli2/CLI.h>

#include <cstdint>
#include <string>

namespace PacBio {
namespace Pancake {

SamplingType ParseSamplingType(const std::string& val);

struct DBFilterSettings
{
    struct Defaults
    {
        static const SamplingType Sampling = SamplingType::None;
        static const int64_t SampleBases = 0;
        static constexpr float BlockSize = 1000.0f;
        static const int64_t RandomSeed = -1;
        static const bool Consolidate = false;
        static const int32_t CompressionLevel = 1;
        static constexpr float BufferSize = 1000.0f;
        static const bool SplitBlocks = false;
    };

    std::string InputPrefix;
    std::string OutputPrefix;
    SamplingType Sampling = Defaults::Sampling;
    int64_t SampleBases = Defaults::SampleBases;
    float BlockSize = Defaults::BlockSize;
    int64_t RandomSeed = Defaults::RandomSeed;
    std::string FilterListPath;
    FilterListType FilterType;

    // Consolidation related options.
    bool Consolidate = Defaults::Consolidate;
    int32_t CompressionLevel = Defaults::CompressionLevel;
    float BufferSize = Defaults::BufferSize;
    bool SplitBlocks = Defaults::SplitBlocks;

    DBFilterSettings();
    DBFilterSettings(const PacBio::CLI_v2::Results& options);
    static PacBio::CLI_v2::Interface CreateCLI();
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_DBFILTER_SETTINGS_H