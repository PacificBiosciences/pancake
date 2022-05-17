// Authors: Ivan Sovic

#ifndef PANCAKE_MAP_CLR_SETTINGS_H
#define PANCAKE_MAP_CLR_SETTINGS_H

#include <cstdint>
#include <string>

#include <pbcopper/cli2/CLI.h>
#include <pancake/AlignmentParameters.hpp>
#include <pancake/MapperCLR.hpp>

namespace PacBio {
namespace Pancake {

struct MapCLRSettings
{
    struct Defaults
    {
        static const size_t NumThreads = 1;

        // Input/output options.
        static const int32_t TargetBlockId = 0;
        static const int32_t QueryBlockStartId = 0;
        static const int32_t QueryBlockEndId = 0;
        static const OverlapWriterFormat OutFormat = OverlapWriterFormat::M4;
        static const bool WriteIds = false;
        static const bool WriteCigar = true;
        static const int32_t CombineBlocks = 1;

        // Seeding options.
        static const int32_t PrimaryKmerSize = 19;
        static const int32_t PrimaryMinimizerWindow = 5;
        static const int32_t PrimarySpacing = 0;
        static const bool PrimaryUseHomopolymerCompression = true;
        static const bool PrimaryUseReverseComplement = true;
        static const int32_t FallbackKmerSize = 19;
        static const int32_t FallbackMinimizerWindow = 5;
        static const int32_t FallbackSpacing = 0;
        static const bool FallbackUseHomopolymerCompression = true;
        static const bool FallbackUseReverseComplement = true;

        static constexpr double FreqPercentile = 0.0002;
        static const int64_t SeedOccurrenceMin = 0;
        static const int64_t SeedOccurrenceMax = 0;
        static const int64_t SeedOccurrenceMaxMemory = 0;
        static const int32_t ChainMaxSkip = 25;
        static const int32_t ChainMaxPredecessors = 500;
        static const int32_t ChainBandwidth = 500;
        static const int32_t SeedJoinDist = 10000;
        static const int32_t LongMergeBandwidth = 10000;
        static const int32_t MinNumSeeds = 3;
        static const int32_t MinCoveredBases = 0;
        static const int32_t MinDPScore = 40;
        static constexpr double SecondaryAllowedOverlapFractionQuery = 0.00;
        static constexpr double SecondaryAllowedOverlapFractionTarget = 0.50;
        static constexpr double SecondaryMinScoreFraction = 0.80;
        static const bool NoLIS = false;
        static const bool Align = false;
        static const int32_t MaxFlankExtensionDistance = 10000;
        static constexpr double FlankExtensionFactor = 1.3;
        static const int32_t MinAlignmentSpan = 200;

        static const bool RefineSeedHits = true;
        static const int32_t RefineMinGap1 = 10;
        static const int32_t RefineDiffThreshold = 40;
        static const int32_t RefineMinGap2 = 30;

        static const int32_t MinQueryLen = 50;
        static const int32_t BestNSecondary = 0;
        static const bool SkipSymmetricOverlaps = false;

        static constexpr AlignmentParameters AlnParams{};
    };

    std::string TargetDBPrefix;
    std::string QueryDBPrefix;
    size_t NumThreads = Defaults::NumThreads;

    int32_t TargetBlockId = Defaults::TargetBlockId;
    int32_t QueryBlockStartId = Defaults::QueryBlockStartId;
    int32_t QueryBlockEndId = Defaults::QueryBlockEndId;
    OverlapWriterFormat OutFormat = Defaults::OutFormat;

    MapperCLRSettings MapperSettings;

    bool WriteIds = Defaults::WriteIds;
    bool WriteCigar = Defaults::WriteCigar;
    int32_t CombineBlocks = Defaults::CombineBlocks;

    MapCLRSettings();
    MapCLRSettings(const PacBio::CLI_v2::Results& options);
    static PacBio::CLI_v2::Interface CreateCLI();
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_MAP_CLR_SETTINGS_H
