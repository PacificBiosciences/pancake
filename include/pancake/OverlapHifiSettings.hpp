// Authors: Ivan Sovic

#ifndef PANCAKE_OVERLAP_HIFI_SETTINGS_H
#define PANCAKE_OVERLAP_HIFI_SETTINGS_H

#include <pancake/OverlapWriterFormat.hpp>

#include <pbcopper/cli2/CLI.h>

#include <cstdint>
#include <string>

namespace PacBio {
namespace Pancake {

struct OverlapHifiSettings
{
    struct Defaults
    {
        static const size_t NumThreads = 1;
        static const int32_t TargetBlockId = 0;
        static const int32_t QueryBlockStartId = 0;
        static const int32_t QueryBlockEndId = 0;

        static constexpr double FreqPercentile = 0.0002;
        static const int64_t MinQueryLen = 50;
        static const int64_t MinTargetLen = 50;
        static const int64_t MaxSeedDistance = 5000;
        static const int64_t MinNumSeeds = 3;
        static const int64_t MinCoveredBases = 30;
        static const int64_t MinChainSpan = 1000;
        static const int64_t ChainBandwidth = 100;
        static constexpr double AlignmentBandwidth = 0.01;
        static constexpr double AlignmentMaxD = 0.03;
        static constexpr double MinIdentity = 98.0;
        static const bool NoSNPsInIdentity = 0.0;
        static const bool NoIndelsInIdentity = 0.0;
        static const int64_t MinMappedLength = 1000;
        static const bool SkipSymmetricOverlaps = false;
        static const bool SkipSelfHits = false;
        static const bool OneHitPerTarget = false;
        static const bool SmartHitPerTarget = false;
        static const bool WriteReverseOverlaps = false;
        static const bool WriteIds = false;
        static const bool WriteCigar = false;
        static const int32_t AllowedDovetailDist = 0;
        static const int32_t AllowedHeuristicExtendDist = 0;
        static const int32_t CombineBlocks = 1;
        static const int32_t BestN = 0;
        static const bool UseHPC = false;
        static const bool UseTraceback = false;
        static const bool MaskHomopolymers = false;
        static const bool MaskSimpleRepeats = false;
        static const bool MaskHomopolymerSNPs = false;
        static const bool MaskHomopolymersArbitrary = false;
        static const OverlapWriterFormat OutFormat = OverlapWriterFormat::M4;
        static const bool MarkSecondary = false;
        static constexpr double SecondaryAllowedOverlapFraction = 0.50;
        static constexpr double SecondaryMinScoreFraction = 0.80;
        static const bool TrimAlignment = false;
        static const int32_t TrimWindowSize = 30;
        static constexpr double TrimWindowMatchFraction = 0.75;
        static const bool TrimToFirstMatch = false;
    };

    std::string TargetDBPrefix;
    std::string QueryDBPrefix;
    size_t NumThreads = Defaults::NumThreads;

    int32_t TargetBlockId = Defaults::TargetBlockId;
    int32_t QueryBlockStartId = Defaults::QueryBlockStartId;
    int32_t QueryBlockEndId = Defaults::QueryBlockEndId;

    double FreqPercentile = Defaults::FreqPercentile;
    int64_t MinQueryLen = Defaults::MinQueryLen;
    int64_t MinTargetLen = Defaults::MinTargetLen;
    int64_t MaxSeedDistance = Defaults::MaxSeedDistance;
    int64_t MinNumSeeds = Defaults::MinNumSeeds;
    int64_t MinCoveredBases = Defaults::MinCoveredBases;
    int64_t MinChainSpan = Defaults::MinChainSpan;
    int64_t ChainBandwidth = Defaults::ChainBandwidth;
    double AlignmentBandwidth = Defaults::AlignmentBandwidth;
    double AlignmentMaxD = Defaults::AlignmentMaxD;
    double MinIdentity = Defaults::MinIdentity;
    bool NoSNPsInIdentity = Defaults::NoSNPsInIdentity;
    bool NoIndelsInIdentity = Defaults::NoIndelsInIdentity;
    int64_t MinMappedLength = Defaults::MinMappedLength;
    bool SkipSymmetricOverlaps = Defaults::SkipSymmetricOverlaps;
    bool SkipSelfHits = Defaults::SkipSelfHits;
    bool OneHitPerTarget = Defaults::OneHitPerTarget;
    bool SmartHitPerTarget = Defaults::SmartHitPerTarget;
    bool WriteReverseOverlaps = Defaults::WriteReverseOverlaps;
    bool WriteIds = Defaults::WriteIds;
    bool WriteCigar = Defaults::WriteCigar;
    int32_t AllowedDovetailDist = Defaults::AllowedDovetailDist;
    int32_t AllowedHeuristicExtendDist = Defaults::AllowedHeuristicExtendDist;
    int32_t CombineBlocks = Defaults::CombineBlocks;
    int32_t BestN = Defaults::BestN;
    bool UseHPC = Defaults::UseHPC;
    bool UseTraceback = Defaults::UseTraceback;
    bool MaskHomopolymers = Defaults::MaskHomopolymers;
    bool MaskSimpleRepeats = Defaults::MaskSimpleRepeats;
    bool MaskHomopolymerSNPs = Defaults::MaskHomopolymerSNPs;
    bool MaskHomopolymersArbitrary = Defaults::MaskHomopolymersArbitrary;
    OverlapWriterFormat OutFormat = Defaults::OutFormat;
    bool MarkSecondary = Defaults::MarkSecondary;
    double SecondaryAllowedOverlapFraction = Defaults::SecondaryAllowedOverlapFraction;
    double SecondaryMinScoreFraction = Defaults::SecondaryMinScoreFraction;
    bool TrimAlignment = Defaults::TrimAlignment;
    int32_t TrimWindowSize = Defaults::TrimWindowSize;
    double TrimWindowMatchFraction = Defaults::TrimWindowMatchFraction;
    bool TrimToFirstMatch = Defaults::TrimToFirstMatch;

    OverlapHifiSettings();
    OverlapHifiSettings(const PacBio::CLI_v2::Results& options);
    static PacBio::CLI_v2::Interface CreateCLI();
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEQDB_SETTINGS_H
