// Author: Ivan Sovic

#ifndef PANCAKE_MAPPER_CLR_SETTINGS_HPP
#define PANCAKE_MAPPER_CLR_SETTINGS_HPP

#include <pancake/AlignerFactory.hpp>
#include <pancake/AlignmentParameters.hpp>
#include <pancake/Range.hpp>
#include <pancake/SeedDBParameters.hpp>

#include <pbcopper/logging/Logging.h>

#include <boost/assign.hpp>
#include <boost/bimap.hpp>

#include <cstdint>
#include <memory>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace PacBio {
namespace Pancake {

enum class MapperSelfHitPolicy
{
    DEFAULT,
    SKIP,
    PERFECT_ALIGNMENT,
};

// clang-format off
static const boost::bimap<std::string, MapperSelfHitPolicy> MapperSelfHitPolicyBimap =
    boost::assign::list_of<boost::bimap<std::string, MapperSelfHitPolicy>::relation>
        ("DEFAULT", Pancake::MapperSelfHitPolicy::DEFAULT)
        ("SKIP", Pancake::MapperSelfHitPolicy::SKIP)
        ("PERFECT_ALIGNMENT", Pancake::MapperSelfHitPolicy::PERFECT_ALIGNMENT);

class MapperCLRMapSettings
{
public:
    // Indexing.
    PacBio::Pancake::SeedDBParameters seedParams{19, 10, 0, false, true, true};
    PacBio::Pancake::SeedDBParameters seedParamsFallback{19, 10, 0, false, true, true};

    // Seed hit occurrence filtering.
    // There are several parameters in play, which are combined according to the following formula:
    //      cutoff = max(seedOccurrenceMin, min(seedOccurrenceMax, seedOccurrenceMemMax, percentileCutoff))
    //  Here, "seedOccurrenceMemMax" is computed from seedOccurrenceMaxMemory and the seed hit histogram, so that
    //  the most abundant seeds are filtered out in order to maximally fill out the specified memory threshold.
    //  The "percentileCutoff" is computed from the top freqPercentile most abundant seeds.
    double freqPercentile = 0.0002;
    int64_t seedOccurrenceMin = 0;          // Minimum value for the occurrence threshold to keep a seed for mapping. If the frequency percentile is smaller than this, the threshold is pinned to this value.
    int64_t seedOccurrenceMax = 0;          // Maximum allowed occurrence of a seed to keep for mapping. Value <= 0 turns off this threshold.
    int64_t seedOccurrenceMaxMemory = 0;    // Maximum allowed memory to be consumed by collected seed hits. This is used to dynamically compute the maximum occurrence cutoff based on the seed hit histogram. Seeds are chosen in the sorted order by their occurrence until the memory threshold is reached. Value <= 0 turns off this threshold.

    // Mapping.
    int32_t chainMaxSkip = 25;
    int32_t chainMaxPredecessors = 500;
    int32_t chainBandwidth = 500;
    int32_t seedJoinDist = 10000;           // Maximum distance in either query or target coordinates to chain neighboring seeds.
    int32_t longMergeBandwidth = 10000;     // Maximum gap distance for chaining (abs(query_dist - target_dist)).
    int32_t minNumSeeds = 3;
    int32_t minCoveredBases = 0;
    int32_t minDPScore = 40;
    double secondaryAllowedOverlapFractionQuery = 0.00;
    double secondaryAllowedOverlapFractionTarget = 0.50;
    double secondaryMinScoreFraction = 0.80;
    bool useLIS = true;

    // Other.
    bool skipSymmetricOverlaps = false;
    MapperSelfHitPolicy selfHitPolicy = MapperSelfHitPolicy::DEFAULT;
    int32_t minQueryLen = 50;
    int32_t bestNSecondary = 0;             // Value < 0 will keep all secondary alignments. Value >= 0 will filter everything but the top N alignments.

    // Alignment region extraction parameters.
    int32_t maxFlankExtensionDist = 10000;  // Maximum length of the query/target sequences to consider when aligning flanks.
    double flankExtensionFactor = 1.3;      // Take this much more of the longer flanking sequence for alignment, to allow for indel errors.
    int32_t minAlignmentSpan = 200;         // If two seeds are closer than this, take the next seed. (Unless there are only 2 seeds left.)

    // Refining seed hits.
    bool refineSeedHits = true;             // Performs seed hit refinement before alignment regions are extracted.
    int32_t refineMinGap1 = 10;             // The minGap parameter of RefineChainedHits.
    int32_t refineDiffThreshold = 40;       // The diffThreshold of RefineChainedHits.
    int32_t refineMinGap2 = 30;             // The minGap parameter of RefineChainedHits2.

    // Reseeding long alignment regions with smaller seeds.
    bool reseedGaps = false;
    int32_t reseedGapMinLength = 500;
    int32_t reseedGapMaxLength = -1;            // Reseeding is applied only on gaps shorter than this value. Value <= 0 deactivates the upper limit.
    PacBio::Pancake::SeedDBParameters reseedSeedParams{10, 5, 0, false, true, true};
    double reseedFreqPercentile = 0.0002;     // Analogous to freqPercentile, but for local reseeded regions.
    int64_t reseedOccurrenceMin = 5;          // Analogous to seedOccurrenceMin, but for local reseeded regions.
    int64_t reseedOccurrenceMax = 100;        // Analogous to seedOccurrenceMax, but for local reseeded regions.
    int64_t reseedOccurrenceMaxMemory = 100'000'000;    // Analogous to seedOccurrenceMaxMemory, but for local reseeded regions.

    // bool oneHitPerTarget = false;
    // int32_t maxSeedDistance = 5000;
    // int32_t minMappedLength = 1000;
    // double minIdentity = 75.0;
};

class MapperCLRAlignSettings
{
public:
    // Alignment.
    bool align = true;
    MapperSelfHitPolicy selfHitPolicy = MapperSelfHitPolicy::DEFAULT;
    AlignerType alignerTypeGlobal = AlignerType::KSW2;
    AlignmentParameters alnParamsGlobal;
    AlignerType alignerTypeExt = AlignerType::KSW2;
    AlignmentParameters alnParamsExt;
};

class MapperCLRSettings
{
public:
    MapperCLRMapSettings map;
    MapperCLRAlignSettings align;
};
// clang-format on

inline std::ostream& operator<<(std::ostream& out, const MapperCLRMapSettings& a)
{
    out << "chainMaxSkip = " << a.chainMaxSkip << "\n"
        << "chainMaxPredecessors = " << a.chainMaxPredecessors << "\n"
        << "chainBandwidth = " << a.chainBandwidth << "\n"
        << "seedJoinDist = " << a.seedJoinDist << "\n"
        << "longMergeBandwidth = " << a.longMergeBandwidth << "\n"
        << "minNumSeeds = " << a.minNumSeeds << "\n"
        << "minCoveredBases = " << a.minCoveredBases << "\n"
        << "minDPScore = " << a.minDPScore << "\n"
        << "useLIS = " << a.useLIS << "\n"

        << "secondaryAllowedOverlapFractionQuery = " << a.secondaryAllowedOverlapFractionQuery
        << "\n"
        << "secondaryAllowedOverlapFractionTarget = " << a.secondaryAllowedOverlapFractionTarget
        << "\n"
        << "secondaryMinScoreFraction = " << a.secondaryMinScoreFraction << "\n"

        << "seedParams.KmerSize = " << a.seedParams.KmerSize << "\n"
        << "seedParams.MinimizerWindow = " << a.seedParams.MinimizerWindow << "\n"
        << "seedParams.Spacing = " << a.seedParams.Spacing << "\n"
        << "seedParams.UseHPC = " << a.seedParams.UseHPC << "\n"
        << "seedParams.UseHPCForSeedsOnly = " << a.seedParams.UseHPCForSeedsOnly << "\n"
        << "seedParams.UseRC = " << a.seedParams.UseRC << "\n"

        << "seedParamsFallback.KmerSize = " << a.seedParamsFallback.KmerSize << "\n"
        << "seedParamsFallback.MinimizerWindow = " << a.seedParamsFallback.MinimizerWindow << "\n"
        << "seedParamsFallback.Spacing = " << a.seedParamsFallback.Spacing << "\n"
        << "seedParamsFallback.UseHPC = " << a.seedParamsFallback.UseHPC << "\n"
        << "seedParamsFallback.UseHPCForSeedsOnly = " << a.seedParamsFallback.UseHPCForSeedsOnly
        << "\n"
        << "seedParamsFallback.UseRC = " << a.seedParamsFallback.UseRC << "\n"

        << "freqPercentile = " << a.freqPercentile << "\n"

        << "maxFlankExtensionDist = " << a.maxFlankExtensionDist << "\n"
        << "flankExtensionFactor = " << a.flankExtensionFactor << "\n"
        << "refineSeedHits = " << a.refineSeedHits << "\n"
        << "minAlignmentSpan = " << a.minAlignmentSpan << "\n"

        << "refineSeedHits = " << a.refineSeedHits << "\n"
        << "refineMinGap1 = " << a.refineMinGap1 << "\n"
        << "refineDiffThreshold = " << a.refineDiffThreshold << "\n"
        << "refineMinGap2 = " << a.refineMinGap2 << "\n"

        << "reseedGaps = " << a.reseedGaps << "\n"
        << "reseedGapMinLength = " << a.reseedGapMinLength << "\n"
        << "reseedGapMaxLength = " << a.reseedGapMaxLength << "\n"
        << "reseedFreqPercentile = " << a.reseedFreqPercentile << "\n"
        << "reseedOccurrenceMin = " << a.reseedOccurrenceMin << "\n"
        << "reseedOccurrenceMax = " << a.reseedOccurrenceMax << "\n"
        << "reseedOccurrenceMaxMemory = " << a.reseedOccurrenceMaxMemory << "\n"
        << "reseedSeedParams.KmerSize = " << a.reseedSeedParams.KmerSize << "\n"
        << "reseedSeedParams.MinimizerWindow = " << a.reseedSeedParams.MinimizerWindow << "\n"
        << "reseedSeedParams.Spacing = " << a.reseedSeedParams.Spacing << "\n"
        << "reseedSeedParams.UseHPC = " << a.reseedSeedParams.UseHPC << "\n"
        << "reseedSeedParams.UseHPCForSeedsOnly = " << a.reseedSeedParams.UseHPCForSeedsOnly << "\n"
        << "reseedSeedParams.UseRC = " << a.reseedSeedParams.UseRC << "\n"

        << "skipSymmetricOverlaps = " << a.skipSymmetricOverlaps << "\n"
        << "minQueryLen = " << a.minQueryLen << "\n"
        << "bestNSecondary = " << a.bestNSecondary << "\n";

    return out;
}

inline std::ostream& operator<<(std::ostream& out, const MapperCLRAlignSettings& a)
{
    out << "align = " << a.align << "\n"
        << "alignerTypeGlobal = " << AlignerTypeToString(a.alignerTypeGlobal) << "\n"
        << "alnParamsGlobal:\n"
        << a.alnParamsGlobal << "alignerTypeExt = " << AlignerTypeToString(a.alignerTypeExt) << "\n"
        << "alnParamsExt:\n"
        << a.alnParamsExt << "\n";
    return out;
}

inline std::ostream& operator<<(std::ostream& out, const MapperCLRSettings& a)
{
    out << a.map << a.align;
    return out;
}

inline std::string MapperSelfHitPolicyToString(const MapperSelfHitPolicy type)
{
    const auto it = MapperSelfHitPolicyBimap.right.find(type);
    if (it != MapperSelfHitPolicyBimap.right.end()) {
        return it->second;
    }
    PBLOG_DEBUG << "Unknown MapperSelfHitPolicy given to the MapperSelfHitPolicyToString function.";
    return "UNKNOWN";
}

inline MapperSelfHitPolicy MapperSelfHitPolicyFromString(const std::string& type)
{
    const auto it = MapperSelfHitPolicyBimap.left.find(type);
    if (it != MapperSelfHitPolicyBimap.left.end()) {
        return it->second;
    }
    throw std::runtime_error("Unknown MapperSelfHitPolicy: '" + type +
                             "' in MapperSelfHitPolicyFromString.");
}

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_MAPPER_CLR_SETTINGS_HPP
