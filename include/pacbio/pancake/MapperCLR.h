// Author: Ivan Sovic

#ifndef PANCAKE_MAPPER_CLR_H
#define PANCAKE_MAPPER_CLR_H

#include <lib/intervaltree/IntervalTree.h>
#include <pacbio/overlaphifi/OverlapHifiSettings.h>
#include <pacbio/pancake/AlignerBase.h>
#include <pacbio/pancake/AlignerFactory.h>
#include <pacbio/pancake/AlignmentParameters.h>
#include <pacbio/pancake/DPChain.h>
#include <pacbio/pancake/FastaSequenceCached.h>
#include <pacbio/pancake/FastaSequenceId.h>
#include <pacbio/pancake/MapperBase.h>
#include <pacbio/pancake/Overlap.h>
#include <pacbio/pancake/Range.h>
#include <pacbio/pancake/Seed.h>
#include <pacbio/pancake/SeedDBIndexCache.h>
#include <pacbio/pancake/SeedIndex.h>
#include <pacbio/pancake/SeqDBReaderCachedBlock.h>
#include <pacbio/pancake/SequenceSeedsCached.h>
#include <pacbio/util/CommonTypes.h>
#include <cstdint>
#include <memory>
#include <unordered_map>
#include <vector>

namespace PacBio {
namespace Pancake {

using IntervalTreeInt32 =
    interval_tree::IntervalTree<int32_t,
                                int32_t>;  // First: interval scalar type, Second: value type.
using IntervalVectorInt32 = IntervalTreeInt32::interval_vector;
using IntervalInt32 = IntervalTreeInt32::interval;

enum class MapperSelfHitPolicy
{
    DEFAULT,
    SKIP,
    PERFECT_ALIGNMENT,
};

// clang-format off
class MapperCLRMapSettings
{
public:
    // Indexing.
    PacBio::Pancake::SeedDB::SeedDBParameters seedParams{19, 10, 0, false, true, true};
    PacBio::Pancake::SeedDB::SeedDBParameters seedParamsFallback{19, 10, 0, false, true, true};
    double freqPercentile = 0.0002;

    // Mapping.
    int32_t chainMaxSkip = 25;
    int32_t chainMaxPredecessors = 500;
    int32_t chainBandwidth = 500;
    int32_t maxGap = 10000;
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
    int32_t bestNSecondary = 0;

    // Alignment region extraction parameters.
    int32_t maxFlankExtensionDist = 10000;  // Maximum length of the query/target sequences to consider when aligning flanks.
    double flankExtensionFactor = 1.3;      // Take this much more of the longer flanking sequence for alignment, to allow for indel errors.
    int32_t minAlignmentSpan = 200;         // If two seeds are closer than this, take the next seed. (Unless there are only 2 seeds left.)

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
        << "maxGap = " << a.maxGap << "\n"
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
        << "minAlignmentSpan = " << a.minAlignmentSpan << "\n"

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

inline std::string MapperSelfHitPolicyToString(MapperSelfHitPolicy mp)
{
    switch (mp) {
        case MapperSelfHitPolicy::DEFAULT:
            return "DEFAULT";
            break;
        case MapperSelfHitPolicy::SKIP:
            return "SKIP";
            break;
        case MapperSelfHitPolicy::PERFECT_ALIGNMENT:
            return "PERFECT_ALIGNMENT";
            break;
        default:
            return "UNKNOWN";
    }
    return "UNKNOWN";
}

class MapperCLR : public MapperBase
{
public:
    MapperCLR(const MapperCLRSettings& settings);
    ~MapperCLR() override;

    /*
     * \brief Runs mapping and alignment of one or more query sequences to one or more target sequences.
     * This is the basic interface for the most simple usage.
     * This interface does not require the SeedIndex or minimizers, because it will compute them internally.
     * The std::string objects are first converted to FastaSequenceCached, and then the MapAndAlign overload
     * is called on these inputs.
    */
    std::vector<MapperBaseResult> MapAndAlign(const std::vector<std::string>& targetSeqs,
                                              const std::vector<std::string>& querySeqs) override;

    /*
     * \brief Runs mapping and alignment of one or more query sequences to one or more target sequences.
     * This interface does not require the SeedIndex or minimizers, because it will compute them internally.
     * This function constructs the FastaSequenceCachedStore from the given FastaSequenceCached objects,
     * and calls the related overload of MapAndAlign with these inputs.
    */
    std::vector<MapperBaseResult> MapAndAlign(
        const std::vector<FastaSequenceCached>& targetSeqs,
        const std::vector<FastaSequenceCached>& querySeqs) override;

    /*
     * \brief Runs mapping and alignment of one or more query sequences to one or more target sequences.
     * This interface does not require the SeedIndex or minimizers, because it will compute them internally.
     * The WrapBuildIndexMapAndAlignWithFallback_ function is called for processing.
     * If there were no alignmentsads produced for a query using the primary seeding parameters,
     * another call to WrapMapAndAlign_ is performed with the seedParamsFallback options.
     * The FastaSequenceCachedStore objects are used for random access to sequences.
    */
    std::vector<MapperBaseResult> MapAndAlign(const FastaSequenceCachedStore& targetSeqs,
                                              const FastaSequenceCachedStore& querySeqs) override;

    /*
     * \brief Runs mapping and alignment of a single query sequence to one or more target sequences.
     * This interface is more suitable for integration into workflow, where target index is provided
     * from the outside, and may be used by other workflow steps too.
     * This function calls WrapMapAndAlign_ to perform mapping and alignment.
     * There is no seed fallback implemented in this function, since the SeedIndex is precomputed
     * and provided from the outside.
    */
    MapperBaseResult MapAndAlignSingleQuery(const FastaSequenceCachedStore& targetSeqs,
                                            const PacBio::Pancake::SeedIndex& index,
                                            const FastaSequenceCached& querySeq,
                                            const std::vector<PacBio::Pancake::Int128t>& querySeeds,
                                            const int32_t queryId, int64_t freqCutoff) override;

    /*
     * \brief Maps the query sequence to the targets, where targets are provided by the SeedIndex.
    */
    MapperBaseResult Map(const FastaSequenceCachedStore& targetSeqs,
                         const PacBio::Pancake::SeedIndex& index,
                         const std::vector<PacBio::Pancake::Int128t>& querySeeds,
                         const int32_t queryLen, const int32_t queryId, int64_t freqCutoff) const;

    /*
     * \brief Aligns a precomputed mapping result.
    */
    MapperBaseResult Align(const FastaSequenceCachedStore& targetSeqs,
                           const FastaSequenceCached& querySeq,
                           const MapperBaseResult& mappingResult);

private:
    MapperCLRSettings settings_;
    AlignerBasePtr alignerGlobal_;
    AlignerBasePtr alignerExt_;

    /*
     * \brief Wraps the entire mapping and alignment process.
     * It calls Map_ and then Align_ if necessary, based on the settings.
     * Also, it filters the final alignments (which is given, since that's part
     * of the mapping process).
    */
    static MapperBaseResult WrapMapAndAlign_(
        const FastaSequenceCachedStore& targetSeqs, const PacBio::Pancake::SeedIndex& index,
        const FastaSequenceCached& querySeq,
        const std::vector<PacBio::Pancake::Int128t>& querySeeds, const int32_t queryId,
        int64_t freqCutoff, const MapperCLRSettings& settings, AlignerBasePtr& alignerGlobal,
        AlignerBasePtr& alignerExt);

    /*
     * This function starts from plain sequences, and constructs the seeds (minimizers),
     * builds the SeedIndex, and runs mapping and alignment using WrapMapAndAlign_.
     * Each query is processed one by one, linearly.
     * If a query did not map, another attempt will be tried with the seedParamsFallback
     * options.
     * If the seedParamsFallback == seedParams, then the second attempt will not be
     * run (nor will the fallback index be generated).
    */
    static std::vector<MapperBaseResult> WrapBuildIndexMapAndAlignWithFallback_(
        const FastaSequenceCachedStore& targetSeqs, const FastaSequenceCachedStore& querySeqs,
        const MapperCLRSettings& settings, AlignerBasePtr& alignerGlobal,
        AlignerBasePtr& alignerExt);

    /*
     * \brief Maps the query sequence to the targets, where targets are provided by the SeedIndex.
    */
    static MapperBaseResult Map_(const FastaSequenceCachedStore& targetSeqs,
                                 const PacBio::Pancake::SeedIndex& index,
                                 const std::vector<PacBio::Pancake::Int128t>& querySeeds,
                                 const int32_t queryLen, const int32_t queryId,
                                 const MapperCLRSettings& settings, int64_t freqCutoff);

    /*
     * \brief Aligns the query sequence to one or more target sequences, based on the mappings and
     * seed hits collected and refined during the mapping process (the Map_ function).
     * This function cannot be const because the member alignerGlobal_ and alignerExt_ can be
     * modified (they have internal memory which gets reused with alignment).
    */
    static MapperBaseResult Align_(const FastaSequenceCachedStore& targetSeqs,
                                   const FastaSequenceCached& querySeq,
                                   const MapperBaseResult& mappingResult,
                                   const MapperCLRSettings& settings, AlignerBasePtr& alignerGlobal,
                                   AlignerBasePtr& alignerExt);

    /*
     * \brief Filters symmetric and self seed hits, based on sequence IDs.
     * Be careful when using this function in case when the query and target DBs are not the same.
    */
    static std::vector<SeedHit> FilterSymmetricAndSelfHits_(const std::vector<SeedHit>& hits,
                                                            const int32_t queryId,
                                                            const bool skipSelfHits,
                                                            const bool skipSymmetricOverlaps);

    /*
     * \brief Performs LIS and DP chaining, then constructs the overlaps from those resulting chains.
    */
    static std::vector<std::unique_ptr<ChainedRegion>> ChainAndMakeOverlap_(
        const FastaSequenceCachedStore& targetSeqs, const std::vector<SeedHit>& hits,
        const std::vector<PacBio::Pancake::Range>& hitGroups, int32_t queryId, int32_t queryLen,
        int32_t chainMaxSkip, int32_t chainMaxPredecessors, int32_t maxGap, int32_t chainBandwidth,
        int32_t minNumSeeds, int32_t minCoveredBases, int32_t minDPScore, bool useLIS);

    /*
     * \brief Takes previously chained regions, collects all remaining seed hits from those regions into
     * a single pile, then runs DP chaining on all of those seed hits.
     * A new set of chained regions is created from the chaining results.
    */
    static std::vector<std::unique_ptr<ChainedRegion>> ReChainSeedHits_(
        const std::vector<std::unique_ptr<ChainedRegion>>& chainedRegions,
        const FastaSequenceCachedStore& targetSeqs, int32_t queryId, int32_t queryLen,
        int32_t chainMaxSkip, int32_t chainMaxPredecessors, int32_t maxGap, int32_t chainBandwidth,
        int32_t minNumSeeds, int32_t minCoveredBases, int32_t minDPScore);

    /*
     * \brief Merges the neighboring chains if they are not overlaping in neither the query nor
     * target coordinates.
     * The right of the two merged chained regions is set to nullptr.
    */
    static void LongMergeChains_(std::vector<std::unique_ptr<ChainedRegion>>& chainedRegions,
                                 int32_t maxGap);

    static std::vector<AlignmentRegion> CollectAlignmentRegions_(const ChainedRegion& singleMapping,
                                                                 int32_t minAlignmentSpan,
                                                                 int32_t maxFlankExtensionDist,
                                                                 double flankExtensionFactor);

    /*
     * \brif A helper function which creates a mocked self-mapping based on the query and target
     * IDs and lengths. The result is an unaligned mapping which spans full length of the
     * query and target.
     * Throws if the query and target lengths are different.
     * This function also does not initialize the Atype and Btype labels (it sets them to Unknown).
    */
    static std::unique_ptr<ChainedRegion> CreateMockedMapping_(const int32_t queryId,
                                                               const int32_t queryLen,
                                                               const int32_t targetId,
                                                               const int32_t targetLen);
};

/*
    * \brief Wraps the labelling of secondary and supplementary chained regions.
*/
void WrapFlagSecondaryAndSupplementary(
    std::vector<std::unique_ptr<ChainedRegion>>& allChainedRegions,
    double secondaryAllowedOverlapFractionQuery, double secondaryAllowedOverlapFractionTarget,
    double secondaryMinScoreFraction);

int32_t CondenseMappings(std::vector<std::unique_ptr<ChainedRegion>>& mappings,
                         int32_t bestNSecondary);

OverlapPtr CreateMockedAlignment(const OverlapPtr& ovl, const int32_t matchScore);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_MAPPER_CLR_H
