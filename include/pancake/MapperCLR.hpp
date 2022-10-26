// Author: Ivan Sovic

#ifndef PANCAKE_MAPPER_CLR_HPP
#define PANCAKE_MAPPER_CLR_HPP

#include <pancake/AlignerBase.hpp>
#include <pancake/AlignerFactory.hpp>
#include <pancake/AlignmentParameters.hpp>
#include <pancake/DPChain.hpp>
#include <pancake/FastaSequenceCached.hpp>
#include <pancake/FastaSequenceId.hpp>
#include <pancake/MapperBase.hpp>
#include <pancake/MapperCLRSettings.hpp>
#include <pancake/Overlap.hpp>
#include <pancake/OverlapHifiSettings.hpp>
#include <pancake/Range.hpp>
#include <pancake/Seed.hpp>
#include <pancake/SeedDBIndexCache.hpp>
#include <pancake/SeedIndex.hpp>
#include <pancake/SeqDBReaderCachedBlock.hpp>
#include <pancake/SequenceSeedsCached.hpp>
#include <pancake/util/CommonTypes.hpp>

#include <cstdint>
#include <memory>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace PacBio {
namespace Pancake {

class MapperCLR : public MapperBase
{
public:
    using MapperBase::MapAndAlign;

    MapperCLR(const MapperCLRSettings& settings);
    ~MapperCLR() override;

    /**
     * \brief Utility function which allows to update the mapping/alignment settings _after_ the MapperCLR
     * object was already constructed.
    */
    void UpdateSettings(const MapperCLRSettings& settings);

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
                                            const SeedIndex& index,
                                            const FastaSequenceCached& querySeq,
                                            const SequenceSeedsCached& querySeeds, int32_t queryId,
                                            int64_t freqCutoff) override;

    /*
     * \brief Maps the query sequence to the targets, where targets are provided by the SeedIndex.
    */
    MapperBaseResult Map(const FastaSequenceCachedStore& targetSeqs, const SeedIndex& index,
                         const FastaSequenceCached& querySeq, const SequenceSeedsCached& querySeeds,
                         int32_t queryLen, int32_t queryId, int64_t freqCutoff);

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
    std::shared_ptr<ChainingScratchSpace> ssChain_;
    std::vector<SeedHit> ssSeedHits_;

    /*
     * \brief Wraps the entire mapping and alignment process.
     * It calls Map_ and then Align_ if necessary, based on the settings.
     * Also, it filters the final alignments (which is given, since that's part
     * of the mapping process).
    */
    static MapperBaseResult WrapMapAndAlign_(
        const FastaSequenceCachedStore& targetSeqs, const SeedIndex& index,
        const FastaSequenceCached& querySeq, const SequenceSeedsCached& querySeeds, int32_t queryId,
        int64_t freqCutoff, const MapperCLRSettings& settings,
        std::shared_ptr<ChainingScratchSpace> ssChain, std::vector<SeedHit>& ssSeedHits,
        AlignerBasePtr& alignerGlobal, AlignerBasePtr& alignerExt);

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
        const MapperCLRSettings& settings, std::shared_ptr<ChainingScratchSpace> ssChain,
        std::vector<SeedHit>& ssSeedHits, AlignerBasePtr& alignerGlobal,
        AlignerBasePtr& alignerExt);

    /*
     * \brief Maps the query sequence to the targets, where targets are provided by the SeedIndex.
    */
    static MapperBaseResult Map_(const FastaSequenceCachedStore& targetSeqs,
                                 const PacBio::Pancake::SeedIndex& index,
                                 const FastaSequenceCached& querySeq,
                                 const SequenceSeedsCached& querySeeds, int32_t queryLen,
                                 int32_t queryId, const MapperCLRSettings& settings,
                                 int64_t freqCutoff, std::shared_ptr<ChainingScratchSpace> ssChain,
                                 std::vector<SeedHit>& ssSeedHits);

    /**
     * @brief This is the core of the Map_ function. It takes a vector of collected seed hits and
     *          performs various filterings and refinement to produce one or more chains of seed hits
     *          which represent the final mappings. This function forms the mapping result which includes
     *          the mapping locations, seed hits, timings, and other information.
     *          NOTE: This function does NOT check the settings.selfHitPolicy. Instead it depends on:
     *              1. Seed hits to have been pre-filtered from self-hits before calling this function, and
     *              2. The addPerfectMapping to be true if a mocked mapping is supposed to be added.
     *
     * @param hits Input seed hits.
     * @param ssChain Scratch space, to reuse the memory for chaining. Nullptr may be provided, and a new scratch space will be allocated internally.
     * @param targetSeqs Target sequences used for mapping.
     * @param queryLen Length of the query sequence.
     * @param queryId ID of the query sequence.
     * @param settings All settings used for mapping.
     * @param addPerfectMapping If true, a mocked perfect mapping will be added to the vector of chained regions. Mocked seed hits
     *                          (consisting of only 2 hits: first and last) will be created.
     *                          NOTE: This function does not check the settings.selfHitPolicy. Instead it uses this parameter.
     * @param timing Map which stores the timing information of various mapping stages that may have occurred before the call to this function. These
     *               will be copied and stored in the returning MapperBaseResult object.
     * @param debugStepId For debug purposes, this is a counter of different mapping steps, used to log information and store debug data. This is a reference
     *                    to allow multiple runs of this function.
     * @return MapperBaseResult All mapping results for this query.
     */
    static MapperBaseResult HitsToMappings_(
        const std::vector<SeedHit>& hits, const std::shared_ptr<ChainingScratchSpace>& ssChain,
        const FastaSequenceCachedStore& targetSeqs, int32_t queryLen, int32_t queryId,
        const MapperCLRMapSettings& settings, bool addPerfectMapping,
        int32_t maxAllowedDistForBadEndRefinement,
        const std::unordered_map<std::string, double>& timing, int32_t& debugStepId);

    static std::vector<SeedHit> ReseedAlignmentRegions_(
        const MapperBaseResult& result, const FastaSequenceCachedStore& targetSeqs,
        const FastaSequenceCached& querySeq, int32_t reseedGapMinLength, int32_t reseedGapMaxLength,
        const PacBio::Pancake::SeedDBParameters& seedParams, double reseedFreqPercentile,
        int64_t reseedOccurrenceMin, int64_t reseedOccurrenceMax,
        int64_t reseedOccurrenceMaxMemory);

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
                                                            int32_t queryId, bool skipSelfHits,
                                                            bool skipSymmetricOverlaps);

    /*
     * \brief Performs LIS and DP chaining, then constructs the overlaps from those resulting chains.
    */
    static std::vector<std::unique_ptr<ChainedRegion>> ChainAndMakeOverlap_(
        const FastaSequenceCachedStore& targetSeqs, const std::vector<SeedHit>& hits,
        const std::vector<PacBio::Pancake::Range>& hitGroups, int32_t queryId, int32_t queryLen,
        int32_t chainMaxSkip, int32_t chainMaxPredecessors, int32_t seedJoinDist,
        int32_t chainBandwidth, int32_t minNumSeeds, int32_t minCoveredBases, int32_t minDPScore,
        bool useLIS, std::shared_ptr<ChainingScratchSpace> ssChain,
        std::unordered_map<std::string, double>& retTimings);

    /*
     * \brief Takes previously chained regions, collects all remaining seed hits from those regions into
     * a single pile, then runs DP chaining on all of those seed hits.
     * A new set of chained regions is created from the chaining results.
    */
    static std::vector<std::unique_ptr<ChainedRegion>> ReChainSeedHits_(
        const std::vector<std::unique_ptr<ChainedRegion>>& chainedRegions,
        const FastaSequenceCachedStore& targetSeqs, int32_t queryId, int32_t queryLen,
        int32_t chainMaxSkip, int32_t chainMaxPredecessors, int32_t seedJoinDist,
        int32_t chainBandwidth, int32_t minNumSeeds, int32_t minCoveredBases, int32_t minDPScore,
        std::shared_ptr<ChainingScratchSpace> ssChain,
        std::unordered_map<std::string, double>& retTimings);

    /*
     * \brief Merges the neighboring chains if they are not overlaping in neither the query nor
     * target coordinates.
     * The right of the two merged chained regions is set to nullptr.
    */
    static void LongMergeChains_(std::vector<std::unique_ptr<ChainedRegion>>& chainedRegions,
                                 int32_t maxBandwidth);

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
    static std::unique_ptr<ChainedRegion> CreateMockedMapping_(int32_t queryId, int32_t queryLen,
                                                               int32_t targetId, int32_t targetLen);
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

OverlapPtr CreateMockedAlignment(const OverlapPtr& ovl, int32_t matchScore);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_MAPPER_CLR_HPP
