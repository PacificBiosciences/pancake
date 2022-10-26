// Authors: Ivan Sovic

#include <pancake/MapperCLR.hpp>

#include <pancake/AlignmentSeeded.hpp>
#include <pancake/AlignmentTools.hpp>
#include <pancake/DPChain.hpp>
#include <pancake/DPChainSimd.hpp>
#include <pancake/DiffCounts.hpp>
#include <pancake/MapperUtility.hpp>
#include <pancake/Minimizers.hpp>
#include <pancake/OverlapWriterBase.hpp>
#include <pancake/Secondary.hpp>
#include <pancake/SeedHitWriter.hpp>
#include <pancake/SesDistanceBanded.hpp>
#include <pancake/third-party/istl/lis.hpp>
#include <pancake/util/DebugTools.hpp>
#include <pancake/util/RunLengthEncoding.hpp>
#include <pancake/util/TicToc.hpp>
#include <pancake/util/Util.hpp>

#include <pancake/third-party/kxsort/kxsort.h>
#include <pancake/third-party/pdqsort/pdqsort.h>
#include <pbcopper/logging/Logging.h>
#include <pbcopper/third-party/edlib.h>
#include <pbcopper/utility/Ssize.h>

#include <algorithm>
#include <array>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <tuple>

// #define PANCAKE_MAP_CLR_DEBUG
// #define PANCAKE_MAP_CLR_DEBUG_2
// #define PANCAKE_MAP_CLR_DEBUG_ALIGN
// #define PANCAKE_MAP_CLR_DEBUG_WRITE_SEED_HITS_TO_FILE
constexpr bool DEBUG_VERBOSE_CHAINS = false;

#if defined(PANCAKE_MAP_CLR_DEBUG) || defined(PANCAKE_MAP_CLR_DEBUG_2)
#include <pbcopper/utility/MemoryConsumption.h>
#endif

#ifdef PANCAKE_USE_SSE41
#define ChainHits ChainHitsSimd
#else
#define ChainHits ChainHitsSisd
#endif

namespace PacBio {
namespace Pancake {

MapperCLR::MapperCLR(const MapperCLRSettings& settings)
    : settings_{settings}
    , alignerGlobal_(nullptr)
    , alignerExt_(nullptr)
    , ssChain_(std::make_shared<ChainingScratchSpace>())
{
    alignerGlobal_ =
        AlignerFactory(settings.align.alignerTypeGlobal, settings.align.alnParamsGlobal);
    alignerExt_ = AlignerFactory(settings.align.alignerTypeExt, settings.align.alnParamsExt);
}

MapperCLR::~MapperCLR() = default;

void MapperCLR::UpdateSettings(const MapperCLRSettings& settings)
{
    settings_ = settings;
    alignerGlobal_ =
        AlignerFactory(settings.align.alignerTypeGlobal, settings.align.alnParamsGlobal);
    alignerExt_ = AlignerFactory(settings.align.alignerTypeExt, settings.align.alnParamsExt);
}

std::vector<MapperBaseResult> MapperCLR::MapAndAlign(const FastaSequenceCachedStore& targetSeqs,
                                                     const FastaSequenceCachedStore& querySeqs)
{
    try {
        return WrapBuildIndexMapAndAlignWithFallback_(targetSeqs, querySeqs, settings_, ssChain_,
                                                      ssSeedHits_, alignerGlobal_, alignerExt_);
    } catch (const std::exception& e) {
        // Log, but do not fail. Important for clients of this class.
        // Return a vector of the size of the input, but with empty results for each query.
        PBLOG_DEBUG << "MapperCLR generated an exception in " << std::string(__FUNCTION__)
                    << ". Message: " << e.what();
        return std::vector<MapperBaseResult>(querySeqs.Size());
    }
}

MapperBaseResult MapperCLR::MapAndAlignSingleQuery(const FastaSequenceCachedStore& targetSeqs,
                                                   const PacBio::Pancake::SeedIndex& index,
                                                   const FastaSequenceCached& querySeq,
                                                   const SequenceSeedsCached& querySeeds,
                                                   const int32_t queryId, int64_t freqCutoff)
{
    try {
        return WrapMapAndAlign_(targetSeqs, index, querySeq, querySeeds, queryId, freqCutoff,
                                settings_, ssChain_, ssSeedHits_, alignerGlobal_, alignerExt_);
    } catch (const std::exception& e) {
        // Log, but do not fail. Important for clients of this class.
        // Return a vector of the size of the input, but with empty results for each query.
        PBLOG_DEBUG << "MapperCLR generated an exception in " << std::string(__FUNCTION__)
                    << ". Message: " << e.what();
        return {};
    }
}

MapperBaseResult MapperCLR::Map(const FastaSequenceCachedStore& targetSeqs,
                                const PacBio::Pancake::SeedIndex& index,
                                const FastaSequenceCached& querySeq,
                                const SequenceSeedsCached& querySeeds, const int32_t queryLen,
                                const int32_t queryId, int64_t freqCutoff)
{
    try {
        return Map_(targetSeqs, index, querySeq, querySeeds, queryLen, queryId, settings_,
                    freqCutoff, ssChain_, ssSeedHits_);
    } catch (const std::exception& e) {
        // Log, but do not fail. Important for clients of this class.
        // Return a vector of the size of the input, but with empty results for each query.
        PBLOG_DEBUG << "MapperCLR generated an exception in " << std::string(__FUNCTION__)
                    << ". Message: " << e.what();
        return {};
    }
}

MapperBaseResult MapperCLR::Align(const FastaSequenceCachedStore& targetSeqs,
                                  const FastaSequenceCached& querySeq,
                                  const MapperBaseResult& mappingResult)
{
    try {
        return Align_(targetSeqs, querySeq, mappingResult, settings_, alignerGlobal_, alignerExt_);
    } catch (const std::exception& e) {
        // Log, but do not fail. Important for clients of this class.
        // Return a vector of the size of the input, but with empty results for each query.
        PBLOG_DEBUG << "MapperCLR generated an exception in " << std::string(__FUNCTION__)
                    << ". Message: " << e.what();
        return {};
    }
}

void DebugWriteSeedHits([[maybe_unused]] const std::vector<SeedHit>& hits,
                        [[maybe_unused]] const int32_t debugStepId,
                        [[maybe_unused]] const std::string& descriptor,
                        [[maybe_unused]] int32_t queryId, [[maybe_unused]] int32_t queryLen)
{
#ifdef PANCAKE_MAP_CLR_DEBUG_WRITE_SEED_HITS_TO_FILE
    std::ostringstream ossFileName;
    ossFileName << "temp-debug/hits-q" << std::to_string(queryId) << "-" << std::setfill('0')
                << std::setw(3) << debugStepId << "-" << descriptor + ".csv";

    WriteSeedHits(ossFileName.str(), hits, 0, hits.size(), 0, "query" + std::to_string(queryId),
                  queryLen, "all_targets", 0, false);
#endif
}

void DebugWriteChainedRegion(const std::vector<std::unique_ptr<ChainedRegion>>& allChainedRegions,
                             [[maybe_unused]] const int32_t debugStepId,
                             [[maybe_unused]] const std::string& descriptor,
                             [[maybe_unused]] const int32_t queryId,
                             [[maybe_unused]] const int32_t queryLen,
                             const bool debugVerboseChains = false)
{

// Special debug output, only if macros are defined.
#if defined(PANCAKE_MAP_CLR_DEBUG) || defined(PANCAKE_MAP_CLR_DEBUG_2)
    const double peakRssGb =
        PacBio::Utility::MemoryConsumption::PeakRss() / 1024.0 / 1024.0 / 1024.0;
    PBLOG_DEBUG << "After: step = " << debugStepId << ", '" << descriptor
                << "'. Chains: " << allChainedRegions.size() << ". Peak RSS: " << std::fixed
                << std::setprecision(3) << peakRssGb;
#endif
#ifdef PANCAKE_MAP_CLR_DEBUG_WRITE_SEED_HITS_TO_FILE
    // Write seed hits.
    for (size_t i = 0; i < allChainedRegions.size(); ++i) {
        auto& region = allChainedRegions[i];
        if (region->priority > 1) {
            continue;
        }
        std::ostringstream ossFileName;
        ossFileName << "temp-debug/hits-q" << std::to_string(queryId) << "-" << std::setfill('0')
                    << std::setw(3) << debugStepId << "-" << descriptor + ".csv";
        WriteSeedHits(ossFileName.str(), region->chain.hits, 0, region->chain.hits.size(), i,
                      "query" + std::to_string(queryId), queryLen,
                      "target" + std::to_string(region->mapping->Bid), region->mapping->Blen,
                      (i > 0));
    }
#endif

    // General debug which can optionally be printed out.
    if (debugVerboseChains) {
        PBLOG_DEBUG << "Chained regions: " << allChainedRegions.size();

        // Write seed hits.
        for (size_t regionId = 0; regionId < allChainedRegions.size(); ++regionId) {
            auto& region = allChainedRegions[regionId];
            // if (region->priority > 1) {
            //     continue;
            // }
            const ChainedRegion& cr = *region;
            PBLOG_DEBUG << "    - [regionId " << regionId
                        << "] chain.hits = " << cr.chain.hits.size()
                        << ", chain.score = " << cr.chain.score
                        << ", chain.covQ = " << cr.chain.coveredBasesQuery
                        << ", chain.covT = " << cr.chain.coveredBasesTarget
                        << ", priority = " << cr.priority
                        << ", isSuppl = " << (cr.isSupplementary ? "true" : "false")
                        << ", ovl: " << *cr.mapping
                        << ", diagStart = " << (cr.mapping->Astart - cr.mapping->Bstart)
                        << ", diagEnd = " << (cr.mapping->Aend - cr.mapping->Bend);
        }
        PBLOG_DEBUG << "";
    }
}

std::vector<MapperBaseResult> MapperCLR::WrapBuildIndexMapAndAlignWithFallback_(
    const FastaSequenceCachedStore& targetSeqs, const FastaSequenceCachedStore& querySeqs,
    const MapperCLRSettings& settings, std::shared_ptr<ChainingScratchSpace> ssChain,
    std::vector<SeedHit>& ssSeedHits, AlignerBasePtr& alignerGlobal, AlignerBasePtr& alignerExt)
{
    // Construct the index.
    std::vector<PacBio::Pancake::Int128t> seeds;
    const auto& seedParams = settings.map.seedParams;
    GenerateMinimizers(seeds, targetSeqs.records(), seedParams.KmerSize, seedParams.MinimizerWindow,
                       seedParams.Spacing, seedParams.UseRC, seedParams.UseHPCForSeedsOnly);
    std::unique_ptr<SeedIndex> seedIndex = std::make_unique<SeedIndex>(std::move(seeds));

    // Calculate the seed frequency statistics, needed for the cutoff.
    int64_t freqMax = 0;
    int64_t freqCutoff = 0;
    double freqAvg = 0.0;
    double freqMedian = 0.0;
    seedIndex->ComputeFrequencyStats(settings.map.freqPercentile, freqMax, freqAvg, freqMedian,
                                     freqCutoff);

    // Construct the fallback index only if required
    int64_t freqCutoffFallback = 0;
    std::unique_ptr<SeedIndex> seedIndexFallback = nullptr;
    if (settings.map.seedParamsFallback != settings.map.seedParams) {
        std::vector<PacBio::Pancake::Int128t> seedsFallback;
        const auto& seedParamsFallback = settings.map.seedParamsFallback;
        GenerateMinimizers(seedsFallback, targetSeqs.records(), seedParamsFallback.KmerSize,
                           seedParamsFallback.MinimizerWindow, seedParamsFallback.Spacing,
                           seedParamsFallback.UseRC, seedParamsFallback.UseHPCForSeedsOnly);
        seedIndexFallback = std::make_unique<SeedIndex>(std::move(seedsFallback));
        seedIndexFallback->ComputeFrequencyStats(settings.map.freqPercentile, freqMax, freqAvg,
                                                 freqMedian, freqCutoffFallback);
    }

    // Run mapping for each query.
    std::vector<MapperBaseResult> results;
    for (int32_t i = 0; i < Utility::Ssize(querySeqs.records()); ++i) {
        const auto& query = querySeqs.records()[i];
        int32_t queryId = query.Id();

#ifdef PANCAKE_MAP_CLR_DEBUG
        PBLOG_INFO << "[i = " << i << "] New query: queryId = " << queryId
                   << ", length = " << query.size();
#endif

        std::vector<PacBio::Pancake::Int128t> querySeeds;
        const std::string_view seq(query.c_str(), query.size());
        int rv = GenerateMinimizers(querySeeds, seq, 0, queryId, settings.map.seedParams.KmerSize,
                                    settings.map.seedParams.MinimizerWindow,
                                    settings.map.seedParams.Spacing, settings.map.seedParams.UseRC,
                                    settings.map.seedParams.UseHPCForSeedsOnly);
        if (rv) {
            throw std::runtime_error("Generating minimizers failed for the query sequence i = " +
                                     std::to_string(i) + ", id = " + std::to_string(queryId));
        }

        // Prepare the view of the seeds.
        const PacBio::Pancake::SequenceSeedsCached querySeedsCached(query.Name(), querySeeds.data(),
                                                                    querySeeds.size(), queryId);

        MapperBaseResult queryResults =
            WrapMapAndAlign_(targetSeqs, *seedIndex, query, querySeedsCached, queryId, freqCutoff,
                             settings, ssChain, ssSeedHits, alignerGlobal, alignerExt);

        if (queryResults.mappings.empty() && seedIndexFallback != nullptr) {
            rv = GenerateMinimizers(
                querySeeds, seq, 0, queryId, settings.map.seedParamsFallback.KmerSize,
                settings.map.seedParamsFallback.MinimizerWindow,
                settings.map.seedParamsFallback.Spacing, settings.map.seedParamsFallback.UseRC,
                settings.map.seedParamsFallback.UseHPCForSeedsOnly);
            if (rv) {
                throw std::runtime_error(
                    "Generating minimizers failed for the query sequence, id = " +
                    std::to_string(queryId));
            }

            // Prepare the view of the fallback seeds.
            const PacBio::Pancake::SequenceSeedsCached querySeedsCachedFallback(
                query.Name(), querySeeds.data(), querySeeds.size(), queryId);

            queryResults = WrapMapAndAlign_(
                targetSeqs, *seedIndexFallback, query, querySeedsCachedFallback, queryId,
                freqCutoffFallback, settings, ssChain, ssSeedHits, alignerGlobal, alignerExt);
        }

        for (const auto& m : queryResults.mappings) {
            m->mapping->Aid = queryId;
        }
        results.emplace_back(std::move(queryResults));

#ifdef PANCAKE_MAP_CLR_DEBUG
        PBLOG_INFO << "\n\n\n";
#endif
    }

    return results;
}

MapperBaseResult MapperCLR::WrapMapAndAlign_(
    const FastaSequenceCachedStore& targetSeqs, const PacBio::Pancake::SeedIndex& index,
    const FastaSequenceCached& querySeq, const SequenceSeedsCached& querySeeds,
    const int32_t queryId, int64_t freqCutoff, const MapperCLRSettings& settings,
    std::shared_ptr<ChainingScratchSpace> ssChain, std::vector<SeedHit>& ssSeedHits,
    AlignerBasePtr& alignerGlobal, AlignerBasePtr& alignerExt)
{
    const int32_t queryLen = querySeq.size();

    // Map the query.
    auto result = Map_(targetSeqs, index, querySeq, querySeeds, queryLen, queryId, settings,
                       freqCutoff, ssChain, ssSeedHits);

    // Align if needed.
    if (settings.align.align) {
        result = Align_(targetSeqs, querySeq, result, settings, alignerGlobal, alignerExt);
    }

    DebugWriteChainedRegion(result.mappings, 300, "result-final", queryId, queryLen,
                            DEBUG_VERBOSE_CHAINS);

    return result;
}

MapperBaseResult MapperCLR::Map_(const FastaSequenceCachedStore& targetSeqs, const SeedIndex& index,
                                 const FastaSequenceCached& querySeq,
                                 const SequenceSeedsCached& querySeeds, const int32_t queryLen,
                                 const int32_t queryId, const MapperCLRSettings& settings,
                                 int64_t freqCutoff, std::shared_ptr<ChainingScratchSpace> ssChain,
                                 std::vector<SeedHit>& ssSeedHits)
{

    // Debug info.
#if defined(PANCAKE_MAP_CLR_DEBUG) || defined(PANCAKE_MAP_CLR_DEBUG_2)
    {
        const double peakRssGb =
            PacBio::Utility::MemoryConsumption::PeakRss() / 1024.0 / 1024.0 / 1024.0;
        PBLOG_DEBUG << "Starting function: " << std::string(__FUNCTION__) << "."
                    << " Peak RSS: " << std::fixed << std::setprecision(3) << peakRssGb;

        const std::vector<std::pair<int64_t, int64_t>> debugSeedHist =
            PacBio::Pancake::ComputeSeedHitHistogram(
                {querySeeds.Seeds(), static_cast<size_t>(querySeeds.Size())}, index.GetHash());

        int64_t totalHits = 0;
        for (size_t i = 0; i < debugSeedHist.size(); ++i) {
            std::cerr << "[hist i = " << i << " / " << debugSeedHist.size()
                      << "] hits = " << debugSeedHist[i].first
                      << ", num_seeds = " << debugSeedHist[i].second
                      << ", multiplied = " << (debugSeedHist[i].first * debugSeedHist[i].second)
                      << "\n";
            totalHits += debugSeedHist[i].first * debugSeedHist[i].second;
        }
        std::cerr << "Total hits: " << totalHits << "\n";
    }
    const bool debugVerboseOccurrenceThreshold = true;
#else
    const bool debugVerboseOccurrenceThreshold = false;
#endif

    std::unordered_map<std::string, double> timings;
    TicToc ttMapAll;
    TicToc ttPartial;

    // Skip short queries.
    if (queryLen < settings.map.minQueryLen) {
        MapperBaseResult result;
        LogTicToc("map-L1-all", ttMapAll, result.time);
        return result;
    }

    ////////////////////////////////////////
    /// Compute the frequency threshold. ///
    ////////////////////////////////////////
    // Compute the seed hit histogram, but only if needed.
    std::vector<std::pair<int64_t, int64_t>> seedHitHistogram;
    if (settings.map.seedOccurrenceMaxMemory > 0) {
        seedHitHistogram = PacBio::Pancake::ComputeSeedHitHistogram(
            {querySeeds.Seeds(), static_cast<size_t>(querySeeds.Size())}, index.GetHash());
    }

    // The occurrence threshold is computed as:
    //      cutoff = max(seedOccurrenceMin, min(seedOccurrenceMax, occThresholdMemMax, seedOccurrenceUserSpecified))
    // where occThresholdMemMax is computed from the histogram if parameter seedOccurrenceMaxMemory > 0, and other parameters are user-provided.
    const int64_t occThreshold = ComputeOccurrenceThreshold(
        seedHitHistogram, settings.map.seedOccurrenceMin, settings.map.seedOccurrenceMax,
        settings.map.seedOccurrenceMaxMemory, freqCutoff, debugVerboseOccurrenceThreshold);

    LogTicToc("map-L1-01-occthresh", ttPartial, timings);

    ////////////////////////////////
    /// Collect seed hits.       ///
    ////////////////////////////////
    auto& hits = ssSeedHits;
    index.CollectHits(querySeeds.Seeds(), querySeeds.Size(), queryLen, hits, occThreshold);
    LogTicToc("map-L1-02-collect", ttPartial, timings);

#if defined(PANCAKE_MAP_CLR_DEBUG) || defined(PANCAKE_MAP_CLR_DEBUG_2)
    {
        const double peakRssGb =
            PacBio::Utility::MemoryConsumption::PeakRss() / 1024.0 / 1024.0 / 1024.0;
        PBLOG_DEBUG << "Collected hits: " << hits.size() << ", occThreshold = " << occThreshold
                    << ", freqCutoff = " << freqCutoff
                    << ", seedOccurrenceMin = " << settings.map.seedOccurrenceMin
                    << ", seedOccurrenceMax = " << settings.map.seedOccurrenceMax << "."
                    << " Peak RSS: " << std::fixed << std::setprecision(3) << peakRssGb;
    }
#endif

    ////////////////////////////////
    /// Mocking check.           ///
    ////////////////////////////////
    // Check if this read has at least one seed hit to itself. If so, self-mapping
    // will be added in the end.
    bool addPerfectMapping = false;
    if (settings.map.selfHitPolicy == MapperSelfHitPolicy::PERFECT_ALIGNMENT) {
        for (const SeedHit& hit : hits) {
            if (hit.targetId == queryId) {
                addPerfectMapping = true;
                break;
            }
        }
    }
    LogTicToc("map-L1-03-perfaln", ttPartial, timings);

    //////////////////////////////////////////
    /// Pre-filter self/symmetric hits.    ///
    //////////////////////////////////////////
    // Filter symmetric and self hits, including self hits for mocked mappings.
    if (settings.map.skipSymmetricOverlaps ||
        settings.map.selfHitPolicy != MapperSelfHitPolicy::DEFAULT) {
        std::vector<SeedHit> newHits = FilterSymmetricAndSelfHits(
            hits, queryId, settings.map.selfHitPolicy != MapperSelfHitPolicy::DEFAULT,
            settings.map.skipSymmetricOverlaps);
        std::swap(hits, newHits);
    }
    LogTicToc("map-L1-04-skipsym", ttPartial, timings);

#if defined(PANCAKE_MAP_CLR_DEBUG) || defined(PANCAKE_MAP_CLR_DEBUG_2)
    {
        const double peakRssGb =
            PacBio::Utility::MemoryConsumption::PeakRss() / 1024.0 / 1024.0 / 1024.0;
        PBLOG_DEBUG << "Sorted hits."
                    << " Peak RSS: " << std::fixed << std::setprecision(3) << peakRssGb;
    }
#endif

    /////////////////////////////////////
    /// Core mapping using seed hits. ///
    /////////////////////////////////////
    // Filter and refine seed hits, and construct mappings.
    int32_t debugStepId = 100;
    MapperBaseResult result = HitsToMappings_(
        hits, ssChain, targetSeqs, queryLen, queryId, settings.map, addPerfectMapping,
        settings.align.alnParamsGlobal.alignBandwidth, timings, debugStepId);

    if (addPerfectMapping) {
        return result;
    }

    /////////////////////////////////////
    /// Reseeding.                    ///
    /////////////////////////////////////

    if (settings.map.reseedGaps) {
        std::vector<SeedHit> newHits = ReseedAlignmentRegions_(
            result, targetSeqs, querySeq, settings.map.reseedGapMinLength,
            settings.map.reseedGapMaxLength, settings.map.reseedSeedParams,
            settings.map.reseedFreqPercentile, settings.map.reseedOccurrenceMin,
            settings.map.reseedOccurrenceMax, settings.map.reseedOccurrenceMaxMemory);

        DebugWriteSeedHits(newHits, debugStepId, "reseeded", querySeq.Id(), querySeq.size());
        ++debugStepId;

        newHits.reserve(hits.size() + newHits.size());

        for (const auto& cr : result.mappings) {
            if (cr == nullptr) {
                continue;
            }
            newHits.insert(newHits.end(), cr->chain.hits.begin(), cr->chain.hits.end());
        }

        result = HitsToMappings_(newHits, ssChain, targetSeqs, queryLen, queryId, settings.map,
                                 addPerfectMapping, settings.align.alnParamsGlobal.alignBandwidth,
                                 timings, debugStepId);
    }

    LogTicToc("map-all", ttMapAll, timings);

    return result;
}

MapperBaseResult MapperCLR::HitsToMappings_(
    const std::vector<SeedHit>& inputHits, const std::shared_ptr<ChainingScratchSpace>& ssChain,
    const FastaSequenceCachedStore& targetSeqs, const int32_t queryLen, const int32_t queryId,
    const MapperCLRMapSettings& settings, const bool addPerfectMapping,
    const int32_t maxAllowedDistForBadEndRefinement,
    const std::unordered_map<std::string, double>& timings, int32_t& debugStepId)
{
    TicToc ttPartial;

    MapperBaseResult result;
    result.time = timings;

    std::vector<SeedHit> sortedHits;

    // Sort the seed hits.
    {
        // Unfortunately, we cannot just pack the values into an 128-bit integer in the order
        // of <targetID, targetRev, diag, targetPos, queryPos>, because the diagonal is internal,
        // and signedness matters here. In that case, negative diagonals would actually get sorted _after_
        // the positive ones. So we have to use pairs/tuples to properly compare the values here.
        // Here, we attempt to optimize the amount of comparisons that need to be performed when sorting:
        //  - First of the pair is the targetID + targetRev.
        //  - Second of the pair stores the diagonal first, then packs the ID of the current seed hit. This way
        //      the diagonal's sign can validly be compared.
        std::vector<std::pair<int32_t, int64_t>> seedHitsPacked(inputHits.size());
        for (size_t i = 0; i < seedHitsPacked.size(); ++i) {
            const auto& sh = inputHits[i];
            seedHitsPacked[i] = std::make_pair(
                ((sh.targetId << 1) | sh.targetRev),
                (static_cast<uint64_t>(sh.targetPos - sh.queryPos) << 32 | (i & 0x0FFFFFFFF)));
        }
        const double timeSortPrepare = ttPartial.GetMicrosecs(true);
        LogTicTocAdd("map-L3-total-sort", timeSortPrepare, result.time);
        LogTicToc("map-L1-05-sort-prepare", ttPartial, result.time);

        pdqsort(seedHitsPacked.begin(), seedHitsPacked.end());

        sortedHits.resize(inputHits.size());
        for (size_t i = 0; i < seedHitsPacked.size(); ++i) {
            sortedHits[i] = inputHits[seedHitsPacked[i].second & 0x0FFFFFFFF];
        }
        const double timeSort = ttPartial.GetMicrosecs(true);
        LogTicTocAdd("map-L3-total-sort", timeSort, result.time);
        LogTicToc("map-L1-05-sort", ttPartial, result.time);
    }

    // Group seed hits by diagonal.
    auto groups = DiagonalGroup(sortedHits, settings.chainBandwidth, true);
    LogTicToc("map-L1-06-diaggroup", ttPartial, result.time);

#if defined(PANCAKE_MAP_CLR_DEBUG) || defined(PANCAKE_MAP_CLR_DEBUG_2)
    {
        const double peakRssGb =
            PacBio::Utility::MemoryConsumption::PeakRss() / 1024.0 / 1024.0 / 1024.0;
        PBLOG_DEBUG << "Diagonal groups: " << groups.size() << "."
                    << " Peak RSS: " << std::fixed << std::setprecision(3) << peakRssGb;
    }
#endif

    // Process each diagonal bin to get the chains.
    std::vector<std::unique_ptr<ChainedRegion>> allChainedRegions = ChainAndMakeOverlap_(
        targetSeqs, sortedHits, groups, queryId, queryLen, settings.chainMaxSkip,
        settings.chainMaxPredecessors, settings.seedJoinDist, settings.chainBandwidth,
        settings.minNumSeeds, settings.minCoveredBases, settings.minDPScore, settings.useLIS,
        ssChain, result.time);
    LogTicToc("map-L1-07-chain", ttPartial, result.time);
    DebugWriteChainedRegion(allChainedRegions, debugStepId, "map-chain-and-make-overlap", queryId,
                            queryLen, DEBUG_VERBOSE_CHAINS);
    ++debugStepId;

    // Take the remaining regions, merge all seed hits, and rechain.
    // Needed because diagonal chaining was greedy and a wide window could have split
    // otherwise good chains. On the other hand, diagonal binning was needed for speed
    // in low-complexity regions.
    allChainedRegions = ReChainSeedHits_(
        allChainedRegions, targetSeqs, queryId, queryLen, settings.chainMaxSkip,
        settings.chainMaxPredecessors, settings.seedJoinDist, settings.chainBandwidth,
        settings.minNumSeeds, settings.minCoveredBases, settings.minDPScore, ssChain, result.time);
    LogTicToc("map-L1-08-rechain", ttPartial, result.time);
    DebugWriteChainedRegion(allChainedRegions, debugStepId, "map-rechain-hits", queryId, queryLen,
                            DEBUG_VERBOSE_CHAINS);
    ++debugStepId;

    // Sort all chains in descending order of the number of hits.
    pdqsort(allChainedRegions.begin(), allChainedRegions.end(),
            [](const auto& a, const auto& b) { return a->chain.score > b->chain.score; });
    LogTicToc("map-L1-09-sortchained", ttPartial, result.time);

    // Secondary/supplementary flagging.
    WrapFlagSecondaryAndSupplementary(
        allChainedRegions, settings.secondaryAllowedOverlapFractionQuery,
        settings.secondaryAllowedOverlapFractionTarget, settings.secondaryMinScoreFraction);
    LogTicToc("map-L1-10-secondary", ttPartial, result.time);
    DebugWriteChainedRegion(allChainedRegions, debugStepId, "map-wrap-flag-secondary-suppl",
                            queryId, queryLen, DEBUG_VERBOSE_CHAINS);
    ++debugStepId;

    // Merge long gaps.
    LongMergeChains(allChainedRegions, settings.longMergeBandwidth);
    LogTicToc("map-L1-11-longmerge", ttPartial, result.time);
    DebugWriteChainedRegion(allChainedRegions, debugStepId, "map-long-merge-chains", queryId,
                            queryLen, DEBUG_VERBOSE_CHAINS);
    ++debugStepId;

    // Add an extra query alignment only if needed (i.e. if the MapperSelfHitPolicy == PERFECT_ALIGNMENT
    // and the query has self-hits ("self" in term of queryId/targetId).
    // NOTE: This function does not check the settings.selfHitPolicy, instead it depends on parameter `addPerfectMapping`
    // to mark this.
    if (addPerfectMapping) {
        const int32_t targetId = queryId;
        const int32_t targetLen = targetSeqs.GetSequence(targetId).size();
        std::unique_ptr<ChainedRegion> newChainedRegion =
            CreateMockedMapping(queryId, queryLen, targetId, targetLen);
        allChainedRegions.emplace_back(std::move(newChainedRegion));
    }
    LogTicToc("map-L1-12-addperfect", ttPartial, result.time);

    // Again relabel, because some chains are longer now.
    WrapFlagSecondaryAndSupplementary(
        allChainedRegions, settings.secondaryAllowedOverlapFractionQuery,
        settings.secondaryAllowedOverlapFractionTarget, settings.secondaryMinScoreFraction);
    LogTicToc("map-L1-13-secondary", ttPartial, result.time);
    DebugWriteChainedRegion(allChainedRegions, debugStepId, "map-wrap-flag-secondary-suppl",
                            queryId, queryLen, DEBUG_VERBOSE_CHAINS);
    ++debugStepId;

    // Sort all chains by priority and then score.
    pdqsort(allChainedRegions.begin(), allChainedRegions.end(),
            [](const std::unique_ptr<ChainedRegion>& a, const std::unique_ptr<ChainedRegion>& b) {
                auto at = std::tuple<int32_t, bool, int32_t>(a->priority, a->isSupplementary,
                                                             a->chain.score);
                auto bt = std::tuple<int32_t, bool, int32_t>(b->priority, b->isSupplementary,
                                                             b->chain.score);
                return at < bt;
            });
    LogTicToc("map-L1-14-sortchained", ttPartial, result.time);

    // Refine seed hits for alignment.
    for (size_t i = 0; i < allChainedRegions.size(); ++i) {
        auto& region = allChainedRegions[i];
        if (region->priority > 1) {
            continue;
        }

        // Always refine seed hits near the ends to avoid edge cases in DP chaining.
        ChainedHits newChain = RefineBadEnds(region->chain, maxAllowedDistForBadEndRefinement,
                                             settings.minDPScore * 2);

        // Optionally refine internal seed hits.
        if (settings.refineSeedHits) {
            newChain =
                RefineChainedHits(newChain, settings.refineMinGap1, settings.refineDiffThreshold,
                                  settings.seedJoinDist / 2, 10);
            newChain =
                RefineChainedHits2(newChain, settings.refineMinGap2, settings.seedJoinDist / 2);
        }

        std::swap(region->chain, newChain);

        if (region->chain.hits.empty()) {
            region->mapping = nullptr;
        } else {
            region->mapping =
                MakeOverlap(region->chain.hits, queryId, queryLen, targetSeqs, 0,
                            region->chain.hits.size(), 0, region->chain.hits.size() - 1);
        }
    }
    LogTicToc("map-L1-15-refineseeds", ttPartial, result.time);
    DebugWriteChainedRegion(allChainedRegions, debugStepId, "map-refining-seed-hits", queryId,
                            queryLen, DEBUG_VERBOSE_CHAINS);
    ++debugStepId;

    // Filter out the mappings.
    result.mappings = std::move(allChainedRegions);
    const int32_t numPrimary = CondenseMappings(result.mappings, settings.bestNSecondary);
    LogTicToc("map-L1-16-condense", ttPartial, result.time);
    DebugWriteChainedRegion(result.mappings, debugStepId, "map-condense", queryId, queryLen,
                            +DEBUG_VERBOSE_CHAINS);
    ++debugStepId;

    // If this occurs, that means that a filtering stage removed the primary alignment for some reason.
    // This can happen in case self-hits are skipped in the overlapping use case (the self-hit is the
    // highest scoring one, and it would get marked for removal).
    // This reruns labeling to produce the next best primary.
    if (numPrimary == 0) {
        WrapFlagSecondaryAndSupplementary(
            allChainedRegions, settings.secondaryAllowedOverlapFractionQuery,
            settings.secondaryAllowedOverlapFractionTarget, settings.secondaryMinScoreFraction);
    }
    LogTicToc("map-L1-17-secondary", ttPartial, result.time);

    // Collect regions for alignment.
    for (size_t i = 0; i < result.mappings.size(); ++i) {
        if (result.mappings[i] == nullptr || result.mappings[i]->mapping == nullptr) {
            continue;
        }
        result.mappings[i]->regionsForAln =
            CollectAlignmentRegions(*result.mappings[i], settings.minAlignmentSpan,
                                    settings.maxFlankExtensionDist, settings.flankExtensionFactor);
    }

    LogTicToc("map-L1-18-collectalnregions", ttPartial, result.time);
    DebugWriteChainedRegion(result.mappings, debugStepId, "map-collect-aln-regions", queryId,
                            queryLen, DEBUG_VERBOSE_CHAINS);
    ++debugStepId;

#ifdef PANCAKE_MAP_CLR_DEBUG_2
    std::cerr << "All hits: sortedHits.size() = " << sortedHits.size() << "\n";
    std::cerr << "Diagonal groups: groups.size() = " << groups.size() << "\n";
    for (size_t i = 0; i < groups.size(); ++i) {
        const int32_t firstDiag = sortedHits[groups[i].start].Diagonal();
        const int32_t lastDiag = sortedHits[groups[i].end - 1].Diagonal();
        std::cerr << "[queryId = " << queryId << ", group " << i << "] start = " << groups[i].start
                  << ", end = " << groups[i].end << ", diagStart = " << firstDiag
                  << ", diagEnd = " << lastDiag << "\n";
    }
#endif

#ifdef PANCAKE_MAP_CLR_DEBUG_WRITE_SEED_HITS_TO_FILE
    // Write ALL seed hits.
    for (size_t i = 0; i < groups.size(); ++i) {
        const auto& g = groups[i];
        const int32_t targetId = sortedHits[g.start].targetId;
        const int32_t targetLen = targetSeqs.GetSequence(targetId).size();
        std::ostringstream ossFile;
        ossFile << "temp-debug/hits-q" << std::to_string(queryId) << "-" << std::setfill('0')
                << std::setw(3) << 0 << "-all-hits-diagonal-groupped.csv";
        WriteSeedHits(ossFile.str(), sortedHits, g.start, g.end, i,
                      "query" + std::to_string(queryId), queryLen,
                      "target" + std::to_string(targetId), targetLen, (i > 0));
    }
#endif

    return result;
}

std::vector<SeedHit> MapperCLR::ReseedAlignmentRegions_(
    const MapperBaseResult& result, const FastaSequenceCachedStore& targetSeqs,
    const FastaSequenceCached& querySeq, const int32_t reseedGapMinLength,
    const int32_t reseedGapMaxLength, const PacBio::Pancake::SeedDBParameters& seedParams,
    const double reseedFreqPercentile, const int64_t reseedOccurrenceMin,
    const int64_t reseedOccurrenceMax, const int64_t reseedOccurrenceMaxMemory)
{

    auto ContainsLargeRegions = [](const MapperBaseResult& inData, const int32_t maxAllowedSpan) {
        for (const auto& cr : inData.mappings) {
            if (cr == nullptr) {
                continue;
            }
            for (const auto& ar : cr->regionsForAln) {
                if (ar.qSpan > maxAllowedSpan || ar.tSpan > maxAllowedSpan) {
                    return true;
                }
            }
        }
        return false;
    };

    std::vector<SeedHit> newHits;

    if (ContainsLargeRegions(result, reseedGapMinLength) == false) {
        return newHits;
    }

    for (const auto& cr : result.mappings) {
        if (cr == nullptr) {
            continue;
        }
        for (size_t regionId = 0; regionId < cr->regionsForAln.size(); ++regionId) {
            const auto& ar = cr->regionsForAln[regionId];

            // Skip short gaps.
            if (ar.qSpan <= reseedGapMinLength && ar.tSpan <= reseedGapMinLength) {
                continue;
            }

            // Optionally skip large gaps.
            if (reseedGapMaxLength > 0 &&
                (ar.qSpan > reseedGapMaxLength || ar.tSpan > reseedGapMaxLength)) {
                continue;
            }

            // We can index query in fwd because indexing will be bidirectional.
            // This is actually preferred, to be in line with other seed hits (computed in the same way).
            const int32_t qStartFwd =
                (ar.queryRev) ? (querySeq.size() - (ar.qStart + ar.qSpan)) : ar.qStart;

            const FastaSequenceCachedStore localQueryStore(
                {FastaSequenceCached("q", querySeq.c_str() + qStartFwd, ar.qSpan, querySeq.Id())});

            // Create the local target sequence.
            const FastaSequenceCached& targetSeq = targetSeqs.GetSequence(cr->mapping->Bid);
            const FastaSequenceCachedStore localTargetStore({FastaSequenceCached(
                "t", targetSeq.c_str() + ar.tStart, ar.tSpan, targetSeq.Id())});

            // Not const because queryPos and targetPos values will be updated below.
            std::vector<std::vector<SeedHit>> localHits = CollectSeedHitsFromSequences(
                localQueryStore, localTargetStore, seedParams.KmerSize, seedParams.MinimizerWindow,
                seedParams.Spacing, seedParams.UseRC, seedParams.UseHPCForSeedsOnly,
                reseedFreqPercentile, reseedOccurrenceMin, reseedOccurrenceMax,
                reseedOccurrenceMaxMemory);

            for (size_t qId = 0; qId < localHits.size(); ++qId) {
                newHits.reserve(newHits.size() + localHits[qId].size());
                for (size_t hitId = 0; hitId < localHits[qId].size(); ++hitId) {
                    auto& hit = localHits[qId][hitId];
                    // Collected seed hits are in the strand of the query, while target is always fwd.
                    // The ar.qStart is also always in the strand of the alignment, so just adding it should line up perfectly.
                    hit.queryPos += ar.qStart;
                    hit.targetPos += ar.tStart;
                    newHits.emplace_back(hit);
                }
            }
        }
    }

    return newHits;
}

MapperBaseResult MapperCLR::Align_(const FastaSequenceCachedStore& targetSeqs,
                                   const FastaSequenceCached& querySeq,
                                   const MapperBaseResult& mappingResult,
                                   const MapperCLRSettings& settings, AlignerBasePtr& alignerGlobal,
                                   AlignerBasePtr& alignerExt)
{
    TicToc ttAlignAll;
    TicToc ttPartial;

    int32_t debugStepId = 200;

#ifdef PANCAKE_MAP_CLR_DEBUG_ALIGN
    std::cerr << "Aligning.\n";
    std::cerr << "alignerTypeGlobal = " << AlignerTypeToString(settings.align.alignerTypeGlobal)
              << "\n";
    std::cerr << "alignerTypeExt = " << AlignerTypeToString(settings.align.alignerTypeExt) << "\n";
#endif

    const int32_t queryLen = querySeq.size();
    const int32_t queryId = querySeq.Id();

    MapperBaseResult alignedResult;

    // Reverse the query sequence, needed for alignment.
    const std::string querySeqRev = PacBio::Pancake::ReverseComplement(
        {querySeq.c_str(), static_cast<size_t>(querySeq.size())}, 0, querySeq.size());

    for (size_t i = 0; i < mappingResult.mappings.size(); ++i) {
        // const auto& chain = result.mappings[i]->chain;
        const auto& chain = mappingResult.mappings[i]->chain;
        const auto& ovl = mappingResult.mappings[i]->mapping;
        const auto& tSeqFwd = targetSeqs.GetSequence(ovl->Bid);
        const bool isIdIdentical = (ovl->Aid == ovl->Bid);

        // Optionally skip, mock or align the overlap.
        if (settings.align.selfHitPolicy == MapperSelfHitPolicy::SKIP && isIdIdentical) {
            // Pass, no need to generate any overlaps.
            PBLOG_TRACE << "(" << __FUNCTION__
                        << ") MapperSelfHitPolicy::SKIP. Skipping self alignment.";

        } else if (settings.align.selfHitPolicy == MapperSelfHitPolicy::PERFECT_ALIGNMENT &&
                   isIdIdentical) {
            // Mock the perfect alignment between the sequence and itself.
            PBLOG_TRACE << "(" << __FUNCTION__
                        << ") MapperSelfHitPolicy::PERFECT_ALIGNMENT. "
                           "Mocking self alignment instead of actually "
                           "aligning.";
            OverlapPtr newOvl =
                CreateMockedAlignment(ovl, settings.align.alnParamsGlobal.matchScore);
            auto newChainedRegion = std::make_unique<ChainedRegion>();
            newChainedRegion->chain = chain;
            newChainedRegion->mapping = std::move(newOvl);
            newChainedRegion->priority = mappingResult.mappings[i]->priority;
            newChainedRegion->isSupplementary = mappingResult.mappings[i]->isSupplementary;
            newChainedRegion->isMockedMapping = mappingResult.mappings[i]->isMockedMapping;
            newChainedRegion->isMockedAlignment = true;
            alignedResult.mappings.emplace_back(std::move(newChainedRegion));

        } else {
            PBLOG_TRACE << "(" << __FUNCTION__ << ") MapperSelfHitPolicy else. Aligning.";

            // Use a custom aligner to align.
            auto newOvl = AlignmentSeeded(ovl, mappingResult.mappings[i]->regionsForAln,
                                          std::string_view(tSeqFwd.data(), tSeqFwd.size()),
                                          std::string_view(querySeq.data(), querySeq.size()),
                                          querySeqRev, alignerGlobal, alignerExt);

            auto newChainedRegion = std::make_unique<ChainedRegion>();
            newChainedRegion->chain = chain;
            newChainedRegion->mapping = std::move(newOvl);
            newChainedRegion->priority = mappingResult.mappings[i]->priority;
            newChainedRegion->isSupplementary = mappingResult.mappings[i]->isSupplementary;
            newChainedRegion->isMockedMapping = mappingResult.mappings[i]->isMockedMapping;
            newChainedRegion->isMockedAlignment = false;
            alignedResult.mappings.emplace_back(std::move(newChainedRegion));
        }

#ifdef PANCAKE_MAP_CLR_DEBUG_ALIGN
        std::cerr << "[mapping i = " << i << ", before alignment] ovl: " << *ovl << "\n";
        const auto& updatedOvl = alignedResult.mappings[i]->mapping;
        std::cerr << "[mapping i = " << i << ", after alignment] ovl: ";
        if (updatedOvl != nullptr) {
            std::cerr << *updatedOvl << "\n";
        } else {
            std::cerr << "nullptr\n";
        }
#endif
    }
    LogTicToc("aln-L1-01-alnseeded", ttPartial, alignedResult.time);

    // Secondary/supplementary flagging.
    WrapFlagSecondaryAndSupplementary(
        alignedResult.mappings, settings.map.secondaryAllowedOverlapFractionQuery,
        settings.map.secondaryAllowedOverlapFractionTarget, settings.map.secondaryMinScoreFraction);

    LogTicToc("aln-L1-02-secondary", ttPartial, alignedResult.time);
    DebugWriteChainedRegion(alignedResult.mappings, debugStepId, "result-after-align", queryId,
                            queryLen, DEBUG_VERBOSE_CHAINS);
    ++debugStepId;

    const int32_t numPrimary =
        CondenseMappings(alignedResult.mappings, settings.map.bestNSecondary);
    LogTicToc("aln-L1-03-condense", ttPartial, alignedResult.time);
    DebugWriteChainedRegion(alignedResult.mappings, debugStepId, "align-condense-mappings", queryId,
                            queryLen, DEBUG_VERBOSE_CHAINS);
    ++debugStepId;

    // If this occurs, that means that a filtering stage removed the primary alignment for some reason.
    // This can happen in case self-hits are skipped in the overlapping use case (the self-hit is the
    // highest scoring one, and it would get marked for removal).
    // However, in this current function this cannot happen right now because flagging happens right
    // before condensing. Still including this here to be future proof.
    // This reruns labeling to produce the next best primary.
    if (numPrimary == 0) {
        WrapFlagSecondaryAndSupplementary(alignedResult.mappings,
                                          settings.map.secondaryAllowedOverlapFractionQuery,
                                          settings.map.secondaryAllowedOverlapFractionTarget,
                                          settings.map.secondaryMinScoreFraction);
    }
    LogTicToc("aln-L1-04-secondary", ttPartial, alignedResult.time);
    DebugWriteChainedRegion(alignedResult.mappings, debugStepId, "align-wrap-flag-secondary-suppl",
                            queryId, queryLen, DEBUG_VERBOSE_CHAINS);
    ++debugStepId;

    LogTicToc("aln-all", ttAlignAll, alignedResult.time);

    return alignedResult;
}

std::vector<std::unique_ptr<ChainedRegion>> MapperCLR::ReChainSeedHits_(
    const std::vector<std::unique_ptr<ChainedRegion>>& chainedRegions,
    const FastaSequenceCachedStore& targetSeqs, int32_t queryId, int32_t queryLen,
    int32_t chainMaxSkip, int32_t chainMaxPredecessors, int32_t seedJoinDist,
    int32_t chainBandwidth, int32_t minNumSeeds, int32_t minCoveredBases, int32_t minDPScore,
    std::shared_ptr<ChainingScratchSpace> ssChain,
    std::unordered_map<std::string, double>& retTimings)
{
    TicToc ttPartial;

#ifdef PANCAKE_MAP_CLR_DEBUG_2
    std::cerr << "(ReChainSeedHits_) Starting to rechain the seed hits.\n";
#endif

    std::vector<SeedHit> hits2;
    {
        // It turns out to be much faster to copy and convert data than to use a custom comparison
        // operator to compare tuples. Also, Pdqsort is more efficient than std::sort.
        // Sort the hits by coordinates.
        // Compute the size for collected hits.
        size_t newSize = 0;
        for (size_t i = 0; i < chainedRegions.size(); ++i) {
            auto& region = chainedRegions[i];
            if (region->priority > 1) {
                continue;
            }
            newSize += region->chain.hits.size();
        }
        // Pack and merge all the hits.
        std::vector<Int128t> seedHitsPacked(newSize);
        for (size_t i = 0, newId = 0; i < chainedRegions.size(); ++i) {
            auto& region = chainedRegions[i];
            if (region->priority > 1) {
                continue;
            }
            for (size_t j = 0; j < region->chain.hits.size(); ++j) {
                seedHitsPacked[newId] = region->chain.hits[j].PackTo128();
                ++newId;
            }
        }
        const double timeSortPrepare = ttPartial.GetMicrosecs(true);
        LogTicTocAdd("map-L3-total-sort", timeSortPrepare, retTimings);
        LogTicToc("map-L2-rechain-01-sort-prepare", ttPartial, retTimings);

        // Sort.
        pdqsort(seedHitsPacked.begin(), seedHitsPacked.end());
        const double timeSort = ttPartial.GetMicrosecs(true);
        LogTicTocAdd("map-L3-total-sort", timeSort, retTimings);
        LogTicToc("map-L2-rechain-01-sort", ttPartial, retTimings);

        // Unpack.
        hits2.resize(seedHitsPacked.size());
        for (size_t i = 0; i < seedHitsPacked.size(); ++i) {
            hits2[i].ParseFrom128(seedHitsPacked[i]);
        }
        const double timeSortPost = ttPartial.GetMicrosecs(true);
        LogTicTocAdd("map-L3-total-sort", timeSortPost, retTimings);
        LogTicToc("map-L2-rechain-01-sort-post", ttPartial, retTimings);
    }

    const auto groups = GroupByTargetAndStrand(hits2);

    LogTicToc("map-L2-rechain-02-group", ttPartial, retTimings);

    std::vector<std::unique_ptr<ChainedRegion>> newChainedRegions;

    for (size_t groupId = 0; groupId < groups.size(); ++groupId) {
        const auto& group = groups[groupId];
#ifdef PANCAKE_MAP_CLR_DEBUG_2
        std::cerr << "(ReChainSeedHits_) group: " << group << ", span = " << group.Span() << "\n";
        std::cerr << "[group " << groupId << "] After GroupByTargetAndStrand:\n";
        for (int32_t i = group.start; i < group.end; ++i) {
            const auto& hit = hits2[i];
            std::cerr << "    [hit i = " << i << "] " << hit << "\n";
        }
        std::cerr << "\n";
#endif

        // DP Chaining of the filtered hits to remove outliers.
        double timeChaining = 0.0, timeBacktrack = 0.0;
        std::vector<ChainedHits> chains =
            ChainHits({&hits2[group.start], static_cast<size_t>(group.end - group.start)},
                      chainMaxSkip, chainMaxPredecessors, seedJoinDist, chainBandwidth, minNumSeeds,
                      minCoveredBases, minDPScore, timeChaining, timeBacktrack, ssChain);

        double timeChainHits = ttPartial.GetMicrosecs(true);
        LogTicTocAdd("map-L2-rechain-03-chainhits", ttPartial, retTimings);
        LogTicTocAdd("map-L3-rechain-01-chaining", timeChaining, retTimings);
        LogTicTocAdd("map-L3-rechain-02-backtrack", timeBacktrack, retTimings);
        LogTicTocAdd("map-L2-total-chainhits-simd", timeChainHits, retTimings);
        LogTicTocAdd("map-L3-total-chaining", timeChaining, retTimings);
        LogTicTocAdd("map-L3-total-backtrack", timeBacktrack, retTimings);

#ifdef PANCAKE_MAP_CLR_DEBUG_2
        std::cerr << "[group " << groupId << "] After ChainHits:\n";
        for (size_t chainId = 0; chainId < chains.size(); ++chainId) {
            const auto& chain = chains[chainId];
            // std::cerr << chain << "\n";
            for (size_t hitId = 0; hitId < chain.hits.size(); ++hitId) {
                const auto& hit = chain.hits[hitId];
                std::cerr << "    [groupId = " << groupId << ", chainId = " << chainId
                          << ", hitId = " << hitId << "] " << hit << "\n";
            }
            std::cerr << "\n";
        }
#endif

        // Accumulate chains and their mapped regions.
        for (size_t i = 0; i < chains.size(); ++i) {
            const auto& chain = chains[i];
            int32_t numHitsInChain = chain.hits.size();

            // Filter.
            if (numHitsInChain < minNumSeeds) {
                continue;
            }

            // Create a new chained region.
            auto ovl = MakeOverlap(chain.hits, queryId, queryLen, targetSeqs, 0, numHitsInChain, 0,
                                   numHitsInChain - 1);
            auto chainedRegion = std::make_unique<ChainedRegion>();
            chainedRegion->chain = std::move(chain);
            chainedRegion->mapping = std::move(ovl);
            chainedRegion->isMockedMapping = false;
            chainedRegion->isMockedAlignment = false;
            newChainedRegions.emplace_back(std::move(chainedRegion));
        }
        LogTicTocAdd("map-L2-rechain-04-makeovl", ttPartial, retTimings);
    }

    return newChainedRegions;
}

std::vector<std::unique_ptr<ChainedRegion>> MapperCLR::ChainAndMakeOverlap_(
    const FastaSequenceCachedStore& targetSeqs, const std::vector<SeedHit>& hits,
    const std::vector<PacBio::Pancake::Range>& hitGroups, int32_t queryId, int32_t queryLen,
    int32_t chainMaxSkip, int32_t chainMaxPredecessors, int32_t seedJoinDist,
    int32_t chainBandwidth, int32_t minNumSeeds, int32_t minCoveredBases, int32_t minDPScore,
    bool useLIS, std::shared_ptr<ChainingScratchSpace> ssChain,
    std::unordered_map<std::string, double>& retTimings)
{
    TicToc ttPartial;

    // Process each diagonal bin to get the final chains.
    std::vector<std::unique_ptr<ChainedRegion>> allChainedRegions;
    for (const auto& range : hitGroups) {
#ifdef PANCAKE_MAP_CLR_DEBUG_2
        std::cerr << "(ChainAndMakeOverlap_) [range] start = " << range.start
                  << ", end = " << range.end << ", span = " << range.Span() << "\n";
#endif
        // Skip diagonals with insufficient hits.
        if (range.Span() < minNumSeeds) {
#ifdef PANCAKE_MAP_CLR_DEBUG_2
            std::cerr << "  -> Skipping, range.Span() < minNumSeeds\n";
#endif
            continue;
        }

        std::vector<SeedHit> groupHits(range.Span());
        {
            // Hits have previously been sorted by diagonals, and not by coordinates, so we need to sort again
            // to get them in proper order.
            // It turns out to be much faster to copy and convert data than to use a custom comparison
            // operator to compare tuples. Also, Pdqsort is more efficient than std::sort.
            std::vector<Int128t> seedHitsPacked(range.Span());
            for (int32_t i = range.start, hitId = 0; i < range.end; ++i, ++hitId) {
                seedHitsPacked[hitId] = hits[i].PackTo128();
            }
            ttPartial.Stop();
            const double timeSortPrepare = ttPartial.GetMicrosecs(true);
            LogTicTocAdd("map-L3-total-sort", timeSortPrepare, retTimings);
            LogTicTocAdd("map-L2-chain-01-sort-prepare", timeSortPrepare, retTimings);
            ttPartial.Start();

            pdqsort(seedHitsPacked.begin(), seedHitsPacked.end());
            ttPartial.Stop();
            const double timeSort = ttPartial.GetMicrosecs(true);
            LogTicTocAdd("map-L3-total-sort", timeSort, retTimings);
            LogTicTocAdd("map-L2-chain-01-sort", timeSort, retTimings);
            ttPartial.Start();

            for (size_t i = 0; i < seedHitsPacked.size(); ++i) {
                groupHits[i].ParseFrom128(seedHitsPacked[i]);
            }
            ttPartial.Stop();
            const double timeSortPost = ttPartial.GetMicrosecs(true);
            LogTicTocAdd("map-L3-total-sort", timeSortPost, retTimings);
            LogTicTocAdd("map-L2-chain-01-sort-post", timeSortPost, retTimings);
            ttPartial.Start();
        }

        // Perform chaining.
        std::vector<ChainedHits> chains;

        if (useLIS) {
            // Longest Increasing Subsequence.
            /**
             * Difference between the sort comparison and the LIS comparison:
             *  - When sorting 2D points we need to restrict ourselves to 1 dimension. In this case, we only sort on A.x < B.x, except in case when A.x == B.x.
             *  - LIS comparison returns true if and only if A is in the upper left corner from B. Meaning, both A.x and A.y need to be smaller.
             * For example:
             *  Point{6, 7} and Point{9, 3} would result in:
             *      - Sort Comparison == true
             *      - LIS Comparison == false
             *      - Packed into uint64_t and using the standard operator< == true (same as the Sort Comparison).
            */
            std::vector<PacBio::Pancake::SeedHit> lisHits =
                istl::LIS(groupHits, [](const SeedHit& a, const SeedHit& b) {
                    return (a.queryPos < b.queryPos && a.targetPos < b.targetPos);
                });
            LogTicTocAdd("map-L2-chain-02-lis", ttPartial, retTimings);

            // DP Chaining of the filtered hits to remove outliers.
            double timeChaining = 0.0, timeBacktrack = 0.0;
            chains = ChainHits(lisHits, chainMaxSkip, chainMaxPredecessors, seedJoinDist,
                               chainBandwidth, minNumSeeds, minCoveredBases, minDPScore,
                               timeChaining, timeBacktrack, ssChain);

            const double timeChainHits = ttPartial.GetMicrosecs(true);
            LogTicTocAdd("map-L2-chain-03-chainhits", ttPartial, retTimings);
            LogTicTocAdd("map-L3-chain-01-chaining", timeChaining, retTimings);
            LogTicTocAdd("map-L3-chain-02-backtrack", timeBacktrack, retTimings);
            LogTicTocAdd("map-L2-total-chainhits-simd", timeChainHits, retTimings);
            LogTicTocAdd("map-L3-total-chaining", timeChaining, retTimings);
            LogTicTocAdd("map-L3-total-backtrack", timeBacktrack, retTimings);

#ifdef PANCAKE_MAP_CLR_DEBUG_2
            std::cerr << "  - Hits before LIS:\n";
            for (size_t ii = 0; ii < groupHits.size(); ++ii) {
                const auto& hit = groupHits[ii];
                std::cerr << "    [groupHits ii = " << ii << "] " << hit << "\n";
            }
            std::cerr << "  - using LIS.\n";
            std::cerr << "  - lisHits.size() = " << lisHits.size() << "\n";
            std::cerr << "  - Hits after LIS:\n";
            for (size_t ii = 0; ii < lisHits.size(); ++ii) {
                const auto& hit = lisHits[ii];
                std::cerr << "    [lisHits ii = " << ii << "] " << hit << "\n";
            }
#endif
        } else {
            double timeChaining = 0.0, timeBacktrack = 0.0;
            chains = ChainHits(groupHits, chainMaxSkip, chainMaxPredecessors, seedJoinDist,
                               chainBandwidth, minNumSeeds, minCoveredBases, minDPScore,
                               timeChaining, timeBacktrack, ssChain);

            const double timeChainHits = ttPartial.GetMicrosecs(true);
            LogTicTocAdd("map-L2-chain-03-chainhits", ttPartial, retTimings);
            LogTicTocAdd("map-L3-chain-01-chaining", timeChaining, retTimings);
            LogTicTocAdd("map-L3-chain-02-backtrack", timeBacktrack, retTimings);
            LogTicTocAdd("map-L2-total-chainhits", timeChainHits, retTimings);
            LogTicTocAdd("map-L3-total-chaining", timeChaining, retTimings);
            LogTicTocAdd("map-L3-total-backtrack", timeBacktrack, retTimings);
#ifdef PANCAKE_MAP_CLR_DEBUG_2
            std::cerr << "  - not using LIS.\n";
#endif
        }

#ifdef PANCAKE_MAP_CLR_DEBUG_2
        std::cerr << "  - Adding chains, chains.size() = " << chains.size() << "\n";
#endif

        // Accumulate chains and their mapped regions.
        for (size_t i = 0; i < chains.size(); ++i) {
#ifdef PANCAKE_MAP_CLR_DEBUG_2
            std::cerr << "    [chain i = " << i << "]\n";
#endif

            const auto& chain = chains[i];
            int32_t numHitsInChain = chain.hits.size();
            // Filter.
            if (numHitsInChain < minNumSeeds) {
#ifdef PANCAKE_MAP_CLR_DEBUG_2
                std::cerr << "      -> Skipping chain, numHitsInChain < minNumSeeds\n";
#endif
                continue;
            }
            // Create a new chained region.
            auto ovl = MakeOverlap(chain.hits, queryId, queryLen, targetSeqs, 0, numHitsInChain, 0,
                                   numHitsInChain - 1);
            auto chainedRegion = std::make_unique<ChainedRegion>();
            chainedRegion->chain = std::move(chain);
            chainedRegion->mapping = std::move(ovl);
            chainedRegion->isMockedMapping = false;
            chainedRegion->isMockedAlignment = false;
            allChainedRegions.emplace_back(std::move(chainedRegion));
        }
        LogTicTocAdd("map-L2-chain-04-makeovl", ttPartial, retTimings);
    }
    return allChainedRegions;
}

}  // namespace Pancake
}  // namespace PacBio
