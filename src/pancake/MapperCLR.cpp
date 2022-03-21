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
#include <pancake/util/RunLengthEncoding.hpp>
#include <pancake/util/TicToc.hpp>
#include <pancake/util/Util.hpp>

#include <pancake/third-party/kxsort/kxsort.h>
#include <pancake/third-party/pdqsort/pdqsort.h>
#include <pbcopper/logging/Logging.h>
#include <pbcopper/third-party/edlib.h>

#include <algorithm>
#include <array>
#include <cstring>
#include <iostream>
#include <sstream>
#include <tuple>

// #define PANCAKE_MAP_CLR_DEBUG
// #define PANCAKE_MAP_CLR_DEBUG_2
// #define PANCAKE_WRITE_SCATTERPLOT
// #define PANCAKE_MAP_CLR_DEBUG_ALIGN

// #define PANCAKE_MAP_CLR_DEBUG_PRINT_CHAINED_REGIONS
// #define PANCAKE_MAP_CLR_DEBUG_WRITE_SEED_HITS_TO_FILE

#include <pancake/util/DebugTools.hpp>

#if defined(PANCAKE_MAP_CLR_DEBUG) || defined(PANCAKE_MAP_CLR_DEBUG_2)
#include <pbcopper/utility/MemoryConsumption.h>
#include <iomanip>
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

std::vector<MapperBaseResult> MapperCLR::MapAndAlign(const std::vector<std::string>& targetSeqs,
                                                     const std::vector<std::string>& querySeqs)
{
    std::vector<FastaSequenceCached> targetSeqsCached;
    for (int32_t i = 0; i < static_cast<int32_t>(targetSeqs.size()); ++i) {
        targetSeqsCached.emplace_back(
            FastaSequenceCached(std::to_string(i), targetSeqs[i].c_str(), targetSeqs[i].size(), i));
    }

    std::vector<FastaSequenceCached> querySeqsCached;
    for (int32_t i = 0; i < static_cast<int32_t>(querySeqs.size()); ++i) {
        querySeqsCached.emplace_back(
            FastaSequenceCached(std::to_string(i), querySeqs[i].c_str(), querySeqs[i].size(), i));
    }

    return MapAndAlign(targetSeqsCached, querySeqsCached);
}

std::vector<MapperBaseResult> MapperCLR::MapAndAlign(
    const std::vector<FastaSequenceCached>& targetSeqs,
    const std::vector<FastaSequenceCached>& querySeqs)
{
    PacBio::Pancake::FastaSequenceCachedStore targetSeqsStore(targetSeqs);
    PacBio::Pancake::FastaSequenceCachedStore querySeqsStore(querySeqs);
    return MapAndAlign(targetSeqsStore, querySeqsStore);
}

std::vector<MapperBaseResult> MapperCLR::MapAndAlign(const FastaSequenceCachedStore& targetSeqs,
                                                     const FastaSequenceCachedStore& querySeqs)
{
    return WrapBuildIndexMapAndAlignWithFallback_(targetSeqs, querySeqs, settings_, ssChain_,
                                                  ssSeedHits_, alignerGlobal_, alignerExt_);
}

MapperBaseResult MapperCLR::MapAndAlignSingleQuery(
    const FastaSequenceCachedStore& targetSeqs, const PacBio::Pancake::SeedIndex& index,
    const FastaSequenceCached& querySeq, const std::vector<PacBio::Pancake::Int128t>& querySeeds,
    const int32_t queryId, int64_t freqCutoff)
{
    return WrapMapAndAlign_(targetSeqs, index, querySeq, querySeeds, queryId, freqCutoff, settings_,
                            ssChain_, ssSeedHits_, alignerGlobal_, alignerExt_);
}

MapperBaseResult MapperCLR::Map(const FastaSequenceCachedStore& targetSeqs,
                                const PacBio::Pancake::SeedIndex& index,
                                const std::vector<PacBio::Pancake::Int128t>& querySeeds,
                                const int32_t queryLen, const int32_t queryId, int64_t freqCutoff)
{
    return Map_(targetSeqs, index, querySeeds, queryLen, queryId, settings_, freqCutoff, ssChain_,
                ssSeedHits_);
}

MapperBaseResult MapperCLR::Align(const FastaSequenceCachedStore& targetSeqs,
                                  const FastaSequenceCached& querySeq,
                                  const MapperBaseResult& mappingResult)
{
    return Align_(targetSeqs, querySeq, mappingResult, settings_, alignerGlobal_, alignerExt_);
}

void DebugWriteChainedRegion(const std::vector<std::unique_ptr<ChainedRegion>>& allChainedRegions,
#if defined(PANCAKE_MAP_CLR_DEBUG) || defined(PANCAKE_MAP_CLR_DEBUG_2) || \
    defined(PANCAKE_MAP_CLR_DEBUG_WRITE_SEED_HITS_TO_FILE)
                             const std::string& descriptor, int32_t queryId, int32_t queryLen,
#else
                             const std::string& /*descriptor*/, int32_t /*queryId*/,
                             int32_t /*queryLen*/,
#endif
                             const bool debugVerboseChains = false)
{

// Special debug output, only if macros are defined.
#if defined(PANCAKE_MAP_CLR_DEBUG) || defined(PANCAKE_MAP_CLR_DEBUG_2)
    const double peakRssGb =
        PacBio::Utility::MemoryConsumption::PeakRss() / 1024.0 / 1024.0 / 1024.0;
    PBLOG_DEBUG << "After: '" << descriptor << "'. Chains: " << allChainedRegions.size()
                << ". Peak RSS: " << std::fixed << std::setprecision(3) << peakRssGb;
#endif
#ifdef PANCAKE_MAP_CLR_DEBUG_WRITE_SEED_HITS_TO_FILE
    // Write seed hits.
    for (size_t i = 0; i < allChainedRegions.size(); ++i) {
        auto& region = allChainedRegions[i];
        if (region->priority > 1) {
            continue;
        }
        WriteSeedHits("temp-debug/hits-q" + std::to_string(queryId) + "-" + descriptor + ".csv",
                      region->chain.hits, 0, region->chain.hits.size(), i,
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
            if (region->priority > 1) {
                continue;
            }
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
    for (int32_t i = 0; i < static_cast<int32_t>(querySeqs.records().size()); ++i) {
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
        if (rv)
            throw std::runtime_error("Generating minimizers failed for the query sequence i = " +
                                     std::to_string(i) + ", id = " + std::to_string(queryId));

        MapperBaseResult queryResults =
            WrapMapAndAlign_(targetSeqs, *seedIndex, query, querySeeds, queryId, freqCutoff,
                             settings, ssChain, ssSeedHits, alignerGlobal, alignerExt);

        if (queryResults.mappings.empty() && seedIndexFallback != nullptr) {
            rv = GenerateMinimizers(
                querySeeds, seq, 0, queryId, settings.map.seedParamsFallback.KmerSize,
                settings.map.seedParamsFallback.MinimizerWindow,
                settings.map.seedParamsFallback.Spacing, settings.map.seedParamsFallback.UseRC,
                settings.map.seedParamsFallback.UseHPCForSeedsOnly);
            if (rv)
                throw std::runtime_error(
                    "Generating minimizers failed for the query sequence, id = " +
                    std::to_string(queryId));

            queryResults = WrapMapAndAlign_(targetSeqs, *seedIndexFallback, query, querySeeds,
                                            queryId, freqCutoffFallback, settings, ssChain,
                                            ssSeedHits, alignerGlobal, alignerExt);
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
    const FastaSequenceCached& querySeq, const std::vector<PacBio::Pancake::Int128t>& querySeeds,
    const int32_t queryId, int64_t freqCutoff, const MapperCLRSettings& settings,
    std::shared_ptr<ChainingScratchSpace> ssChain, std::vector<SeedHit>& ssSeedHits,
    AlignerBasePtr& alignerGlobal, AlignerBasePtr& alignerExt)
{
    const int32_t queryLen = querySeq.size();

    // Map the query.
    auto result = Map_(targetSeqs, index, querySeeds, queryLen, queryId, settings, freqCutoff,
                       ssChain, ssSeedHits);

    // Align if needed.
    if (settings.align.align) {
        result = Align_(targetSeqs, querySeq, result, settings, alignerGlobal, alignerExt);
    }

#ifdef PANCAKE_MAP_CLR_DEBUG_PRINT_CHAINED_REGIONS
    DebugWriteChainedRegion(result.mappings, "9-result-final", queryId, queryLen, true);
#else
    DebugWriteChainedRegion(result.mappings, "9-result-final", queryId, queryLen, false);
#endif

    return result;
}

MapperBaseResult MapperCLR::Map_(const FastaSequenceCachedStore& targetSeqs,
                                 const PacBio::Pancake::SeedIndex& index,
                                 const std::vector<PacBio::Pancake::Int128t>& querySeeds,
                                 const int32_t queryLen, const int32_t queryId,
                                 const MapperCLRSettings& settings, int64_t freqCutoff,
                                 std::shared_ptr<ChainingScratchSpace> ssChain,
                                 std::vector<SeedHit>& ssSeedHits)
{
    TicToc ttMapAll;
    TicToc ttPartial;

    MapperBaseResult result;

    // Skip short queries.
    if (queryLen < settings.map.minQueryLen) {
        LogTicToc("map-L1-all", ttMapAll, result.time);
        return result;
    }

#if defined(PANCAKE_MAP_CLR_DEBUG) || defined(PANCAKE_MAP_CLR_DEBUG_2)
    {
        const double peakRssGb =
            PacBio::Utility::MemoryConsumption::PeakRss() / 1024.0 / 1024.0 / 1024.0;
        PBLOG_DEBUG << "Starting function: " << std::string(__FUNCTION__) << "."
                    << " Peak RSS: " << std::fixed << std::setprecision(3) << peakRssGb;

        const std::vector<std::pair<int64_t, int64_t>> debugSeedHist =
            PacBio::Pancake::ComputeSeedHitHistogram(querySeeds, index.GetHash());

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

    // Compute the seed hit histogram, but only if needed.
    std::vector<std::pair<int64_t, int64_t>> seedHitHistogram;
    if (settings.map.seedOccurrenceMaxMemory > 0) {
        seedHitHistogram = PacBio::Pancake::ComputeSeedHitHistogram(querySeeds, index.GetHash());
    }

    const int64_t occThreshold = ComputeOccurrenceThreshold(
        seedHitHistogram, settings.map.seedOccurrenceMin, settings.map.seedOccurrenceMax,
        settings.map.seedOccurrenceMaxMemory, freqCutoff, debugVerboseOccurrenceThreshold);

    LogTicToc("map-L1-01-occthresh", ttPartial, result.time);

    // Collect seed hits.
    auto& hits = ssSeedHits;
    index.CollectHits(querySeeds.data(), querySeeds.size(), queryLen, hits, occThreshold);
    LogTicToc("map-L1-02-collect", ttPartial, result.time);

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
    LogTicToc("map-L1-03-perfaln", ttPartial, result.time);

    // Filter symmetric and self hits.
    if (settings.map.skipSymmetricOverlaps ||
        settings.map.selfHitPolicy != MapperSelfHitPolicy::DEFAULT) {
        std::vector<SeedHit> newHits = FilterSymmetricAndSelfHits_(
            hits, queryId, settings.map.selfHitPolicy != MapperSelfHitPolicy::DEFAULT,
            settings.map.skipSymmetricOverlaps);
        std::swap(hits, newHits);
    }
    LogTicToc("map-L1-04-skipsym", ttPartial, result.time);

#if defined(PANCAKE_MAP_CLR_DEBUG) || defined(PANCAKE_MAP_CLR_DEBUG_2)
    {
        const double peakRssGb =
            PacBio::Utility::MemoryConsumption::PeakRss() / 1024.0 / 1024.0 / 1024.0;
        PBLOG_DEBUG << "Sorted hits."
                    << " Peak RSS: " << std::fixed << std::setprecision(3) << peakRssGb;
    }
#endif

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
        std::vector<std::pair<int32_t, int64_t>> seedHitsPacked(hits.size());
        for (size_t i = 0; i < hits.size(); ++i) {
            const auto& sh = hits[i];
            seedHitsPacked[i] = std::make_pair(
                ((sh.targetId << 1) | sh.targetRev),
                (static_cast<uint64_t>(sh.targetPos - sh.queryPos) << 32 | (i & 0x0FFFFFFFF)));
        }
        const double timeSortPrepare = ttPartial.GetMicrosecs(true);
        LogTicTocAdd("map-L3-total-sort", timeSortPrepare, result.time);
        LogTicToc("map-L1-05-sort-prepare", ttPartial, result.time);

        pdqsort(seedHitsPacked.begin(), seedHitsPacked.end());

        std::vector<PacBio::Pancake::SeedHit> sortedHits(hits.size());
        for (size_t i = 0; i < seedHitsPacked.size(); ++i) {
            sortedHits[i] = hits[seedHitsPacked[i].second & 0x0FFFFFFFF];
        }
        std::swap(hits, sortedHits);
        const double timeSort = ttPartial.GetMicrosecs(true);
        LogTicTocAdd("map-L3-total-sort", timeSort, result.time);
        LogTicToc("map-L1-05-sort", ttPartial, result.time);
    }

    // Group seed hits by diagonal.
    auto groups = DiagonalGroup(hits, settings.map.chainBandwidth, true);
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
        targetSeqs, hits, groups, queryId, queryLen, settings.map.chainMaxSkip,
        settings.map.chainMaxPredecessors, settings.map.seedJoinDist, settings.map.chainBandwidth,
        settings.map.minNumSeeds, settings.map.minCoveredBases, settings.map.minDPScore,
        settings.map.useLIS, ssChain, result.time);
    DebugWriteChainedRegion(allChainedRegions, "1-chain-and-make-overlap", queryId, queryLen);
    LogTicToc("map-L1-07-chain", ttPartial, result.time);

    // Take the remaining regions, merge all seed hits, and rechain.
    // Needed because diagonal chaining was greedy and a wide window could have split
    // otherwise good chains. On the other hand, diagonal binning was needed for speed
    // in low-complexity regions.
    allChainedRegions = ReChainSeedHits_(
        allChainedRegions, targetSeqs, queryId, queryLen, settings.map.chainMaxSkip,
        settings.map.chainMaxPredecessors, settings.map.seedJoinDist, settings.map.chainBandwidth,
        settings.map.minNumSeeds, settings.map.minCoveredBases, settings.map.minDPScore, ssChain,
        result.time);
    DebugWriteChainedRegion(allChainedRegions, "2-rechain-hits", queryId, queryLen);
    LogTicToc("map-L1-08-rechain", ttPartial, result.time);

    // Sort all chains in descending order of the number of hits.
    pdqsort(allChainedRegions.begin(), allChainedRegions.end(),
            [](const auto& a, const auto& b) { return a->chain.score > b->chain.score; });
    LogTicToc("map-L1-09-sortchained", ttPartial, result.time);

    // Secondary/supplementary flagging.
    WrapFlagSecondaryAndSupplementary(
        allChainedRegions, settings.map.secondaryAllowedOverlapFractionQuery,
        settings.map.secondaryAllowedOverlapFractionTarget, settings.map.secondaryMinScoreFraction);
    DebugWriteChainedRegion(allChainedRegions, "3-wrap-flag-secondary-suppl-1", queryId, queryLen);
    LogTicToc("map-L1-10-secondary", ttPartial, result.time);

    // Merge long gaps.
    LongMergeChains_(allChainedRegions, settings.map.longMergeBandwidth);
    DebugWriteChainedRegion(allChainedRegions, "4-long-merge-chains", queryId, queryLen);
    LogTicToc("map-L1-11-longmerge", ttPartial, result.time);

    // Add an extra query alignment only if needed (i.e. if the MapperSelfHitPolicy == PERFECT_ALIGNMENT
    // and the query has self-hits ("self" in term of queryId/targetId).
    if (addPerfectMapping) {
        const int32_t targetId = queryId;
        const int32_t targetLen = targetSeqs.GetSequence(targetId).size();
        std::unique_ptr<ChainedRegion> newChainedRegion =
            CreateMockedMapping_(queryId, queryLen, targetId, targetLen);
        allChainedRegions.emplace_back(std::move(newChainedRegion));
    }
    LogTicToc("map-L1-12-addperfect", ttPartial, result.time);

    // Again relabel, because some chains are longer now.
    WrapFlagSecondaryAndSupplementary(
        allChainedRegions, settings.map.secondaryAllowedOverlapFractionQuery,
        settings.map.secondaryAllowedOverlapFractionTarget, settings.map.secondaryMinScoreFraction);
    DebugWriteChainedRegion(allChainedRegions, "5-wrap-flag-secondary-suppl-2", queryId, queryLen);
    LogTicToc("map-L1-13-secondary", ttPartial, result.time);

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
        ChainedHits newChain =
            RefineBadEnds(region->chain, settings.align.alnParamsGlobal.alignBandwidth,
                          settings.map.minDPScore * 2);
        newChain = RefineChainedHits(newChain, 10, 40, settings.map.seedJoinDist / 2, 10);
        newChain = RefineChainedHits2(newChain, 30, settings.map.seedJoinDist / 2);

        std::swap(region->chain, newChain);

        if (region->chain.hits.empty()) {
            region->mapping = nullptr;
        } else {
            region->mapping =
                MakeOverlap(region->chain.hits, queryId, queryLen, targetSeqs, 0,
                            region->chain.hits.size(), 0, region->chain.hits.size() - 1);
        }
    }
    DebugWriteChainedRegion(allChainedRegions, "6-refining-seed-hits", queryId, queryLen);
    LogTicToc("map-L1-15-refineseeds", ttPartial, result.time);

    // Filter out the mappings.
    result.mappings = std::move(allChainedRegions);
    const int32_t numPrimary = CondenseMappings(result.mappings, settings.map.bestNSecondary);
    LogTicToc("map-L1-16-condense", ttPartial, result.time);

    // If this occurs, that means that a filtering stage removed the primary alignment for some reason.
    // This can happen in case self-hits are skipped in the overlapping use case (the self-hit is the
    // highest scoring one, and it would get marked for removal).
    // This reruns labeling to produce the next best primary.
    if (numPrimary == 0) {
        WrapFlagSecondaryAndSupplementary(allChainedRegions,
                                          settings.map.secondaryAllowedOverlapFractionQuery,
                                          settings.map.secondaryAllowedOverlapFractionTarget,
                                          settings.map.secondaryMinScoreFraction);
    }
    LogTicToc("map-L1-17-secondary", ttPartial, result.time);

    // Collect regions for alignment.
    for (size_t i = 0; i < result.mappings.size(); ++i) {
        if (result.mappings[i] == nullptr || result.mappings[i]->mapping == nullptr) {
            continue;
        }
        result.mappings[i]->regionsForAln = CollectAlignmentRegions_(
            *result.mappings[i], settings.map.minAlignmentSpan, settings.map.maxFlankExtensionDist,
            settings.map.flankExtensionFactor);
    }

    DebugWriteChainedRegion(result.mappings, "7-result-mappings", queryId, queryLen);
    LogTicToc("map-L1-18-collectalnregions", ttPartial, result.time);

#ifdef PANCAKE_MAP_CLR_DEBUG_2
    std::cerr << "All hits: hits.size() = " << hits.size() << "\n";
    std::cerr << "Diagonal groups: groups.size() = " << groups.size() << "\n";
    for (size_t i = 0; i < groups.size(); ++i) {
        const int32_t firstDiag = hits[groups[i].start].Diagonal();
        const int32_t lastDiag = hits[groups[i].end - 1].Diagonal();
        std::cerr << "[queryId = " << queryId << ", group " << i << "] start = " << groups[i].start
                  << ", end = " << groups[i].end << ", diagStart = " << firstDiag
                  << ", diagEnd = " << lastDiag << "\n";
    }

    // Write ALL seed hits.
    for (size_t i = 0; i < groups.size(); ++i) {
        const auto& g = groups[i];
        const int32_t targetId = hits[g.start].targetId;
        const int32_t targetLen = targetSeqs.GetSequence(targetId).size();
        WriteSeedHits(
            "temp-debug/hits-q" + std::to_string(queryId) + "-0-all-hits-diagonal-groupped.csv",
            hits, g.start, g.end, i, "query" + std::to_string(queryId), queryLen,
            "target" + std::to_string(targetId), targetLen, (i > 0));
    }
#endif

    LogTicToc("map-all", ttMapAll, result.time);

    return result;
}

MapperBaseResult MapperCLR::Align_(const FastaSequenceCachedStore& targetSeqs,
                                   const FastaSequenceCached& querySeq,
                                   const MapperBaseResult& mappingResult,
                                   const MapperCLRSettings& settings, AlignerBasePtr& alignerGlobal,
                                   AlignerBasePtr& alignerExt)
{
    TicToc ttAlignAll;
    TicToc ttPartial;

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
    const std::string querySeqRev =
        PacBio::Pancake::ReverseComplement(querySeq.c_str(), 0, querySeq.size());

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

    DebugWriteChainedRegion(alignedResult.mappings, "8-result-after-align", queryId, queryLen);
    LogTicToc("aln-L1-02-secondary", ttPartial, alignedResult.time);

    const int32_t numPrimary =
        CondenseMappings(alignedResult.mappings, settings.map.bestNSecondary);
    LogTicToc("aln-L1-03-condense", ttPartial, alignedResult.time);

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

    LogTicToc("aln-all", ttAlignAll, alignedResult.time);

    return alignedResult;
}

std::vector<SeedHit> MapperCLR::FilterSymmetricAndSelfHits_(const std::vector<SeedHit>& hits,
                                                            const int32_t queryId,
                                                            const bool skipSelfHits,
                                                            const bool skipSymmetricOverlaps)
{
    std::vector<SeedHit> newHits(hits.size());
    size_t newHitId = 0;
    for (const SeedHit& hit : hits) {
        if ((skipSelfHits && hit.targetId == queryId) ||
            (skipSymmetricOverlaps && hit.targetId > queryId)) {
            continue;
        }
        newHits[newHitId] = hit;
        ++newHitId;
    }
    newHits.resize(newHitId);
    return newHits;
}

std::vector<AlignmentRegion> MapperCLR::CollectAlignmentRegions_(const ChainedRegion& singleMapping,
                                                                 int32_t minAlignmentSpan,
                                                                 int32_t maxFlankExtensionDist,
                                                                 double flankExtensionFactor)
{
    if (singleMapping.mapping == nullptr) {
        return {};
    }

    const auto& ovl = singleMapping.mapping;
    const auto& chain = singleMapping.chain;
    const auto& sortedHits = chain.hits;

    // Sanity checks.
    if (ovl->Arev) {
        throw std::runtime_error("(CollectAlignmentRegions) The ovl->Arev should always be false!");
    }
    const int32_t Astart = ovl->Brev ? (ovl->Alen - ovl->Aend) : ovl->Astart;
    const int32_t Aend = ovl->Brev ? (ovl->Alen - ovl->Astart) : ovl->Aend;
    const int32_t Bstart = ovl->Brev ? (ovl->Blen - ovl->Bend) : ovl->Bstart;
    const int32_t Bend = ovl->Brev ? (ovl->Blen - ovl->Bstart) : ovl->Bend;

    if (Astart != sortedHits.front().queryPos ||
        Aend != (sortedHits.back().queryPos + sortedHits.back().querySpan) ||
        Bstart != sortedHits.front().targetPos ||
        Bend != (sortedHits.back().targetPos + sortedHits.back().targetSpan)) {
        std::ostringstream oss;
        oss << "(" << __FUNCTION__
            << ") Provided overlap coordinates do not match the first/last seed hit!"
            << " ovl: " << *ovl << "; sortedHits.front() = " << sortedHits.front()
            << "; sortedHits.back() = " << sortedHits.back() << "; Astart = " << Astart
            << ", Aend = " << Aend << ", Bstart = " << Bstart << ", Bend = " << Bend;
        throw std::runtime_error(oss.str());
    }
    if (ovl->Brev != sortedHits.front().targetRev || ovl->Brev != sortedHits.back().targetRev) {
        std::ostringstream oss;
        oss << "(" << __FUNCTION__
            << ") Strand of the provided overlap does not match the first/last "
               "seed hit."
            << " ovl->Brev = " << (ovl->Brev ? "true" : "false")
            << ", sortedHits.front().targetRev = "
            << (sortedHits.front().targetRev ? "true" : "false")
            << ", sortedHits.back().targetRev = "
            << (sortedHits.back().targetRev ? "true" : "false");
        throw std::runtime_error(oss.str());
    }

    // Extract the regions.
    std::vector<AlignmentRegion> regions =
        ExtractAlignmentRegions(chain.hits, ovl->Alen, ovl->Blen, ovl->Brev, minAlignmentSpan,
                                maxFlankExtensionDist, flankExtensionFactor);

    return regions;
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
            allChainedRegions.emplace_back(std::move(chainedRegion));
        }
        LogTicTocAdd("map-L2-chain-04-makeovl", ttPartial, retTimings);
    }
    return allChainedRegions;
}

void MapperCLR::LongMergeChains_(std::vector<std::unique_ptr<ChainedRegion>>& chainedRegions,
                                 const int32_t maxBandwidth)
{
    if (maxBandwidth < 0) {
        throw std::runtime_error(
            "(LongMergeChains_) maxBandwidth cannot be negative. maxBandwidth = " +
            std::to_string(maxBandwidth));
    }

    if (chainedRegions.empty()) {
        return;
    }

    std::vector<int32_t> candidates;
    for (size_t i = 0; i < chainedRegions.size(); ++i) {
        if (chainedRegions[i] == nullptr) {
            continue;
        }
        auto& cr = chainedRegions[i];
        if (cr->priority < 0 || cr->priority > 1) {
            continue;
        }
        candidates.emplace_back(i);
    }

    pdqsort(candidates.begin(), candidates.end(), [&](const auto& a, const auto& b) {
        const auto& ar = chainedRegions[a];
        const auto& br = chainedRegions[b];
        return std::tuple(ar->mapping->Bid, ar->mapping->Brev, ar->mapping->Astart,
                          ar->mapping->Bstart) < std::tuple(br->mapping->Bid, br->mapping->Brev,
                                                            br->mapping->Astart,
                                                            br->mapping->Bstart);
    });

    std::unordered_set<int32_t> doSort;

    int32_t lastId = candidates[0];
    for (int32_t i = 1; i < static_cast<int32_t>(candidates.size()); ++i) {
        const int32_t currId = candidates[i];
        auto& last = chainedRegions[lastId];
        const auto& curr = chainedRegions[currId];

        // Skip if there is an overlap (or if the target is wrong), to preserve colinearity.
        if (curr->mapping->Bid != last->mapping->Bid ||
            curr->mapping->Brev != last->mapping->Brev ||
            curr->mapping->Astart <= last->mapping->Aend ||
            curr->mapping->Bstart <= last->mapping->Bend) {
            lastId = currId;
            continue;
        }

        const int32_t gap = std::abs((curr->mapping->Bstart - last->mapping->Bend) -
                                     (curr->mapping->Astart - last->mapping->Aend));

        if (gap > maxBandwidth) {
            lastId = currId;
            continue;
        }

        last->chain.hits.insert(last->chain.hits.end(), curr->chain.hits.begin(),
                                curr->chain.hits.end());
        last->chain.score += curr->chain.score;
        last->chain.coveredBasesQuery += curr->chain.coveredBasesQuery;
        last->chain.coveredBasesTarget += curr->chain.coveredBasesTarget;

        last->mapping->Aend = curr->mapping->Aend;
        last->mapping->Bend = curr->mapping->Bend;
        last->mapping->Score += curr->mapping->Score;
        last->mapping->NumSeeds += curr->mapping->NumSeeds;

        chainedRegions[currId] = nullptr;

        doSort.emplace(lastId);
    }

    // Sort only the extended chains.
    for (const auto& i : doSort) {
        if (chainedRegions[i] == nullptr) {
            continue;
        }
        pdqsort(chainedRegions[i]->chain.hits.begin(), chainedRegions[i]->chain.hits.end());
    }

    // Remove the merged ones.
    std::vector<std::unique_ptr<ChainedRegion>> ret;
    for (size_t i = 0; i < chainedRegions.size(); ++i) {
        if (chainedRegions[i] == nullptr) {
            continue;
        }
        if (chainedRegions[i]->mapping == nullptr) {
            continue;
        }
        ret.emplace_back(std::move(chainedRegions[i]));
    }

    std::swap(chainedRegions, ret);
}

void WrapFlagSecondaryAndSupplementary(
    std::vector<std::unique_ptr<ChainedRegion>>& allChainedRegions,
    double secondaryAllowedOverlapFractionQuery, double secondaryAllowedOverlapFractionTarget,
    double secondaryMinScoreFraction)
{
    /*
     * Edits the objects in place.
    */
    // Copy the overlaps so we can satisfy the FlagSecondaryAndSupplementary API.
    std::vector<OverlapPtr> tmpOverlaps;
    std::vector<int32_t> inputOverlapToTmpOverlap(allChainedRegions.size(), -1);
    for (size_t i = 0; i < allChainedRegions.size(); ++i) {
        if (allChainedRegions[i] == nullptr || allChainedRegions[i]->mapping == nullptr) {
            continue;
        }
        inputOverlapToTmpOverlap[i] = tmpOverlaps.size();
        tmpOverlaps.emplace_back(CreateOverlap(allChainedRegions[i]->mapping));
    }
    // Flag the secondary and supplementary overlaps.
    std::vector<OverlapPriority> overlapPriorities = FlagSecondaryAndSupplementary(
        tmpOverlaps, secondaryAllowedOverlapFractionQuery, secondaryAllowedOverlapFractionTarget,
        secondaryMinScoreFraction);

    // Set the flags.
    for (size_t i = 0; i < allChainedRegions.size(); ++i) {
        if (allChainedRegions[i] == nullptr || allChainedRegions[i]->mapping == nullptr) {
            continue;
        }
        const int32_t tmpOverlapId = inputOverlapToTmpOverlap[i];
        if (tmpOverlapId < 0) {
            continue;
        }
        auto& cr = allChainedRegions[i];
        cr->mapping->IsSecondary = (overlapPriorities[tmpOverlapId].priority > 0);
        cr->mapping->IsSupplementary = overlapPriorities[tmpOverlapId].isSupplementary;
        cr->priority = overlapPriorities[tmpOverlapId].priority;
        cr->isSupplementary = overlapPriorities[tmpOverlapId].isSupplementary;
    }
}

int32_t CondenseMappings(std::vector<std::unique_ptr<ChainedRegion>>& mappings,
                         int32_t bestNSecondary)
{
    // Filter mappings.
    size_t numValid = 0;
    int32_t numPrimary = 0;
    int32_t numSelectedSecondary = 0;
    for (size_t i = 0; i < mappings.size(); ++i) {
        if (mappings[i] == nullptr) {
            continue;
        }
        const auto& region = mappings[i];
        if (region->mapping == nullptr) {
            continue;
        }
        if (region->priority > 1) {
            continue;
        }
        if (bestNSecondary >= 0 && region->priority == 1 &&
            numSelectedSecondary >= bestNSecondary) {
            continue;
        }
        if (region->priority == 0) {
            ++numPrimary;
        }

        if (bestNSecondary < 0 ||
            (region->priority == 1 && numSelectedSecondary < bestNSecondary)) {
            ++numSelectedSecondary;
        }
        if (i != numValid) {
            std::swap(mappings[i], mappings[numValid]);
        }
        ++numValid;
    }
    mappings.resize(numValid);
    return numPrimary;
}

std::unique_ptr<ChainedRegion> MapperCLR::CreateMockedMapping_(const int32_t queryId,
                                                               const int32_t queryLen,
                                                               const int32_t targetId,
                                                               const int32_t targetLen)
{
    if (queryLen != targetLen) {
        throw std::runtime_error(
            "Cannot mock a perfect mapping between the two sequences with same ID, the "
            "lengths are different. The sequences might be mislabelled. queryLen = " +
            std::to_string(queryLen) + ", targetLen = " + std::to_string(targetLen));
    }

    const int32_t numSeeds = queryLen;
    const int32_t alnScore = queryLen;
    const float identity = 0.0;
    const int32_t editDist = -1;
    const bool isRev = false;

    ChainedHits newChain{targetId,
                         false,
                         {
                             SeedHit(targetId, false, 0, 0, 0, 0, 0),
                             SeedHit(targetId, false, targetLen, queryLen, 0, 0, 0),
                         },
                         alnScore,
                         queryLen,
                         targetLen};

    OverlapPtr newOvl = CreateOverlap(queryId, targetId, alnScore, identity, isRev, 0, queryLen,
                                      queryLen, false, 0, targetLen, targetLen, editDist, numSeeds,
                                      OverlapType::Unknown, OverlapType::Unknown);

    auto newChainedRegion = std::make_unique<ChainedRegion>();
    newChainedRegion->chain = std::move(newChain);
    newChainedRegion->regionsForAln = {};  // This will be generated later.
    newChainedRegion->mapping = std::move(newOvl);
    newChainedRegion->priority = 0;
    newChainedRegion->isSupplementary = false;

    return newChainedRegion;
}

OverlapPtr CreateMockedAlignment(const OverlapPtr& ovl, const int32_t matchScore)
{

    if (ovl->Alen != ovl->Blen) {
        std::ostringstream oss;
        oss << "Cannot mock a perfect alignment between the two sequences with same ID, the "
               "lengths are different. The sequences might be mislabelled. Overlap: "
            << *ovl;
        throw std::runtime_error(oss.str());
    }

    const int32_t score = ovl->Alen * matchScore;
    const float identity = 1.0;
    const int32_t editDist = 0;
    const int32_t numSeeds = ovl->Alen;
    OverlapPtr newOvl = CreateOverlap(
        ovl->Aid, ovl->Bid, score, identity, false, 0, ovl->Alen, ovl->Alen, false, 0, ovl->Blen,
        ovl->Blen, editDist, numSeeds, OverlapType::Unknown, OverlapType::Unknown,
        // OverlapType::Contains, OverlapType::Contains,
        PacBio::BAM::Cigar(std::to_string(ovl->Alen) + "="), "", "", false, false, false);
    return newOvl;
}

}  // namespace Pancake
}  // namespace PacBio
