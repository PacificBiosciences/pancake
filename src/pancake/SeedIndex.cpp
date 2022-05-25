// Authors: Ivan Sovic

#include <pancake/SeedIndex.hpp>

#include <pancake/Minimizers.hpp>
#include <pancake/Seed.hpp>
#include <pancake/util/TicToc.hpp>

#include <pancake/third-party/kxsort/kxsort.h>
#include <pbcopper/logging/Logging.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>

namespace PacBio {
namespace Pancake {

SeedIndex::SeedIndex(std::vector<PacBio::Pancake::SeedRaw>&& seeds)
    : seeds_(std::move(seeds)), minSeedSpan_(0), maxSeedSpan_(0), avgSeedSpan_(0.0)
{
#ifdef SEED_INDEX_USING_DENSEHASH
    hash_.set_empty_key(
        SEED_INDEX_EMPTY_HASH_KEY);  // Densehash requires this to be defined on top.
#endif

    BuildHash_(seeds_, hash_, minSeedSpan_, maxSeedSpan_, avgSeedSpan_);
}

SeedIndex::~SeedIndex() = default;

void SeedIndex::BuildHash_(std::vector<PacBio::Pancake::SeedRaw>& seeds, SeedHashType& retHash,
                           int32_t& retMinSeedSpan, int32_t& retMaxSeedSpan, double& retAvgSeedSpan)
{
    // Clear the return values.
    retMinSeedSpan = 0;
    retMaxSeedSpan = 0;
    retAvgSeedSpan = 0.0;
    retHash.clear();

    // Stop early.
    if (seeds.empty()) {
        return;
    }

    // Sort by key first.
    kx::radix_sort(seeds.begin(), seeds.end());

#ifdef SEED_INDEX_USING_DENSEHASH
    retHash.resize(seeds.size());
#elif defined SEED_INDEX_USING_SPARSEHASH
    retHash.resize(seeds.size());
#elif defined SEED_INDEX_USING_UNORDERED_MAP
    retHash.reserve(seeds.size());
#endif

    // Fill out the hash table.
    int64_t start = 0;
    int64_t end = 0;
    uint64_t prevKey = PacBio::Pancake::Seed::DecodeKey(seeds[0]);
    retMinSeedSpan = seeds.empty() ? 0 : PacBio::Pancake::Seed::DecodeSpan(seeds[0]);
    retMaxSeedSpan = retMinSeedSpan;
    retAvgSeedSpan = 0.0;
    for (size_t i = 0; i < seeds.size(); ++i) {
        const uint64_t key = PacBio::Pancake::Seed::DecodeKey(seeds[i]);
        const int32_t span = PacBio::Pancake::Seed::DecodeSpan(seeds[i]);
        retMinSeedSpan = (span < retMinSeedSpan) ? span : retMinSeedSpan;
        retMaxSeedSpan = (span > retMaxSeedSpan) ? span : retMaxSeedSpan;
        retAvgSeedSpan += static_cast<double>(span);
        if (key != prevKey) {
            retHash[prevKey] = std::make_pair(start, end);
            start = i;
            end = i;
        }
        ++end;
        prevKey = key;
    }
    if (end > start) {
        retHash[prevKey] = std::make_pair(start, end);
    }
    retAvgSeedSpan = seeds.empty() ? 0.0 : (retAvgSeedSpan / static_cast<double>(seeds.size()));
}

void SeedIndex::ComputeFrequencyStats(double percentileCutoff, int64_t& retFreqMax,
                                      double& retFreqAvg, double& retFreqMedian,
                                      int64_t& retFreqCutoff) const
{
    retFreqMax = 0;
    retFreqAvg = 0.0;
    retFreqMedian = 0.0;
    retFreqCutoff = 0;

    // Sanity check.
    if (percentileCutoff < 0.0 || percentileCutoff > 1.0) {
        std::ostringstream oss;
        oss << "Invalid percentileCutoff value, should be in range [0.0, 1.0] but provided value = "
            << percentileCutoff;
        throw std::runtime_error(oss.str());
    }

    // Empty input.
    if (hash_.empty()) {
        return;
    }

    // Collect all frequencies as a vector.
    std::vector<int64_t> freqs;
    double sumFreqs = 0.0;
    int64_t numValidKeys = 0;
    freqs.reserve(hash_.size());
    for (auto it = hash_.begin(); it != hash_.end(); ++it) {
        int64_t start = std::get<0>(it->second);
        int64_t end = std::get<1>(it->second);
        int64_t span = end - start;
        if (span == 0) {
            continue;
        }
        freqs.emplace_back(span);
        sumFreqs += static_cast<double>(span);
        ++numValidKeys;
    }

    // Sanity check that there actually are enough valid keys in the hash.
    if (numValidKeys <= 0) {
        throw std::runtime_error("Invalid number of valid keys! numValidKeys = " +
                                 std::to_string(numValidKeys));
    }

    // Sort the vector for percentile calculation.
    kx::radix_sort(freqs.begin(), freqs.end());

    // Find the percentile cutoff ID.
    const double numKeysDouble = numValidKeys;
    const int64_t cutoffId = std::max(0.0, std::floor(numKeysDouble * (1.0 - percentileCutoff)));

    // Return values.
    retFreqMax = freqs.back();
    retFreqCutoff = (cutoffId >= numValidKeys) ? (freqs[numValidKeys - 1] + 1) : freqs[cutoffId];
    retFreqAvg = sumFreqs / numKeysDouble;
    retFreqMedian = (static_cast<double>(freqs[numValidKeys / 2]) +
                     static_cast<double>(freqs[(numValidKeys - 1) / 2])) /
                    2.0;
}

int64_t SeedIndex::GetSeeds(uint64_t key, std::vector<PacBio::Pancake::SeedRaw>& seeds) const
{
    seeds.clear();
    const auto it = hash_.find(key);
    if (it == hash_.end()) {
        return 0;
    }
    const int64_t start = std::get<0>(it->second);
    const int64_t end = std::get<1>(it->second);
    seeds.insert(seeds.end(), seeds_.begin() + start, seeds_.begin() + end);
    return (end - start);
}

void SeedIndex::CollectHits(const std::vector<PacBio::Pancake::SeedRaw>& querySeeds,
                            int32_t queryLen, std::vector<SeedHit>& hits, int64_t freqCutoff) const
{
    CollectHits(querySeeds.data(), querySeeds.size(), queryLen, hits, freqCutoff);
}

void SeedIndex::CollectHits(const PacBio::Pancake::SeedRaw* querySeeds, int64_t querySeedsSize,
                            int32_t queryLen, std::vector<SeedHit>& hits, int64_t freqCutoff) const
{
    PacBio::Pancake::CollectSeedHits(hits, {querySeeds, static_cast<size_t>(querySeedsSize)},
                                     queryLen, hash_, seeds_, freqCutoff);
}

std::vector<std::pair<int64_t, int64_t>> ComputeSeedHitHistogram(
    const std::span<const PacBio::Pancake::SeedRaw> querySeeds, const SeedHashType& hash)
{
    std::unordered_map<int64_t, int64_t> hist;
    for (int64_t seedId = 0; seedId < std::ssize(querySeeds); ++seedId) {
        const auto& querySeed = querySeeds[seedId];
        auto decodedQuery = PacBio::Pancake::Seed(querySeed);
        auto it = hash.find(decodedQuery.key);
        if (it != hash.end()) {
            int64_t start = std::get<0>(it->second);
            int64_t end = std::get<1>(it->second);
            hist[end - start] += 1;
        } else {
            hist[0] += 1;
        }
    }

    std::vector<std::pair<int64_t, int64_t>> histVec(hist.begin(), hist.end());
    std::sort(histVec.begin(), histVec.end());

    return histVec;
}

int64_t ComputeOccurrenceThreshold(const std::vector<std::pair<int64_t, int64_t>>& seedHitHistogram,
                                   const int64_t seedOccurrenceMin, const int64_t seedOccurrenceMax,
                                   const int64_t seedOccurrenceMaxMemory,
                                   const int64_t seedOccurrenceUserSpecified,
                                   const bool debugVerbose)
{
    // Compute the memory-bound maximum occurrence from the histogram, if requested.
    int64_t occThresholdMemMax = std::numeric_limits<int64_t>::max();
    if (seedOccurrenceMaxMemory > 0) {
        // We can fit at most this many seed hits to satisfy the memory requirements.
        const int64_t maxHitsToFit = std::ceil(static_cast<double>(seedOccurrenceMaxMemory) /
                                               static_cast<double>(sizeof(SeedHit)));

        // Loop through all seed hits until we fill the bucket.
        int64_t totalHits = 0;
        for (size_t i = 0; i < seedHitHistogram.size(); ++i) {
            const int64_t currentBinSize = seedHitHistogram[i].first * seedHitHistogram[i].second;
            // std::cerr << "[i = " << i << " / " << seedHitHistogram.size() << "] hits = " << seedHitHistogram[i].first << ", numSeeds = "
            // << seedHitHistogram[i].second << ", currentBinSize = " << currentBinSize << ", totalHits = " << totalHits
            // << ", maxHitsToFit = " << maxHitsToFit << "\n";
            if ((totalHits + currentBinSize) > maxHitsToFit) {
                break;
            }
            totalHits += currentBinSize;
            occThresholdMemMax = seedHitHistogram[i].first + 1;
        }
    }

    const int64_t occThresholdMax =
        seedOccurrenceMax > 0 ? seedOccurrenceMax : std::numeric_limits<int64_t>::max();
    const int64_t occPercentileCutoff = seedOccurrenceUserSpecified > 0
                                            ? seedOccurrenceUserSpecified
                                            : std::numeric_limits<int64_t>::max();

    const int64_t occThreshold =
        std::max(seedOccurrenceMin,
                 std::min(std::min(occThresholdMax, occThresholdMemMax), occPercentileCutoff));

    if (debugVerbose) {
        std::cerr << "Seed hit occurrence threshold: occThreshold = " << occThreshold
                  << " (occUserSpec = " << seedOccurrenceUserSpecified
                  << ", occMin = " << seedOccurrenceMin << ", occMax = " << seedOccurrenceMax
                  << ", occMemMax = " << occThresholdMemMax
                  << ", maxMemoryInBytes = " << seedOccurrenceMaxMemory << ")"
                  << "\n";
    }

    return occThreshold;
}

void CollectSeedHits(std::vector<SeedHit>& hits,
                     const std::span<const PacBio::Pancake::SeedRaw> querySeeds,
                     const int32_t queryLen, const SeedHashType& hash,
                     const std::span<const PacBio::Pancake::SeedRaw> targetSeeds,
                     const int64_t freqCutoff)
{
    hits.clear();

    for (int64_t seedId = 0; seedId < std::ssize(querySeeds); ++seedId) {
        const auto& querySeed = querySeeds[seedId];
        const auto decodedQuery = PacBio::Pancake::Seed(querySeed);
        const auto it = hash.find(decodedQuery.key);

        if (it != hash.end()) {
            const int64_t start = std::get<0>(it->second);
            const int64_t end = std::get<1>(it->second);

            // Skip very frequent seeds.
            if (freqCutoff > 0 && (end - start) > freqCutoff) {
                continue;
            }

            for (int64_t i = start; i < end; ++i) {
                const auto decodedTarget = PacBio::Pancake::Seed(targetSeeds[i]);
                const int32_t targetPos = decodedTarget.pos;
                const int32_t querySpan = decodedQuery.span;
                const int32_t targetSpan = decodedTarget.span;
                bool isRev = false;
                int32_t queryPos = decodedQuery.pos;

                if (decodedQuery.IsRev() != decodedTarget.IsRev()) {
                    isRev = true;
                    // End pos in fwd is start pos in rev.
                    queryPos = queryLen - (decodedQuery.pos + querySpan);
                }

                const SeedHit hit{static_cast<int32_t>(decodedTarget.seqID),
                                  isRev,
                                  targetPos,
                                  queryPos,
                                  targetSpan,
                                  querySpan,
                                  0};

                hits.emplace_back(hit);
            }
        }
    }
}

}  // namespace Pancake
}  // namespace PacBio
