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

bool SeedIndex::CollectHits(const std::vector<PacBio::Pancake::SeedRaw>& querySeeds,
                            int32_t queryLen, std::vector<SeedHit>& hits, int64_t freqCutoff) const
{
    return CollectHits(querySeeds.data(), querySeeds.size(), queryLen, hits, freqCutoff);
}

bool SeedIndex::CollectHits(const PacBio::Pancake::SeedRaw* querySeeds, int64_t querySeedsSize,
                            int32_t queryLen, std::vector<SeedHit>& hits, int64_t freqCutoff) const
{
    return PacBio::Pancake::CollectSeedHits<SeedHashType>(hits, querySeeds, querySeedsSize,
                                                          queryLen, hash_, seeds_.data(),
                                                          seeds_.size(), freqCutoff);
}

}  // namespace Pancake
}  // namespace PacBio
