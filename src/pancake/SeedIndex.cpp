// Authors: Ivan Sovic

#include <lib/kxsort/kxsort.h>
#include <pacbio/pancake/Minimizers.h>
#include <pacbio/pancake/Seed.h>
#include <pacbio/pancake/SeedIndex.h>
#include <pacbio/util/TicToc.h>
#include <pbcopper/logging/Logging.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>

namespace PacBio {
namespace Pancake {

SeedIndex::SeedIndex(std::shared_ptr<PacBio::Pancake::SeedDBIndexCache>& seedDBCache,
                     std::vector<PacBio::Pancake::SeedDB::SeedRaw>&& seeds)
    : seeds_(std::move(seeds)), seedParams_(seedDBCache->seedParams)
{
#ifdef SEED_INDEX_USING_DENSEHASH
    hash_.set_empty_key(
        SEED_INDEX_EMPTY_HASH_KEY);  // Densehash requires this to be defined on top.
#endif

    sequenceLengths_.resize(seedDBCache->seedLines.size());
    for (size_t i = 0; i < seedDBCache->seedLines.size(); ++i) {
        sequenceLengths_[i] = seedDBCache->seedLines[i].numBases;
    }

    BuildHash_();
}

SeedIndex::SeedIndex(const PacBio::Pancake::SeedDB::SeedDBParameters& seedParams,
                     const std::vector<int32_t>& sequenceLengths,
                     std::vector<PacBio::Pancake::SeedDB::SeedRaw>&& seeds)
    : seeds_(std::move(seeds)), seedParams_(seedParams), sequenceLengths_(sequenceLengths)
{
#ifdef SEED_INDEX_USING_DENSEHASH
    hash_.set_empty_key(
        SEED_INDEX_EMPTY_HASH_KEY);  // Densehash requires this to be defined on top.
#endif

    BuildHash_();
}

SeedIndex::~SeedIndex() = default;

void SeedIndex::BuildHash_()
{
    // Stop early.
    if (seeds_.empty()) {
        return;
    }

    // Sort by key first.
    TicToc ttSort;
    kx::radix_sort(seeds_.begin(), seeds_.end());
    ttSort.Stop();

    // Clear the storage for the hash.
    hash_.clear();

#ifdef SEED_INDEX_USING_DENSEHASH
    hash_.resize(seeds_.size());
#elif defined SEED_INDEX_USING_SPARSEHASH
    hash_.resize(seeds_.size());
#elif defined SEED_INDEX_USING_UNORDERED_MAP
    hash_.reserve(seeds_.size());
#endif

    // Fill out the hash table.
    int64_t start = 0;
    int64_t end = 0;
    uint64_t prevKey = PacBio::Pancake::SeedDB::Seed::DecodeKey(seeds_[0]);
    for (size_t i = 0; i < seeds_.size(); ++i) {
        uint64_t key = PacBio::Pancake::SeedDB::Seed::DecodeKey(seeds_[i]);
        if (key == prevKey) {
            ++end;
        } else {
            hash_[prevKey] = std::make_pair(start, end);
            start = i;
            end = i + 1;
        }
        prevKey = key;
    }
    if (end > start) {
        hash_[prevKey] = std::make_pair(start, end);
    }
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

int64_t SeedIndex::GetSeeds(uint64_t key,
                            std::vector<PacBio::Pancake::SeedDB::SeedRaw>& seeds) const
{
    seeds.clear();
    auto it = hash_.find(key);
    if (it == hash_.end()) {
        return 0;
    }
    int64_t start = std::get<0>(it->second);
    int64_t end = std::get<1>(it->second);
    seeds.insert(seeds.end(), seeds_.begin() + start, seeds_.begin() + end);
    return (end - start);
}

bool SeedIndex::CollectHits(const std::vector<PacBio::Pancake::SeedDB::SeedRaw>& querySeeds,
                            int32_t queryLen, std::vector<SeedHit>& hits, int64_t freqCutoff) const
{
    return CollectHits(&querySeeds[0], querySeeds.size(), queryLen, hits, freqCutoff);
}

bool SeedIndex::CollectHits(const PacBio::Pancake::SeedDB::SeedRaw* querySeeds,
                            int64_t querySeedsSize, int32_t queryLen, std::vector<SeedHit>& hits,
                            int64_t freqCutoff) const
{
    return PacBio::Pancake::SeedDB::CollectSeedHits<SeedHashType>(
        hits, querySeeds, querySeedsSize, queryLen, hash_, &seeds_[0], seeds_.size(),
        sequenceLengths_, seedParams_.KmerSize, seedParams_.Spacing, freqCutoff);
}

}  // namespace Pancake
}  // namespace PacBio
