// Author: Ivan Sovic

#ifndef PANCAKE_OVERLAPHIFI_SEEDINDEX_H
#define PANCAKE_OVERLAPHIFI_SEEDINDEX_H

#include <cstdint>
#include <memory>
#include <pancake/Seed.hpp>
#include <pancake/SeedDBIndexCache.hpp>
#include <pancake/SeedHit.hpp>
#include <unordered_map>
#include <vector>

// struct SeedHasher {
//     size_t operator()(const uint64_t& x) const {
//         return x;
//     }
// };

/*
 * Experimental code test performance with different hash maps.
 */
// #define SEED_INDEX_USING_DENSEHASH
// #define SEED_INDEX_USING_SPARSEHASH
// #define SEED_INDEX_USING_UNORDERED_MAP
#define SEED_INDEX_USING_FLATHASHMAP

const uint64_t SEED_INDEX_EMPTY_HASH_KEY = 0xFFFFFFFFFFFFFFFF;

#ifdef SEED_INDEX_USING_UNORDERED_MAP
#include <unordered_map>
// Key: kmer hash, Value: pair of <startId, endId> in the seeds_ vector.
typedef std::unordered_map<uint64_t, std::pair<int64_t, int64_t>, std::hash<uint64_t>> SeedHashType;
#endif

#ifdef SEED_INDEX_USING_DENSEHASH
#include <sparsehash/dense_hash_map>
using google::dense_hash_map;  // namespace where class lives by default
// Key: kmer hash, Value: pair of <startId, endId> in the seeds_ vector.
typedef dense_hash_map<uint64_t, std::pair<int64_t, int64_t>, std::hash<uint64_t>> SeedHashType;
#endif

#ifdef SEED_INDEX_USING_SPARSEHASH
#include <sparsehash/sparse_hash_map>
using google::sparse_hash_map;  // namespace where class lives by default
// Key: kmer hash, Value: pair of <startId, endId> in the seeds_ vector.
typedef sparse_hash_map<uint64_t, std::pair<int64_t, int64_t>, std::hash<uint64_t>> SeedHashType;
#endif

#ifdef SEED_INDEX_USING_FLATHASHMAP
#include <flat_hash_map/flat_hash_map.hpp>
// Key: kmer hash, Value: pair of <startId, endId> in the seeds_ vector.
typedef ska::flat_hash_map<uint64_t, std::pair<int64_t, int64_t>> SeedHashType;
#endif

namespace PacBio {
namespace Pancake {

class SeedIndex
{
public:
    SeedIndex(std::vector<PacBio::Pancake::SeedDB::SeedRaw>&& seeds);
    ~SeedIndex();

    void ComputeFrequencyStats(double percentileCutoff, int64_t& retFreqMax, double& retFreqAvg,
                               double& retFreqMedian, int64_t& retFreqCutoff) const;
    int64_t GetSeeds(uint64_t key, std::vector<PacBio::Pancake::SeedDB::SeedRaw>& seeds) const;
    bool CollectHits(const std::vector<PacBio::Pancake::SeedDB::SeedRaw>& querySeeds,
                     int32_t queryLen, std::vector<SeedHit>& hits, int64_t freqCutoff) const;
    bool CollectHits(const PacBio::Pancake::SeedDB::SeedRaw* querySeeds, int64_t querySeedsSize,
                     int32_t queryLen, std::vector<SeedHit>& hits, int64_t freqCutoff) const;

    int32_t GetMinSeedSpan() const { return minSeedSpan_; }
    int32_t GetMaxSeedSpan() const { return maxSeedSpan_; }
    double GetAvgSeedSpan() const { return avgSeedSpan_; }

    const SeedHashType& GetHash() const { return hash_; }

private:
    std::vector<PacBio::Pancake::SeedDB::SeedRaw> seeds_;
    SeedHashType hash_;
    int32_t minSeedSpan_;
    int32_t maxSeedSpan_;
    double avgSeedSpan_;

    static void BuildHash_(std::vector<PacBio::Pancake::SeedDB::SeedRaw>& seeds,
                           SeedHashType& retHash, int32_t& retMinSeedSpan, int32_t& retMaxSeedSpan,
                           double& retAvgSeedSpan);
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEQDB_READER_H
