// Author: Ivan Sovic

#ifndef PANCAKE_SEED_INDEX_HPP
#define PANCAKE_SEED_INDEX_HPP

#include <pancake/Seed.hpp>
#include <pancake/SeedDBIndexCache.hpp>
#include <pancake/SeedHit.hpp>

#include <cstdint>
#include <memory>
#include <span>
#include <unordered_map>
#include <vector>

// #define SEED_INDEX_USING_DENSEHASH
// #define SEED_INDEX_USING_SPARSEHASH
// #define SEED_INDEX_USING_UNORDERED_MAP
#define SEED_INDEX_USING_FLATHASHMAP

#ifdef SEED_INDEX_USING_UNORDERED_MAP
#include <unordered_map>
// Key: kmer hash, Value: pair of <startId, endId> in the seeds_ vector.
typedef std::unordered_map<uint64_t, std::pair<int64_t, int64_t>, std::hash<uint64_t>> SeedHashType;
#endif

#ifdef SEED_INDEX_USING_DENSEHASH
#include <pancake/third-party/sparsehash/dense_hash_map>
using google::dense_hash_map;  // namespace where class lives by default
// Key: kmer hash, Value: pair of <startId, endId> in the seeds_ vector.
const uint64_t SEED_INDEX_EMPTY_HASH_KEY = 0xFFFFFFFFFFFFFFFF;
typedef dense_hash_map<uint64_t, std::pair<int64_t, int64_t>, std::hash<uint64_t>> SeedHashType;
#endif

#ifdef SEED_INDEX_USING_SPARSEHASH
#include <pancake/third-party/sparsehash/sparse_hash_map>
using google::sparse_hash_map;  // namespace where class lives by default
// Key: kmer hash, Value: pair of <startId, endId> in the seeds_ vector.
typedef sparse_hash_map<uint64_t, std::pair<int64_t, int64_t>, std::hash<uint64_t>> SeedHashType;
#endif

#ifdef SEED_INDEX_USING_FLATHASHMAP
#include <pancake/third-party/flat_hash_map/flat_hash_map.hpp>
// Key: kmer hash, Value: pair of <startId, endId> in the seeds_ vector.
typedef ska::flat_hash_map<uint64_t, std::pair<int64_t, int64_t>> SeedHashType;
#endif

namespace PacBio {
namespace Pancake {

/**
 * @brief For a given set of query seeds, fetches the count of seed hits and produces
 *          a histogram of seed hits.
 * @return The histogram is in a vector form, where each element of the vector consists of
 *          <occurrence_count, num query seeds that have this count>
*/
std::vector<std::pair<int64_t, int64_t>> ComputeSeedHitHistogram(
    const std::span<const PacBio::Pancake::SeedRaw> querySeeds, const SeedHashType& hash);

/**
 * \brief Given several types of occurrence thresholds, computes the final occurrence threshold that can be used
 *          for collecting seed hits. Uses the following formula to compute the return value:
 *              cutoff = max(seedOccurrenceMin, min(seedOccurrenceMax, occThresholdMemMax, seedOccurrenceUserSpecified))
 *          where occThresholdMemMax is computed from the histogram if parameter seedOccurrenceMaxMemory > 0, and other parameters are user-provided via the API.
 * \param seedHitHistogram Histogram of seed hit occurrences, given as a vector of pairs: <seed occurrence, number of query seeds with this occurrence>.
 *                          Can be an empty vector. Used only when seedOccurrenceMaxMemory > 0. The vector should be sorted in the ascending order of seed occurrence.
 * \param seedOccurrenceMin Minimum value for the occurrence threshold. If the other cutoff values result in a smaller value than the minimum, then the return value is pinned to this.
 * \param seedOccurrenceMax Maximum allowed occurrence of a seed to keep for mapping. Value <= 0 turns off this threshold.
 * \param seedOccurrenceMaxMemory Maximum allowed memory to be consumed by collected seed hits. This is used to dynamically compute the maximum
 *                                  occurrence cutoff based on the seed hit histogram. Seeds are chosen in the sorted order by their occurrence
 *                                  until the memory threshold is reached. Value <= 0 turns off this threshold.
 * \param seedOccurrenceUserSpecified User-provided cutoff threshold, computed, for example, from the frequency percentile threshold
 *                                          (e.g. top 0.002% of most abundant seeds in the index should be skipped.).
 * \return Single value that consolidates the occurrence cutoff. If the cutoff is not applied (e.g. all or some provided values are zero), it returns std::numeric_limits<int64_t>::max().
*/
int64_t ComputeOccurrenceThreshold(const std::vector<std::pair<int64_t, int64_t>>& seedHitHistogram,
                                   int64_t seedOccurrenceMin, int64_t seedOccurrenceMax,
                                   int64_t seedOccurrenceMaxMemory,
                                   int64_t seedOccurrenceUserSpecified, bool debugVerbose);

/**
 * @brief Given a vector of query seeds and a set of hashed target seeds, finds all seed hits and
 *          returns them.
 * @param hits Return value, all collected hits. Collected hits are returned via parameter to allow memory reuse.
 * @param querySeeds Input query seeds.
 * @param queryLen Length of the query sequence.
 * @param hash Hash of the target seeds for random lookup.
 * @param targetSeeds Input target seeds.
 * @param freqCutoff Maximum allowed occurrence of a seed hit to collect it.
 */
void CollectSeedHits(std::vector<SeedHit>& hits,
                     std::span<const PacBio::Pancake::SeedRaw> querySeeds, int32_t queryLen,
                     const SeedHashType& hash,
                     std::span<const PacBio::Pancake::SeedRaw> targetSeeds, int64_t freqCutoff);

class SeedIndex
{
public:
    SeedIndex(std::vector<PacBio::Pancake::SeedRaw>&& seeds);
    ~SeedIndex();

    void ComputeFrequencyStats(double percentileCutoff, int64_t& retFreqMax, double& retFreqAvg,
                               double& retFreqMedian, int64_t& retFreqCutoff) const;
    int64_t GetSeeds(uint64_t key, std::vector<PacBio::Pancake::SeedRaw>& seeds) const;
    void CollectHits(const std::vector<PacBio::Pancake::SeedRaw>& querySeeds, int32_t queryLen,
                     std::vector<SeedHit>& hits, int64_t freqCutoff) const;
    void CollectHits(const PacBio::Pancake::SeedRaw* querySeeds, int64_t querySeedsSize,
                     int32_t queryLen, std::vector<SeedHit>& hits, int64_t freqCutoff) const;

    int32_t GetMinSeedSpan() const { return minSeedSpan_; }
    int32_t GetMaxSeedSpan() const { return maxSeedSpan_; }
    double GetAvgSeedSpan() const { return avgSeedSpan_; }

    const SeedHashType& GetHash() const { return hash_; }

private:
    std::vector<PacBio::Pancake::SeedRaw> seeds_;
    SeedHashType hash_;
    int32_t minSeedSpan_;
    int32_t maxSeedSpan_;
    double avgSeedSpan_;

    static void BuildHash_(std::vector<PacBio::Pancake::SeedRaw>& seeds, SeedHashType& retHash,
                           int32_t& retMinSeedSpan, int32_t& retMaxSeedSpan,
                           double& retAvgSeedSpan);
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEED_INDEX_HPP
