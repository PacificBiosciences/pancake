// Author: Ivan Sovic

#ifndef PANCAKE_MINIMIZERS_HPP
#define PANCAKE_MINIMIZERS_HPP

#include <pancake/FastaSequenceCached.hpp>
#include <pancake/Seed.hpp>
#include <pancake/SeedHit.hpp>
#include <pancake/util/CommonTypes.hpp>

#include <array>
#include <cstdint>
#include <deque>
#include <unordered_map>
#include <vector>

namespace PacBio {
namespace Pancake {
namespace SeedDB {

static inline uint64_t InvertibleHash(uint64_t key, uint64_t mask)
{
    /*
    Credit: Heng Li, Minimap2.
    */
    key = (~key + (key << 21)) & mask;  // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask;  // key * 265
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask;  // key * 21
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;
    return key;
}

static inline uint64_t ComputeKmerMask(int32_t kmerSize)
{
    const uint64_t mask =
        (kmerSize < 32) ? ((((uint64_t)1) << (2 * kmerSize)) - 1)
                        : 0xFFFFFFFFFFFFFFFF;  // Mask the number of required bits for the kmer.
    return mask;
}

/*
 * \brief Computes minimizers for a single input sequence.
*/
int GenerateMinimizers(std::vector<PacBio::Pancake::Int128t>& minimizers, const uint8_t* seq,
                       const int32_t seqLen, const int32_t seqOffset, const int32_t seqId,
                       const int32_t kmerSize, const int32_t winSize, const int32_t spacing,
                       const bool useReverseComplement, const bool useHPC);

/*
 * \brief Computes minimizers for a set of input sequences, given as a vector of FastaSequenceCached objects.
*/
void GenerateMinimizers(std::vector<PacBio::Pancake::Int128t>& retSeeds,
                        const std::vector<FastaSequenceCached>& targetSeqs, const int32_t kmerSize,
                        const int32_t winSize, const int32_t spacing,
                        const bool useReverseComplement, const bool useHPC);

/*
 * \brief Computes minimizers for a set of input sequences, given as a vector of std::string objects.
*/
void GenerateMinimizers(std::vector<PacBio::Pancake::Int128t>& retSeeds,
                        const std::vector<std::string>& targetSeqs, const int32_t kmerSize,
                        const int32_t winSize, const int32_t spacing,
                        const bool useReverseComplement, const bool useHPC);

/*
 * \brief Computes minimizers for a set of input sequences, given as a vector of FastaSequenceCached objects.
 *        Also collects sequence lengths for all given input sequences.
*/
void GenerateMinimizers(std::vector<PacBio::Pancake::Int128t>& retSeeds,
                        std::vector<int32_t>& retSequenceLengths,
                        const std::vector<FastaSequenceCached>& targetSeqs, const int32_t kmerSize,
                        const int32_t winSize, const int32_t spacing,
                        const bool useReverseComplement, const bool useHPC);

/*
 * \brief Computes minimizers for a set of input sequences, given as a vector of std::string objects.
 *        Also collects sequence lengths for all given input sequences.
*/
void GenerateMinimizers(std::vector<PacBio::Pancake::Int128t>& retSeeds,
                        std::vector<int32_t>& retSequenceLengths,
                        const std::vector<std::string>& targetSeqs, const int32_t kmerSize,
                        const int32_t winSize, const int32_t spacing,
                        const bool useReverseComplement, const bool useHPC);

/**
 * \brief For a given set of query seeds, fetches the count of seed hits and produces
 *          a histogram of seed hits.
 * \return The histogram is in a vector form, where each element of the vector consists of
 *          <occurrence_count, num query seeds that have this count>
*/
template <class TargetHashType>
std::vector<std::pair<int64_t, int64_t>> ComputeSeedHitHistogram(
    const PacBio::Pancake::SeedDB::SeedRaw* querySeeds, const int64_t querySeedsSize,
    const TargetHashType& hash)
{
    std::unordered_map<int64_t, int64_t> hist;
    for (int64_t seedId = 0; seedId < querySeedsSize; ++seedId) {
        const auto& querySeed = querySeeds[seedId];
        auto decodedQuery = PacBio::Pancake::SeedDB::Seed(querySeed);
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

template <class TargetHashType>
bool CollectSeedHits(std::vector<SeedHit>& hits, const PacBio::Pancake::SeedDB::SeedRaw* querySeeds,
                     const int64_t querySeedsSize, const int32_t queryLen,
                     const TargetHashType& hash,
                     const PacBio::Pancake::SeedDB::SeedRaw* targetSeeds,
                     const int64_t /*targetSeedsSize*/, const int64_t freqCutoff)
{
    hits.clear();

    // The +1 is because for every seed base there are Spacing spaces, and the subtraction
    // is because after the last seed base the spaces shouldn't be counted.
    //    const int32_t seedSize = kmerSize * (spacing + 1) - spacing;

    for (int64_t seedId = 0; seedId < querySeedsSize; ++seedId) {
        const auto& querySeed = querySeeds[seedId];
        const auto decodedQuery = PacBio::Pancake::SeedDB::Seed(querySeed);
        const auto it = hash.find(decodedQuery.key);
        if (it != hash.end()) {
            const int64_t start = std::get<0>(it->second);
            const int64_t end = std::get<1>(it->second);
            // Skip very frequent seeds.
            if (freqCutoff > 0 && (end - start) > freqCutoff) {
                continue;
            }
            for (int64_t i = start; i < end; ++i) {
                const auto decodedTarget = PacBio::Pancake::SeedDB::Seed(targetSeeds[i]);
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

                SeedHit hit{static_cast<int32_t>(decodedTarget.seqID),
                            isRev,
                            targetPos,
                            queryPos,
                            targetSpan,
                            querySpan,
                            0};

                hits.emplace_back(std::move(hit));
            }
        }
    }

    return hits.size() > 0;
}

}  // namespace SeedDB
}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_MINIMIZERS_HPP
