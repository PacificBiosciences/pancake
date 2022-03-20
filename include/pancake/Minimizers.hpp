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
#include <span>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace PacBio {
namespace Pancake {

static inline uint64_t ComputeKmerMask(const int32_t kmerSize)
{
    const uint64_t mask =
        (kmerSize < 32) ? ((((uint64_t)1) << (2 * kmerSize)) - 1)
                        : 0xFFFFFFFFFFFFFFFF;  // Mask the number of required bits for the kmer.
    return mask;
}

/*
 * \brief Computes minimizers for a single input sequence.
*/
int GenerateMinimizers(std::vector<PacBio::Pancake::Int128t>& minimizers, std::string_view seq,
                       int32_t seqOffset, int32_t seqId, int32_t kmerSize, int32_t winSize,
                       int32_t spacing, bool useReverseComplement, bool useHPC);

/*
 * \brief Computes minimizers for a set of input sequences, given as a vector of FastaSequenceCached objects.
*/
void GenerateMinimizers(std::vector<PacBio::Pancake::Int128t>& retSeeds,
                        const std::vector<FastaSequenceCached>& seqs, int32_t kmerSize,
                        int32_t winSize, int32_t spacing, bool useReverseComplement, bool useHPC);

/*
 * \brief Computes minimizers for a set of input sequences, given as a vector of std::string objects.
*/
void GenerateMinimizers(std::vector<PacBio::Pancake::Int128t>& retSeeds,
                        const std::vector<std::string>& seqs, int32_t kmerSize, int32_t winSize,
                        int32_t spacing, bool useReverseComplement, bool useHPC);

/*
 * \brief Computes minimizers for a set of input sequences, given as a vector of FastaSequenceCached objects.
 *        Also collects sequence lengths for all given input sequences.
*/
void GenerateMinimizers(std::vector<PacBio::Pancake::Int128t>& retSeeds,
                        std::vector<int32_t>& retSequenceLengths,
                        const std::vector<FastaSequenceCached>& seqs, int32_t kmerSize,
                        int32_t winSize, int32_t spacing, bool useReverseComplement, bool useHPC);

/*
 * \brief Computes minimizers for a set of input sequences, given as a vector of std::string objects.
 *        Also collects sequence lengths for all given input sequences.
*/
void GenerateMinimizers(std::vector<PacBio::Pancake::Int128t>& retSeeds,
                        std::vector<int32_t>& retSequenceLengths,
                        const std::vector<std::string>& seqs, int32_t kmerSize, int32_t winSize,
                        int32_t spacing, bool useReverseComplement, bool useHPC);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_MINIMIZERS_HPP
