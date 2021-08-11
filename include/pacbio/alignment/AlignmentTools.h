// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNMENT_TOOLS_H
#define PANCAKE_ALIGNMENT_TOOLS_H

#include <pacbio/alignment/DiffCounts.h>
#include <pbbam/Cigar.h>
#include <pbbam/CigarOperation.h>

namespace PacBio {
namespace Pancake {

struct TrimmingInfo
{
    int32_t queryFront = 0;
    int32_t targetFront = 0;
    int32_t queryBack = 0;
    int32_t targetBack = 0;
};
inline bool operator==(const TrimmingInfo& lhs, const TrimmingInfo& rhs)
{
    return lhs.queryFront == rhs.queryFront && lhs.queryBack == rhs.queryBack &&
           lhs.targetFront == rhs.targetFront && lhs.targetBack == rhs.targetBack;
}

PacBio::BAM::Cigar EdlibAlignmentToCigar(const unsigned char* aln, int32_t alnLen,
                                         Alignment::DiffCounts& retDiffs);

void EdlibAlignmentDiffCounts(const unsigned char* aln, int32_t alnLen, int32_t& numEq,
                              int32_t& numX, int32_t& numI, int32_t& numD);

void CigarDiffCounts(const PacBio::BAM::Cigar& cigar, int32_t& numEq, int32_t& numX, int32_t& numI,
                     int32_t& numD);

Alignment::DiffCounts CigarDiffCounts(const PacBio::BAM::Cigar& cigar);

void AppendToCigar(PacBio::BAM::Cigar& cigar, PacBio::BAM::CigarOperationType newOp,
                   int32_t newLen);

PacBio::BAM::Cigar ExpandMismatches(const char* query, int64_t queryLen, const char* target,
                                    int64_t targetLen, const PacBio::BAM::Cigar& cigar);

void ValidateCigar(const char* query, int64_t queryLen, const char* target, int64_t targetLen,
                   const PacBio::BAM::Cigar& cigar, const std::string& label);

void ExtractVariantString(const char* query, int64_t queryLen, const char* target,
                          int64_t targetLen, const PacBio::BAM::Cigar& cigar, bool maskHomopolymers,
                          bool maskSimpleRepeats, bool maskHomopolymerSNPs,
                          bool maskHomopolymersArbitrary, std::string& retQueryVariants,
                          std::string& retTargetVariants, Alignment::DiffCounts& retDiffsPerBase,
                          Alignment::DiffCounts& retDiffsPerEvent);

Alignment::DiffCounts ComputeDiffCounts(const PacBio::BAM::Cigar& cigar,
                                        const std::string& queryVariants,
                                        const std::string& targetVariants,
                                        bool throwOnPartiallyMaskedIndels);

/// \brief For a given query position finds the corresponding target position based
///         on a provided CIGAR string.
///         Lineraly scans through all CIGAR operations to perform the mapping.
int32_t FindTargetPosFromCigar(const BAM::Cigar& cigar, int32_t queryPos);

void NormalizeAlignmentInPlace(std::string& queryAln, std::string& targetAln);

void ConvertCigarToM5(const char* query, int64_t queryLen, const char* target, int64_t targetLen,
                      const Data::Cigar& cigar, std::string& retQueryAln,
                      std::string& retTargetAln);

Data::Cigar ConvertM5ToCigar(const std::string& queryAln, const std::string& targetAln);

Data::Cigar NormalizeCigar(const char* query, int64_t queryLen, const char* target,
                           int64_t targetLen, const Data::Cigar& cigar);

bool TrimCigar(const PacBio::BAM::Cigar& cigar, int32_t windowSize, int32_t minMatches,
               bool clipOnFirstMatch, PacBio::BAM::Cigar& retTrimmedCigar,
               TrimmingInfo& retTrimming);

int32_t ScoreCigarAlignment(const PacBio::BAM::Cigar& cigar, int32_t match, int32_t mismatch,
                            int32_t gapOpen, int32_t gapExt);

void MergeCigars(PacBio::Data::Cigar& dest, const PacBio::Data::Cigar& src);

/**
 * \brief Computes a vector of the length of the input sequence, where each position has an
 *          8-bit unsigned int indicating whether the base is masked or not.
 *          Value 0 means there is no masking. Multiple levels of simple repeats can be marked
 *          in the same element of the vector: HPs have a value of (1 << 0), dinucs a value of (1 << 1),
 *          trinucs (1 << 2), etc. So if a base is marked as a homopolymer, the corresponding position
 *          in the return vector would have a value of 1. If the base is both a part of a HP and a dinuc
 *          repeat, it would have a value of (1 + 2 = 3), and so on.
 * \param seq C-style string of the sequence. Not null-terminated.
 * \param seqLen Length of the input sequence.
 * \param maxWindowSize The maximum level of simple repeats for masking: 0 means no masking, 1 will mask homopolymers,
 *          2 will mask homopolymers and dinucleotide repeats, 3 will mask HPs + dinucs + trinucs, etc.
 *          Complexity of computation is O(seqLen * maxWindowSize).
 * \return Vector with a mask for each sequence base indicating whether the base is masked (value > 0) or not.
*/
std::vector<uint8_t> ComputeSimpleRepeatMask(const char* seq, int32_t seqLen,
                                             int32_t maxWindowSize);

/** \brief Walks through the CIGAR vector and for every position computes the diagonal.
 *          If (diag > (bandwidth - 1) || diag < -(bandwidth - 1)), function returns true
 *          to indicate suboptimal alignments.
*/
bool CheckAlignmentOutOfBand(const PacBio::Data::Cigar& cigar, const int32_t bandwidth);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ALIGNMENT_TOOLS_H