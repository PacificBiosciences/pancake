// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNMENT_TOOLS_HPP
#define PANCAKE_ALIGNMENT_TOOLS_HPP

#include <pancake/DiffCounts.hpp>

#include <pbbam/Cigar.h>
#include <pbbam/CigarOperation.h>

#include <cstdint>
#include <span>
#include <string>
#include <string_view>

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

/**
 * @brief Converts the Edlib-style alignment to a PacBio::Cigar type.
 *
 * @param aln Input Edlib alignment.
 * @param retDiffs Return parameter, counts of CIGAR operations.
 * @return PacBio::BAM::Cigar The CIGAR version of the input alignment.
 */
PacBio::BAM::Cigar EdlibAlignmentToCigar(std::span<const unsigned char> aln, DiffCounts& retDiffs);

/**
 * @brief Computes the diff counts (matches, mismatches, insertions and deletions) from a given
 *          Edlib-style alignment.
 *
 * @param aln Input alignment, Edlib-style.
 * @return DiffCounts Counts of Cigar operations.
 */
DiffCounts EdlibAlignmentDiffCounts(std::span<const unsigned char> aln);

/**
 * @brief Computes the diff counts from a given CIGAR: matches, mismatches, insertions and deletions.
 *
 * @param cigar Input CIGAR.
 * @return DiffCounts Counts of Cigar operations.
 */
DiffCounts CigarDiffCounts(const PacBio::BAM::Cigar& cigar);

/**
 * @brief Adds a single CIGAR operation to an existing Cigar object. Takes care to check if
 *          the op matches the last op in the existing CIGAR, and in that case it just increments
 *          the count.
 *
 * @param cigar Cigar alignment, modified in place.
 * @param newOp Cigar operation to add.
 * @param newLen Length of the cigar operation to add.
 */
void AppendToCigar(PacBio::BAM::Cigar& cigar, PacBio::BAM::CigarOperationType newOp,
                   int32_t newLen);

/**
 * @brief Finds successive INS/DEL operations in the CIGAR, and converts them to a combination of
 *          match/mismatch operations + a remaining gap.
 * @param query Query sequence in alignment.
 * @param target Target sequence in alignment.
 * @param cigar. Alignment.
 * @return PacBio::BAM::Cigar New CIGAR with resolved neighboring INS/DEL operations.
 */
PacBio::BAM::Cigar ExpandMismatches(std::string_view query, std::string_view target,
                                    const PacBio::BAM::Cigar& cigar);

/**
 * @brief Validates the correctness of the CIGAR with respect to the query/target pair. THROWS.
 *          Throws if the CIGAR is not valid. Examples where it may throw:
 *          - Query/target coordinats out of bounds of the query/target sequences when walking
 *              down the CIGAR vector.
 *          - CIGAR match operation given, but sequences do not match.
 *          - CIGAR mismatch operation given, but some of the bases in the span are equal in query/target.
 *          - CIGAR insertion/soft clip operation given, but the span of this op steps out of bounds of the query.
 *          - CIGAR deletion/reference skip operation given, but the span of this op steps out of bounds of the target.
 *          - Unsupported CIGAR operation detected.
 *          - Computed query/target lengths from the given CIGAR does not match the input query/target sequences lengths.
 *
 * @param query Query sequence used to produce the CIGAR.
 * @param target Target sequence used to produce the CIGAR.
 * @param cigar Alignment between the query and target to validate.
 * @param label Label for the exception message, for debug purposes.
 */
void ValidateCigar(std::string_view query, std::string_view target, const PacBio::BAM::Cigar& cigar,
                   std::string_view label);

/**
 * @brief For a given set of query sequence, target sequence and CIGAR, it extracts all
 *          non-matching bases covered by the alignment into two strings:
 *          - queryVariants - mismatches and insertions.
 *          - targetVariants - mismatches and deletions.
 *          These describe all variants between the query and the target sequence.
 *          Considering that valid alignments are composed of >50% of matching bases,
 *          this can reduce space in storing only the variant positions instead all three
 *          components (query, target, CIGAR).
 *          This function also allows to mask homopolymers and simple repeats.
 *          Simple repeats are tandem expansions of short kmers.
 *          Masking the variants means that their bases will be reported in lower case. Unmasked
 *          variants (default) are upper case.
 *
 * @param query Query sequence.
 * @param target Target sequence.
 * @param cigar CIGAR alignment.
 * @param maskHomopolymers Masks homopolymer bases.
 * @param maskSimpleRepeats Checks if an indel is exactly the same as the preceding or following bases in either the
 *                          query or the target, and masks the events if so.
 * @param maskHomopolymerSNPs If there is a SNP in what appears to be a homopolymer (in either query or target) this will
 *                              mask the SNP.
 * @param maskHomopolymersArbitrary Allows masking "homopolymer" events where there is a non-HP base in the middle of the event.
 * @param retQueryVariants Return string containing query variants.
 * @param retTargetVariants Return string containing target vriants.
 * @param retDiffsPerBase Return count of CIGAR operation differences, per base.
 * @param retDiffsPerEvent Return count of CIGAR operation differences, indels are computed per event (mismatches per base).
 */
void ExtractVariantString(std::string_view query, std::string_view target,
                          const PacBio::BAM::Cigar& cigar, bool maskHomopolymers,
                          bool maskSimpleRepeats, bool maskHomopolymerSNPs,
                          bool maskHomopolymersArbitrary, std::string& retQueryVariants,
                          std::string& retTargetVariants, DiffCounts& retDiffsPerBase,
                          DiffCounts& retDiffsPerEvent);

/**
 * @brief Computes the CIGAR diff counts, but ignores masked variants.
 *
 * @param cigar Input alignment.
 * @param queryVariants Query variant string containing mismatches and insertions. Masked bases are lowercase.
 * @param targetVariants Target variant string containing mismatches and deletions. Masked bases are lowercase.
 * @param throwOnPartiallyMaskedIndels Throws if an indel event has a mix of masked/unmasked bases.
 * @return DiffCounts Counts of differences, sans masked bases.
 */
DiffCounts ComputeMaskedDiffCounts(const PacBio::BAM::Cigar& cigar, std::string_view queryVariants,
                                   std::string_view targetVariants,
                                   bool throwOnPartiallyMaskedIndels);

/**
 * @brief For a given query position finds the corresponding target position based in the
 *          input CIGAR alignment. Lineraly scans through all CIGAR operations.
 *
 * @param cigar CIGAR alignment.
 * @param queryPos Query position to search for.
 * @return int32_t
 */
int32_t FindTargetPosFromCigar(const BAM::Cigar& cigar, int32_t queryPos);

/**
 * @brief This function normalizes gaps by pushing them towards the ends of the
 *          query and target sequences. It also takes care of mismatches, shifting gaps through them.
 *          Example of what this function does.
 *          TTGACACT       TTGACACT
 *          ||| X|||   ->  |||X |||
 *          TTG-TACT       TTGT-ACT
 *
 *          Throws if input alignment lengths differ.
 *
 * @param queryAln M5-style query alignment.
 * @param targetAln M5-style target alignment.
 */
void NormalizeM5AlignmentInPlace(std::string& queryAln, std::string& targetAln);

/**
 * @brief Converts a given alignment from CIGAR style to the M5-style consisting
 *          of two strings (query and target) with all bases in it, interspersed with
 *          '-' ops for gaps. Example:
 *          CIGAR:     3=1D1X1=1I2=
 *          queryAln:  ACT-AGGAT
 *          targetAln: ACTACG-AT
 *
 * @param query Input query sequence.
 * @param target Input target sequence.
 * @param cigar Input CIGAR alignment.
 * @param retQueryAln Aligned query string, M5-style.
 * @param retTargetAln Aligned target string, M5-style.
 */
void ConvertCigarToM5(std::string_view query, std::string_view target, const Data::Cigar& cigar,
                      std::string& retQueryAln, std::string& retTargetAln);

/**
 * @brief Converts the M5-formatted alignment into the CIGAR format.
 *
 * @param queryAln Query portion of the M5 alignment.
 * @param targetAln Target portion of the M5 alignment.
 * @return Data::Cigar
 */
Data::Cigar ConvertM5ToCigar(std::string_view queryAln, std::string_view targetAln);

/**
 * @brief Normalizes the gaps in a CIGAR alignment.
 *
 * @param query Query sequence.
 * @param target Target sequence.
 * @param cigar CIGAR alignment.
 * @return Data::Cigar Gap-normalized CIGAR.
 */
Data::Cigar NormalizeCigar(std::string_view query, std::string_view target,
                           const Data::Cigar& cigar);

/**
 * @brief Trims the CIGAR alignment on the 5' and 3' ends with a sliding window.
 *          A window is slid from the left (or right) and the number of matches and diffs
 *          computed. A window is valid if it has at least minMatches matches.
 *          The alignment is clipped from the beginning to the first base of a first valid window
 *          (analogously, the 3' end is processed in a similar way).
 *          If clipOnFirstMatch is true, then the first base of a valid window also needs to be
 *          a match event.
 *
 * @param cigar Input alignment.
 * @param windowSize Window size for the sliding window analysis.
 * @param minMatches Minimum number of matches in a window to call it valid.
 * @param clipOnFirstMatch Clipping will only happen if the first CIGAR op in a valid window is a match.
 * @param retTrimmedCigar Return value, the trimmed CIGAR.
 * @param retTrimming Return value, structure with the number of bases that were clipped from the target/query front/back.
 * @return True if everything went fine.
 */
bool TrimCigar(const PacBio::BAM::Cigar& cigar, int32_t windowSize, int32_t minMatches,
               bool clipOnFirstMatch, PacBio::BAM::Cigar& retTrimmedCigar,
               TrimmingInfo& retTrimming);

/**
 * @brief Computes the alignment score from a given CIGAR vector.
 *
 * @param cigar Input alignment in CIGAR format.
 * @param match Match score.
 * @param mismatch Mismatch score (positive value).
 * @param gapOpen Gap open score (positive value).
 * @param gapExt Gap extend score (positive value).
 * @return int32_t Alignment score.
 */
int32_t ScoreCigarAlignment(const PacBio::BAM::Cigar& cigar, int32_t match, int32_t mismatch,
                            int32_t gapOpen, int32_t gapExt);

/**
 * @brief Computes the alignment score from a given CIGAR vector, using double affine gap penalties.
 *
 * @param cigar Input alignment in CIGAR format.
 * @param match Match score.
 * @param mismatch Mismatch score (positive value).
 * @param gapOpen1 Gap open score for the first affine function (positive value).
 * @param gapExt1 Gap extend score for the first affine function (positive value).
 * @param gapOpen2 Gap open score for the second affine function (positive value).
 * @param gapExt2 Gap extend score for the second affine function (positive value).
 * @return std::pair<int32_t, PacBio::Pancake::DiffCounts> Pair: (alignment score, diff counts).
 */
std::pair<int32_t, PacBio::Pancake::DiffCounts> ScoreCigarAlignment(
    const PacBio::BAM::Cigar& cigar, int32_t match, int32_t mismatch, int32_t gapOpen1,
    int32_t gapExt1, int32_t gapOpen2, int32_t gapExt2);

/**
 * @brief Merges the src CIGAR vector into the existing dest vector.
 *
 * @param dest Destination of the merge. Operations from src will be appended to the back of dst.
 * @param src Source for merging.
 */
void MergeCigars(PacBio::Data::Cigar& dest, const PacBio::Data::Cigar& src);

/**
 * @brief Computes a vector of the length of the input sequence, where each position has an
 *          8-bit unsigned int indicating whether the base is masked or not.
 *          Value 0 means there is no masking. Multiple levels of simple repeats can be marked
 *          in the same element of the vector: HPs have a value of (1 << 0), dinucs a value of (1 << 1),
 *          trinucs (1 << 2), etc. So if a base is marked as a homopolymer, the corresponding position
 *          in the return vector would have a value of 1. If the base is both a part of a HP and a dinuc
 *          repeat, it would have a value of (1 + 2 = 3), and so on.
 * @param seq C-style string of the sequence. Not null-terminated.
 * @param seqLen Length of the input sequence.
 * @param maxWindowSize The maximum level of simple repeats for masking: 0 means no masking, 1 will mask homopolymers,
 *          2 will mask homopolymers and dinucleotide repeats, 3 will mask HPs + dinucs + trinucs, etc.
 *          Complexity of computation is O(seqLen * maxWindowSize).
 * @return Vector with a mask for each sequence base indicating whether the base is masked (value > 0) or not.
*/
std::vector<uint8_t> ComputeSimpleRepeatMask(std::string_view seq, int32_t maxWindowSize);

/**
 * @brief Walks through the CIGAR vector and for every position computes the diagonal.
 *          If (diag > (bandwidth - 1) || diag < -(bandwidth - 1)), function returns true
 *          to indicate suboptimal alignments.
 * @param cigar Input alignment in CIGAR format.
 * @param bandwidth Maximum allowed bandwidth in alignment.
 * @return true if the alignment path touches or exceeds the diagonal defined by bandwidth, otherwise false.
 */
bool CheckAlignmentOutOfBand(const PacBio::Data::Cigar& cigar, int32_t bandwidth);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ALIGNMENT_TOOLS_HPP
