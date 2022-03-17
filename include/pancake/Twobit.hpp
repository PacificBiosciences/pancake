// Author: Ivan Sovic

#ifndef PANCAKE_TWOBIT_HPP
#define PANCAKE_TWOBIT_HPP

#include <pancake/Range.hpp>

#include <string>
#include <vector>

namespace PacBio {
namespace Pancake {

/// \brief Converts a sequence to a 2-bit packed array.
///
/// \note  Any contiguous sequence of non-ACTG bases is called a "range".
///        If a non-ACTG base is found, a new range is constructed.
///        This is a lossy compression, because any non-ACTGN base will be
///        considered as an N.
///
///        To understand the purpose of a Range better, imagine an input sequence
///        "ACTGNNNTTT". If we were to split this sequence on N-bases,
///        we would get two subsequences: "ACTG" and "TTT".
///        The "ACTG" subsequence begins at position 0 and ends at position 4, so
///        it's range is (0, 4).
///        The "TTT" subsequence begins at position 7 and ends at position 10, so
///        it's range is (7, 10).
///        Since N-bases cannot be encoded in 2 bits, we remove them
///        from the final compressed sequence. All other bases are just
///        "condensed" into: "ACTGTTT" during compression, and we can decode
///        where the N-bases need to go through the vector of Ranges.
///
/// \param[in]  bases          Input sequence in ASCII encoding.
/// \param[out] twobit         Resulting 2-bit packed vector.
/// \param[out] ranges         A vector of contiguous non-ACTG ranges in the
///                            packed 2-bit vector.
///
/// \returns Number of bases that were compressed to twobit representation.
///
int32_t CompressSequence(const std::string& bases, std::vector<uint8_t>& twobit,
                         std::vector<PacBio::Pancake::Range>& ranges);

/// \brief Decompresses a 2-bit compressed sequence into a string.
///
/// \note  Any contiguous sequence of non-ACTG bases is called a "range".
///        Ranges are decompressed by converting the 2-bit representation of
///        any base into the ASCII form, and the space in-between ranges is
///        translated into a series of N-bases.
///        The total number of bases of the decompressed sequence is needed so
///        that the trailing N-bases can be added (they are not encoded in the
///        ranges array).
///
/// \param[in] twobit          The 2-bit packed vector.
/// \param[in] numBases        The length of the uncompressed sequence.
/// \param[in] ranges          A vector of contiguous non-ACTG ranges in the
///                            packed 2-bit vector.
/// \param[out]  bases         Decompressed sequence in ASCII encoding.
///
/// \throws std::runtime_error if ranges is empty or if numBases is invalid.
///
void DecompressSequence(const std::vector<uint8_t>& twobit, int32_t numBases,
                        const std::vector<PacBio::Pancake::Range>& ranges, std::string& bases);

/// \brief Decompresses a 2-bit compressed sequence into a string.
///        Same as the previous overloaded function, but the difference is that
///        this provides C-type interface, and does not perform internal
///        allocation for the return bases.
///
///        The space for outBases needs to be allocated prior to calling this function!
///
///        This is a specialized function, intended for use in e.g. fast loading
///        of compressed sequences, without having to pay for the high overhead
///        of std::string in cases when we actually want to store the data
///        into a preallocated vector of uint8_t.
///
void DecompressSequence(const uint8_t* twobit, int64_t twobitLen, int32_t numBases,
                        const std::vector<PacBio::Pancake::Range>& ranges, uint8_t* outBases);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_TWOBIT_HPP
