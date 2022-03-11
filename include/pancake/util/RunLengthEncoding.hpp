// Author: Ivan Sovic

#ifndef PANCAKE_RLE_H
#define PANCAKE_RLE_H

#include <cstdint>
#include <string>
#include <vector>

namespace PacBio {
namespace Pancake {

/// \brief Compresses the homopolymers from the input sequence. Returns (via
///         parameters) the compressed sequence, and the run lengths for each base.
void RunLengthEncoding(const std::string& seq, std::string& encodedSeq,
                       std::vector<int32_t>& runLengths);

/// \brief Compresses the homopolymers from the input sequence. Returns (via
///         parameters) the compressed sequence, and the run lengths for each base.
void RunLengthEncoding(const char* seq, int64_t seqLen, std::string& encodedSeq,
                       std::vector<int32_t>& runLengths);

/// \brief Compresses the homopolymers from the input sequence in-place.
int64_t RunLengthEncoding(char* seq, int64_t seqLen, std::vector<int32_t>& runLengths);

/// \brief Compresses the homopolymers from the input sequence. Returns (via
///         parameters) the compressed sequence and two vectors for converting the
///         positions from the compressed space back to the original sequence
///         space (hpcToSeqCoords) and vice versa (seqToHPCCoords).
void RunLengthEncoding(const std::string& seq, std::string& encodedSeq,
                       std::vector<int32_t>& seqToHPCCoords, std::vector<int32_t>& hpcToSeqCoords);

/// \brief Compresses the homopolymers from the input sequence. Returns (via
///         parameters) the compressed sequence and two vectors for converting the
///         positions from the compressed space back to the original sequence
///         space (hpcToSeqCoords) and vice versa (seqToHPCCoords).
void RunLengthEncoding(const char* seq, int64_t seqLen, std::string& encodedSeq,
                       std::vector<int32_t>& seqToHPCCoords, std::vector<int32_t>& hpcToSeqCoords);

/// \brief Compresses the homopolymers from the input sequence in place. Returns (via
///         parameters) the compressed sequence and two vectors for converting the
///         positions from the compressed space back to the original sequence
///         space (hpcToSeqCoords) and vice versa (seqToHPCCoords).
int64_t RunLengthEncoding(char* seq, int64_t seqLen, std::vector<int32_t>& seqToHPCCoords,
                          std::vector<int32_t>& hpcToSeqCoords);

/// \brief Specialized function which should be used carefuly! The output buffers:
///         destHPC, seqToHPCCoords and hpcToSeqCoords will always be larger than the
///         length of the HPC sequence to optimize the number of reallocations.
//          That's why the "hpcLen" output value should be used to determine how
//          many elements are valid in destHPC and hpcToSeqCoords.
int64_t RunLengthEncoding(const char* seq, int64_t seqLen, std::vector<char>& destHPC,
                          int32_t& hpcLen, std::vector<int32_t>& seqToHPCCoords,
                          std::vector<int32_t>& hpcToSeqCoords);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_RLE_H