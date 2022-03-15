// Author: Ivan Sovic

#ifndef PANCAKE_COMPRESSED_SEQUENCE_HPP
#define PANCAKE_COMPRESSED_SEQUENCE_HPP

#include <pancake/Range.hpp>
#include <pancake/Twobit.hpp>

#include <cstdint>
#include <memory>
#include <vector>

namespace PacBio {
namespace Pancake {

class CompressedSequence
{
public:
    CompressedSequence();
    CompressedSequence(std::string_view bases);
    ~CompressedSequence();

    void SetFromBases(std::string_view bases);
    const std::vector<uint8_t>& GetTwobit() const { return twobit_; }
    const std::vector<PacBio::Pancake::Range>& GetRanges() const { return ranges_; }
    int64_t GetNumUncompressedBases() const { return numUncompressedBases_; }
    int64_t GetNumCompressedBases() const { return numCompressedBases_; }

private:
    std::vector<uint8_t> twobit_;
    std::vector<PacBio::Pancake::Range> ranges_;
    int64_t numUncompressedBases_ = 0;
    int64_t numCompressedBases_ = 0;
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_COMPRESSED_SEQUENCE_HPP
