// Author: Ivan Sovic

#ifndef PANCAKE_ADAPTER_FINDER_HPP
#define PANCAKE_ADAPTER_FINDER_HPP

#include <pancake/FastaSequenceCached.hpp>
#include <pancake/Range.hpp>

#include <cstdint>
#include <memory>
#include <string_view>
#include <vector>

namespace PacBio {
namespace Pancake {
namespace Adapter {

enum class AdapterDetectionType
{
    PALINDROMIC_ROUGH,
    PALINDROMIC_FINE,
    ALIGNED,
    UNKNOWN,
};

class AdapterPosition
{
public:
    int32_t Start = 0;
    int32_t End = 0;
    AdapterDetectionType type = AdapterDetectionType::UNKNOWN;

    bool operator==(const AdapterPosition&) const noexcept = default;
};

/**
 * @brief
 *
 * @param seq Input read sequence for analysis.
 * @param adapterSeq Adapter sequence for alignment.
 * @param minNumSeeds Ignore mappings with fewer seed hits than this.
 * @param minIdentity Ignore alignments with identity less than this.
 * @param minDiagonalSpanFraction Mapping/Alignment must be at least this big compared to the diagonal where it lies.
 * @param maxAllowedFlank Mapping/Alignment can have at most this many unaligned bases at the beginning/end of the diagonal.
 * @param maxAdapterEditDist Maximum edit distsance of an aligned adapter sequence.
 * @return std::vector<AdapterPosition> A vector of candidate adapter positions (or positions where an adapter is missing), sorted in ascending order.
 */
std::vector<AdapterPosition> FindPotentialAdapterLocations(
    const FastaSequenceCached& seq, const std::string_view adapterSeq, int32_t k, int32_t w,
    int32_t minNumSeeds, double minIdentity, double minDiagonalSpanFraction,
    int32_t maxAllowedFlank, int32_t maxAdapterEditDist);

inline std::ostream& operator<<(std::ostream& os, const AdapterDetectionType& b)
{
    os << ((b == AdapterDetectionType::PALINDROMIC_ROUGH)  ? "PALINDROMIC_ROUGH"
           : (b == AdapterDetectionType::PALINDROMIC_FINE) ? "PALINDROMIC_FINE"
                                                           : "ALIGNED");
    return os;
}

inline std::ostream& operator<<(std::ostream& os, const AdapterPosition& b)
{
    os << "Start = " << b.Start << ", End = " << b.End << ", type = " << b.type;
    return os;
}

}  // namespace Adapter
}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ADAPTER_FINDER_HPP
