// Author: Ivan Sovic

#ifndef PANCAKE_MAPPER_BASE_HPP
#define PANCAKE_MAPPER_BASE_HPP

#include <pancake/AlignmentSeeded.hpp>
#include <pancake/DPChain.hpp>
#include <pancake/FastaSequenceCached.hpp>
#include <pancake/FastaSequenceId.hpp>
#include <pancake/Overlap.hpp>
#include <pancake/OverlapWriterBase.hpp>
#include <pancake/SeedIndex.hpp>

#include <cstdint>
#include <memory>
#include <unordered_map>
#include <vector>

namespace PacBio {
namespace Pancake {

struct ChainedRegion
{
    ChainedHits chain;
    std::vector<AlignmentRegion> regionsForAln;
    OverlapPtr mapping;
    // Priority 0 means primary alignment, 1 secondary, < 0 not set, and > 1 filtered.
    int32_t priority = 0;
    bool isSupplementary = false;
};

inline bool operator==(const ChainedRegion& lhs, const ChainedRegion& rhs) noexcept
{
    return lhs.chain == rhs.chain && lhs.regionsForAln == rhs.regionsForAln &&
           ((lhs.mapping == nullptr && rhs.mapping == nullptr) ||
            (lhs.mapping != nullptr && rhs.mapping != nullptr &&
             *(lhs.mapping) == *(rhs.mapping))) &&
           lhs.priority == rhs.priority && lhs.isSupplementary == rhs.isSupplementary;
}

inline void VerboseChainedRegion(std::ostream& os, const ChainedRegion& b, bool detailed)
{
    if (b.mapping == nullptr) {
        os << "Mapping: nullptr!";
    } else {
        os << "Mapping: \"" << *b.mapping << "\"";
    }
    os << ", priority = " << b.priority << ", isSuppl = " << (b.isSupplementary ? "true" : "false");
    if (detailed) {
        os << "\n";
        os << "Chain: " << b.chain;
        os << "RegionsForAln:\n";
        for (size_t i = 0; i < b.regionsForAln.size(); ++i) {
            os << "  [region " << i << "] " << b.regionsForAln[i] << "\n";
        }
    }
}

inline std::ostream& operator<<(std::ostream& os, const ChainedRegion& b) noexcept
{
    VerboseChainedRegion(os, b, false);
    return os;
}

class MapperBaseResult
{
public:
    std::vector<std::unique_ptr<ChainedRegion>> mappings;
    std::unordered_map<std::string, double> time;
};

inline bool operator==(const MapperBaseResult& lhs, const MapperBaseResult& rhs) noexcept
{
    if (lhs.mappings.size() != rhs.mappings.size()) {
        return false;
    }
    for (size_t i = 0; i < lhs.mappings.size(); ++i) {
        if ((lhs.mappings[i] == nullptr && rhs.mappings[i] == nullptr) ||
            (lhs.mappings[i] != nullptr && rhs.mappings[i] != nullptr &&
             *lhs.mappings[i] == *rhs.mappings[i])) {
            continue;
        }
        return false;
    }
    return true;
}

inline void VerboseMapperBaseResult(std::ostream& os, const MapperBaseResult& b,
                                    const bool detailed)
{
    for (size_t i = 0; i < b.mappings.size(); ++i) {
        os << "[mapping " << i << "] ";
        if (b.mappings[i] == nullptr) {
            os << "nullptr\n";
        } else {
            VerboseChainedRegion(os, *b.mappings[i], detailed);
            os << "\n";
        }
    }
}

inline std::ostream& operator<<(std::ostream& os, const MapperBaseResult& b) noexcept
{
    VerboseMapperBaseResult(os, b, false);
    return os;
}

class MapperBase
{
public:
    virtual ~MapperBase() = default;

    virtual std::vector<MapperBaseResult> MapAndAlign(
        const std::vector<std::string>& targetSeqs, const std::vector<std::string>& querySeqs) = 0;

    virtual std::vector<MapperBaseResult> MapAndAlign(
        const std::vector<FastaSequenceCached>& targetSeqs,
        const std::vector<FastaSequenceCached>& querySeqs) = 0;

    virtual std::vector<MapperBaseResult> MapAndAlign(
        const FastaSequenceCachedStore& targetSeqs, const FastaSequenceCachedStore& querySeqs) = 0;

    virtual MapperBaseResult MapAndAlignSingleQuery(
        const FastaSequenceCachedStore& targetSeqs, const PacBio::Pancake::SeedIndex& index,
        const FastaSequenceCached& querySeq,
        const std::vector<PacBio::Pancake::Int128t>& querySeeds, int32_t queryId,
        int64_t freqCutoff) = 0;
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_MAPPER_BASE_HPP
