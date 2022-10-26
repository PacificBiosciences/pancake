// Author: Ivan Sovic

#ifndef PANCAKE_MAPPER_BASE_HPP
#define PANCAKE_MAPPER_BASE_HPP

#include <pbcopper/utility/Ssize.h>
#include <pancake/AlignmentSeeded.hpp>
#include <pancake/DPChain.hpp>
#include <pancake/FastaSequenceCached.hpp>
#include <pancake/FastaSequenceId.hpp>
#include <pancake/Overlap.hpp>
#include <pancake/OverlapWriterBase.hpp>
#include <pancake/SeedIndex.hpp>
#include <pancake/SequenceSeedsCached.hpp>

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
    bool isMockedMapping = false;    // True if the mapping was mocked.
    bool isMockedAlignment = false;  // True if the alignment was mocked.
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

    virtual MapperBaseResult MapAndAlignSingleQuery(const FastaSequenceCachedStore& targetSeqs,
                                                    const PacBio::Pancake::SeedIndex& index,
                                                    const FastaSequenceCached& querySeq,
                                                    const SequenceSeedsCached& querySeeds,
                                                    const int32_t queryId, int64_t freqCutoff) = 0;

    virtual std::vector<MapperBaseResult> MapAndAlign(
        const FastaSequenceCachedStore& targetSeqs, const FastaSequenceCachedStore& querySeqs) = 0;

    /*
     * \brief Runs mapping and alignment of one or more query sequences to one or more target sequences.
     * This is the basic interface for the most simple usage.
     * This interface does not require the SeedIndex or minimizers, because it will compute them internally.
     * The std::string objects are first converted to FastaSequenceCached, and then the MapAndAlign overload
     * is called on these inputs.
    */
    virtual std::vector<MapperBaseResult> MapAndAlign(const std::vector<std::string>& targetSeqs,
                                                      const std::vector<std::string>& querySeqs);

    /*
     * \brief Runs mapping and alignment of one or more query sequences to one or more target sequences.
     * This interface does not require the SeedIndex or minimizers, because it will compute them internally.
     * This function constructs the FastaSequenceCachedStore from the given FastaSequenceCached objects,
     * and calls the related overload of MapAndAlign with these inputs.
    */
    virtual std::vector<MapperBaseResult> MapAndAlign(
        const std::vector<FastaSequenceCached>& targetSeqs,
        const std::vector<FastaSequenceCached>& querySeqs);
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_MAPPER_BASE_HPP
