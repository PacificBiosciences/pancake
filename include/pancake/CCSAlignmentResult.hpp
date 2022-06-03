#ifndef PANCAKE_CCS_ALIGNMENT_RESULT_HPP
#define PANCAKE_CCS_ALIGNMENT_RESULT_HPP

#include <pancake/MapperBase.hpp>

#include <pbcopper/data/Cigar.h>

#include <cstdint>
#include <iosfwd>
#include <memory>

namespace PacBio {
namespace Pancake {

struct CCSAlignmentResult
{
    CCSAlignmentResult(const ChainedRegion& chainedRegion);

    CCSAlignmentResult(int32_t rId, bool isReversed, int64_t rStart, int64_t rEnd, int64_t qStart,
                       int64_t qEnd, int64_t qLen, Data::Cigar cigar, uint8_t mapq, int32_t score,
                       bool isAligned, bool isSupplementary, bool isSecondary);

    CCSAlignmentResult(bool isReversed, int64_t rStart, int64_t rEnd, int64_t qStart, int64_t qEnd,
                       int64_t qLen, Data::Cigar cigar);

    bool operator==(const CCSAlignmentResult& op) const noexcept;

    Data::Cigar Cigar;
    int64_t RefStart = 0;
    int64_t RefEnd = 0;
    int64_t QueryStart = 0;
    int64_t QueryEnd = 0;
    int64_t QueryLen = 0;
    int32_t RefId = 0;
    int32_t Score = 0;
    uint8_t MapQ = 0;
    bool IsReversed = false;
    bool IsAligned = false;
    bool IsSupplementary = false;
    bool IsSecondary = false;
};

std::ostream& operator<<(std::ostream& out, const CCSAlignmentResult& a);

using CCSAlnResults = std::vector<std::unique_ptr<Pancake::CCSAlignmentResult>>;
}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_CCS_ALIGNMENT_RESULT_HPP
