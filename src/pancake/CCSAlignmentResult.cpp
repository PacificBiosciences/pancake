#include <pancake/CCSAlignmentResult.hpp>

#include <cstddef>
#include <ostream>
#include <stdexcept>
#include <tuple>

namespace PacBio {
namespace Pancake {

CCSAlignmentResult::CCSAlignmentResult(int32_t rId, bool isReversed, int64_t rStart, int64_t rEnd,
                                       int64_t qStart, int64_t qEnd, int64_t qLen,
                                       Data::Cigar cigar, uint8_t mapq, int32_t score,
                                       bool isAligned, bool isSupplementary, bool isSecondary)
    : Cigar{std::move(cigar)}
    , RefStart{rStart}
    , RefEnd{rEnd}
    , QueryStart{qStart}
    , QueryEnd{qEnd}
    , QueryLen{qLen}
    , RefId{rId}
    , Score{score}
    , MapQ{mapq}
    , IsReversed{isReversed}
    , IsAligned{isAligned}
    , IsSupplementary{isSupplementary}
    , IsSecondary{isSecondary}
{}

CCSAlignmentResult::CCSAlignmentResult(bool isReversed, int64_t rStart, int64_t rEnd,
                                       int64_t qStart, int64_t qEnd, int64_t qLen,
                                       Data::Cigar cigar)
    : CCSAlignmentResult{0,  isReversed, rStart, rEnd,  qStart, qEnd, qLen, std::move(cigar),
                         60, 0,          true,   false, false}
{}

bool CCSAlignmentResult::operator==(const CCSAlignmentResult& op) const noexcept
{
    return std::tie(RefId, IsReversed, RefStart, RefEnd, QueryStart, QueryEnd, QueryLen, Cigar,
                    MapQ, Score, IsAligned, IsSupplementary, IsSecondary) ==
           std::tie(op.RefId, op.IsReversed, op.RefStart, op.RefEnd, op.QueryStart, op.QueryEnd,
                    op.QueryLen, op.Cigar, op.MapQ, op.Score, op.IsAligned, op.IsSupplementary,
                    op.IsSecondary);
}

std::ostream& operator<<(std::ostream& out, const CCSAlignmentResult& a)
{
    out << "query" << '\t' << a.QueryLen << '\t' << a.QueryStart << '\t' << a.QueryEnd << '\t'
        << (a.IsReversed ? '-' : '+') << '\t' << a.RefId << '\t' << a.RefStart << '\t' << a.RefEnd
        << '\t' << static_cast<int>(a.MapQ) << '\t' << a.Score << '\t' << a.IsAligned << '\t'
        << a.IsSupplementary << '\t' << a.IsSecondary << '\t' << a.Cigar.ToStdString();
    return out;
}

}  // namespace Pancake
}  // namespace PacBio
