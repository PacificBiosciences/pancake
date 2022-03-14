// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNMENT_SES_RESULTS_H
#define PANCAKE_ALIGNMENT_SES_RESULTS_H

#include <pancake/DiffCounts.hpp>

#include <pbbam/Cigar.h>
#include <pbbam/CigarOperation.h>

#include <cstdint>
#include <limits>
#include <memory>
#include <sstream>

namespace PacBio {
namespace Pancake {
namespace Alignment {

const int32_t MINUS_INF =
    std::numeric_limits<int32_t>::min() + 10000;  // We need a margin to avoid overflows.

class SESPathPoint
{
public:
    int32_t x = 0;
    int32_t y = 0;

    SESPathPoint() = default;
    SESPathPoint(int32_t _x, int32_t _y) : x(_x), y(_y) {}
    bool operator==(const SESPathPoint& b) const { return x == b.x && y == b.y; }
};

class SESTracebackPoint
{
public:
    // int32_t x1 = MINUS_INF;
    int32_t x2 = MINUS_INF;
    int32_t prevK = MINUS_INF;
};

class SESScratchSpace
{
public:
    std::vector<int32_t> v;             // Working row.
    std::vector<int32_t> u;             // Banding.
    std::vector<SESTracebackPoint> v2;  // Traceback matrix.
    std::vector<SESPathPoint> alnPath;
    std::vector<std::pair<int32_t, int32_t>> dStart;  // <row start, minK>
};

class SesResults
{
public:
    int32_t lastQueryPos = 0;
    int32_t lastTargetPos = 0;
    DiffCounts diffCounts;
    int32_t numDiffs = 0;
    bool valid = false;
    PacBio::BAM::Cigar cigar;

    SesResults() = default;
    SesResults(int32_t _lastQueryPos, int32_t _lastTargetPos, int32_t _diffs, bool _valid)
        : lastQueryPos(_lastQueryPos)
        , lastTargetPos(_lastTargetPos)
        , numDiffs(_diffs)
        , valid(_valid)
    {}
    SesResults(int32_t _lastQueryPos, int32_t _lastTargetPos, int32_t _diffs, int32_t _numEq,
               int32_t _numX, int32_t _numI, int32_t _numD, bool _valid,
               const PacBio::BAM::Cigar& _cigar)
        : lastQueryPos(_lastQueryPos)
        , lastTargetPos(_lastTargetPos)
        , diffCounts(DiffCounts(_numEq, _numX, _numI, _numD))
        , numDiffs(_diffs)
        , valid(_valid)
        , cigar(_cigar)
    {}
    SesResults(int32_t _lastQueryPos, int32_t _lastTargetPos, DiffCounts _diffCounts, bool _valid,
               const PacBio::BAM::Cigar& _cigar)
        : lastQueryPos(_lastQueryPos)
        , lastTargetPos(_lastTargetPos)
        , diffCounts(_diffCounts)
        , numDiffs(_diffCounts.NumDiffs())
        , valid(_valid)
        , cigar(_cigar)
    {}

    bool operator==(const SesResults& b) const
    {
        return lastQueryPos == b.lastQueryPos && lastTargetPos == b.lastTargetPos &&
               diffCounts == b.diffCounts && numDiffs == b.numDiffs && valid == b.valid &&
               cigar == b.cigar;
    }
    friend std::ostream& operator<<(std::ostream& os, const SesResults& r);
};
inline std::ostream& operator<<(std::ostream& os, const SesResults& a)
{
    os << "lastQueryPos = " << a.lastQueryPos << ", lastTargetPos = " << a.lastTargetPos
       << ", diffs = " << a.numDiffs << ", numEq = " << a.diffCounts.numEq
       << ", numX = " << a.diffCounts.numX << ", numI = " << a.diffCounts.numI
       << ", numD = " << a.diffCounts.numD << ", valid = " << a.valid << ", cigar = '"
       << a.cigar.ToStdString() << "'";
    return os;
}
}  // namespace Alignment
}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ALIGNMENT_SES_RESULTS_H
