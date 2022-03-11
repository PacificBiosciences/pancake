// Author: Ivan Sovic

#ifndef PANCAKE_DIFF_COUNTS_H
#define PANCAKE_DIFF_COUNTS_H

#include <algorithm>
#include <cstdint>
#include <ostream>

namespace PacBio {
namespace Pancake {
namespace Alignment {

class DiffCounts
{
public:
    int32_t numEq = 0;
    int32_t numX = 0;
    int32_t numI = 0;
    int32_t numD = 0;

    DiffCounts() = default;
    DiffCounts(int32_t _numEq, int32_t _numX, int32_t _numI, int32_t _numD)
        : numEq(_numEq), numX(_numX), numI(_numI), numD(_numD)
    {}

    void Clear() { numEq = numX = numI = numD = 0; }

    int32_t NumDiffs() const { return numX + numI + numD; }

    DiffCounts operator+(const DiffCounts& b) const
    {
        DiffCounts ret;
        ret.numEq = numEq + b.numEq;
        ret.numX = numX + b.numX;
        ret.numI = numI + b.numI;
        ret.numD = numD + b.numD;
        return ret;
    }

    DiffCounts& operator+=(const DiffCounts& b)
    {
        numEq += b.numEq;
        numX += b.numX;
        numI += b.numI;
        numD += b.numD;
        return *this;
    }

    bool operator==(const DiffCounts& rhs) const
    {
        return numEq == rhs.numEq && numX == rhs.numX && numI == rhs.numI && numD == rhs.numD;
    }

    void Identity(bool noSNPs, bool noIndels, float& retIdentity, int32_t& retEditDist) const
    {
        retIdentity = 0.0;
        retEditDist = 0;
        const float span = (numEq + numX) + std::max(numI, numD);
        if (span == 0.0f) {
            return;
        }
        retEditDist = EditDistance(noSNPs, noIndels);
        retIdentity = (span - static_cast<float>(retEditDist)) / span;
    }

    int32_t EditDistance(bool noSNPs, bool noIndels) const
    {
        return (noSNPs ? 0 : numX) + (noIndels ? 0 : (numI + numD));
    }

    friend std::ostream& operator<<(std::ostream& os, const DiffCounts& r);
};
inline std::ostream& operator<<(std::ostream& os, const DiffCounts& a)
{
    os << "numEq = " << a.numEq << ", numX = " << a.numX << ", numI = " << a.numI
       << ", numD = " << a.numD;
    return os;
}
}  // namespace Alignment
}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_DIFF_COUNTS_H
