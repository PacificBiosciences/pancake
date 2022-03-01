// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNMENT_SES2_ALIGN_BANDED_H
#define PANCAKE_ALIGNMENT_SES2_ALIGN_BANDED_H

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <vector>

#include <pacbio/alignment/SesOptions.h>
#include <pacbio/alignment/SesResults.h>

// #define SES2_DEBUG

namespace PacBio {
namespace Pancake {
namespace Alignment {

template <SESAlignMode ALIGN_MODE, SESTrimmingMode TRIM_MODE, SESTracebackMode TRACEBACK>
SesResults SES2AlignBanded(const char* query, size_t queryLen, const char* target, size_t targetLen,
                           int32_t maxDiffs, int32_t bandwidth,
                           std::shared_ptr<SESScratchSpace> ss = nullptr)
{
    SesResults ret;

    if (queryLen == 0 || targetLen == 0) {
        ret.valid = true;
        return ret;
    }

    // Allocate scratch space memory if required.
    if (ss == nullptr) {
        ss = std::make_shared<SESScratchSpace>();
    }

    bandwidth = std::min(bandwidth, maxDiffs);

    // Define the required variables.
    const int32_t maxAllowedDiffs = std::max(maxDiffs, bandwidth);
    const int32_t qlen = queryLen;
    const int32_t tlen = targetLen;
    const int32_t zero_offset = maxAllowedDiffs + 1;
    const int32_t bandTolerance = bandwidth / 2 + 1;
    const int32_t rowLen = (2 * maxAllowedDiffs + 3);

    // Working space for regular alignment (without traceback).
    auto& W = ss->v;  // Y for a diagonal k. 'W' is taken from the pseudocode. Working row.

    // Trimming related options.
    std::vector<uint64_t> B;  // Bitmask for trimming.
    std::vector<int32_t> M;   // Match count.
    uint64_t b = 0;
    int32_t m = 0;
    const uint64_t C = 60;
    const uint64_t MASKC = static_cast<uint64_t>(1) << (C - 1);

    // Banding info.
    auto& u = ss->u;
    int32_t minK = 0;
    int32_t maxK = 0;
    int32_t best_u = 0;

    // Traceback info.
    int32_t lastK = 0;
    int32_t lastD = 0;
    int32_t prevK = -1;
    auto& WMatrix =
        ss->v2;  // Traceback matrix, implemented as a flat vector. We track the start of each row with dStart.
    auto& dStart =
        ss->dStart;  // Start of each diff's row in the WMatrix vector. dStart[d] = <WMatrixPos, minK>,
    // where WMatrixPos is the index of the element in the WMatrix's flat vector where the row begins,
    // and minK is the banding related minimum K for the inner loop.
    auto& alnPath = ss->alnPath;  // Alignment path during traceback.
    int32_t WMatrixPos =
        0;  // Tracks the current location in the WMatrix (which is implemented as a flat vector).

    // A useless void cast to prevent the compiler from complaining
    // about unused variables when the constexpr if condition is not met.
    (void)lastK;
    (void)lastD;
    (void)prevK;

    // Allocate memory for basic alignment.
    if (rowLen > static_cast<int32_t>(W.capacity())) {
        W.resize(rowLen, -1);
        u.resize(rowLen, MINUS_INF);
    }
    // Allocate the memory for trimming.
    if constexpr (TRIM_MODE == SESTrimmingMode::Enabled) {
        B.resize(rowLen, MINUS_INF);
        M.resize(rowLen, MINUS_INF);
    }
    // Allocate memory for traceback.
    // clang-format off
    if constexpr (TRACEBACK == SESTracebackMode::Enabled) {
        if (rowLen > static_cast<int32_t>(dStart.capacity())) {
            dStart.resize(rowLen, {0, 0});
        }
        if (rowLen > static_cast<int32_t>(alnPath.capacity())) {
            alnPath.resize(rowLen);
        }
        if ((rowLen * maxDiffs) > static_cast<int32_t>(WMatrix.capacity())) {
            WMatrix.resize(rowLen * maxDiffs);
        }
    }
    // clang-format on

    // Initialize the alignment vectors.
    u[zero_offset] = u[zero_offset + 1] = u[zero_offset - 1] = MINUS_INF;
    W[zero_offset] = W[zero_offset + 1] = W[zero_offset - 1] = -1;

    for (int32_t d = 0; d < maxDiffs; ++d) {
        ret.numDiffs = d;
        if ((maxK - minK) > bandwidth) {
            ret.valid = false;
            break;
        }

        // clang-format off
        if constexpr (TRACEBACK == SESTracebackMode::Enabled) {
            // Location where to store the traceback info.
            // Each row is wide at most as the number of diffs.
            // Because of banding, we start filling up the row from the beginning which
            // corresponds to minK, and we need to keep track of what minK is.
            dStart[d] = {WMatrixPos, minK};
        }
        // clang-format on

        W[zero_offset + minK - 1] = W[zero_offset + maxK + 1] = W[zero_offset + maxK + 0] = -1;
        u[zero_offset + minK - 1] = u[zero_offset + maxK + 1] = MINUS_INF;

        int32_t ym = MINUS_INF;
        int32_t yc = MINUS_INF;
        int32_t yp = -1;
        int32_t y = MINUS_INF;

#ifdef SES2_DEBUG
        std::cerr << "\n";
        std::cerr << "W[" << minK << " - 1, " << maxK << " + 1]:";
        for (int32_t k = (minK - 1); k <= (maxK + 1); ++k) {
            std::cerr << " " << W[zero_offset + k];
        }
        std::cerr << "\n";
#endif

        for (int32_t k = (minK); k <= (maxK); ++k) {
            int32_t kz = k + zero_offset;

            ym = yc;
            yc = yp;
            yp = W[kz + 1];

            int32_t maxY = std::max(yc, std::max(ym, yp));

#ifdef SES2_DEBUG
            std::cerr << "[d = " << d << ", k = " << k << "] ym = " << ym << ", yc = " << yc
                      << ", yp = " << yp << ", maxY = " << maxY;
#endif

            if (yc == maxY && yc < tlen) {
                y = yc + 1;
                // clang-format off
                if constexpr (TRACEBACK == SESTracebackMode::Enabled) {
                    prevK = k;
#ifdef SES2_DEBUG
                    std::cerr << ": (else) y = yc + 1 = " << y << ", prevK = " << prevK;
#endif
                }
                if constexpr (TRIM_MODE == SESTrimmingMode::Enabled) {
                    m = M[kz];
                    b = B[kz];
                }
                // clang-format on
            } else if (k == minK || (k != maxK && yp == maxY) || yc >= tlen) {
                y = yp +
                    1;  // Unlike 1986 paper, here we update y instead of x, so the +1 goes to the move to right (yp) instead of down (ym).
                // clang-format off
                if constexpr (TRACEBACK == SESTracebackMode::Enabled) {
                    prevK = k + 1;
#ifdef SES2_DEBUG
                    std::cerr << ": (yp) y = yp + 1 = " << y << ", prevK = " << prevK;
#endif
                }
                if constexpr (TRIM_MODE == SESTrimmingMode::Enabled) {
                    m = M[kz + 1];
                    b = B[kz + 1];
                }
                // clang-format on
            } else {
                y = ym;
                // clang-format off
                if constexpr (TRACEBACK == SESTracebackMode::Enabled) {
                    prevK = k - 1;
#ifdef SES2_DEBUG
                    std::cerr << ": (ym) y = ym = " << y << ", prevK = " << prevK;
#endif
                }
                if constexpr (TRIM_MODE == SESTrimmingMode::Enabled) {
                    m = M[kz - 1];
                    b = B[kz - 1];
                }
                // clang-format on
            }

            // clang-format off
            if constexpr (TRIM_MODE == SESTrimmingMode::Enabled) {
                if ((b & MASKC) != 0) {
                    --m;
                }
                b <<= 1;
            }
            // clang-format on

            int32_t x = y + k;
            int32_t minLeft = std::min(qlen - x, tlen - y);
            const char* querySub = query + x;
            const char* targetSub = target + y;
            int32_t moves = 0;

            while (moves < minLeft && querySub[moves] == targetSub[moves]) {
                ++moves;
                if constexpr (TRIM_MODE == SESTrimmingMode::Enabled) {
                    if ((b & MASKC) == 0) {
                        ++m;
                    }
                    b = (b << 1) | 1;
                }
            }
            y += moves;
            x += moves;
            W[kz] = y;
            u[kz] = y + k + y;  // x + y = 2*y + k

            // clang-format off
            if constexpr (TRACEBACK == SESTracebackMode::Enabled) {
                WMatrix[WMatrixPos] = {x, prevK};
                ++WMatrixPos;
            }
            if constexpr (TRIM_MODE == SESTrimmingMode::Enabled) {
                M[kz] = m;
                B[kz] = b;
            }
            // clang-format on

            lastK = k;
            lastD = d;
            ret.lastQueryPos = x;
            ret.lastTargetPos = y;

            if (best_u <= u[kz]) {
                best_u = u[kz];
            }

#ifdef SES2_DEBUG
            std::cerr << "; x2 = " << x << ", y2 = " << y << ", u[kz] = " << u[kz]
                      << ", lastK = " << lastK << ", lastD = " << lastD;
            std::cerr << "\n";
#endif

            // clang-format off
            if constexpr(ALIGN_MODE == SESAlignMode::Global) {
                if (x >= qlen && y >= tlen) {
                    ret.valid = true;
                    break;
                }

            } else {
                if (x >= qlen || y >= tlen) {
                    ret.valid = true;
                    break;
                }
            }
            // clang-format on
        }

        if (ret.valid) {
            break;
        }

        int32_t newMinK = maxK;
        int32_t newMaxK = minK;
        for (int32_t k = (minK - 1); k <= (maxK + 1); ++k) {
            // Is there a bug here? Should this also have '&& u[k + zero_offset] <= (best_u + bandTolerance'?
            if (u[k + zero_offset] >= (best_u - bandTolerance)) {
                newMinK = std::min(k, newMinK);
                newMaxK = std::max(k, newMaxK);
            }
        }
        minK = newMinK - 1;
        maxK = newMaxK + 1;
    }

    // clang-format off
    if constexpr (TRACEBACK == SESTracebackMode::Enabled) {

#ifdef SES2_DEBUG
        for (int32_t d = 1; d <= lastD; ++d) {
            int32_t b = dStart[d-1].first;
            int32_t e = dStart[d].first;
            int32_t minK = dStart[d-1].second;
            int32_t maxK = minK + (e - b);
            std::cerr << "[d = " << (d - 1) << "] b = " << b << ", e = " << e << ", minK = " << minK << ", maxK = " << maxK << "\n";
            std::cerr << "    ";
            for (int32_t k = minK; k < maxK; ++k) {
                const auto& w = WMatrix[b + k - minK];
                std::cerr << "\t{k = " << k << ", w.x2 = " << w.x2 << ", y2 = " << (w.x2 - k) << ", w.prevK = " << w.prevK << "}\n";
            }
            // std::cerr << "\n";
            std::cerr << "\n";
        }
        std::cerr << "[d = " << lastD << "] b = " << dStart[lastD].first << ", minK = " << minK << ", lastK = " << lastK << "\n";
        std::cerr << "";
        for (int32_t k = dStart[lastD].second; k <= lastK; ++k) {
            const auto& w = WMatrix[dStart[lastD].first + k - dStart[lastD].second];
            std::cerr << "\t{k = " << k << ", w.x2 = " << w.x2 << ", y2 = " << (w.x2 - k) << ", w.prevK = " << w.prevK << "}\n";
        }
        std::cerr << "lastD = " << lastD << ", lastK = " << lastK << "\n";
        std::cerr << "\n";
#endif

        auto AppendToCigar = [](PacBio::BAM::Cigar& cigar, PacBio::BAM::CigarOperationType newOp, int32_t newLen)
        {
            if (cigar.empty() || newOp != cigar.back().Type()) {
                cigar.emplace_back(PacBio::BAM::CigarOperation(newOp, newLen));
            } else {
                cigar.back().Length(cigar.back().Length() + newLen);
            }
        };

        int32_t currD = lastD;
        int32_t currK = lastK;
        ret.cigar.clear();
        ret.cigar.reserve(currD);
        while (currD > 0) {
            int32_t currRowStart = dStart[currD].first;
            int32_t currMinK = dStart[currD].second;
            const auto& currW = WMatrix[currRowStart + currK - currMinK];
            int32_t x2 = currW.x2;
            int32_t y2 = currW.x2 - currK;
            int32_t trPrevK = currW.prevK;

            int32_t prevRowStart = dStart[currD - 1].first;
            int32_t prevMinK = dStart[currD - 1].second;
            const auto& prevW = WMatrix[prevRowStart + trPrevK - prevMinK];
            int32_t prevX2 = prevW.x2;
            int32_t prevY2 = prevW.x2 - trPrevK;

            if (currK > trPrevK) {
                int32_t matches = std::min(x2 - prevX2, y2 - prevY2);
                if (matches > 0) {
                    AppendToCigar(ret.cigar, PacBio::BAM::CigarOperationType::SEQUENCE_MATCH, matches);
                }
                AppendToCigar(ret.cigar, PacBio::BAM::CigarOperationType::INSERTION, 1);
                ret.diffCounts.numEq += matches;
                ++ret.diffCounts.numI;
            } else if (currK < trPrevK) {
                int32_t matches = std::min(x2 - prevX2, y2 - prevY2);
                if (matches > 0) {
                    AppendToCigar(ret.cigar, PacBio::BAM::CigarOperationType::SEQUENCE_MATCH, matches);
                }
                AppendToCigar(ret.cigar, PacBio::BAM::CigarOperationType::DELETION, 1);
                ret.diffCounts.numEq += matches;
                ++ret.diffCounts.numD;
            } else {
                int32_t matches = std::min(x2 - prevX2, y2 - prevY2) - 1;
                if (matches > 0) {
                    AppendToCigar(ret.cigar, PacBio::BAM::CigarOperationType::SEQUENCE_MATCH, matches);
                }
                AppendToCigar(ret.cigar, PacBio::BAM::CigarOperationType::SEQUENCE_MISMATCH, 1);
                ret.diffCounts.numEq += matches;
                ++ret.diffCounts.numX;
            }
            currK = trPrevK;
            --currD;
        }
        {
            int32_t currRowStart = dStart[currD].first;
            int32_t currMinK = dStart[currD].second;
            const auto& currW = WMatrix[currRowStart + currK - currMinK];
            int32_t x2 = currW.x2;
            int32_t y2 = currW.x2 - currK;
            int32_t prevX2 = 0;
            int32_t prevY2 = 0;

            int32_t matches = std::min(x2 - prevX2, y2 - prevY2);
            if (matches > 0) {
                AppendToCigar(ret.cigar, PacBio::BAM::CigarOperationType::SEQUENCE_MATCH, matches);
            }
            ret.diffCounts.numEq += matches;
        }
        ret.numDiffs = ret.diffCounts.NumDiffs();

        std::reverse(ret.cigar.begin(), ret.cigar.end());
    }
    // clang-format on

    return ret;
}
}  // namespace Alignment
}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ALIGNMENT_SES2_ALIGN_BANDED_H
