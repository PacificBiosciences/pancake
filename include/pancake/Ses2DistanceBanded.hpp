// Author: Ivan Sovic

#ifndef PANCAKE_SES2DISTANCE_BANDED_HPP
#define PANCAKE_SES2DISTANCE_BANDED_HPP

#include <pancake/SesOptions.hpp>
#include <pancake/SesResults.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <vector>

// #define SES2_DEBUG

namespace PacBio {
namespace Pancake {

template <SESAlignMode ALIGN_MODE, SESTrimmingMode TRIM_MODE>
SesResults SES2DistanceBanded(std::string_view query, std::string_view target,
                              const int32_t maxDiffs, int32_t bandwidth)
{
    SesResults ret;

    const char* queryData = query.data();
    const char* targetData = target.data();
    const int32_t queryLen = query.size();
    const int32_t targetLen = target.size();

    if (queryLen == 0 || targetLen == 0) {
        ret.valid = true;
        return ret;
    }

    bandwidth = std::min(bandwidth, maxDiffs);

    const int32_t maxAllowedDiffs = std::max(maxDiffs, bandwidth);
    const int32_t qlen = queryLen;
    const int32_t tlen = targetLen;
    const int32_t zero_offset = maxAllowedDiffs + 1;
    const int32_t rowLen = (2 * maxAllowedDiffs + 3);
    std::vector<int32_t> W(rowLen, MINUS_INF);  // Y for a diagonal k.
    std::vector<uint64_t> B;                    // Bitmask for trimming.
    std::vector<int32_t> M;                     // Match count.

    if constexpr (TRIM_MODE == SESTrimmingMode::Enabled) {
        B.resize(rowLen, MINUS_INF);
        M.resize(rowLen, MINUS_INF);
    }

    // Banding info.
    std::vector<int32_t> u(rowLen, MINUS_INF);
    const int32_t bandTolerance = bandwidth / 2 + 1;
    int32_t minK = 0;
    int32_t maxK = 0;
    int32_t best_u = 0;

    // Trimming related options.
    uint64_t b = 0;
    int32_t m = 0;
    const uint64_t C = 60;
    const uint64_t MASKC = static_cast<uint64_t>(1) << (C - 1);

    for (int32_t d = 0; d < maxDiffs; ++d) {
        ret.numDiffs = d;
        if ((maxK - minK) > bandwidth) {
            ret.valid = false;
            break;
        }

        W[zero_offset + minK - 1] = W[zero_offset + maxK + 1] = W[zero_offset + maxK + 0] = -1;

        int32_t ym = MINUS_INF;
        int32_t yc = MINUS_INF;
        int32_t yp = -1;
        int32_t y = MINUS_INF;

        for (int32_t k = (minK); k <= (maxK); ++k) {
            int32_t kz = k + zero_offset;

            ym = yc;
            yc = yp;
            yp = W[kz + 1];

            int32_t maxY = std::max(yc, std::max(ym, yp));

            if (yc == maxY && yc < tlen) {
                y = yc + 1;
                // clang-format off
                if constexpr (TRIM_MODE == SESTrimmingMode::Enabled) {
                    m = M[kz];
                    b = B[kz];
                }
                // clang-format on
            } else if (k == minK || (k != maxK && yp == maxY) || yc >= tlen) {
                // Unlike 1986 paper, here we update y instead of x, so the +1 goes to the move to right (yp) instead of down (ym).
                y = yp + 1;
                // clang-format off
                if constexpr (TRIM_MODE == SESTrimmingMode::Enabled) {
                    m = M[kz + 1];
                    b = B[kz + 1];
                }
                // clang-format on
            } else {
                y = ym;
                // clang-format off
                if constexpr (TRIM_MODE == SESTrimmingMode::Enabled) {
                    m = M[kz - 1];
                    b = B[kz - 1];
                }
                // clang-format on
            }

            int32_t x = y + k;
            int32_t minLeft = std::min(qlen - x, tlen - y);
            const char* querySub = queryData + x;
            const char* targetSub = targetData + y;
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
            if constexpr (TRIM_MODE == SESTrimmingMode::Enabled) {
                M[kz] = m;
                B[kz] = b;
            }

            u[kz] = y + k + y;  // x + y = 2*y + k
            if (best_u <= u[kz]) {
                best_u = u[kz];
                ret.lastQueryPos = x;
                ret.lastTargetPos = y;
            }

            if constexpr (ALIGN_MODE == SESAlignMode::Global) {
                if (x >= qlen && y >= tlen) {
                    ret.valid = true;
                    ret.lastQueryPos = x;
                    ret.lastTargetPos = y;
                    break;
                }

            } else {
                if (x >= qlen || y >= tlen) {
                    ret.valid = true;
                    ret.lastQueryPos = x;
                    ret.lastTargetPos = y;
                    break;
                }
            }
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

    return ret;
}
}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SES2DISTANCE_BANDED_HPP
