// Author: Ivan Sovic

#include <pancake/SesDistanceBanded.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#define DEBUG_SES_ALIGN_GLOBAL

namespace PacBio {
namespace Pancake {
namespace Alignment {

SesResults SESDistanceBanded(const char* query, size_t queryLen, const char* target,
                             size_t targetLen, int32_t maxDiffs, int32_t bandwidth)
{
    SesResults ret;

    int32_t N = queryLen;
    int32_t M = targetLen;

    bandwidth = std::min(bandwidth, maxDiffs);

    const int32_t maxAllowedDiffs = std::max(maxDiffs, bandwidth);
    int32_t zero_offset = maxAllowedDiffs + 1;
    const int32_t rowLen = 2 * maxAllowedDiffs + 3;
    std::vector<int32_t> v(rowLen, 0);
    std::vector<int32_t> u(rowLen, MINUS_INF);

    int32_t bandTolerance = bandwidth / 2 + 1;
    int32_t minK = 0;
    int32_t maxK = 0;
    int32_t best_u = 0;

    for (int32_t d = 0; d < maxDiffs; ++d) {
        ret.numDiffs = d;
        if ((maxK - minK) > bandwidth) {
            ret.valid = false;
            break;
        }

        int32_t x = 0;
        for (int32_t k = minK; k <= maxK; k += 2) {
            int32_t kz = k + zero_offset;
            if (k == minK || (k != maxK && v[kz - 1] < v[kz + 1])) {
                x = v[kz + 1];
            } else {
                x = v[kz - 1] + 1;
            }

            int32_t y = x - k;
            while (x < N && y < M && query[x] == target[y]) {
                ++x;
                ++y;
            }
            v[kz] = x;
            u[kz] = x + y;

            if (best_u <= u[kz]) {
                best_u = u[kz];
                ret.lastQueryPos = x;
                ret.lastTargetPos = y;
            }

            if (x >= N || y >= M) {
                ret.valid = true;
                ret.lastQueryPos = x;
                ret.lastTargetPos = y;
                break;
            }
        }

        if (ret.valid) {
            break;
        }

        int32_t new_minK = maxK;
        int32_t new_maxK = minK;
        for (int32_t k = minK; k <= maxK; k += 2) {
            // Is there a bug here? Should this also have '&& u[k + zero_offset] <= (best_u + bandTolerance'?
            if (u[k + zero_offset] >= (best_u - bandTolerance)) {
                new_minK = std::min(k, new_minK);
                new_maxK = std::max(k, new_maxK);
            }
        }
        minK = new_minK - 1;
        maxK = new_maxK + 1;
    }

    return ret;
}
}  // namespace Alignment
}  // namespace Pancake
}  // namespace PacBio
