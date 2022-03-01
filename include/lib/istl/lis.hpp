/*
 * lis.hpp
 *
 *      Author: Ivan Sovic
 *      GitHub: @isovic
 *      Copyright: Ivan Sovic, 2017, 2020, 2021
 *      Licence: MIT
 *
 * A generic implementation of the Longest Increasing Subsequence algorithm
 * which allows for custom data types, provided that a suitable comparison
 * function is given.
 */

#ifndef ISTL_LIS_H_
#define ISTL_LIS_H_

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <functional>
#include <vector>

namespace istl {

template <class T, class COMP>
std::vector<T> LIS(const T* points, const int64_t nPoints, const COMP& compLessThan)
{
    // Sanity check.
    assert(nPoints >= 0);
    if (nPoints < 0) {
        return {};
    }
    if (nPoints == 0) {
        return {};
    }

    // Prepare the DP storage.
    std::vector<int64_t> dp(nPoints + 1, 0);
    std::vector<int64_t> pred(nPoints + 1, 0);
    int64_t len = 0;

    // Compute the LIS.
    for (int64_t i = 0; i < nPoints; ++i) {
        int32_t low = 1;
        int32_t high = len;
        while (low <= high) {
            const int32_t mid = (low + high) >> 1;
            if (compLessThan(points[dp[mid]], points[i])) {
                low = mid + 1;
            } else {
                high = mid - 1;
            }
        }
        const int32_t newLen = low;
        pred[i] = dp[newLen - 1];
        dp[newLen] = i;
        if (newLen > len) {
            len = newLen;
        }
    }

    // Backtrack.
    std::vector<T> lis(len);
    int64_t k = dp[len];
    for (int64_t i = (len - 1); i >= 0; --i) {
        lis[i] = points[k];
        k = pred[k];
    }

    return lis;
}

template <class T>
std::vector<T> LIS(const std::vector<T>& points, const int64_t begin, const int64_t end)
{
    assert(end >= begin);
    if (end < begin) {
        return {};
    }

    return LIS(points.data() + begin, (end - begin), std::less<T>{});
}

template <class T, class COMP>
std::vector<T> LIS(const std::vector<T>& points, const int64_t begin, const int64_t end,
                   const COMP& compLessThan)
{
    assert(end >= begin);
    if (end < begin) {
        return {};
    }

    return LIS(points.data() + begin, (end - begin), compLessThan);
}

template <class T, class COMP>
std::vector<T> LIS(const std::vector<T>& points, const COMP& compLessThan)
{
    return LIS(points.data(), points.size(), compLessThan);
}

}  // namespace istl

#endif
