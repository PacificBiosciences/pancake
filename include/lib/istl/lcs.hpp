/*
 * lcs.hpp
 *
 *      Author: Ivan Sovic
 *      GitHub: @isovic
 *      Copyright: Ivan Sovic, 2017, 2020, 2021, 2022
 *      Licence: MIT
 *
 * A generic implementation of the Longest Common Subsequence algorithm
 * which allows for custom numeric data types.
 *
 * Based on the algorithm described here:
 * https://arxiv.org/abs/1407.2407
 * https://github.com/fpavetic/lcskpp
 */

#ifndef ISTL_LCS_H_
#define ISTL_LCS_H_

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <functional>
#include <vector>
#include "fenwick.hpp"

namespace istl {

template <typename T>
std::vector<T> LCSkCore(const std::vector<std::pair<T, T>>& points, const int k)
{
    if (points.empty()) {
        return {};
    }

    const int32_t numMatches = points.size();
    std::vector<std::tuple<int32_t, int32_t, int32_t>> events;
    std::vector<int64_t> pointIndices(points.size(), -1);
    events.reserve(2 * numMatches);

    int32_t n = 0;
    for (int32_t idx = 0; idx < numMatches; ++idx) {
        const auto& point = points[idx];
        events.push_back(std::make_tuple(point.first, point.second, idx + numMatches));
        events.push_back(std::make_tuple(point.first + k, point.second + k, idx));
        n = std::max(n, point.first + k);
        n = std::max(n, point.second + k);
    }
    std::sort(events.begin(), events.end());

    // Fenwick tree: <first = dp value, second = index in points>.
    FenwickMax<std::pair<int32_t, int32_t>> dpColumnMax(n);
    std::vector<int32_t> dp(points.size());
    std::vector<int32_t> pred(points.size());
    std::vector<int32_t> continues(points.size(), -1);

    // For k > 1 kmers can overlap. This algorithm can extend overlapping matches.
    if (k > 1) {
        for (size_t i = 0; i < points.size(); ++i) {
            const auto& p = points[i];
            auto G = std::make_pair(p.first - 1, p.second - 1);
            const auto prev = std::lower_bound(points.begin(), points.end(), G);
            if (*prev == G) {
                continues[i] = prev - points.begin();
            }
        }
    }

    int32_t bestIdx = 0;
    int32_t lcsLen = 0;
    for (const auto& event : events) {
        const auto& [i, j, rawIdx] = event;
        const int32_t idx = (rawIdx >= numMatches) ? (rawIdx - numMatches) : rawIdx;
        const bool isBeginning = rawIdx >= numMatches;

        if (isBeginning) {
            const std::pair<int, int> prevDP = dpColumnMax.get(j);
            dp[idx] = k;
            pred[idx] = -1;
            if (prevDP.first > 0) {
                dp[idx] = prevDP.first + k;
                pred[idx] = prevDP.second;
            }

        } else {
            if (continues[idx] != -1 && (dp[continues[idx]] + 1) > dp[idx]) {
                dp[idx] = dp[continues[idx]] + 1;
                pred[idx] = continues[idx];
            }

            dpColumnMax.update(j, std::make_pair(dp[idx], idx));

            if (dp[idx] > lcsLen) {
                lcsLen = dp[idx];
                bestIdx = idx;
            }
        }
    }

    // Traceback.
    std::vector<int32_t> retIndices;
    retIndices.reserve(points.size());
    int32_t tbIndex = bestIdx;
    while (tbIndex >= 0) {
        retIndices.emplace_back(tbIndex);
        tbIndex = pred[tbIndex];
    }
    std::reverse(retIndices.begin(), retIndices.end());

    return retIndices;
}

template <typename T>
std::vector<std::pair<T, T>> LCSk(const std::vector<std::pair<T, T>>& points, const T k)
{
    const std::vector<int32_t> lcsIndices = LCSkCore(points, k);
    std::vector<std::pair<T, T>> lcs;
    lcs.reserve(points.size());
    for (const auto idx : lcsIndices) {
        lcs.emplace_back(points[idx]);
    }
    return lcs;
}

template <typename T>
std::vector<int32_t> LCSkIndices(const std::vector<std::pair<T, T>>& points, const T k)
{
    return LCSkCore(points, k);
}
}  // namespace istl

#endif