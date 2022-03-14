// Author: Ivan Sovic

#ifndef PANCAKE_OVERLAPHIFI_SECONDARY_H
#define PANCAKE_OVERLAPHIFI_SECONDARY_H

#include <pancake/Overlap.hpp>

#include <pancake/third-party/intervaltree/IntervalTree.h>

#include <algorithm>
#include <cstdint>
#include <functional>
#include <memory>
#include <unordered_map>
#include <vector>

namespace PacBio {
namespace Pancake {

using IntervalTreeInt32 =
    interval_tree::IntervalTree<int32_t,
                                int32_t>;  // First: interval scalar type, Second: value type.
using IntervalVectorInt32 = IntervalTreeInt32::interval_vector;
using IntervalInt32 = IntervalTreeInt32::interval;

struct OverlapPriority
{
    int32_t priority = 0;
    bool isSupplementary = false;

    bool operator==(const OverlapPriority& b) const
    {
        return priority == b.priority && isSupplementary == b.isSupplementary;
    }
};

std::vector<OverlapPriority> FlagSecondaryAndSupplementary(std::vector<OverlapPtr>& overlaps,
                                                           double allowedOverlapFractionQuery,
                                                           double allowedOverlapFractionTarget,
                                                           double minSecondaryScoreFraction);

void CreateRegionIntervalTrees(const std::vector<OverlapPtr>& overlaps,
                               const std::vector<OverlapPriority>& priorities,
                               std::function<bool(int32_t a)> CompPriority,
                               IntervalVectorInt32& queryIntervals, IntervalTreeInt32& queryTrees,
                               std::unordered_map<int32_t, IntervalVectorInt32>& targetIntervals,
                               std::unordered_map<int32_t, IntervalTreeInt32>& targetTrees);

bool CheckRegionSupplementary(const std::vector<OverlapPtr>& overlaps, const OverlapPtr& currentOvl,
                              IntervalTreeInt32& queryTrees,
                              std::unordered_map<int32_t, IntervalTreeInt32>& targetTrees,
                              double allowedOverlapFractionQuery,
                              double allowedOverlapFractionTarget);

int32_t CalcIntervalOverlap(int32_t s1, int32_t e1, int32_t s2, int32_t e2);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_OVERLAPHIFI_SECONDARY_H
