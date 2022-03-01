// Authors: Ivan Sovic

#include <pacbio/pancake/Secondary.h>

#include <sstream>

namespace PacBio {
namespace Pancake {

std::vector<OverlapPriority> FlagSecondaryAndSupplementary(std::vector<OverlapPtr>& overlaps,
                                                           double allowedOverlapFractionQuery,
                                                           double allowedOverlapFractionTarget,
                                                           double minSecondaryScoreFraction)
{
    if (overlaps.empty()) {
        return {};
    }

    std::vector<OverlapPriority> priorities(overlaps.size());

    // Find the maximum score.
    size_t maxScoreId = 0;
    for (size_t i = 0; i < overlaps.size(); ++i) {
        auto& ovl = overlaps[i];
        if (ovl == nullptr) {
            throw std::runtime_error("Overlap is nullptr in FlagSecondaryAndSupplementary!");
        }
        if (std::abs(ovl->Score) > std::abs(overlaps[maxScoreId]->Score)) {
            maxScoreId = i;
        }
    }
    const double topScore = std::abs(static_cast<double>(overlaps[maxScoreId]->Score));
    const double minSecondaryScore = topScore * minSecondaryScoreFraction;

    // Mark every region except the maximum as secondary.
    priorities[maxScoreId].priority = 0;
    priorities[maxScoreId].isSupplementary = false;
    for (size_t i = 0; i < priorities.size(); ++i) {
        if (i == maxScoreId) {
            continue;
        }
        priorities[i].priority = 1;
        priorities[i].isSupplementary = false;
    }

    IntervalVectorInt32 queryIntervals;
    std::unordered_map<int32_t, IntervalVectorInt32> targetIntervals;
    IntervalTreeInt32 queryTrees;
    std::unordered_map<int32_t, IntervalTreeInt32> targetTrees;

    // Initial construction of the interval trees.
    CreateRegionIntervalTrees(
        overlaps, priorities, [](int32_t a) { return a == 0; }, queryIntervals, queryTrees,
        targetIntervals, targetTrees);

    // Loop through all mappings, and check if they overlap the primary mappings and by how much.
    // If the overlap amount is below the specified fraction, mark the region as supplementary and
    // add it to the interval trees.
    for (size_t i = 0; i < overlaps.size(); ++i) {
        const auto& aln = overlaps[i];
        auto& alnPriority = priorities[i];
        if (alnPriority.priority == 0) {
            continue;
        }
        const double score = std::abs(static_cast<double>(aln->Score));

        if (CheckRegionSupplementary(overlaps, aln, queryTrees, targetTrees,
                                     allowedOverlapFractionQuery, allowedOverlapFractionTarget)) {
            // This is a supplementary alignment.
            alnPriority.isSupplementary = true;
            alnPriority.priority = 0;

            // Add the new supplementary mapping to the query tree.
            queryIntervals.emplace_back(IntervalInt32(aln->AstartFwd(), aln->AendFwd() - 1, i));
            auto queryIntervalsCopy = queryIntervals;
            queryTrees = IntervalTreeInt32(std::move(queryIntervalsCopy));

            // Add the new supplementary mapping to the target trees.
            targetIntervals[aln->Bid].emplace_back(
                IntervalInt32(aln->BstartFwd(), aln->BendFwd() - 1, i));
            auto targetIntervalsCopy = targetIntervals[aln->Bid];
            targetTrees[aln->Bid] = IntervalTreeInt32(std::move(targetIntervalsCopy));

        } else if (score >= minSecondaryScore) {
            // This is a secondary mapping.
            alnPriority.isSupplementary = false;
            alnPriority.priority = 1;

        } else {
            // Score is too low to be marked as secondary. Flag it with a special
            // priority value.
            alnPriority.isSupplementary = false;
            alnPriority.priority = 2;
        }
    }

    return priorities;
}

int32_t CalcIntervalOverlap(int32_t s1, int32_t e1, int32_t s2, int32_t e2)
{
    return std::max(0, std::min(e1, e2) - std::max(s1, s2));
}

void CreateRegionIntervalTrees(const std::vector<OverlapPtr>& overlaps,
                               const std::vector<OverlapPriority>& priorities,
                               std::function<bool(int32_t a)> CompPriority,
                               IntervalVectorInt32& queryIntervals, IntervalTreeInt32& queryTrees,
                               std::unordered_map<int32_t, IntervalVectorInt32>& targetIntervals,
                               std::unordered_map<int32_t, IntervalTreeInt32>& targetTrees)
{
    if (overlaps.size() != priorities.size()) {
        std::ostringstream oss;
        oss << "Size of the overlaps and of the priorities is not the same! overlaps.size() = "
            << overlaps.size() << ", priorities.size() = " << priorities.size();
        throw std::runtime_error(oss.str());
    }

    queryIntervals.clear();
    targetIntervals.clear();
    targetTrees.clear();

    // IntervalVectorInt64 query_intervals_secondary;
    // std::unordered_map<int64_t, IntervalVectorInt64> target_intervals_secondary;
    for (int32_t i = 0; i < static_cast<int32_t>(overlaps.size()); ++i) {
        const auto& aln = overlaps[i];
        // Skip regions that are not of interest.
        if (CompPriority(priorities[i].priority) == false) {
            continue;
        }
        // Accumulate the intervals.
        queryIntervals.emplace_back(IntervalInt32(aln->AstartFwd(), aln->AendFwd() - 1, i));

        targetIntervals[aln->Bid].emplace_back(
            IntervalInt32(aln->BstartFwd(), aln->BendFwd() - 1, i));
    }

    // Construct the query interval trees from the intervals.
    auto queryIntervalsCopy = queryIntervals;
    queryTrees = IntervalTreeInt32(std::move(queryIntervalsCopy));

    // Construct the target interval trees from the intervals.
    for (auto it = targetIntervals.begin(); it != targetIntervals.end(); ++it) {
        auto targetIntervalsCopy = targetIntervals[it->first];
        targetTrees[it->first] = IntervalTreeInt32(std::move(targetIntervalsCopy));
    }
}

bool CheckRegionSupplementary(const std::vector<OverlapPtr>& overlaps, const OverlapPtr& currentOvl,
                              IntervalTreeInt32& queryTrees,
                              std::unordered_map<int32_t, IntervalTreeInt32>& targetTrees,
                              double allowedOverlapFractionQuery,
                              double allowedOverlapFractionTarget)
{
    /////////////////////////////////////////////////////
    /// Check the query coordinate space for overlaps ///
    /////////////////////////////////////////////////////
    // Find if the current alignment has any overlaps in the query coordinate space.
    auto foundQueryIntervals =
        queryTrees.findOverlapping(currentOvl->AstartFwd(), currentOvl->AendFwd() - 1);
    if (allowedOverlapFractionQuery <= 0.0 && foundQueryIntervals.size()) {
        // For speed, if no overlap between mappings is allowed but some are found, return.
        return false;
    }
    for (const auto& interval : foundQueryIntervals) {
        auto& otherAln = overlaps[interval.value];

        // Amount of overlap in the query coordinates, as a fraction.
        int32_t ovlQuery = CalcIntervalOverlap(currentOvl->AstartFwd(), currentOvl->AendFwd(),
                                               otherAln->AstartFwd(), otherAln->AendFwd());
        double ovlQueryFrac = static_cast<double>(ovlQuery) /
                              static_cast<double>(std::min(currentOvl->ASpan(), otherAln->ASpan()));

        // Amount of overlap in the target coordinates, as a fraction.
        int32_t ovlTarget = CalcIntervalOverlap(currentOvl->BstartFwd(), currentOvl->BendFwd(),
                                                otherAln->BstartFwd(), otherAln->BendFwd());
        double ovlTargetFrac =
            static_cast<double>(ovlTarget) /
            static_cast<double>(std::min(currentOvl->BSpan(), otherAln->BSpan()));

        // The two mappings shouldn't overlap by more than the allowed fraction in either the
        // query and target coordinates if one of them is supplementary.
        if (ovlQueryFrac > allowedOverlapFractionQuery ||
            ovlTargetFrac > allowedOverlapFractionTarget) {
            return false;
        }
    }

    //////////////////////////////////////////////////////
    /// Check the target coordinate space for overlaps ///
    //////////////////////////////////////////////////////
    auto itTrees = targetTrees.find(currentOvl->Bid);
    if (itTrees == targetTrees.end()) {
        // No overlaps in query and no overlaps in target becuase the target
        // was not previously added to the tree.
        return true;
    }
    auto foundTargetIntervals =
        itTrees->second.findOverlapping(currentOvl->BstartFwd(), currentOvl->BendFwd() - 1);
    if (allowedOverlapFractionTarget <= 0.0 && foundTargetIntervals.size()) {
        // For speed, if no overlap between mappings is allowed but some are found, return.
        return false;
    }
    for (const auto& interval : foundTargetIntervals) {
        auto& otherAln = overlaps[interval.value];

        // Amount of overlap in the query coordinates, as a fraction.
        int32_t ovlQuery = CalcIntervalOverlap(currentOvl->AstartFwd(), currentOvl->AendFwd(),
                                               otherAln->AstartFwd(), otherAln->AendFwd());
        double ovlQueryFrac = static_cast<double>(ovlQuery) /
                              static_cast<double>(std::min(currentOvl->ASpan(), otherAln->ASpan()));

        // Amount of overlap in the target coordinates, as a fraction.
        int32_t ovlTarget = CalcIntervalOverlap(currentOvl->BstartFwd(), currentOvl->BendFwd(),
                                                otherAln->BstartFwd(), otherAln->BendFwd());
        double ovlTargetFrac =
            static_cast<double>(ovlTarget) /
            static_cast<double>(std::min(currentOvl->BSpan(), otherAln->BSpan()));

        // The two mappings shouldn't overlap by more than the allowed fraction in either the
        // query and target coordinates if one of them is supplementary.
        if (ovlQueryFrac > allowedOverlapFractionQuery ||
            ovlTargetFrac > allowedOverlapFractionTarget) {
            return false;
        }
    }

    return true;
}

}  // namespace Pancake
}  // namespace PacBio
