// Authors: Ivan Sovic

#include <pancake/SeedHit.hpp>

#include <pbcopper/utility/Ssize.h>

#include <sstream>

namespace PacBio {
namespace Pancake {

std::pair<int32_t, int32_t> CalcHitCoverage(const std::vector<SeedHit>& hits,
                                            const int32_t hitsBegin, int32_t hitsEnd)
{
    /*
      Expects the seed hits to be sorted!
    */
    int32_t coveredBasesQuery = 0;
    int32_t coveredBasesTarget = 0;

    if (hits.empty() || hitsBegin >= static_cast<int32_t>(Utility::Ssize(hits)) ||
        (hitsEnd - hitsBegin) <= 0) {
        return {coveredBasesQuery, coveredBasesTarget};
    }

    // Add the left part of the tile to the covered bases.
    coveredBasesQuery = hits.front().querySpan;
    coveredBasesTarget = hits.front().targetSpan;
    hitsEnd = std::min<int32_t>(hitsEnd, Utility::Ssize(hits));

    for (int32_t i = (hitsBegin + 1); i < hitsEnd; i++) {
        if (hits[i].queryPos < hits[i - 1].queryPos || hits[i].targetPos < hits[i - 1].targetPos) {
            std::ostringstream oss;
            oss << "Invalid seed hit ordering, hits are either not sorted properly or not "
                   "monotonically increasing in both query and target coordinates. "
                << "hits[i - 1].queryPos = " << hits[i - 1].queryPos
                << ", hits[i - 1].targetPos = " << hits[i - 1].targetPos
                << ", hits[i].queryPos = " << hits[i].queryPos
                << ", hits[i].targetPos = " << hits[i].targetPos << ", i = " << i;
            oss << "\n";
            for (int32_t j = hitsBegin; j < hitsEnd; j++) {
                oss << "[hit " << j << "] " << hits[j] << "\n";
            }
            throw std::runtime_error(oss.str());
        }
        coveredBasesQuery += std::min(static_cast<int32_t>(hits[i].querySpan),
                                      (hits[i].queryPos - hits[i - 1].queryPos));
        coveredBasesTarget += std::min(static_cast<int32_t>(hits[i].targetSpan),
                                       (hits[i].targetPos - hits[i - 1].targetPos));
    }

    return {coveredBasesQuery, coveredBasesTarget};
}

}  // namespace Pancake
}  // namespace PacBio
