// Authors: Robert Vaser

#include <pancake/BurstExcision.hpp>

#include <pbcopper/logging/Logging.h>
#include <pbcopper/utility/Ssize.h>
#include <pancake/FastaSequenceCachedStore.hpp>
#include <pancake/MapperCLR.hpp>
#include <pancake/SeedHit.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/stats.hpp>

#include <algorithm>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace PacBio {
namespace Pancake {

namespace {

struct Burst
{
    int32_t TargetBegin = 0;
    int32_t TargetEnd = 0;
    int32_t QueryId = 0;
    int32_t QueryBegin = 0;
    int32_t QueryEnd = 0;
    bool IsReversed = false;
    bool IsInTarget = false;
    bool IsCandidate = true;

    Burst(const int32_t targetBegin, const int32_t targetEnd, const int32_t queryId,
          const int32_t queryBegin, const int32_t queryEnd, const int32_t queryLen,
          const bool isReversed)
        : TargetBegin(targetBegin)
        , TargetEnd(targetEnd)
        , QueryId(queryId)
        , QueryBegin(isReversed ? queryLen - queryEnd : queryBegin)
        , QueryEnd(isReversed ? queryLen - queryBegin : queryEnd)
        , IsReversed(isReversed)
        , IsInTarget((targetEnd - targetBegin) > (queryEnd - queryBegin))
    {}

    bool Overlaps(const Burst& other) const
    {
        return (TargetEnd >= other.TargetBegin) && (other.TargetEnd >= TargetBegin);
    }
};

}  // namespace

std::vector<std::vector<std::pair<int32_t, int32_t>>> BurstExcision(
    const FastaSequenceCachedStore& zmw, const MapperCLRSettings& mapperSettings,
    const int32_t minBurstSize, const double maxHitDistanceRatio,
    const double lowerInterPulseDistanceRatio, const double upperInterPulseDistanceRatio,
    const double lowerPulseWidthRatio, const double upperPulseWidthRatio)
{
    std::vector<std::vector<std::pair<int32_t, int32_t>>> dst(zmw.Size());

    if (zmw.Size() < 2) {
        PBLOG_DEBUG << "Empty input subread store";
        return dst;
    }
    for (const auto& it : zmw.records()) {
        if ((it.IPD() == nullptr) || (it.PW() == nullptr)) {
            PBLOG_DEBUG << "Subread " << it.Name() << " has missing kinetics";
            return dst;
        }
    }

    PacBio::Pancake::MapperCLRSettings settings{mapperSettings};
    settings.map.refineSeedHits = false;  // critical!
    settings.align.align = false;         // no need for alignment

    PacBio::Pancake::MapperCLR mapper{settings};

    FastaSequenceCachedStore backbone;
    backbone.AddRecord(zmw.records()[0]);

    const auto mappingResult = mapper.MapAndAlign(backbone, zmw);

    // find candidate regions
    std::vector<Burst> bursts;
    for (int32_t i = 1; i < Utility::Ssize(mappingResult); ++i) {
        const auto& subreadResult = mappingResult[i];
        for (const auto& mapping : subreadResult.mappings) {
            if ((mapping->priority > 1) || mapping->chain.hits.empty()) {
                continue;
            }

            const auto isBurstCandidate = [minBurstSize, maxHitDistanceRatio](const SeedHit& lhs,
                                                                              const SeedHit& rhs) {
                const int32_t queryGap = rhs.queryPos - lhs.queryPos;
                const int32_t targetGap = rhs.targetPos - lhs.targetPos;
                return (std::abs(queryGap - targetGap) > minBurstSize) &&
                       (std::min(queryGap, targetGap) <
                        std::max(queryGap, targetGap) * maxHitDistanceRatio);
            };

            auto it = std::cbegin(mapping->chain.hits);
            const auto last = std::prev(std::cend(mapping->chain.hits));
            for (; it != last; ++it) {
                if (isBurstCandidate(*it, *(it + 1))) {
                    bursts.emplace_back(it->targetPos, (it + 1)->targetPos, i, it->queryPos,
                                        (it + 1)->queryPos, zmw.records()[i].Size(), it->targetRev);
                }
            }
        }
    }

    std::sort(std::begin(bursts), std::end(bursts),
              [](const Burst& lhs, const Burst& rhs) { return lhs.TargetBegin < rhs.TargetBegin; });

    // check occurrence of all burst regions, given their backbone coordinates
    // - each backbone burst will update its boundaries and nullify all subsequent
    //   overlapping bursts
    // - each subread burst which occurs too many times will nullify itself and all
    //   other overlapping bursts
    std::vector<std::vector<Burst>> burstsPerSubread(zmw.Size());
    {
        const int32_t majorityThreshold = (zmw.Size() / 2ULL) - 1;

        int32_t numCandidates = 0;
        for (int32_t i = 0; i < Utility::Ssize(bursts); ++i) {
            if (!bursts[i].IsCandidate) {
                continue;
            }

            if (bursts[i].IsInTarget) {
                int32_t numOcc = 0;
                for (int32_t j = (i + 1); j < Utility::Ssize(bursts); ++j) {
                    if (bursts[j].IsInTarget && bursts[i].Overlaps(bursts[j])) {
                        // only accept the first backbone burst and update its coordinates
                        bursts[i].TargetBegin =
                            std::min(bursts[i].TargetBegin, bursts[j].TargetBegin);
                        bursts[i].TargetEnd = std::max(bursts[i].TargetEnd, bursts[j].TargetEnd);
                        bursts[j].IsCandidate = false;  // nullify overlapping backbone burst
                        ++numOcc;
                    }
                }
                if (numOcc > majorityThreshold) {
                    std::swap(bursts[i].QueryBegin, bursts[i].TargetBegin);
                    std::swap(bursts[i].QueryEnd, bursts[i].TargetEnd);
                    burstsPerSubread[0].emplace_back(bursts[i]);
                    ++numCandidates;
                }

            } else {
                std::vector<int32_t> inserts = {bursts[i].QueryEnd - bursts[i].QueryBegin};
                for (int32_t j = (i + 1); j < Utility::Ssize(bursts); ++j) {
                    if ((bursts[i].IsReversed == bursts[j].IsReversed) &&
                        bursts[i].Overlaps(bursts[j])) {
                        inserts.emplace_back(bursts[j].QueryEnd - bursts[j].QueryBegin);
                    }
                }
                bool evilContext = false;
                if (Utility::Ssize(inserts) > majorityThreshold) {
                    double mean = std::accumulate(std::cbegin(inserts), std::cend(inserts), 0.0) /
                                  Utility::Ssize(inserts);
                    for (const auto& it : inserts) {
                        const double deviationFromMean = std::abs(it - mean) / mean;
                        if (deviationFromMean > 0.5) {
                            evilContext = true;
                            break;
                        }
                    }
                }

                if (!evilContext && (Utility::Ssize(inserts) > majorityThreshold)) {
                    // nullify all overlapping subread bursts
                    for (int32_t j = (i + 1); j < Utility::Ssize(bursts); ++j) {
                        if ((bursts[i].IsReversed == bursts[j].IsReversed) &&
                            bursts[i].Overlaps(bursts[j])) {
                            bursts[j].IsCandidate = false;
                        }
                    }
                } else {
                    burstsPerSubread[bursts[i].QueryId].emplace_back(bursts[i]);
                    ++numCandidates;
                }
            }
        }

        if (numCandidates == 0) {
            return dst;
        }
    }

    // find zmw kinetics median
    constexpr uint32_t WINDOW_SIZE = 32;
    constexpr int32_t WINDOW_HALF_SIZE = 16;

    int32_t ipdMedian = 0;
    int32_t pwMedian = 0;
    {
        using KineticsAccumulator = boost::accumulators::accumulator_set<
            int64_t, boost::accumulators::stats<boost::accumulators::tag::median>>;

        KineticsAccumulator ipdAccumulator;
        KineticsAccumulator pwAccumulator;

        const auto accumulateKinetics = [WINDOW_HALF_SIZE](const std::vector<uint16_t>& kinetics,
                                                           KineticsAccumulator& accumulator) {
            int32_t window =
                std::accumulate(std::cbegin(kinetics),
                                std::cbegin(kinetics) +
                                    std::min<int32_t>(Utility::Ssize(kinetics), WINDOW_HALF_SIZE),
                                0);
            for (int32_t i = 0; i < Utility::Ssize(kinetics); ++i) {
                if (i >= WINDOW_HALF_SIZE) {
                    window -= kinetics[i - WINDOW_HALF_SIZE];
                }
                if (i < (Utility::Ssize(kinetics) - WINDOW_HALF_SIZE)) {
                    window += kinetics[i + WINDOW_HALF_SIZE];
                }
                accumulator(window / WINDOW_SIZE);
            }
        };

        for (const auto& subread : zmw.records()) {
            accumulateKinetics(subread.IPD()->Data(), ipdAccumulator);
            accumulateKinetics(subread.PW()->Data(), pwAccumulator);
        }

        ipdMedian = std::round(boost::accumulators::median(ipdAccumulator));
        pwMedian = std::round(boost::accumulators::median(pwAccumulator));
    }
    const uint32_t ipdLowerBound =
        std::max<uint32_t>(0, ipdMedian * (1 - lowerInterPulseDistanceRatio));
    const uint32_t ipdUpperBound = ipdMedian * (1 + upperInterPulseDistanceRatio);

    const uint32_t pwLowerBound = std::max<uint32_t>(0, pwMedian * (1 - lowerPulseWidthRatio));
    const uint32_t pwUpperBound = pwMedian * (1 + upperPulseWidthRatio);

    // find burst begin/end points and remove them
    {
        for (int32_t i = 0; i < Utility::Ssize(burstsPerSubread); ++i) {
            if (burstsPerSubread[i].empty()) {
                continue;
            }

            const auto& subread = zmw.records()[i];
            const int32_t subreadSize = subread.Size();

            const auto ipd = subread.IPD()->Data();
            const auto pw = subread.PW()->Data();

            int32_t numBursts = 0;
            for (auto& it : burstsPerSubread[i]) {
                int32_t ipdWindow = std::accumulate(
                    std::cbegin(ipd) + it.QueryBegin - std::min(it.QueryBegin, WINDOW_HALF_SIZE),
                    std::cbegin(ipd) + it.QueryBegin +
                        std::min(subreadSize - it.QueryBegin, WINDOW_HALF_SIZE),
                    0);
                int32_t pwWindow = std::accumulate(
                    std::cbegin(pw) + it.QueryBegin - std::min(it.QueryBegin, WINDOW_HALF_SIZE),
                    std::cbegin(pw) + it.QueryBegin +
                        std::min(subreadSize - it.QueryBegin, WINDOW_HALF_SIZE),
                    0);

                int32_t burstBegin = 0;
                int32_t burstEnd = 0;

                for (int32_t j = it.QueryBegin, b = 0; j < it.QueryEnd; ++j) {
                    if (j >= WINDOW_HALF_SIZE) {
                        ipdWindow -= ipd[j - WINDOW_HALF_SIZE];
                        pwWindow -= pw[j - WINDOW_HALF_SIZE];
                    }
                    if (j < (subreadSize - WINDOW_HALF_SIZE)) {
                        ipdWindow += ipd[j + WINDOW_HALF_SIZE];
                        pwWindow += pw[j + WINDOW_HALF_SIZE];
                    }
                    if (((ipdWindow / WINDOW_SIZE) < ipdLowerBound) ||
                        ((ipdWindow / WINDOW_SIZE) > ipdUpperBound) ||
                        ((pwWindow / WINDOW_SIZE) < pwLowerBound) ||
                        ((pwWindow / WINDOW_SIZE) > pwUpperBound)) {
                        if (b == 0) {
                            b = 1;
                            burstBegin = j;
                        }
                        burstEnd = j;
                    }
                }

                if (burstBegin && burstEnd) {
                    it.QueryBegin = burstBegin;
                    it.QueryEnd = burstEnd;
                    ++numBursts;
                } else {
                    it.IsCandidate = false;
                }
            }

            if (numBursts > 0) {
                std::sort(std::begin(burstsPerSubread[i]), std::end(burstsPerSubread[i]),
                          [](const Burst& lhs, const Burst& rhs) {
                              return lhs.QueryBegin < rhs.QueryBegin;
                          });

                dst[i].reserve(numBursts);
                for (const auto& it : burstsPerSubread[i]) {
                    if (it.IsCandidate) {
                        dst[i].emplace_back(it.QueryBegin, it.QueryEnd);
                    }
                }
            }
        }
    }

    return dst;
}

}  // namespace Pancake
}  // namespace PacBio
