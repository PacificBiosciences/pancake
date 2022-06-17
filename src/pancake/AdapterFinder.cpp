// Authors: Ivan Sovic

#include <pancake/AdapterFinder.hpp>

#include <pancake/AlignerEdlib.hpp>
#include <pancake/MapperCLR.hpp>

#include <pbcopper/utility/Ssize.h>

#include <algorithm>
#include <deque>

// #define PANCAKE_ADAPTER_FINDER_DEBUG

namespace PacBio {
namespace Pancake {

namespace Adapter {

std::vector<AdapterPosition> FindPotentialAdapterLocations(
    const FastaSequenceCached& seq, const std::string_view adapterSeq, const int32_t k,
    const int32_t w, const int32_t minNumSeeds, const double minIdentity,
    const double minDiagonalSpanFraction, const int32_t maxAllowedFlank,
    const int32_t maxAdapterEditDist)
{
    const auto InitPancakeSettings = [](const int32_t userK, const int32_t userW) {
        PacBio::Pancake::MapperCLRSettings settings;

        // Indexing.
        settings.map.seedParams.KmerSize = userK;
        settings.map.seedParams.MinimizerWindow = userW;
        settings.map.seedParams.Spacing = 0;
        settings.map.seedParams.UseHPCForSeedsOnly = true;
        settings.map.seedParams.UseRC = false;  // <- IMPORTANT. Do not map to revcmp.
        settings.map.seedParamsFallback = settings.map.seedParams;

        // Seed hit occurrence filtering.
        settings.map.freqPercentile = 0.000;
        settings.map.seedOccurrenceMaxMemory = 100'000'000;
        settings.map.seedOccurrenceMax = 100;
        settings.map.seedOccurrenceMin = 5;

        // Mapping.
        settings.map.chainMaxSkip = 25;
        settings.map.chainMaxPredecessors = 500;
        settings.map.chainBandwidth = 500;
        settings.map.seedJoinDist = 10000;
        settings.map.longMergeBandwidth = 10000;
        settings.map.minNumSeeds = 3;
        settings.map.minCoveredBases = 0;
        settings.map.minDPScore = 50;
        settings.map.secondaryAllowedOverlapFractionQuery = 0.0;
        settings.map.secondaryAllowedOverlapFractionTarget = 0.0;
        settings.map.secondaryMinScoreFraction = 0.00;
        settings.map.useLIS = true;

        // Other.
        settings.map.skipSymmetricOverlaps = false;
        settings.map.selfHitPolicy = MapperSelfHitPolicy::DEFAULT;
        settings.map.minQueryLen = 50;
        settings.map.bestNSecondary = -1;

        // Alignment region extraction parameters.
        settings.map.maxFlankExtensionDist = 10000;
        settings.map.flankExtensionFactor = 1.3;
        settings.map.minAlignmentSpan = 200;

        // Refining seed hits.
        settings.map.refineSeedHits = false;
        settings.map.refineMinGap1 = 50;
        settings.map.refineDiffThreshold = 40;
        settings.map.refineMinGap2 = 30;

        // Reseeding long alignment regions with smaller seeds.
        settings.map.reseedGaps = false;
        settings.map.reseedGapMinLength = 500;
        settings.map.reseedGapMaxLength = -1;
        settings.map.reseedSeedParams =
            PacBio::Pancake::SeedDBParameters{10, 5, 0, false, true, true};
        settings.map.reseedFreqPercentile = 0.0002;
        settings.map.reseedOccurrenceMin = 5;
        settings.map.reseedOccurrenceMax = 100;
        settings.map.reseedOccurrenceMaxMemory = 100'000'000;

        //Alignment parameters.
        settings.align.align = true;
        settings.align.selfHitPolicy = MapperSelfHitPolicy::DEFAULT;
        settings.align.alignerTypeGlobal = PacBio::Pancake::AlignerType::EDLIB;
        settings.align.alnParamsGlobal.zdrop = 400;
        settings.align.alnParamsGlobal.zdrop2 = 200;
        settings.align.alnParamsGlobal.alignBandwidth = 500;
        settings.align.alnParamsGlobal.endBonus = 1000;
        settings.align.alnParamsGlobal.matchScore = 2;
        settings.align.alnParamsGlobal.mismatchPenalty = 4;
        settings.align.alnParamsGlobal.gapOpen1 = 4;
        settings.align.alnParamsGlobal.gapExtend1 = 2;
        settings.align.alnParamsGlobal.gapOpen2 = 24;
        settings.align.alnParamsGlobal.gapExtend2 = 1;
        settings.align.alignerTypeExt = PacBio::Pancake::AlignerType::KSW2;
        settings.align.alnParamsExt = settings.align.alnParamsGlobal;

        return settings;
    };

    auto PointPositionVsAntidiag = [](const int32_t seqLen, const int32_t x,
                                      const int32_t y) -> int32_t {
        /**
         * This function returns one of [-1, 0, 1] if the point is below, on or above the antidiagonal.
         * NOTE: A point is above/below the antidiagonal if:
         *          c = L - x - y
         *       where "L" is the sequence length (both query and target are assumed to have the same length), and "x" and "y" are the point.
         *       The "c" is < 0 if the point is above the line, > 0 if it's below the line and == 0 if the point is on the line.
         *       This equation can be obtained from the cross-product of the line vs the point.
         */
        const int32_t val = (seqLen - x - y);
        return std::clamp(val, -1, 1);
    };

    struct Line
    {
        int32_t x1 = 0;
        int32_t y1 = 0;
        int32_t x2 = 0;
        int32_t y2 = 0;
    };
    auto LineLineIntersection = [](const Line& line1, const Line& line2) {
        /**
         * @brief Returns the intersection point between two lines.
         * \returns Tuple consisting of: <result_valid, intersection_x, intersection_y>.
         *              If two lines are parallel, the first value in the tuple is "false".
         *
         * Note: This could be implemented with a cross-product alternatively too.
         */
        if ((line1.x2 != line1.x1) && (line2.x2 != line2.x1)) {
            const double k1 = (line1.y2 - line1.y1) / static_cast<double>(line1.x2 - line1.x1);
            const double l1 = -k1 * line1.x1 + line1.y1;

            const double k2 = (line2.y2 - line2.y1) / static_cast<double>(line2.x2 - line2.x1);
            const double l2 = -k2 * line2.x1 + line2.y1;

            if (k1 == k2) {
                return std::tuple<bool, int32_t, int32_t>(false,
                                                          std::numeric_limits<int32_t>::max(),
                                                          std::numeric_limits<int32_t>::max());
            }

            const double x = (l2 - l1) / (k1 - k2);
            const double y = k1 * x + l1;

            return std::tuple<bool, int32_t, int32_t>(true, x, y);

        } else if ((line1.x2 == line1.x1) && (line2.x2 != line2.x1)) {
            const double k2 = (line2.y2 - line2.y1) / static_cast<double>(line2.x2 - line2.x1);
            const double l2 = -k2 * line2.x1 + line2.y1;

            const double x = line1.x1;
            const double y = k2 * x + l2;

            return std::tuple<bool, int32_t, int32_t>(true, x, y);

        } else if ((line1.x2 != line1.x1) && (line2.x2 == line2.x1)) {
            const double k1 = (line1.y2 - line1.y1) / static_cast<double>(line1.x2 - line1.x1);
            const double l1 = -k1 * line1.x1 + line1.y1;

            const double x = line2.x1;
            const double y = k1 * x + l1;

            return std::tuple<bool, int32_t, int32_t>(true, x, y);
        }

        return std::tuple<bool, int32_t, int32_t>(false, std::numeric_limits<int32_t>::max(),
                                                  std::numeric_limits<int32_t>::max());
    };

    /// Reverse the query sequence to only align against that.
    const std::string seqRev = PacBio::Pancake::ReverseComplement(
        {seq.c_str(), static_cast<size_t>(seq.size())}, 0, seq.size());
    const FastaSequenceCached seqRevCached("rev", seqRev.c_str(), seqRev.size(), seq.Id());

    /// Construct adapter sequence views.
    const std::string adapterSeqRev =
        PacBio::Pancake::ReverseComplement(adapterSeq, 0, adapterSeq.size());
    const std::string_view adapterSeqView(adapterSeq);
    const std::string_view adapterSeqRevView(adapterSeqRev);

    /// Map the revcmp sequence (query) to its forward self (target). Target coordinates will be in the space of the original seq.
    const PacBio::Pancake::MapperCLRSettings settings = InitPancakeSettings(k, w);
    PacBio::Pancake::MapperCLR mapper(settings);
    const std::vector<PacBio::Pancake::MapperBaseResult> mapResults =
        mapper.MapAndAlign({seq}, {seqRevCached});

    /// Sanity check.
    if (mapResults.size() != 1) {
        PBLOG_DEBUG << "Wrong number of queries aligned! Expected exactly one query, but got "
                    << mapResults.size() << ". Skipping.\n";
        return {};
    }

#ifdef PANCAKE_ADAPTER_FINDER_DEBUG
    std::cerr.flush();
#endif

    /// Define the Edlib aligner.
    const AlignmentParameters alnParams;
    PacBio::Pancake::AlignerEdlib aligner(alnParams);
    auto AlignAdapter = [&aligner](const std::string_view tseq, const std::string_view adapterFwd,
                                   const std::string_view adapterRev, const int32_t maxEditDist) {
        const auto alnFwd = aligner.Infix(adapterFwd, tseq);
        const auto alnRev = aligner.Infix(adapterRev, tseq);
        const int32_t edFwd = alnFwd.diffs.EditDistance(false, false);
        const int32_t edRev = alnRev.diffs.EditDistance(false, false);
        if (edFwd < edRev) {
            const bool isValid = (edFwd <= maxEditDist);
            return std::tuple(isValid, alnFwd, false);
        }
        const bool isValid = (edRev <= maxEditDist);
        return std::tuple(isValid, alnRev, true);
    };

    /// Process each self-mapping to find potential adapter locations.
    std::vector<AdapterPosition> ret;
    const auto& mapResult = mapResults[0];
    for (size_t alnId = 0; alnId < mapResult.mappings.size(); ++alnId) {
        if (mapResult.mappings[alnId] == nullptr) {
            continue;
        }

        /// Chained region.
        const auto& cr = *(mapResult.mappings[alnId]);
        const auto& aln = cr.mapping;

        /// Compute the rough maximum span of the diagonal where the current mapping lies, allowing
        /// differences between the start and end diagonals. This consists of 3 parts:
        ///     <unaligned_left_flank> + aligned_middle + <unaligned_right_flank>
        const int32_t mappedSpan = std::max(aln->ASpan(), aln->BSpan());
        const int32_t diagMaxSpan = std::min(aln->Astart, aln->Bstart) + mappedSpan +
                                    std::min((aln->Alen - aln->Aend), (aln->Blen - aln->Bend));
        const double diagCoverage =
            diagMaxSpan == 0 ? 0.0 : (static_cast<double>(mappedSpan) / diagMaxSpan);
        const int32_t minFlankLength = (std::min(aln->Astart, aln->Bstart),
                                        std::min((aln->Alen - aln->Aend), (aln->Blen - aln->Bend)));

#ifdef PANCAKE_ADAPTER_FINDER_DEBUG
        // Debug verbose.
        const int32_t diagStart = aln->Bstart - aln->Astart;
        const int32_t diagEnd = aln->Bend - aln->Aend;
        std::cerr << "[alnId = " << alnId << "] "
                  << PacBio::Pancake::OverlapWriterBase::PrintOverlapAsM4(*aln, "", "", true, false)
                  << ", NumSeeds = " << cr.mapping->NumSeeds << ", diagStart = " << diagStart
                  << ", diagEnd = " << diagEnd << ", mappedSpan = " << mappedSpan
                  << ", diagMaxSpan = " << diagMaxSpan << ", diagCoverage = " << diagCoverage
                  << "\n";
#endif

        /// Only fwd alignments allowed, becauase the query is reverse complemented.
        assert((cr.mapping->Arev == false) && (cr.mapping->Brev == false));
        if (cr.mapping->Arev || cr.mapping->Brev) {
            PBLOG_DEBUG << "Strange: cr.mapping->Arev || cr.mapping->Brev. This shouldn't be "
                           "happening due to mapping parameters.\n";
            continue;
        }

        /// Skip poor mappings. If alignment was applied then filter on identity, otherwise filter
        /// on number of chained seed hits.
        const bool isAligned = cr.mapping->Identity > 0.0f;
        if ((isAligned && (static_cast<double>(cr.mapping->Identity) < minIdentity)) ||
            ((isAligned == false) && (cr.mapping->NumSeeds < minNumSeeds))) {
#ifdef PANCAKE_ADAPTER_FINDER_DEBUG
            std::cerr << "-> Skipped. isAligned = " << isAligned
                      << ", identity = " << cr.mapping->Identity
                      << ", numSeeds = " << cr.mapping->NumSeeds
                      << ", minIdentity = " << minIdentity << ", minNumSeeds = " << minNumSeeds
                      << "\n";
#endif
            continue;
        }

        if ((minDiagonalSpanFraction > 0) && (diagCoverage < minDiagonalSpanFraction)) {
#ifdef PANCAKE_ADAPTER_FINDER_DEBUG
            std::cerr << "-> Skipped. diagCoverage = " << diagCoverage
                      << ", minDiagonalSpanFraction = " << minDiagonalSpanFraction << "\n";
#endif
            continue;
        }

        if ((maxAllowedFlank > 0) && (minFlankLength > maxAllowedFlank)) {
#ifdef PANCAKE_ADAPTER_FINDER_DEBUG
            std::cerr << "-> Skipped. minFlankLength = " << minFlankLength
                      << ", maxAllowedFlank = " << maxAllowedFlank << "\n";
#endif
            continue;
        }

        /// Find candidate position for an adapter.
        int32_t xAdapter = -1;  // Not valid because coords are [0, seqLen].
        int32_t yAdapter = -1;
        int32_t searchDist = 300;  // Distance from rough adapter position for alignment search.
        bool isPalindromic = false;
        AdapterPosition adapter;

        /// Stage 1: Rough estimate of the intersection point between the diagonal/antidiagonal.
        {
            /// Compute the intersection point of this mapping vs the antidiagonal.
            const int32_t startRelativePos =
                PointPositionVsAntidiag(aln->Blen, aln->Astart, aln->Bstart);
            const int32_t endRelativePos = PointPositionVsAntidiag(aln->Blen, aln->Aend, aln->Bend);
            bool rv = false;
            std::tie(rv, xAdapter, yAdapter) = LineLineIntersection(
                {0, aln->Blen, aln->Alen, 0}, {aln->Astart, aln->Bstart, aln->Aend, aln->Bend});
            isPalindromic = rv & (startRelativePos != endRelativePos);

            /// Store the results.
            adapter.Start = yAdapter;
            adapter.End = yAdapter;
            adapter.type = AdapterDetectionType::PALINDROMIC_ROUGH;

#ifdef PANCAKE_ADAPTER_FINDER_DEBUG
            std::cerr << "[Stage 1: rough] Found intersection:\n";
            std::cerr << "  - xAdapter = " << xAdapter << ", yAdapter = " << yAdapter
                      << ", seq.size() = " << seq.size() << ", isPalindromic = " << isPalindromic
                      << "\n"
                      << "  - startRelativePos = " << startRelativePos
                      << ", endRelativePos = " << endRelativePos << "\n";
            std::cerr << "  - Adapter: " << adapter << "\n";
#endif
        }

        /// If there is no intersection with the antidiagonal, skip this alignment.
        if ((isPalindromic == false) || (xAdapter < 0) || (yAdapter < 0)) {
            continue;
        }

        /// Stage 2: Find a more accurate intersection point based on seed hits, if seed hits are available.
        /// This is important if there is a missing adapter so that the breakpoint is as accurate as possible.
        /// Looks for the last seed hit before the antidiagonal and first seed hit after the antidiagonal.
        /// No need to pay special attention to reverse hits because only fwd are allowed.
        {
            for (size_t i = 1; i < cr.chain.hits.size(); ++i) {
                const auto& hit = cr.chain.hits[i];
                const auto& hitPrev = cr.chain.hits[i - 1];

                /// Check that the two points are on the opposite sides of the antidiagonal.
                const int32_t startRelativePos =
                    PointPositionVsAntidiag(aln->Blen, hitPrev.queryPos, hitPrev.targetPos);
                const int32_t endRelativePos =
                    PointPositionVsAntidiag(aln->Blen, hit.queryPos, hit.targetPos);

                if (startRelativePos == endRelativePos) {
                    continue;
                }

                const auto [rv, xAdapterNew, yAdapterNew] = LineLineIntersection(
                    {0, aln->Blen, aln->Alen, 0},
                    {hitPrev.queryPos, hitPrev.targetPos, hit.queryPos, hit.targetPos});

                /// Update the values.
                xAdapter = xAdapterNew;
                yAdapter = yAdapterNew;
                adapter.Start = yAdapterNew;
                adapter.End = yAdapterNew;
                adapter.type = AdapterDetectionType::PALINDROMIC_FINE;

#ifdef PANCAKE_ADAPTER_FINDER_DEBUG
                std::cerr << "[Stage 2: detailed] Found intersection:\n";
                std::cerr << "  - hitPrev: " << hitPrev << "\n";
                std::cerr << "  - hit:     " << hit << "\n";
                std::cerr << "  - xAdapterNew = " << xAdapterNew
                          << ", yAdapterNew = " << yAdapterNew << ", seq.size() = " << seq.size()
                          << "\n";
                std::cerr << "  - Adapter: " << adapter << "\n";
                std::cerr << "\n";
#endif
                break;
            }
        }

        /// Stage 3: Align the adapter in a region around the computed y-coordinate.
        {
            const int32_t startPos = std::max(0, yAdapter - searchDist);
            const int32_t span = std::min(aln->Blen, yAdapter + searchDist) - startPos;
            const std::string_view tseq(seq.c_str() + startPos, span);
            const auto [alnValid, alnAdapterResult, isAlnRev] =
                AlignAdapter(tseq, adapterSeqView, adapterSeqRevView, maxAdapterEditDist);

            const int32_t alnAdapterEnd = startPos + alnAdapterResult.lastTargetPos + 1;
            const int32_t alnAdapterStart =
                alnAdapterEnd - (alnAdapterResult.diffs.numEq + alnAdapterResult.diffs.numX +
                                 alnAdapterResult.diffs.numD);

            /// Update results only if alignment of the adapter was succesful.
            if (alnValid) {
                adapter.End = alnAdapterEnd;
                adapter.Start = alnAdapterStart;
                adapter.type = AdapterDetectionType::ALIGNED;
            }

#ifdef PANCAKE_ADAPTER_FINDER_DEBUG
            std::cerr << "[Stage 3: edlib] Alignment results:\n";
            std::cerr << "  - startPos = " << startPos << ", span = " << span << "\n";
            std::cerr << "  - adapterSeq.size() = " << adapterSeq.size()
                      << ", tseq.size() = " << tseq.size() << "\n";
            std::cerr << "  - " << alnAdapterResult << "; alnValid = " << alnValid
                      << ", isAlnRev = " << isAlnRev << "\n";
            std::cerr << "  - alnAdapterStart = " << alnAdapterStart
                      << ", alnAdapterEnd = " << alnAdapterEnd << "\n";
            std::cerr << "  - Adapter: " << adapter << "\n";
#endif
        }

        /// Store the potential adapter location.
        ret.emplace_back(adapter);

#ifdef PANCAKE_ADAPTER_FINDER_DEBUG
        std::cerr << "[final] Adapter position:\n";
        std::cerr << "  - Adapter: " << adapter << "\n";
        std::cerr << "\n";
#endif
    }

    /// Always return sorted adapter locations.
    std::sort(std::begin(ret), std::end(ret),
              [](const auto& a, const auto& b) { return a.Start < b.Start; });

    return ret;
}

}  // namespace Adapter
}  // namespace Pancake
}  // namespace PacBio
