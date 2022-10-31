// Authors: Ivan Sovic

#include <pancake/AlignmentSeeded.hpp>

#include <pancake/AlignmentTools.hpp>
#include <pancake/OverlapWriterBase.hpp>

#include <pbcopper/logging/Logging.h>

#include <iostream>
#include <string_view>

namespace PacBio {
namespace Pancake {

// #define DEBUG_ALIGNMENT_SEEDED

AlignmentResult AlignSingleRegion(const std::string_view targetSeq,
                                  const std::string_view querySeqFwd,
                                  const std::string_view querySeqRev, AlignerBasePtr& alignerGlobal,
                                  AlignerBasePtr& alignerExt, const AlignmentRegion& region)
{
    const char* querySeqInStrand = region.queryRev ? querySeqRev.data() : querySeqFwd.data();
    const char* targetSeqInStrand = targetSeq.data();
    int32_t qStart = region.qStart;
    int32_t tStart = region.tStart;
    const int32_t qSpan = region.qSpan;
    const int32_t tSpan = region.tSpan;
    const int32_t queryLen = querySeqFwd.size();
    const int32_t targetLen = targetSeq.size();

    if (qSpan == 0 && tSpan == 0) {
        return {};
    }

    if (qStart >= queryLen || (qStart + qSpan) > queryLen || tStart >= targetLen ||
        (tStart + tSpan) > targetLen || qSpan < 0 || tSpan < 0) {
        std::ostringstream oss;
        oss << "(AlignmentRegion) Coordinates out of bounds in AlignRegionsGeneric!"
            << " qStart = " << qStart << ", qSpan = " << qSpan << ", queryLen = " << queryLen
            << ", tStart = " << tStart << ", tSpan = " << tSpan << ", targetLen = " << targetLen
            << ", regionId = " << region.regionId;
        throw std::runtime_error(oss.str());
    }
    if (targetSeq.data() == NULL || querySeqFwd.data() == NULL || querySeqRev.data() == NULL) {
        std::ostringstream oss;
        oss << "(AlignmentRegion) NULL sequence passed to AlignmentRegion!"
            << " qStart = " << qStart << ", qSpan = " << qSpan << ", queryLen = " << queryLen
            << ", tStart = " << tStart << ", tSpan = " << tSpan << ", targetLen = " << targetLen
            << ", regionId = " << region.regionId;
        throw std::runtime_error(oss.str());
    }
    if (querySeqFwd.size() != querySeqRev.size()) {
        std::ostringstream oss;
        oss << "(AlignmentRegion) Forward and reverse query sequences are not of the same length!";
        throw std::runtime_error(oss.str());
    }

    // Prepare the reversed front sequence if required.
    std::string qSubSeq;
    std::string tSubSeq;
    if (region.type == RegionType::FRONT) {
        qSubSeq = std::string(querySeqInStrand + region.qStart, region.qSpan);
        tSubSeq = std::string(targetSeqInStrand + region.tStart, region.tSpan);
        std::reverse(qSubSeq.begin(), qSubSeq.end());
        std::reverse(tSubSeq.begin(), tSubSeq.end());
        querySeqInStrand = qSubSeq.c_str();
        targetSeqInStrand = tSubSeq.c_str();
        qStart = 0;
        tStart = 0;
    }

    // Align.
    AlignmentResult alnRes;
    if (region.type == RegionType::FRONT || region.type == RegionType::BACK) {
        alnRes = alignerExt->Extend(std::string_view(querySeqInStrand + qStart, qSpan),
                                    std::string_view(targetSeqInStrand + tStart, tSpan));
    } else {
        alnRes = alignerGlobal->Global(std::string_view(querySeqInStrand + qStart, qSpan),
                                       std::string_view(targetSeqInStrand + tStart, tSpan));
    }

    if (region.type == RegionType::FRONT) {
        std::reverse(alnRes.cigar.begin(), alnRes.cigar.end());
    }

    return alnRes;
}

AlignRegionsGenericResult AlignRegionsGeneric(const std::string_view targetSeq,
                                              const std::string_view querySeqFwd,
                                              const std::string_view querySeqRev,
                                              const std::vector<AlignmentRegion>& regions,
                                              AlignerBasePtr& alignerGlobal,
                                              AlignerBasePtr& alignerExt)
{
    AlignRegionsGenericResult ret;

    std::vector<AlignmentResult> alignedRegions;

    for (size_t i = 0; i < regions.size(); ++i) {
        const auto& region = regions[i];

        auto alnRes = AlignSingleRegion(targetSeq, querySeqFwd, querySeqRev, alignerGlobal,
                                        alignerExt, region);

        if (region.type == RegionType::FRONT) {
            ret.offsetFrontQuery = alnRes.lastQueryPos;
            ret.offsetFrontTarget = alnRes.lastTargetPos;
        } else if (region.type == RegionType::BACK) {
            ret.offsetBackQuery = alnRes.lastQueryPos;
            ret.offsetBackTarget = alnRes.lastTargetPos;
        }

        // Store the results.
        alignedRegions.emplace_back(std::move(alnRes));

#ifdef DEBUG_ALIGNMENT_SEEDED
        std::cerr << "[aln region i = " << i << " / " << regions.size() << "] " << region
                  << ", CIGAR: " << alignedRegions.back().cigar.ToStdString() << "\n"
                  << alnRes << "\n\n";
#endif
    }

    // Merge the CIGAR chunks.
    int32_t score = 0;
    for (const auto& alnRegion : alignedRegions) {
        const auto& currCigar = alnRegion.cigar;
        if (currCigar.empty()) {
            continue;
        }
        if (ret.cigar.empty() || ret.cigar.back().Type() != currCigar.front().Type()) {
            ret.cigar.emplace_back(currCigar.front());
        } else {
            ret.cigar.back().Length(ret.cigar.back().Length() + currCigar.front().Length());
        }
        ret.cigar.insert(ret.cigar.end(), currCigar.begin() + 1, currCigar.end());
        score += alnRegion.score;
    }
    ret.score = score;

    return ret;
}

OverlapPtr AlignmentSeeded(const OverlapPtr& ovl, const std::vector<AlignmentRegion>& alnRegions,
                           const std::string_view targetSeq, const std::string_view querySeqFwd,
                           const std::string_view querySeqRev, AlignerBasePtr& alignerGlobal,
                           AlignerBasePtr& alignerExt)
{
    const int32_t queryLen = querySeqFwd.size();
    const int32_t targetLen = targetSeq.size();

    // Sanity checks.
    if (ovl->Arev) {
        throw std::runtime_error("(AlignmentSeeded) The ovl->Arev should always be false!");
    }
    if (ovl->Alen != queryLen) {
        std::ostringstream oss;
        oss << "(AlignmentSeeded) The query length in the overlap is not the same as the provided "
               "sequence! ovl->Alen = "
            << ovl->Alen << ", queryLen = " << queryLen;
        throw std::runtime_error(oss.str());
    }
    if (ovl->Blen != targetLen) {
        std::ostringstream oss;
        oss << "(AlignmentSeeded) The target length in the overlap is not the same as the provided "
               "sequence! ovl->Blen = "
            << ovl->Blen << ", targetLen = " << targetLen;
        throw std::runtime_error(oss.str());
    }
    if (alnRegions.empty()) {
        std::ostringstream oss;
        oss << "(AlignmentSeeded) There needs to be at least one region to align.";
        throw std::runtime_error(oss.str());
    }

    // Run the alignment.
    AlignRegionsGenericResult alns = AlignRegionsGeneric(targetSeq, querySeqFwd, querySeqRev,
                                                         alnRegions, alignerGlobal, alignerExt);

    // Process the alignment results and make a new overlap.
    int32_t globalAlnQueryStart = 0;
    int32_t globalAlnTargetStart = 0;
    int32_t globalAlnQueryEnd = 0;
    int32_t globalAlnTargetEnd = 0;
    // Find the leftmost coordinate for global alignment.
    for (int32_t i = 0; i < static_cast<int32_t>(alnRegions.size()); ++i) {
        if (alnRegions[i].type != RegionType::GLOBAL) {
            continue;
        }
        globalAlnQueryStart = alnRegions[i].qStart;
        globalAlnTargetStart = alnRegions[i].tStart;
        break;
    }
    // Find the rightmost coordinate for global alignment.
    for (int32_t i = static_cast<int32_t>(alnRegions.size()) - 1; i >= 0; --i) {
        if (alnRegions[i].type != RegionType::GLOBAL) {
            continue;
        }
        globalAlnQueryEnd = alnRegions[i].qStart + alnRegions[i].qSpan;
        globalAlnTargetEnd = alnRegions[i].tStart + alnRegions[i].tSpan;
        break;
    }
    // Construct the new overlap.
    OverlapPtr ret = CreateOverlap(ovl);
    ret->Astart = globalAlnQueryStart - alns.offsetFrontQuery;
    ret->Aend = globalAlnQueryEnd + alns.offsetBackQuery;
    ret->Bstart = globalAlnTargetStart - alns.offsetFrontTarget;
    ret->Bend = globalAlnTargetEnd + alns.offsetBackTarget;
    ret->Cigar = std::move(alns.cigar);
    ret->Score = alns.score;

    // Validate here, before we reverse the coordinates and the CIGAR to satisfy the internal convention.
    {
        const char* querySeqForValidation = (ret->Brev) ? querySeqRev.data() : querySeqFwd.data();
        try {
            ValidateCigar(std::string_view(querySeqForValidation + ret->Astart, ret->ASpan()),
                          std::string_view(targetSeq.data() + ret->Bstart, ret->BSpan()),
                          ret->Cigar, "Full length validation.");
        } catch (const std::exception& e) {
            PBLOG_DEBUG << "[Note: Exception when aligning!] " << e.what() << "\n";
            PBLOG_DEBUG << "Q: " << std::string_view(querySeqForValidation, ret->ASpan()) << "\n";
            PBLOG_DEBUG << "T: " << std::string_view(targetSeq.data() + ret->Bstart, ret->BSpan())
                        << "\n";
            return nullptr;
        }
    }

    /**
     * IMPORTANT:
     * Internal convention is as follows:
     *      - Query is ALWAYS FWD in the internal representation.
     *      - Target can be fwd or rev - in the strand of the alignment.
     *      - CIGAR is oriented so that the query is always FWD and the target is complemented.
     *      - To avoid copying the target and reversing it, we need to reorient the coordinates here.
     *
     * Opposed to that, this function performs alignment in the following manner:
     *      - Query coordinate is in strand of alignment (fwd or rev).
     *      - Target coordinate is ALWAYS fwd.
     *      - CIGAR is oriented so that the target is always FWD and the query is copmlemented.
     * This is done so that we can reuse the query sequence for which we already have the reverse complement here.
     * Also, there is a general assumption that the query sequence is smaller than the target sequence.
     *
     * The following block is then intended to convert the conventions.
     *
     * TODO: Change the internal representation. It could be a substantial effort though.
     */
    // Reverse the CIGAR and the coordinates if needed.
    if (ovl->Brev) {
        // CIGAR reversal.
        std::reverse(ret->Cigar.begin(), ret->Cigar.end());

        // Reverse the query coordinates.
        std::swap(ret->Astart, ret->Aend);
        ret->Astart = ret->Alen - ret->Astart;
        ret->Aend = ret->Alen - ret->Aend;

        // Get the forward-oriented target coordinates.
        std::swap(ret->Bstart, ret->Bend);
        ret->Bstart = ret->Blen - ret->Bstart;
        ret->Bend = ret->Blen - ret->Bend;
    }

    // Set the alignment identity and edit distance.
    DiffCounts diffs = CigarDiffCounts(ret->Cigar);
    diffs.Identity(false, false, ret->Identity, ret->EditDistance);

    return ret;
}

}  // namespace Pancake
}  // namespace PacBio
