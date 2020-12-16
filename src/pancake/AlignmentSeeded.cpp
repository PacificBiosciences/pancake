// Authors: Ivan Sovic

#include <pacbio/alignment/AlignmentTools.h>
#include <pacbio/pancake/AlignmentSeeded.h>
#include <pacbio/pancake/OverlapWriterBase.h>
#include <pbcopper/logging/Logging.h>
#include <iostream>

namespace PacBio {
namespace Pancake {

// #define DEBUG_ALIGNMENT_SEEDED

std::vector<AlignmentRegion> ExtractAlignmentRegions(const std::vector<SeedHit>& inSortedHits,
                                                     int32_t qLen, int32_t tLen, bool isRev,
                                                     int32_t minAlignmentSpan,
                                                     int32_t maxFlankExtensionDist,
                                                     double flankExtensionFactor)
{
    if (qLen < 0) {
        throw std::runtime_error("Invalid function parameter in ExtractAlignmentRegions! qLen = " +
                                 std::to_string(qLen) + ", should be >= 0.");
    }
    if (tLen < 0) {
        throw std::runtime_error("Invalid function parameter in ExtractAlignmentRegions! tLen = " +
                                 std::to_string(tLen) + ", should be >= 0.");
    }
    if (flankExtensionFactor < 1.0) {
        throw std::runtime_error(
            "Invalid function parameter in ExtractAlignmentRegions! flankExtensionFactor = " +
            std::to_string(flankExtensionFactor) + ", should be >= 1.0.");
    }

    if (inSortedHits.empty()) {
        return {};
    }

    if (maxFlankExtensionDist < 0) {
        maxFlankExtensionDist = std::max(qLen, tLen);
    }

    // NOTE: This is only required if the hit coordinates are always
    // fwd in the query sequence. If they are not, this can be removed.
    //
    // Reverse the hit coordinates if required, to make the
    // alignment simpler. This is required because the original coordinate
    // system keeps query in fwd and target in strand, whereas here we
    // use query in strand and target in fwd.
    std::vector<SeedHit> hits = inSortedHits;
    if (isRev) {
        for (size_t i = 0; i < hits.size(); ++i) {
            auto& hit = hits[i];
            // std::swap(hit.queryPos, hit.targetPos
            hit.queryPos = qLen - hit.queryPos;
            hit.targetPos = tLen - hit.targetPos;

            if (i > 0) {
                bool isLongCurr = hits[i].CheckFlagLongJoin();
                bool isLongPrev = hits[i - 1].CheckFlagLongJoin();
                hits[i - 1].SetFlagLongJoin(isLongCurr);
                hits[i].SetFlagLongJoin(isLongPrev);
            }
        }
        std::reverse(hits.begin(), hits.end());
    }

    std::vector<AlignmentRegion> ret;
    int32_t globalAlnQueryStart = hits.front().queryPos;
    int32_t globalAlnTargetStart = hits.front().targetPos;
    int32_t globalAlnQueryEnd = hits.back().queryPos;
    int32_t globalAlnTargetEnd = hits.back().targetPos;
    int32_t numRegions = 0;

    // Extract the front.
    if (globalAlnQueryStart > 0 && globalAlnTargetStart > 0) {
        // Determine the maximum flank we want to align, in both target and query.
        const int32_t projTStart =
            std::max(0.0, globalAlnTargetStart - flankExtensionFactor * globalAlnQueryStart);
        const int32_t projQStart = 0;
        const int32_t qExtLen = std::min(globalAlnQueryStart - projQStart, maxFlankExtensionDist);
        const int32_t tExtLen = std::min(globalAlnTargetStart - projTStart, maxFlankExtensionDist);

        // Find the frame of the sequences to compare.
        int32_t qStart = globalAlnQueryStart - qExtLen;
        int32_t tStart = globalAlnTargetStart - tExtLen;
        const int32_t qSpan = globalAlnQueryStart - qStart;
        const int32_t tSpan = globalAlnTargetStart - tStart;

        AlignmentRegion region;
        region.qStart = qStart;
        region.qSpan = qSpan;
        region.tStart = tStart;
        region.tSpan = tSpan;
        region.queryRev = isRev;
        region.type = RegionType::FRONT;
        region.regionId = numRegions;
        ++numRegions;
        ret.emplace_back(std::move(region));
    }

    // Align between seeds.
    const int32_t nHits = hits.size();
    std::vector<Data::Cigar> cigarChunks;
    int32_t startId = 0;
    for (int32_t i = 1; i < nHits; ++i) {
        // Shorthands.
        const auto& h1 = hits[startId];
        const auto& h2 = hits[i];

        // Compute the new region.
        AlignmentRegion region;
        region.qStart = h1.queryPos;
        region.qSpan = h2.queryPos - h1.queryPos;
        region.tStart = h1.targetPos;
        region.tSpan = h2.targetPos - h1.targetPos;
        region.type = RegionType::GLOBAL;
        region.queryRev = isRev;
        region.regionId = numRegions;

        // Sanity check.
        if (region.qSpan < 0 || region.tSpan < 0) {
            std::ostringstream oss;
            oss << "Region span not valid, in ExtractAlignmentRegions! qStart = " << region.qStart
                << ", qSpan = " << region.qSpan << ", tStart = " << region.tStart
                << ", tSpan = " << region.tSpan << "\n";
            throw std::runtime_error(oss.str());
        }

        // Skip short alignment portions.
        // if ((i + 1) < nHits && h2.CheckFlagLongJoin() == false &&
        if ((i + 1) < nHits &&
            (region.qSpan < minAlignmentSpan || region.tSpan < minAlignmentSpan)) {
            continue;
        }

        // Update the start ID for the next iteration.
        startId = i;
        ++numRegions;

        // Add the new region.
        ret.emplace_back(std::move(region));
    }

    // Back chunk.
    if (globalAlnQueryEnd < qLen && globalAlnTargetEnd < tLen) {
        // Determine the maximum flank we want to align, in both target and query.
        const int32_t qFlankLen = qLen - globalAlnQueryEnd;
        const int32_t projTEnd = std::min(
            tLen, static_cast<int32_t>(globalAlnTargetEnd + flankExtensionFactor * qFlankLen));
        const int32_t projQEnd = qLen;
        const int32_t qExtLen = std::min(projQEnd - globalAlnQueryEnd, maxFlankExtensionDist);
        const int32_t tExtLen = std::min(projTEnd - globalAlnTargetEnd, maxFlankExtensionDist);

        // Compute the new region.
        AlignmentRegion region;
        region.qStart = globalAlnQueryEnd;
        region.qSpan = qExtLen;
        region.tStart = globalAlnTargetEnd;
        region.tSpan = tExtLen;
        region.type = RegionType::BACK;
        region.queryRev = isRev;
        region.regionId = numRegions;
        ++numRegions;
        ret.emplace_back(std::move(region));
    }

    return ret;
}

AlignmentResult AlignSingleRegion(const char* targetSeq, int32_t targetLen, const char* querySeqFwd,
                                  const char* querySeqRev, int32_t queryLen,
                                  AlignerBasePtr& alignerGlobal, AlignerBasePtr& alignerExt,
                                  const AlignmentRegion& region)
{
    const char* querySeqInStrand = region.queryRev ? querySeqRev : querySeqFwd;
    const char* targetSeqInStrand = targetSeq;
    int32_t qStart = region.qStart;
    int32_t tStart = region.tStart;
    const int32_t qSpan = region.qSpan;
    const int32_t tSpan = region.tSpan;

    if (qSpan == 0 && tSpan == 0) {
        return {};
    }

    if (qStart >= queryLen || (qStart + qSpan) > queryLen || tStart >= targetLen ||
        (tStart + tSpan) > targetLen || qSpan < 0 || tSpan < 0) {
        std::ostringstream oss;
        oss << "AlignmentRegion coordinates out of bounds in AlignRegionsGeneric!";
        oss << " qStart = " << qStart << ", qSpan = " << qSpan << ", queryLen = " << queryLen
            << ", tStart = " << tStart << ", tSpan = " << tSpan << ", targetLen = " << targetLen
            << ", regionId = " << region.regionId;
        throw std::runtime_error(oss.str());
    }
    if (targetSeq == NULL || querySeqFwd == NULL || querySeqRev == NULL) {
        std::ostringstream oss;
        oss << "NULL sequence passed to AlignmentRegion!";
        oss << " qStart = " << qStart << ", qSpan = " << qSpan << ", queryLen = " << queryLen
            << ", tStart = " << tStart << ", tSpan = " << tSpan << ", targetLen = " << targetLen
            << ", regionId = " << region.regionId;
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
        alnRes =
            alignerExt->Extend(querySeqInStrand + qStart, qSpan, targetSeqInStrand + tStart, tSpan);
    } else {
        alnRes = alignerGlobal->Global(querySeqInStrand + qStart, qSpan, targetSeqInStrand + tStart,
                                       tSpan);
    }

    if (region.type == RegionType::FRONT) {
        std::reverse(alnRes.cigar.begin(), alnRes.cigar.end());
    }

    return alnRes;
}

AlignRegionsGenericResult AlignRegionsGeneric(const char* targetSeq, const int32_t targetLen,
                                              const char* queryFwd, const char* queryRev,
                                              const int32_t queryLen,
                                              const std::vector<AlignmentRegion>& regions,
                                              AlignerBasePtr& alignerGlobal,
                                              AlignerBasePtr& alignerExt)
{
    AlignRegionsGenericResult ret;

    std::vector<AlignmentResult> alignedRegions;

    for (size_t i = 0; i < regions.size(); ++i) {
        const auto& region = regions[i];

        auto alnRes = AlignSingleRegion(targetSeq, targetLen, queryFwd, queryRev, queryLen,
                                        alignerGlobal, alignerExt, region);

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
    }

    return ret;
}

OverlapPtr AlignmentSeeded(const OverlapPtr& ovl, const std::vector<SeedHit>& sortedHits,
                           const char* targetSeq, const int32_t targetLen, const char* queryFwd,
                           const char* queryRev, const int32_t queryLen, int32_t minAlignmentSpan,
                           int32_t maxFlankExtensionDist, AlignerBasePtr& alignerGlobal,
                           AlignerBasePtr& alignerExt)
{
    // Sanity checks.
    if (ovl->Arev) {
        throw std::runtime_error("The ovl->Arev should always be false! In Align_.");
    }
    if (sortedHits.empty()) {
        throw std::runtime_error("There are no seed hits provided to Align_ for alignment.");
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
               "sequence! ovl->Alen = "
            << ovl->Blen << ", queryLen = " << targetLen;
        throw std::runtime_error(oss.str());
    }
    if (ovl->Astart != sortedHits.front().queryPos || ovl->Aend != sortedHits.back().queryPos ||
        ovl->Bstart != sortedHits.front().targetPos || ovl->Bend != sortedHits.back().targetPos) {
        std::ostringstream oss;
        oss << "(AlignmentSeeded) Provided overlap coordinates do not match the first/last seed "
               "hit!"
            << " ovl: " << OverlapWriterBase::PrintOverlapAsM4(ovl, "", "", true, false)
            << "; sortedHits.front() = " << sortedHits.front()
            << "; sortedHits.back() = " << sortedHits.back();
        throw std::runtime_error(oss.str());
    }
    if (ovl->Brev != sortedHits.front().targetRev || ovl->Brev != sortedHits.back().targetRev) {
        std::ostringstream oss;
        oss << "(AlignmentSeeded) Strand of the provided overlap does not match the first/last "
               "seed hit."
            << " ovl->Brev = " << (ovl->Brev ? "true" : "false")
            << ", sortedHits.front().targetRev = "
            << (sortedHits.front().targetRev ? "true" : "false")
            << ", sortedHits.back().targetRev = "
            << (sortedHits.back().targetRev ? "true" : "false");
        throw std::runtime_error(oss.str());
    }

    // Prepare the regions for alignment.
    std::vector<AlignmentRegion> regions = ExtractAlignmentRegions(
        sortedHits, ovl->Alen, ovl->Blen, ovl->Brev, minAlignmentSpan, maxFlankExtensionDist, 1.3);

    // Run the alignment.
    AlignRegionsGenericResult alns = AlignRegionsGeneric(
        targetSeq, targetLen, queryFwd, queryRev, queryLen, regions, alignerGlobal, alignerExt);

    // Process the alignment results and make a new overlap.
    int32_t globalAlnQueryStart = 0;
    int32_t globalAlnTargetStart = 0;
    int32_t globalAlnQueryEnd = 0;
    int32_t globalAlnTargetEnd = 0;
    // Find the leftmost coordinate for global alignment.
    for (int32_t i = 0; i < static_cast<int32_t>(regions.size()); ++i) {
        if (regions[i].type != RegionType::GLOBAL) {
            continue;
        }
        globalAlnQueryStart = regions[i].qStart;
        globalAlnTargetStart = regions[i].tStart;
        break;
    }
    // Find the rightmost coordinate for global alignment.
    for (int32_t i = static_cast<int32_t>(regions.size()) - 1; i >= 0; --i) {
        if (regions[i].type != RegionType::GLOBAL) {
            continue;
        }
        globalAlnQueryEnd = regions[i].qStart + regions[i].qSpan;
        globalAlnTargetEnd = regions[i].tStart + regions[i].tSpan;
        break;
    }
    // Construct the new overlap.
    OverlapPtr ret = createOverlap(ovl);
    ret->Cigar.clear();
    ret->Astart = globalAlnQueryStart - alns.offsetFrontQuery;
    ret->Aend = globalAlnQueryEnd + alns.offsetBackQuery;
    ret->Bstart = globalAlnTargetStart - alns.offsetFrontTarget;
    ret->Bend = globalAlnTargetEnd + alns.offsetBackTarget;
    ret->Cigar = std::move(alns.cigar);

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
    Alignment::DiffCounts diffs = CigarDiffCounts(ret->Cigar);
    diffs.Identity(false, false, ret->Identity, ret->EditDistance);
    ret->Score = -diffs.numEq;

    const int32_t qstart = ret->Astart;
    const int32_t tstart = ret->Bstart;
    const char* querySeqFwd = queryFwd;
    const std::string targetSeqForValidation =
        (ret->Brev) ? PacBio::Pancake::ReverseComplement(targetSeq, 0, targetLen) : targetSeq;

    // Validation. In case an alignment was dropped.
    try {
        ValidateCigar(querySeqFwd + qstart, ret->ASpan(), targetSeqForValidation.c_str() + tstart,
                      ret->BSpan(), ret->Cigar, "Full length validation.");
    } catch (std::exception& e) {
        // std::cerr << "[Note: Exception when aligning!] " << e.what() << "\n";
        // std::cerr << "Input overlap: "
        //           << OverlapWriterBase::PrintOverlapAsM4(ovl, "", "", true, true) << "\n";
        // std::cerr << "ASpan = " << ret->ASpan() << ", BSpan = " << ret->BSpan() << "\n";
        // std::cerr << "Q: " << std::string(querySeqFwd, queryLen) << "\n";
        // std::cerr << "T: " << targetSeqForValidation << "\n";
        // std::cerr << "\n";
        // std::cerr << "regions.size() = " << regions.size() << "\n";
        // for (size_t i = 0; i < regions.size(); ++i) {
        //     std::cerr << "[region i = " << i << "] " << regions[i] << "\n";
        // }
        // std::cerr.flush();
        PBLOG_DEBUG << "[Note: Exception when aligning!] " << e.what() << "\n";
        PBLOG_DEBUG << "Q: " << std::string(querySeqFwd, ret->ASpan()) << "\n";
        PBLOG_DEBUG << "T: " << targetSeqForValidation << "\n";
        ret = nullptr;
    }

    return ret;
}

}  // namespace Pancake
}  // namespace PacBio
