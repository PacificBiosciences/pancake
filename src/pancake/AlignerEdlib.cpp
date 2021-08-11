// Authors: Ivan Sovic

#include <pacbio/alignment/AlignmentTools.h>
#include <pacbio/pancake/AlignerEdlib.h>
#include <pbcopper/third-party/edlib.h>

namespace PacBio {
namespace Pancake {

std::shared_ptr<AlignerBase> CreateAlignerEdlib(const AlignmentParameters& opt)
{
    return std::shared_ptr<AlignerBase>(new AlignerEdlib(opt));
}

AlignerEdlib::AlignerEdlib(const AlignmentParameters& opt) : opt_(opt) {}

AlignerEdlib::~AlignerEdlib() {}

AlignmentResult AlignerEdlib::Global(const char* qseq, int64_t qlen, const char* tseq, int64_t tlen)
{
    if (qlen == 0 || tlen == 0) {
        AlignmentResult ret = EdgeCaseAlignmentResult(
            qlen, tlen, opt_.matchScore, opt_.mismatchPenalty, opt_.gapOpen1, opt_.gapExtend1);
        return ret;
    }

    EdlibAlignResult edlibResult = edlibAlign(
        qseq, qlen, tseq, tlen, edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));

    if (edlibResult.numLocations == 0) {
        edlibFreeAlignResult(edlibResult);
        return {};
    }

    AlignmentResult ret;

    Data::Cigar cigar =
        EdlibAlignmentToCigar(edlibResult.alignment, edlibResult.alignmentLength, ret.diffs);
    bool valid = true;

    try {
        cigar = NormalizeCigar(qseq, qlen, tseq, tlen, cigar);
    } catch (std::exception& e) {
        valid = false;
        cigar.clear();
    }

    ret.cigar = std::move(cigar);
    ret.score = ScoreCigarAlignment(ret.cigar, opt_.matchScore, opt_.mismatchPenalty, opt_.gapOpen1,
                                    opt_.gapExtend1);
    ret.valid = valid;
    ret.maxScore = ret.score;
    ret.zdropped = false;
    ret.lastQueryPos = qlen;
    ret.lastTargetPos = tlen;
    ret.maxQueryPos = qlen;
    ret.maxTargetPos = tlen;

    edlibFreeAlignResult(edlibResult);

    return ret;
}

AlignmentResult AlignerEdlib::Extend(const char* /*qseq*/, int64_t /*qlen*/, const char* /*tseq*/,
                                     int64_t /*tlen*/)
{
    AlignmentResult ret;
    ret.valid = false;
    ret.lastQueryPos = 0;
    ret.lastTargetPos = 0;
    return ret;
}

}  // namespace Pancake
}  // namespace PacBio
