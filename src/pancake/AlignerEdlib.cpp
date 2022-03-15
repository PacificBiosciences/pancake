// Authors: Ivan Sovic

#include <pancake/AlignerEdlib.hpp>

#include <pancake/AlignmentTools.hpp>

#include <pbcopper/third-party/edlib.h>

namespace PacBio {
namespace Pancake {

std::shared_ptr<AlignerBase> CreateAlignerEdlib(const AlignmentParameters& opt)
{
    return std::shared_ptr<AlignerBase>(new AlignerEdlib(opt));
}

AlignerEdlib::AlignerEdlib(const AlignmentParameters& opt) : opt_(opt) {}

AlignerEdlib::~AlignerEdlib() {}

AlignmentResult AlignerEdlib::Global(const std::string_view qseq, const std::string_view tseq)
{
    if (qseq.empty() || tseq.empty()) {
        AlignmentResult ret =
            EdgeCaseAlignmentResult(qseq.size(), tseq.size(), opt_.matchScore, opt_.mismatchPenalty,
                                    opt_.gapOpen1, opt_.gapExtend1);
        return ret;
    }

    EdlibAlignResult edlibResult =
        edlibAlign(qseq.data(), qseq.size(), tseq.data(), tseq.size(),
                   edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));

    if (edlibResult.numLocations == 0) {
        edlibFreeAlignResult(edlibResult);
        return {};
    }

    AlignmentResult ret;

    Data::Cigar cigar = EdlibAlignmentToCigar(
        {edlibResult.alignment, static_cast<size_t>(edlibResult.alignmentLength)}, ret.diffs);
    bool valid = true;

    try {
        cigar = NormalizeCigar(qseq, tseq, cigar);
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
    ret.lastQueryPos = qseq.size();
    ret.lastTargetPos = tseq.size();
    ret.maxQueryPos = qseq.size();
    ret.maxTargetPos = tseq.size();

    edlibFreeAlignResult(edlibResult);

    return ret;
}

AlignmentResult AlignerEdlib::Extend([[maybe_unused]] const std::string_view qseq,
                                     [[maybe_unused]] const std::string_view tseq)
{
    AlignmentResult ret;
    ret.valid = false;
    ret.lastQueryPos = 0;
    ret.lastTargetPos = 0;
    return ret;
}

}  // namespace Pancake
}  // namespace PacBio
