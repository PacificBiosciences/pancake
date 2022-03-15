// Authors: Ivan Sovic

#include <pancake/AlignerSES2.hpp>

#include <pancake/AlignmentTools.hpp>
#include <pancake/Ses2AlignBanded.hpp>

#include <pbcopper/third-party/edlib.h>

namespace PacBio {
namespace Pancake {

std::shared_ptr<AlignerBase> CreateAlignerSES2(const AlignmentParameters& opt)
{
    return std::shared_ptr<AlignerBase>(new AlignerSES2(opt));
}

AlignerSES2::AlignerSES2(const AlignmentParameters& opt)
    : opt_(opt), sesScratch_{std::make_shared<Pancake::SESScratchSpace>()}
{}

AlignerSES2::~AlignerSES2() {}

AlignmentResult AlignerSES2::Global(const std::string_view qseq, const std::string_view tseq)
{
    if (qseq.empty() || tseq.empty()) {
        AlignmentResult ret =
            EdgeCaseAlignmentResult(qseq.size(), tseq.size(), opt_.matchScore, opt_.mismatchPenalty,
                                    opt_.gapOpen1, opt_.gapExtend1);
        return ret;
    }

    const double alignMaxDiff = 1.0;  // 0.30;
    // const double alignBandwidth = 0.30;
    const int32_t maxDiffs = std::max(10, static_cast<int32_t>(qseq.size() * alignMaxDiff));
    // const int32_t bandwidth =
    //     std::max(10, static_cast<int32_t>(std::min(tseq.size(), qseq.size()) * alignBandwidth));

    // Compute the actual bandwidth. If this was a long join, we need to allow more room.
    // const int32_t bw = opt_.alignBandwidth;  // * 1.5 + 1.;
    // const int32_t longestSpan = std::max(qseq.size(), tseq.size());
    // const int32_t spanDiff = std::abs(tseq.size() - qseq.size());
    // const int32_t actualBandwidth = ((static_cast<double>(spanDiff) / static_cast<double>(longestSpan)) > 0.05) ? longestSpan : bw;
    const int32_t actualBandwidth = qseq.size() + tseq.size();

    auto aln =
        SES2AlignBanded<SESAlignMode::Global, SESTrimmingMode::Disabled, SESTracebackMode::Enabled>(
            qseq, tseq, maxDiffs, actualBandwidth, sesScratch_);

    Data::Cigar cigar = NormalizeCigar(qseq, tseq, aln.cigar);

    AlignmentResult ret;
    ret.cigar = std::move(cigar);
    ret.score = ScoreCigarAlignment(ret.cigar, opt_.matchScore, opt_.mismatchPenalty, opt_.gapOpen1,
                                    opt_.gapExtend1);
    ret.valid = aln.valid;
    ret.maxScore = ret.score;
    ret.zdropped = false;
    ret.lastQueryPos = qseq.size();
    ret.lastTargetPos = tseq.size();
    ret.maxQueryPos = qseq.size();
    ret.maxTargetPos = tseq.size();
    ret.diffs = CigarDiffCounts(ret.cigar);

    if (ret.valid == false) {
        ret.cigar.clear();
    }

    return ret;
}

AlignmentResult AlignerSES2::Extend(const std::string_view /*qseq*/,
                                    const std::string_view /*tseq*/)
{
    AlignmentResult ret;
    ret.valid = false;
    ret.lastQueryPos = 0;
    ret.lastTargetPos = 0;
    return ret;
}

}  // namespace Pancake
}  // namespace PacBio
