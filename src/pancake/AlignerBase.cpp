// Authors: Ivan Sovic

#include <pacbio/alignment/AlignmentTools.h>
#include <pacbio/pancake/AlignerBase.h>

namespace PacBio {
namespace Pancake {

AlignmentResult EdgeCaseAlignmentResult(int32_t qlen, int32_t tlen, int32_t matchScore,
                                        int32_t mismatchPenalty, int32_t gapOpen, int32_t gapExtend)
{
    AlignmentResult ret;
    if (qlen == 0 && tlen > 0) {
        ret.cigar = PacBio::Data::Cigar(std::to_string(tlen) + "D");
        ret.valid = true;
        ret.lastQueryPos = 0;
        ret.lastTargetPos = tlen;
        ret.maxQueryPos = 0;
        ret.maxTargetPos = tlen;
        ret.score = ScoreCigarAlignment(ret.cigar, matchScore, mismatchPenalty, gapOpen, gapExtend);
        ret.maxScore = ret.score;
        ret.zdropped = false;
        ret.diffs = Alignment::DiffCounts(0, 0, 0, tlen);
    } else if (qlen > 0 && tlen == 0) {
        ret.cigar = PacBio::Data::Cigar(std::to_string(qlen) + "I");
        ret.valid = true;
        ret.lastQueryPos = qlen;
        ret.lastTargetPos = 0;
        ret.maxQueryPos = qlen;
        ret.maxTargetPos = 0;
        ret.score = ScoreCigarAlignment(ret.cigar, matchScore, mismatchPenalty, gapOpen, gapExtend);
        ret.maxScore = ret.score;
        ret.zdropped = false;
        ret.diffs = Alignment::DiffCounts(0, 0, qlen, 0);
    } else {
        ret.cigar = PacBio::Data::Cigar();
        ret.valid = true;
        ret.lastQueryPos = qlen;
        ret.lastTargetPos = tlen;
        ret.maxQueryPos = qlen;
        ret.maxTargetPos = tlen;
        ret.score = ScoreCigarAlignment(ret.cigar, matchScore, mismatchPenalty, gapOpen, gapExtend);
        ret.maxScore = ret.score;
        ret.zdropped = false;
        ret.diffs = Alignment::DiffCounts(0, 0, 0, 0);
    }
    return ret;
}

}  // namespace Pancake
}  // namespace PacBio
