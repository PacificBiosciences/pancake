// Authors: Ivan Sovic

#include <pacbio/alignment/AlignmentTools.h>
#include <pacbio/pancake/AlignerKSW2.h>
#include <pacbio/pancake/Lookups.h>
#include <array>
#include <cstring>
#include <iostream>

namespace PacBio {
namespace Pancake {

mm_tbuf_t* mm_tbuf_init(void)
{
    mm_tbuf_t* b;
    b = (mm_tbuf_t*)calloc(1, sizeof(mm_tbuf_t));
    // if (!(mm_dbg_flag & 1)) b->km = km_init();              // This flag means: #define MM_DBG_NO_KALLOC     0x1
    return b;
}

void mm_tbuf_destroy(mm_tbuf_t* b)
{
    if (b == 0) return;
    km_destroy(b->km);
    free(b);
}

std::shared_ptr<AlignerBase> CreateAlignerKSW2(const AlignmentParameters& opt)
{
    return std::shared_ptr<AlignerBase>(new AlignerKSW2(opt));
}

AlignerKSW2::AlignerKSW2(const AlignmentParameters& opt)
    : opt_(opt), buffer_(mm_tbuf_init(), mm_tbuf_destroy)
{
    // Minimap2 alignment matrix.
    // From Minimap2: int sc_ambi; // score when one or both bases are "N"
    const int32_t scAmbi = 1;
    GenerateSimpleMatrix_(5, mat_, opt.matchScore, opt.mismatchPenalty, scAmbi);
}

AlignerKSW2::~AlignerKSW2() = default;

AlignmentResult AlignerKSW2::Global(const char* qseq, int64_t qlen, const char* tseq, int64_t tlen)
{
    if (qlen == 0 || tlen == 0) {
        AlignmentResult ret = EdgeCaseAlignmentResult(
            qlen, tlen, opt_.matchScore, opt_.mismatchPenalty, opt_.gapOpen1, opt_.gapExtend1);
        return ret;
    }

    const int32_t extra_flag = 0;

    // Bandwidth heuristic, as it is in Minimap2.
    const int32_t bw = (int)(opt_.alignBandwidth * 1.5 + 1.);

    // Memory allocations required for KSW2.

    // Convert the subsequence's alphabet from ACTG to [0123].
    const std::vector<uint8_t> qseqInt = ConvertSeqAlphabet_(qseq, qlen, &BaseToTwobit[0]);
    const std::vector<uint8_t> tseqInt = ConvertSeqAlphabet_(tseq, tlen, &BaseToTwobit[0]);

    // Compute the actual bandwidth. If this was a long join, we need to allow more room.
    const int32_t longestSpan = std::max(qlen, tlen);
    const int32_t spanDiff = std::abs(tlen - qlen);

    // This is slighlty modified compared to how it is in Minimap2:
    // const int32_t actualBandwidth =
    //     (h2.CheckFlagLongJoin() || spanDiff > (bw / 2)) ? longestSpan : bw;
    const int32_t actualBandwidth = (spanDiff > (bw / 2)) ? longestSpan : bw;

    AlignmentResult ret;

    std::vector<int32_t> bws;
    if (opt_.dynamicBandwidth) {
        bws.reserve(5);
        for (const int32_t curBw : (longestSpan < 500) ? std::vector<int32_t>{50, 100, 250}
                                                       : std::vector<int32_t>{250, 500}) {
            if (curBw < actualBandwidth) {
                bws.emplace_back(curBw);
            }
        }
    }
    bws.emplace_back(actualBandwidth);
    bws.emplace_back(longestSpan);
    for (int32_t curBw : bws) {
        ret = AlignmentResult{};
        ksw_extz_t ez;
        memset(&ez, 0, sizeof(ksw_extz_t));
        AlignPair_(buffer_->km, qlen, &qseqInt[0], tlen, &tseqInt[0], mat_, curBw, -1, -1,
                   extra_flag | KSW_EZ_APPROX_MAX, &ez, opt_.gapOpen1, opt_.gapExtend1,
                   opt_.gapOpen2, opt_.gapExtend2);

        PacBio::Data::Cigar currCigar;
        int32_t qAlnLen = 0;
        int32_t tAlnLen = 0;
        ConvertMinimap2CigarToPbbam_(ez.cigar, ez.n_cigar, qseqInt, tseqInt, currCigar, qAlnLen,
                                     tAlnLen, ret.diffs);

        // If full-width bandwidth is used, don't run the check.
        const bool isSuboptimal =
            (curBw != longestSpan) ? CheckAlignmentOutOfBand(currCigar, curBw) : false;

        ret.valid = !isSuboptimal && !ez.zdropped && (ez.n_cigar != 0) && (qAlnLen == qlen) &&
                    (tAlnLen == tlen);
        if (ret.valid) {
            ret.cigar = std::move(currCigar);
        }
        ret.lastQueryPos = qlen;
        ret.lastTargetPos = tlen;
        ret.maxQueryPos = ez.max_q;
        ret.maxTargetPos = ez.max_t;
        ret.score = ez.score;
        ret.maxScore = ez.max;
        ret.zdropped = ez.zdropped;

        // Free KSW2 memory.
        kfree(buffer_->km, ez.cigar);

        // Early return if bandwidth lead to an optimal alignment
        if (ret.valid) {
            break;
        }
    }
    return ret;
}

AlignmentResult AlignerKSW2::Extend(const char* qseq, int64_t qlen, const char* tseq, int64_t tlen)
{
    if (qlen == 0 || tlen == 0) {
        AlignmentResult ret = EdgeCaseAlignmentResult(
            qlen, tlen, opt_.matchScore, opt_.mismatchPenalty, opt_.gapOpen1, opt_.gapExtend1);
        return ret;
    }

    const int32_t extra_flag = 0;

    // Bandwidth heuristic, as it is in Minimap2.
    const int32_t bw = (int)(opt_.alignBandwidth * 1.5 + 1.);

    // Memory allocations required for KSW2.
    ksw_extz_t ez;
    memset(&ez, 0, sizeof(ksw_extz_t));

    // Convert the subsequence's alphabet from ACTG to [0123].
    const std::vector<uint8_t> qseqInt = ConvertSeqAlphabet_(qseq, qlen, &BaseToTwobit[0]);
    const std::vector<uint8_t> tseqInt = ConvertSeqAlphabet_(tseq, tlen, &BaseToTwobit[0]);

    AlignPair_(buffer_->km, qlen, &qseqInt[0], tlen, &tseqInt[0], mat_, bw, opt_.endBonus,
               opt_.zdrop, extra_flag | KSW_EZ_EXTZ_ONLY | KSW_EZ_RIGHT, &ez, opt_.gapOpen1,
               opt_.gapExtend1, opt_.gapOpen2, opt_.gapExtend2);

    AlignmentResult ret;

    PacBio::Data::Cigar currCigar;
    int32_t qAlnLen = 0;
    int32_t tAlnLen = 0;
    ConvertMinimap2CigarToPbbam_(ez.cigar, ez.n_cigar, qseqInt, tseqInt, currCigar, qAlnLen,
                                 tAlnLen, ret.diffs);

    ret.cigar = std::move(currCigar);
    ret.valid = true;
    ret.lastQueryPos = (ez.reach_end ? qlen : ez.max_q + 1);
    ret.lastTargetPos = (ez.reach_end ? ez.mqe_t + 1 : ez.max_t + 1);
    ret.maxQueryPos = ez.max_q;
    ret.maxTargetPos = ez.max_t;
    ret.zdropped = ez.zdropped;
    ret.maxScore = ez.max;
    ret.score = std::max(ret.maxScore, ez.score);

    // Free KSW2 memory.
    kfree(buffer_->km, ez.cigar);

    return ret;
}

std::vector<uint8_t> AlignerKSW2::ConvertSeqAlphabet_(const char* seq, size_t seqlen,
                                                      const int8_t* conv_table)
{
    std::vector<uint8_t> ret(seqlen);
    for (size_t i = 0; i < seqlen; i++) {
        ret[i] = (int8_t)conv_table[(uint8_t)seq[i]];
    }
    return ret;
}

void AlignerKSW2::ConvertMinimap2CigarToPbbam_(
    uint32_t* mm2Cigar, int32_t cigarLen, const std::vector<uint8_t>& qseq,
    const std::vector<uint8_t>& tseq, PacBio::Data::Cigar& retCigar, int32_t& retQueryAlignmentLen,
    int32_t& retTargetAlignmentLen, Alignment::DiffCounts& retDiffs)
{
    retCigar.clear();
    retDiffs.Clear();
    retQueryAlignmentLen = 0;
    retTargetAlignmentLen = 0;

    std::array<int32_t, 9> counts{0, 0, 0, 0, 0, 0, 0, 0, 0};

    int32_t qPos = 0;
    int32_t tPos = 0;
    for (int32_t opId = 0; opId < cigarLen; ++opId) {
        char op = "MIDNSH"[mm2Cigar[opId] & 0xf];
        int32_t count = mm2Cigar[opId] >> 4;
        // cigar.emplace_back(Data::CigarOperation(op, count));

        if (op == 'M') {
            // std::cerr << "[opId = " << opId << "] op: " << count << op << "\n";

            char prevOp = 'u';
            int32_t span = 0;
            for (int32_t m = 0; m < count; ++m) {
                char currOp = 'X';
                if (qseq[qPos + m] == tseq[tPos + m]) {
                    currOp = '=';
                }
                if (currOp == prevOp) {
                    ++span;
                } else {
                    if (span > 0) {
                        retCigar.emplace_back(Data::CigarOperation(prevOp, span));
                        counts[CigarCharToNum[static_cast<int32_t>(prevOp)]] += span;
                        // std::cerr << "  -> Added (mid):  " << span << prevOp << "\n";
                    }
                    span = 1;
                }
                prevOp = currOp;
                // ++qPos;
                // ++tPos;
            }
            if (span > 0) {
                retCigar.emplace_back(Data::CigarOperation(prevOp, span));
                counts[CigarCharToNum[static_cast<int32_t>(prevOp)]] += span;
                // std::cerr << "  -> Added (last): " << span << prevOp << "\n";
            }
            qPos += count;
            tPos += count;
        } else if (op == '=' || op == 'X') {
            retCigar.emplace_back(Data::CigarOperation(op, count));
            counts[CigarCharToNum[static_cast<int32_t>(op)]] += count;
            qPos += count;
            tPos += count;
        } else if (op == 'I' || op == 'S') {
            retCigar.emplace_back(Data::CigarOperation(op, count));
            counts[CigarCharToNum[static_cast<int32_t>(op)]] += count;
            qPos += count;
        } else if (op == 'D' || op == 'N') {
            retCigar.emplace_back(Data::CigarOperation(op, count));
            counts[CigarCharToNum[static_cast<int32_t>(op)]] += count;
            tPos += count;
        }
    }
    retQueryAlignmentLen = qPos;
    retTargetAlignmentLen = tPos;

    retDiffs.numEq = counts[CIGAR_OP_EQ];
    retDiffs.numX = counts[CIGAR_OP_X];
    retDiffs.numI = counts[CIGAR_OP_I];
    retDiffs.numD = counts[CIGAR_OP_D];
}

void AlignerKSW2::AlignPair_(void* km, int qlen, const uint8_t* qseq, int tlen, const uint8_t* tseq,
                             const int8_t* mat, int w, int endBonus, int zdrop, int flag,
                             ksw_extz_t* ez, int q, int e, int q2, int e2)
{
    if (q == q2 && e == e2)
        ksw_extz2_simde(km, qlen, qseq, tlen, tseq, 5, mat, q, e, w, zdrop, endBonus, flag, ez);
    else
        ksw_extd2_simde(km, qlen, qseq, tlen, tseq, 5, mat, q, e, q2, e2, w, zdrop, endBonus, flag,
                        ez);
}

void AlignerKSW2::GenerateSimpleMatrix_(int m, int8_t* mat, int8_t a, int8_t b, int8_t scAmbi)
{
    int i, j;
    a = a < 0 ? -a : a;
    b = b > 0 ? -b : b;
    scAmbi = scAmbi > 0 ? -scAmbi : scAmbi;
    for (i = 0; i < m - 1; ++i) {
        for (j = 0; j < m - 1; ++j)
            mat[i * m + j] = i == j ? a : b;
        mat[i * m + m - 1] = scAmbi;
    }
    for (j = 0; j < m; ++j)
        mat[(m - 1) * m + j] = scAmbi;
}

}  // namespace Pancake
}  // namespace PacBio
