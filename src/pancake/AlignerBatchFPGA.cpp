// Copyright (c) 2019, Pacific Biosciences of California, Inc.
// All rights reserved.
// See LICENSE.txt.
//
// Contributions from NVIDIA are Copyright (c) 2021, NVIDIA Corporation.
// All rights reserved.
// SPDX-License-Identifier: BSD-3-Clause-Clear
//
// Authors: Ivan Sovic

#include <pacbio/alignment/AlignmentTools.h>
#include <pacbio/pancake/AlignerBatchFPGA.h>
#include <pacbio/pancake/Lookups.h>
#include <pacbio/util/Util.h>
#include <pbcopper/logging/Logging.h>
#include <pbcopper/utility/Stopwatch.h>
#include <cstdint>
#include <iostream>
#include <sstream>
#include "../third-party/huxlerate/cores/ksw2_cigar.hpp"
#include "../third-party/huxlerate/hugenomic.hpp"

namespace PacBio {
namespace Pancake {
namespace {
std::vector<uint8_t> ConvertSeqAlphabet(const char* seq, size_t seqlen, const int8_t* conv_table)
{
    std::vector<uint8_t> ret(seqlen);
    for (size_t i = 0; i < seqlen; i++) {
        ret[i] = (int8_t)conv_table[(uint8_t)seq[i]];
    }
    return ret;
}

void ConvertKSW2CigarToPbbam(const std::vector<hugenomic::ksw2c_cigar_op>& cigar,
                             const std::vector<uint8_t>& qseq, const std::vector<uint8_t>& tseq,
                             PacBio::Data::Cigar& retCigar, int32_t& retQueryAlignmentLen,
                             int32_t& retTargetAlignmentLen, Alignment::DiffCounts& retDiffs)
{
    retCigar.clear();
    retDiffs.Clear();
    retQueryAlignmentLen = 0;
    retTargetAlignmentLen = 0;

    std::array<int32_t, 9> counts{0, 0, 0, 0, 0, 0, 0, 0, 0};

    int32_t qPos = 0;
    int32_t tPos = 0;
    for (const hugenomic::ksw2c_cigar_op singleCigar : cigar) {
        char op = hugenomic::ksw2c_get_cigar_op_type(singleCigar);
        int32_t count = hugenomic::ksw2c_get_cigar_op_length(singleCigar);

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
}  // namespace

AlignerBatchFPGA::AlignerBatchFPGA(const AlignmentParameters& alnParams, uint32_t deviceId)
    : alnParams_(alnParams), deviceId_{deviceId}
{
    hugenomic::fpga_configure(hugenomic::hugenomic_core::ksw2_cigar, deviceId);
}

AlignerBatchFPGA::~AlignerBatchFPGA() { hugenomic::fpga_release(deviceId_); }

void AlignerBatchFPGA::Clear()
{
    alnResults_.clear();
    fpgaBatch_.clear();
}

void AlignerBatchFPGA::AddSequencePair(const char* query, int32_t queryLen, const char* target,
                                       int32_t targetLen)
{
    const std::vector<uint8_t> qseqInt = ConvertSeqAlphabet(query, queryLen, &BaseToTwobit[0]);
    const std::vector<uint8_t> tseqInt = ConvertSeqAlphabet(target, targetLen, &BaseToTwobit[0]);
    fpgaBatch_.emplace_back(std::move(qseqInt), std::move(tseqInt));
}

std::pair<int64_t, int64_t> AlignerBatchFPGA::AlignAll(const bool extension)
{
#ifdef PANCAKE_TIMINGS
    int64_t cpuTime = 0;
    int64_t fpgaTime = 0;
    PacBio::Utility::Stopwatch timer;
#endif
    // clang-format off
    int8_t match = alnParams_.matchScore;
    int8_t mismatch = -1*alnParams_.mismatchPenalty;
    std::vector<int8_t> SCORING_MATRIX = {
        match,    mismatch, mismatch, mismatch,  0,
        mismatch,   match,  mismatch, mismatch,  0,
        mismatch, mismatch,    match, mismatch,  0,
        mismatch, mismatch, mismatch,    match,  0,
               0,        0,        0,        0,  0
    };
    // clang-format on
    std::vector<hugenomic::ksw2c_response> huResults = hugenomic::ksw2c_batch(
        deviceId_, fpgaBatch_, alnParams_.gapOpen1, alnParams_.gapExtend1, alnParams_.gapOpen2,
        alnParams_.gapExtend2, SCORING_MATRIX, alnParams_.alignBandwidth,
        extension ? alnParams_.zdrop : -1, extension ? alnParams_.endBonus : -1, extension);

#ifdef PANCAKE_TIMINGS
    fpgaTime += timer.ElapsedNanoseconds();
    timer.Reset();
#endif

    for (int32_t i = 0; i < huResults.size(); ++i) {
        AlignmentResult ret;
        const auto& qseqInt = fpgaBatch_[i].first;
        const auto& tseqInt = fpgaBatch_[i].second;
        const int32_t qlen = qseqInt.size();
        const int32_t tlen = tseqInt.size();
        const auto& huResult = huResults[i];

        PacBio::Data::Cigar currCigar;
        int32_t qAlnLen = 0;
        int32_t tAlnLen = 0;
        ConvertKSW2CigarToPbbam(huResult.cigar, qseqInt, tseqInt, currCigar, qAlnLen, tAlnLen,
                                ret.diffs);
        // PBLOG_INFO << qAlnLen << ' ' << qlen << " | " << tAlnLen << ' ' << tlen << '\t'
        //            << currCigar.ToStdString();

        ret.cigar = std::move(currCigar);
        ret.valid =
            extension ||
            ((huResult.zdropped || ret.cigar.empty() || qAlnLen != qlen || tAlnLen != tlen) ? false
                                                                                            : true);
        // PBLOG_INFO << ret.valid;
        ret.lastQueryPos = qlen;
        ret.lastTargetPos = tlen;
        if (!extension) {
            ret.maxQueryPos = huResult.global_target_query_end_position;
            ret.maxTargetPos = huResult.global_query_target_end_position;
            ret.score = huResult.global_query_max_score;
            ret.maxScore = huResult.extension_max_score;
        } else {
            ret.maxQueryPos = huResult.extension_query_end_position;
            ret.maxTargetPos = huResult.extension_target_end_position;
            ret.maxScore = huResult.extension_max_score;
            ret.score = std::max(ret.maxScore, huResult.global_query_max_score);
        }
        ret.zdropped = huResult.zdropped;
        alnResults_.emplace_back(std::move(ret));
    }
#ifdef PANCAKE_TIMINGS
    cpuTime += timer.ElapsedNanoseconds();
    timer.Reset();
#endif

#ifdef PANCAKE_TIMINGS
    return std::make_pair(cpuTime, fpgaTime);
#else
    return std::make_pair(-1, -1);
#endif
}

}  // namespace Pancake
}  // namespace PacBio
