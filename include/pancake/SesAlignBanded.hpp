// Author: Ivan Sovic

#ifndef PANCAKE_SES_ALIGN_BANDED_HPP
#define PANCAKE_SES_ALIGN_BANDED_HPP

#include <pancake/AlignmentTools.hpp>
#include <pancake/SesOptions.hpp>
#include <pancake/SesResults.hpp>

#include <pbbam/Cigar.h>
#include <pbbam/CigarOperation.h>

#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>

namespace PacBio {
namespace Pancake {

/// \brief Alignment function based on the O(nd) algorithm with modifications for banded alignment.
///         this implementation supports both global and semiglobal alignment modes, where semiglobal
///         doesn't penalize gaps at the end of either query or target.
///         Also, this implementation implements a traceback, so the CIGAR vector can be obtained.
///         Traceback and alignment modes are specified at compile time via the template parameters
///         for efficiency.
///
/// \param ALIGN_MODE Alignment mode: SESAlignMode::Global or SESAlignMode::Semiglobal.
/// \param TRACEBACK Specifies whether the traceback should be computed.
///                     Either: SESTracebackMode::Enabled or SESTracebackMode::Disabled.
/// \param query Query sequence as a C-type string consisting of [ACTG] characters.
/// \param target Target sequence as a C-type string consisting of [ACTG] characters.
/// \param maxDiffs Maximum allowed number of indel diffs for the alignment. This is a parameter
///         specified by the original algorithm ("D"). If you want the maximum number of diffs.
/// \param bandwidth Bandwidth for banded alignment.
/// \param ss Scratch space for computation. Can be nullptr if no scratch space is provided. In this
///             case, scratch space will be generated for one time use internally.
///             It's much more efficient to provide the scratch space to this function to prevent
///             unnecessary allocations when high efficiency is needed.
///             To create a scratch space externally exactly the same approach can be used as it
///             is done internally:
///                 auto ss = std::make_shared<SESScratchSpace>();
///             and then provide that to the function call.
///
/// \returns A SesResults object containing alignment coordinates, scores and CIGAR vector (if computed).
template <SESAlignMode ALIGN_MODE, SESTracebackMode TRACEBACK>
SesResults SESAlignBanded(std::string_view query, std::string_view target, const int32_t maxDiffs,
                          int32_t bandwidth, std::shared_ptr<SESScratchSpace> ss = nullptr)
{
    SesResults ret;

    if (ss == nullptr) {
        ss = std::make_shared<SESScratchSpace>();
    }

    const int32_t queryLen = query.size();
    const int32_t targetLen = target.size();

    bandwidth = std::min(bandwidth, maxDiffs);

    const int32_t maxAllowedDiffs = std::max(maxDiffs, bandwidth);
    const int32_t N = queryLen;                       // ss->N;
    const int32_t M = targetLen;                      // ss->M;
    const int32_t zero_offset = maxAllowedDiffs + 1;  // ss->zero_offset;
    const int32_t bandTolerance = bandwidth / 2 + 1;  // ss->bandTolerance;
    int32_t lastK = 0;                                // ss->lastK;
    int32_t lastD = 0;
    int32_t minK = 0;             // ss->minK;
    int32_t maxK = 0;             // ss->maxK;
    int32_t best_u = 0;           // ss->best_u;
    auto& v = ss->v;              // Working row.
    auto& u = ss->u;              // Banding.
    auto& v2 = ss->v2;            // Traceback matrix.
    auto& alnPath = ss->alnPath;  // Alignment path during traceback.
    auto& dStart = ss->dStart;    // Start of each diff's row in the v2 vector.

    const int32_t rowLen = (2 * maxAllowedDiffs + 3);

    if (rowLen > static_cast<int32_t>(v.capacity())) {
        v.resize(rowLen, 0);
        u.resize(rowLen, MINUS_INF);
        dStart.resize(rowLen, {0, 0});
        ss->alnPath.resize(rowLen);
    }
    // clang-format off
    if constexpr (TRACEBACK == SESTracebackMode::Enabled) {
        if ((rowLen * maxAllowedDiffs) > static_cast<int32_t>(v2.capacity())) {
            v2.resize(rowLen * maxAllowedDiffs);
        }
    }
    // clang-format on

    int32_t v2Pos = 0;
    int32_t prevK = -1;
    int32_t x = 0;

    v[zero_offset + 1] = 0;

    for (int32_t d = 0; d < maxDiffs; ++d) {
        ret.numDiffs = d;
        if ((maxK - minK) > bandwidth) {
            ret.valid = false;
            break;
        }

        // clang-format off
        if constexpr (TRACEBACK == SESTracebackMode::Enabled) {
            // Location where to store the traceback info.
            dStart[d] = {v2Pos, minK};
        }
        // clang-format on

        for (int32_t k = minK; k <= maxK; k += 2) {
            int32_t kz = k + zero_offset;
            if (k == minK || (k != maxK && v[kz - 1] < v[kz + 1])) {
                x = v[kz + 1];
                // clang-format off
                if constexpr (TRACEBACK == SESTracebackMode::Enabled) {
                    prevK = k + 1;
                }
                // clang-format on
            } else {
                x = v[kz - 1] + 1;
                // clang-format off
                if constexpr (TRACEBACK == SESTracebackMode::Enabled) {
                    prevK = k - 1;
                }
                // clang-format on
            }

            int32_t y = x - k;
            while (x < N && y < M && query[x] == target[y]) {
                ++x;
                ++y;
            }
            v[kz] = x;
            u[kz] = x + y;

            // clang-format off
            if constexpr (TRACEBACK == SESTracebackMode::Enabled) {
                v2[v2Pos] = {x, prevK};
                ++v2Pos;
            }
            // clang-format on

            ret.lastQueryPos = x;
            ret.lastTargetPos = y;
            lastK = k;
            lastD = d;

            if (best_u <= u[kz]) {
                best_u = u[kz];
            }

            // clang-format off
            if constexpr(ALIGN_MODE == SESAlignMode::Global) {
                if (x >= N && y >= M) {
                    ret.valid = true;
                    break;
                }
            } else {
                if (x >= N || y >= M) {
                    ret.valid = true;
                    break;
                }
            }
            // clang-format on
        }

        if (ret.valid) {
            break;
        }

        int32_t new_minK = maxK;
        int32_t new_maxK = minK;
        for (int32_t k1 = minK; k1 <= maxK; k1 += 2) {
            // Is there a bug here? Should this also have '&& u[k + zero_offset] <= (best_u + bandTolerance'?
            if (u[k1 + zero_offset] >= (best_u - bandTolerance)) {
                new_minK = std::min(k1, new_minK);
                new_maxK = std::max(k1, new_maxK);
            }
        }
        minK = new_minK - 1;
        maxK = new_maxK + 1;
    }

    // clang-format off
    if constexpr (TRACEBACK == SESTracebackMode::Enabled) {
        int32_t currD = lastD;
        int32_t currK = lastK;
        const int32_t numPoints = (currD + 1) * 2;
        int32_t currPoint = numPoints - 1;

        while (currD > 0) {
            const auto& elCurr = v2[dStart[currD].first + ((currK - dStart[currD].second) >> 1)];
            const auto& elPrev = v2[dStart[currD - 1].first + ((elCurr.prevK - dStart[currD - 1].second) >> 1)];

            int32_t y2 = elCurr.x2 - currK;
            alnPath[currPoint--] = {elCurr.x2, y2};

            int32_t yPrev = elPrev.x2 - elCurr.prevK;
            if (currK > elCurr.prevK) {
                alnPath[currPoint--] = {elPrev.x2 + 1, yPrev};
            } else {
                alnPath[currPoint--] = {elPrev.x2, yPrev + 1};
            }

            currK = elCurr.prevK;
            --currD;
        }
        // Handle the first two points separately so that the above while loop
        // can be iterated witohut one more branching.
        const auto& elCurr = ss->v2[dStart[currD].first + ((currK - dStart[currD].second) >> 2)];
        alnPath[1] = {elCurr.x2, elCurr.x2 - currK};
        alnPath[0] = {0, 0};

        // Convert the trace to CIGAR.
        ret.cigar.reserve(numPoints / 2);
        ret.diffCounts.Clear();
        PacBio::BAM::CigarOperationType op = PacBio::BAM::CigarOperationType::UNKNOWN_OP;
        PacBio::BAM::CigarOperationType prevOp = PacBio::BAM::CigarOperationType::UNKNOWN_OP;
        uint32_t count = 0;
        uint32_t prevCount = 0;

        auto ConvertMismatchesAndAppend = [&]() {
            if (ret.cigar.size() > 0 &&
                ((prevOp == PacBio::BAM::CigarOperationType::DELETION &&
                    ret.cigar.back().Type() == PacBio::BAM::CigarOperationType::INSERTION) ||
                    (prevOp == PacBio::BAM::CigarOperationType::INSERTION &&
                    ret.cigar.back().Type() == PacBio::BAM::CigarOperationType::DELETION))) {

                auto& last = ret.cigar.back();
                uint32_t lastCount = last.Length();
                int32_t minLen = std::min(prevCount, lastCount);  // Number of mismatches.
                int32_t leftHang = static_cast<int32_t>(lastCount) - minLen;   // Remaining indels to the left.
                int32_t rightHang = static_cast<int32_t>(prevCount) - minLen;  // Remaining indels to the right.
                if (leftHang == 0) {
                    last = PacBio::BAM::CigarOperation(
                        PacBio::BAM::CigarOperationType::SEQUENCE_MISMATCH, minLen);
                } else {
                    last.Length(leftHang);
                    ret.cigar.emplace_back(PacBio::BAM::CigarOperation(
                        PacBio::BAM::CigarOperationType::SEQUENCE_MISMATCH, minLen));
                }
                prevCount = rightHang;
                ret.diffCounts.numX += minLen;
                ret.diffCounts.numD -= minLen;
                ret.diffCounts.numI -= minLen;
            }
            if (prevCount > 0) {
                ret.cigar.emplace_back(PacBio::BAM::CigarOperation(prevOp, prevCount));
            }
        };

        for (int32_t i = 1; i < numPoints; ++i) {
            const auto& prevMove = ss->alnPath[i - 1];
            const auto& alnMove = ss->alnPath[i];
            if (prevMove.x == alnMove.x && prevMove.y == alnMove.y) {
                continue;
            }
            // Determine the CIGAR op.
            op = PacBio::BAM::CigarOperationType::UNKNOWN_OP;
            count = 0;
            if (prevMove.x == alnMove.x && prevMove.y != alnMove.y) {
                op = PacBio::BAM::CigarOperationType::DELETION;
                count = abs(alnMove.y - prevMove.y);
                ret.diffCounts.numD += count;
            } else if (prevMove.x != alnMove.x && prevMove.y == alnMove.y) {
                op = PacBio::BAM::CigarOperationType::INSERTION;
                count = abs(alnMove.x - prevMove.x);
                ret.diffCounts.numI += count;
            } else {
                op = PacBio::BAM::CigarOperationType::SEQUENCE_MATCH;
                count = abs(alnMove.x - prevMove.x);
                ret.diffCounts.numEq += count;
            }
            if (op != prevOp && prevOp != PacBio::BAM::CigarOperationType::UNKNOWN_OP) {
                ConvertMismatchesAndAppend();
                prevCount = 0;
            }
            prevCount += count;
            prevOp = op;
        }
        // Add the final CIGAR operation.
        ConvertMismatchesAndAppend();
        ret.numDiffs = ret.diffCounts.NumDiffs();
    }
    // clang-format on

    return ret;
}
}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SES_ALIGN_BANDED_HPP
