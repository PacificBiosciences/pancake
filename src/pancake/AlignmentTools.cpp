// Authors: Ivan Sovic

#include <pancake/AlignmentTools.hpp>

#include <pancake/Lookups.hpp>

#include <pbcopper/third-party/edlib.h>
#include <pbcopper/utility/Ssize.h>

#include <array>
#include <cassert>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <string>

namespace PacBio {
namespace Pancake {

Data::Cigar EdlibAlignmentToCigar(const std::span<const unsigned char> aln, DiffCounts& retDiffs)
{
    retDiffs.Clear();

    const int32_t alnLen = aln.size();
    if (alnLen == 0) {
        return {};
    }

    // Edlib move codes: 0: '=', 1: 'I', 2: 'D', 3: 'X'
    const std::array<Data::CigarOperationType, 4> opToCigar = {
        Data::CigarOperationType::SEQUENCE_MATCH, Data::CigarOperationType::INSERTION,
        Data::CigarOperationType::DELETION, Data::CigarOperationType::SEQUENCE_MISMATCH};

    std::array<int32_t, 4> counts{0, 0, 0, 0};

    Data::CigarOperationType prevOp = Data::CigarOperationType::UNKNOWN_OP;
    unsigned char prevOpRaw = 0;
    int32_t count = 0;
    Data::Cigar ret;

    for (int32_t i = 0; i <= alnLen; i++) {
        if (i == alnLen ||
            (opToCigar[aln[i]] != prevOp && prevOp != Data::CigarOperationType::UNKNOWN_OP)) {
            ret.emplace_back(Data::CigarOperation(prevOp, count));
            counts[prevOpRaw] += count;
            count = 0;
        }
        if (i < alnLen) {
            prevOp = opToCigar[aln[i]];
            prevOpRaw = aln[i];
            count += 1;
        }
    }
    retDiffs.numEq = counts[EDLIB_EDOP_MATCH];
    retDiffs.numX = counts[EDLIB_EDOP_MISMATCH];
    retDiffs.numI = counts[EDLIB_EDOP_INSERT];
    retDiffs.numD = counts[EDLIB_EDOP_DELETE];
    return ret;
}

DiffCounts EdlibAlignmentDiffCounts(const std::span<const unsigned char> aln)
{
    const int32_t alnLen = aln.size();

    if (alnLen <= 0) {
        return {};
    }

    DiffCounts ret;
    for (int32_t i = 0; i < alnLen; i++) {
        switch (aln[i]) {
            case EDLIB_EDOP_MATCH:
                ++ret.numEq;
                break;
            case EDLIB_EDOP_MISMATCH:
                ++ret.numX;
                break;
            case EDLIB_EDOP_INSERT:
                ++ret.numI;
                break;
            case EDLIB_EDOP_DELETE:
                ++ret.numD;
                break;
            default:
                throw std::runtime_error("Unknown Edlib operation: " +
                                         std::to_string(static_cast<int32_t>(aln[i])));
                break;
        }
    }
    return ret;
}

DiffCounts CigarDiffCounts(const Data::Cigar& cigar)
{
    DiffCounts ret;
    for (const auto& op : cigar) {
        if (op.Type() == Data::CigarOperationType::SEQUENCE_MATCH) {
            ret.numEq += op.Length();
        } else if (op.Type() == Data::CigarOperationType::SEQUENCE_MISMATCH) {
            ret.numX += op.Length();
        } else if (op.Type() == Data::CigarOperationType::INSERTION) {
            ret.numI += op.Length();
        } else if (op.Type() == Data::CigarOperationType::DELETION) {
            ret.numD += op.Length();
        }
    }
    return ret;
}

void AppendToCigar(Data::Cigar& cigar, const Data::CigarOperationType newOp, const int32_t newLen)
{
    if (newLen == 0) {
        return;
    }
    if (cigar.empty() || newOp != cigar.back().Type()) {
        cigar.emplace_back(Data::CigarOperation(newOp, newLen));
    } else {
        cigar.back().Length(cigar.back().Length() + newLen);
    }
}

Data::Cigar ExpandMismatches(const std::string_view query, const std::string_view target,
                             const Data::Cigar& cigar)
{
    if (cigar.size() <= 1) {
        return cigar;
    }

    const int64_t queryLen = query.size();
    const int64_t targetLen = target.size();

    int64_t queryPos = 0;
    int64_t targetPos = 0;
    int32_t lastAddedOp = -1;
    int32_t numCigarOps = cigar.size();
    Data::Cigar ret;

    for (int32_t i = 1; i < numCigarOps; ++i) {
        const auto& prevOp = cigar[i - 1];
        const auto& currOp = cigar[i];

        if (queryPos >= queryLen || targetPos >= targetLen) {
            std::ostringstream oss;
            oss << "Invalid CIGAR string: queryPos = " << queryPos << ", targetPos = " << targetPos
                << ", queryLen = " << queryLen << ", targetLen = " << targetLen;
            throw std::runtime_error(oss.str());
        }

        // Check if we have an INS+DEL or DEL+INS pair. If so, we'll convert them
        // into a single diagonal set of MATCH/MISMATCH operations, plus the left
        // hang and right hang indel operations.
        if ((prevOp.Type() == Data::CigarOperationType::INSERTION &&
             currOp.Type() == Data::CigarOperationType::DELETION) ||
            (prevOp.Type() == Data::CigarOperationType::DELETION &&
             currOp.Type() == Data::CigarOperationType::INSERTION)) {

            const int32_t minLen = std::min(prevOp.Length(), currOp.Length());
            const int32_t leftHang = static_cast<int32_t>(prevOp.Length()) - minLen;
            const int32_t rightHang = static_cast<int32_t>(currOp.Length()) - minLen;

            AppendToCigar(ret, prevOp.Type(), leftHang);
            if (prevOp.Type() == Data::CigarOperationType::DELETION) {
                targetPos += leftHang;
            } else {
                queryPos += leftHang;
            }

            for (int32_t pos = 0; pos < minLen; ++pos) {
                if (query[queryPos] == target[targetPos]) {
                    AppendToCigar(ret, Data::CigarOperationType::SEQUENCE_MATCH, 1);
                } else {
                    AppendToCigar(ret, Data::CigarOperationType::SEQUENCE_MISMATCH, 1);
                }
                ++queryPos;
                ++targetPos;
            }

            AppendToCigar(ret, currOp.Type(), rightHang);
            if (currOp.Type() == Data::CigarOperationType::DELETION) {
                targetPos += rightHang;
            } else {
                queryPos += rightHang;
            }
            lastAddedOp = i;
            ++i;
            continue;
        } else {
            AppendToCigar(ret, prevOp.Type(), prevOp.Length());
            lastAddedOp = i - 1;
            if (prevOp.Type() != Data::CigarOperationType::DELETION) {
                queryPos += prevOp.Length();
            }
            if (prevOp.Type() != Data::CigarOperationType::INSERTION) {
                targetPos += prevOp.Length();
            }
        }
    }
    // Any remaining ops just get passed in.
    for (int32_t i = lastAddedOp + 1; i < numCigarOps; ++i) {
        AppendToCigar(ret, cigar[i].Type(), cigar[i].Length());
    }
    return ret;
}

void ValidateCigar(const std::string_view query, const std::string_view target,
                   const Data::Cigar& cigar, const std::string_view label)
{
    const int64_t queryLen = query.size();
    const int64_t targetLen = target.size();
    const char* queryData = query.data();
    const char* targetData = target.data();
    const int32_t numCigarOps = cigar.size();
    int64_t queryPos = 0;
    int64_t targetPos = 0;

    for (int32_t i = 0; i < numCigarOps; ++i) {
        const auto& op = cigar[i];

        if (queryPos > queryLen || targetPos > targetLen) {
            std::ostringstream oss;
            oss << "Invalid CIGAR string (global): "
                << "coordinates out of bounds! "
                << "queryPos = " << queryPos << ", targetPos = " << targetPos
                << ", offending CIGAR op: " << op.Length() << op.TypeToChar(op.Type())
                << ", queryLen = " << queryLen << ", targetLen = " << targetLen
                << ", CIGAR: " << cigar.ToStdString() << ", label: '" << label << "'";
            throw std::runtime_error(oss.str());
        }

        if (op.Type() == Data::CigarOperationType::SEQUENCE_MATCH) {
            if ((queryPos + op.Length()) > queryLen || (targetPos + op.Length()) > targetLen) {
                std::ostringstream oss;
                oss << "Invalid CIGAR string (SEQUENCE_MATCH): "
                    << "coordinates out of bounds! "
                    << "queryPos = " << queryPos << ", targetPos = " << targetPos
                    << ", offending CIGAR op: " << op.Length() << op.TypeToChar(op.Type())
                    << ", queryLen = " << queryLen << ", targetLen = " << targetLen
                    << ", CIGAR: " << cigar.ToStdString() << ", label: '" << label << "'";
                throw std::runtime_error(oss.str());
            }
            if (std::strncmp(queryData + queryPos, targetData + targetPos, op.Length()) != 0) {
                std::ostringstream oss;
                oss << "Invalid CIGAR string (SEQUENCE_MATCH): "
                    << "sequences are not equal even though they are delimited by a SEQUENCE_MATCH "
                       "operation! "
                    << "queryPos = " << queryPos << ", targetPos = " << targetPos
                    << ", offending CIGAR op: " << op.Length() << op.TypeToChar(op.Type())
                    << ", queryLen = " << queryLen << ", targetLen = " << targetLen
                    << ", CIGAR: " << cigar.ToStdString() << ", label: '" << label << "'";
                throw std::runtime_error(oss.str());
            }
            queryPos += op.Length();
            targetPos += op.Length();
        } else if (op.Type() == Data::CigarOperationType::SEQUENCE_MISMATCH) {
            if ((queryPos + op.Length()) > queryLen || (targetPos + op.Length()) > targetLen) {
                std::ostringstream oss;
                oss << "Invalid CIGAR string (SEQUENCE_MISMATCH): "
                    << "coordinates out of bounds! "
                    << "queryPos = " << queryPos << ", targetPos = " << targetPos
                    << ", offending CIGAR op: " << op.Length() << op.TypeToChar(op.Type())
                    << ", queryLen = " << queryLen << ", targetLen = " << targetLen
                    << ", CIGAR: " << cigar.ToStdString() << ", label: '" << label << "'";
                throw std::runtime_error(oss.str());
            }

            for (int64_t pos = 0; pos < static_cast<int64_t>(op.Length()); ++pos) {
                if (queryData[queryPos + pos] == targetData[targetPos + pos]) {
                    std::ostringstream oss;
                    oss << "Invalid CIGAR string (SEQUENCE_MISMATCH): "
                        << "sequences are equal even though they are delimited by a "
                           "SEQUENCE_MISMATCH operation! "
                        << "queryPos = " << queryPos << ", targetPos = " << targetPos
                        << ", offending CIGAR op: " << op.Length() << op.TypeToChar(op.Type())
                        << ", queryLen = " << queryLen << ", targetLen = " << targetLen
                        << ", CIGAR: " << cigar.ToStdString() << ", label: '" << label << "'";
                    throw std::runtime_error(oss.str());
                }
            }
            queryPos += op.Length();
            targetPos += op.Length();
        } else if (op.Type() == Data::CigarOperationType::INSERTION ||
                   op.Type() == Data::CigarOperationType::SOFT_CLIP) {
            if ((queryPos + op.Length()) > queryLen) {
                std::ostringstream oss;
                oss << "Invalid CIGAR string (INSERTION): "
                    << "coordinates out of bounds! "
                    << "queryPos = " << queryPos << ", targetPos = " << targetPos
                    << ", offending CIGAR op: " << op.Length() << op.TypeToChar(op.Type())
                    << ", queryLen = " << queryLen << ", targetLen = " << targetLen
                    << ", CIGAR: " << cigar.ToStdString() << ", label: '" << label << "'";
                throw std::runtime_error(oss.str());
            }
            queryPos += op.Length();
        } else if (op.Type() == Data::CigarOperationType::DELETION ||
                   op.Type() == Data::CigarOperationType::REFERENCE_SKIP) {
            if ((targetPos + op.Length()) > targetLen) {
                std::ostringstream oss;
                oss << "Invalid CIGAR string (DELETION): "
                    << "coordinates out of bounds! "
                    << "queryPos = " << queryPos << ", targetPos = " << targetPos
                    << ", offending CIGAR op: " << op.Length() << op.TypeToChar(op.Type())
                    << ", queryLen = " << queryLen << ", targetLen = " << targetLen
                    << ", CIGAR: " << cigar.ToStdString() << ", label: '" << label << "'";
                throw std::runtime_error(oss.str());
            }
            targetPos += op.Length();
        } else if (op.Type() == Data::CigarOperationType::HARD_CLIP) {
            // Do nothing.
        } else {
            std::ostringstream oss;
            oss << "CIGAR operation '" << op.TypeToChar(op.Type())
                << "' not supported by the validation function.";
            throw std::runtime_error(oss.str());
        }
    }

    if (queryPos != queryLen || targetPos != targetLen) {
        const DiffCounts diffs = CigarDiffCounts(cigar);
        std::ostringstream oss;
        oss << "Invalid CIGAR string (length): "
            << "Computed query or target length does not match the sequence length. "
            << "queryPos = " << queryPos << ", targetPos = " << targetPos
            << ", queryLen = " << queryLen << ", targetLen = " << targetLen
            << ", alnQuerySpan = " << (diffs.numEq + diffs.numX + diffs.numI)
            << ", alnTargetSpan = " << (diffs.numEq + diffs.numX + diffs.numD) << ", " << diffs
            << ", edit_dist = " << diffs.NumDiffs() << ", CIGAR: " << cigar.ToStdString()
            << ", label: '" << label << "'";
        throw std::runtime_error(oss.str());
    }
}

void ExtractVariantString(const std::string_view query, const std::string_view target,
                          const Data::Cigar& cigar, const bool maskHomopolymers,
                          const bool maskSimpleRepeats, const bool maskHomopolymerSNPs,
                          const bool maskHomopolymersArbitrary, std::string& retQueryVariants,
                          std::string& retTargetVariants, DiffCounts& retDiffsPerBase,
                          DiffCounts& retDiffsPerEvent)
{
    const int64_t queryLen = query.size();
    const int64_t targetLen = target.size();
    const char* queryData = query.data();
    const char* targetData = target.data();
    const int32_t numCigarOps = cigar.size();
    int64_t queryPos = 0;
    int64_t targetPos = 0;

    int32_t varStrQuerySize = 0;
    int32_t varStrTargetSize = 0;
    for (int32_t i = 0; i < numCigarOps; ++i) {
        const auto& op = cigar[i];
        if (op.Type() == Data::CigarOperationType::SEQUENCE_MISMATCH) {
            varStrQuerySize += op.Length();
            varStrTargetSize += op.Length();
        } else if (op.Type() == Data::CigarOperationType::INSERTION) {
            varStrQuerySize += op.Length();
        } else if (op.Type() == Data::CigarOperationType::DELETION) {
            varStrTargetSize += op.Length();
        }
    }

    std::string varStrQuery(varStrQuerySize, '0');
    std::string varStrTarget(varStrTargetSize, '0');
    int32_t varStrQueryPos = 0;
    int32_t varStrTargetPos = 0;
    DiffCounts diffsPerBase;
    DiffCounts diffsPerEvent;

    for (int32_t i = 0; i < numCigarOps; ++i) {
        const auto& op = cigar[i];
        const int32_t opLen = op.Length();

        if (queryPos > queryLen || targetPos > targetLen) {
            std::ostringstream oss;
            oss << "Invalid CIGAR string (global): "
                << "coordinates out of bounds! "
                << "queryPos = " << queryPos << ", targetPos = " << targetPos
                << ", offending CIGAR op: " << op.Length() << op.TypeToChar(op.Type())
                << ", queryLen = " << queryLen << ", targetLen = " << targetLen
                << ", CIGAR: " << cigar.ToStdString();
            throw std::runtime_error(oss.str());
        }

        if (op.Type() == Data::CigarOperationType::SEQUENCE_MATCH) {
            // If it's a match, just move down the sequences.
            // Sanity check.
            if ((queryPos + op.Length()) > queryLen || (targetPos + op.Length()) > targetLen) {
                std::ostringstream oss;
                oss << "Invalid CIGAR string (SEQUENCE_MATCH): "
                    << "coordinates out of bounds! "
                    << "queryPos = " << queryPos << ", targetPos = " << targetPos
                    << ", offending CIGAR op: " << op.Length() << op.TypeToChar(op.Type())
                    << ", queryLen = " << queryLen << ", targetLen = " << targetLen
                    << ", CIGAR: " << cigar.ToStdString();
                throw std::runtime_error(oss.str());
            }
            // Compute diffs.
            diffsPerBase.numEq += op.Length();
            diffsPerEvent.numEq += op.Length();
            // Move down.
            queryPos += op.Length();
            targetPos += op.Length();
        } else if (op.Type() == Data::CigarOperationType::SEQUENCE_MISMATCH) {
            // For a mismatch, include both alleles.
            // Sanity check.
            if ((queryPos + op.Length()) > queryLen || (targetPos + op.Length()) > targetLen) {
                std::ostringstream oss;
                oss << "Invalid CIGAR string (SEQUENCE_MISMATCH): "
                    << "coordinates out of bounds! "
                    << "queryPos = " << queryPos << ", targetPos = " << targetPos
                    << ", offending CIGAR op: " << op.Length() << op.TypeToChar(op.Type())
                    << ", queryLen = " << queryLen << ", targetLen = " << targetLen
                    << ", CIGAR: " << cigar.ToStdString();
                throw std::runtime_error(oss.str());
            }

            bool isMasked = false;
            if (maskHomopolymerSNPs) {
                auto IsSequenceHP = [](const char* seq, const int32_t len) {
                    int32_t baseSwitches = 0;
                    char prevBase = seq[0];
                    for (int64_t pos = 0; pos < len; ++pos) {
                        if (seq[pos] != prevBase) {
                            ++baseSwitches;
                            prevBase = seq[pos];
                            break;
                        }
                    }
                    return baseSwitches == 0;
                };
                const bool isQueryHP = IsSequenceHP(queryData + queryPos, op.Length());
                const bool isTargetHP = IsSequenceHP(targetData + targetPos, op.Length());

                // If the query variant sequence is a homopolymer, then check if it matches
                // one base before or after it:
                //  Q: TTTTT
                //     |||X|
                //  T: TTTGT
                if (isQueryHP && ((queryPos > 0 && query[queryPos - 1] == query[queryPos]) ||
                                  ((queryPos + opLen) < queryLen &&
                                   query[queryPos + opLen - 1] == query[queryPos + opLen]))) {
                    isMasked = true;
                }
                // Same for the target sequence.
                //  Q: TTTGT
                //     |||X|
                //  T: TTTTT
                if (isTargetHP && ((targetPos > 0 && target[targetPos - 1] == target[targetPos]) ||
                                   ((targetPos + opLen) < targetLen &&
                                    target[targetPos + opLen - 1] == target[targetPos + opLen]))) {
                    isMasked = true;
                }
            }

            if (isMasked) {
                for (int64_t pos = 0; pos < static_cast<int64_t>(op.Length()); ++pos) {
                    varStrTarget[varStrTargetPos + pos] = std::tolower(target[targetPos + pos]);
                    varStrQuery[varStrQueryPos + pos] = std::tolower(query[queryPos + pos]);
                }
            } else {
                for (int64_t pos = 0; pos < static_cast<int64_t>(op.Length()); ++pos) {
                    varStrTarget[varStrTargetPos + pos] = target[targetPos + pos];
                    varStrQuery[varStrQueryPos + pos] = query[queryPos + pos];
                }
                // Compute diffs.
                diffsPerBase.numX += op.Length();
                diffsPerEvent.numX += op.Length();
            }

            varStrQueryPos += op.Length();
            varStrTargetPos += op.Length();

            // Move down.
            queryPos += op.Length();
            targetPos += op.Length();
        } else if (op.Type() == Data::CigarOperationType::INSERTION) {
            // Sanity check.
            if ((queryPos + op.Length()) > queryLen) {
                std::ostringstream oss;
                oss << "Invalid CIGAR string (INSERTION): "
                    << "coordinates out of bounds! "
                    << "queryPos = " << queryPos << ", targetPos = " << targetPos
                    << ", offending CIGAR op: " << op.Length() << op.TypeToChar(op.Type())
                    << ", queryLen = " << queryLen << ", targetLen = " << targetLen
                    << ", CIGAR: " << cigar.ToStdString();
                throw std::runtime_error(oss.str());
            }

            bool isMasked = false;

            if (maskHomopolymers) {
                // All bases in the event need to be the same to be a homopolymer event.
                int32_t baseSwitches = 0;
                char prevBase = query[queryPos + 0];
                for (int64_t pos = 0; pos < static_cast<int64_t>(op.Length()); ++pos) {
                    if (query[queryPos + pos] != prevBase) {
                        ++baseSwitches;
                        prevBase = query[queryPos + pos];
                    }
                }
                // Check if the current event bases are the same as the previous/next base
                // in either query or target to call it homopolymer.
                if (baseSwitches == 0 &&
                    ((queryPos > 0 && query[queryPos - 1] == prevBase) ||
                     ((queryPos + 1) < queryLen && query[queryPos + 1] == prevBase) ||
                     (maskHomopolymersArbitrary && queryPos > 0 && (queryPos + 1) < queryLen &&
                      query[queryPos - 1] ==
                          query[queryPos + 1]) ||  // Insertion of different base into a HP.
                     (targetPos < targetLen && target[targetPos] == prevBase) ||
                     (targetPos > 0 && target[targetPos - 1] == prevBase))) {
                    isMasked = true;
                }
            }

            // Check if the indel is exactly the same as preceding or following bases in
            // either query or target.
            if (maskSimpleRepeats && isMasked == false && op.Length() > 1) {
                if (queryPos >= op.Length() &&
                    std::strncmp(queryData + queryPos - op.Length(), queryData + queryPos,
                                 op.Length()) == 0) {
                    isMasked = true;
                } else if ((queryPos + 2 * op.Length()) <= queryLen &&
                           std::strncmp(queryData + queryPos, queryData + queryPos + op.Length(),
                                        op.Length()) == 0) {
                    isMasked = true;
                } else if (targetPos >= op.Length() &&
                           std::strncmp(targetData + targetPos - op.Length(), queryData + queryPos,
                                        op.Length()) == 0) {
                    isMasked = true;
                } else if ((targetPos + op.Length()) <= targetLen &&
                           std::strncmp(queryData + queryPos, targetData + targetPos,
                                        op.Length()) == 0) {
                    // Note: using "(targetPos + op.Length()) <= targetLen" instead of "(targetPos + 2 * op.Length()) <= targetLen"
                    // because the bases don't exist in the target so we need to start at the current position.
                    isMasked = true;
                }
            }

            // Add the query (insertion) bases.
            if (isMasked) {
                for (int64_t pos = 0; pos < static_cast<int64_t>(op.Length()); ++pos) {
                    varStrQuery[varStrQueryPos + pos] = std::tolower(query[queryPos + pos]);
                }
            } else {
                for (int64_t pos = 0; pos < static_cast<int64_t>(op.Length()); ++pos) {
                    varStrQuery[varStrQueryPos + pos] = query[queryPos + pos];
                }
                // Compute diffs.
                diffsPerBase.numI += op.Length();
                ++diffsPerEvent.numI;
            }
            varStrQueryPos += op.Length();

            // Move down.
            queryPos += op.Length();
        } else if (op.Type() == Data::CigarOperationType::DELETION) {
            // Sanity check.
            if ((targetPos + op.Length()) > targetLen) {
                std::ostringstream oss;
                oss << "Invalid CIGAR string (DELETION): "
                    << "coordinates out of bounds! "
                    << "queryPos = " << queryPos << ", targetPos = " << targetPos
                    << ", offending CIGAR op: " << op.Length() << op.TypeToChar(op.Type())
                    << ", queryLen = " << queryLen << ", targetLen = " << targetLen
                    << ", CIGAR: " << cigar.ToStdString();
                throw std::runtime_error(oss.str());
            }

            bool isMasked = false;

            if (maskHomopolymers) {
                // All bases in the event need to be the same to be a homopolymer event.
                int32_t baseSwitches = 0;
                char prevBase = target[targetPos + 0];
                for (int64_t pos = 0; pos < static_cast<int64_t>(op.Length()); ++pos) {
                    if (target[targetPos + pos] != prevBase) {
                        ++baseSwitches;
                        prevBase = target[targetPos + pos];
                    }
                }
                // Check if the current event bases are the same as the previous/next base
                // in either query or target to call it homopolymer.
                if (baseSwitches == 0 &&
                    ((targetPos > 0 && target[targetPos - 1] == prevBase) ||
                     ((targetPos + 1) < targetLen && target[targetPos + 1] == prevBase) ||
                     (maskHomopolymersArbitrary && targetPos > 0 && (targetPos + 1) < targetLen &&
                      target[targetPos - 1] ==
                          target[targetPos + 1]) ||  // Insertion of different base into a HP.
                     (queryPos < queryLen && query[queryPos] == prevBase) ||
                     (queryPos > 0 && query[queryPos - 1] == prevBase))) {
                    isMasked = true;
                }
            }

            // Check if the indel is exactly the same as preceding or following bases in
            // either query or target.
            if (maskSimpleRepeats && isMasked == false && op.Length() > 1) {
                if (targetPos >= op.Length() &&
                    std::strncmp(targetData + targetPos - op.Length(), targetData + targetPos,
                                 op.Length()) == 0) {
                    isMasked = true;
                } else if ((targetPos + 2 * op.Length()) <= targetLen &&
                           std::strncmp(targetData + targetPos,
                                        targetData + targetPos + op.Length(), op.Length()) == 0) {
                    isMasked = true;
                } else if (queryPos >= op.Length() &&
                           std::strncmp(queryData + queryPos - op.Length(), targetData + targetPos,
                                        op.Length()) == 0) {
                    isMasked = true;
                } else if ((queryPos + op.Length()) <= queryLen &&
                           std::strncmp(queryData + queryPos, targetData + targetPos,
                                        op.Length()) == 0) {
                    // Note: using "(queryPos + op.Length()) <= queryLen" instead of "(queryPos + 2 * op.Length()) <= queryLen"
                    // because the bases don't exist in the query so we need to start at the current position.
                    isMasked = true;
                }
            }

            // Add the target (deletion) bases.
            if (isMasked) {
                for (int64_t pos = 0; pos < static_cast<int64_t>(op.Length()); ++pos) {
                    varStrTarget[varStrTargetPos + pos] = std::tolower(target[targetPos + pos]);
                }
            } else {
                for (int64_t pos = 0; pos < static_cast<int64_t>(op.Length()); ++pos) {
                    varStrTarget[varStrTargetPos + pos] = target[targetPos + pos];
                }
                // Compute diffs.
                diffsPerBase.numD += op.Length();
                ++diffsPerEvent.numD;
            }
            varStrTargetPos += op.Length();

            // Move down.
            targetPos += op.Length();

        } else if (op.Type() == Data::CigarOperationType::SOFT_CLIP) {
            // Sanity check.
            if ((queryPos + op.Length()) > queryLen) {
                std::ostringstream oss;
                oss << "Invalid CIGAR string (SOFT_CLIP): "
                    << "coordinates out of bounds! "
                    << "queryPos = " << queryPos << ", targetPos = " << targetPos
                    << ", offending CIGAR op: " << op.Length() << op.TypeToChar(op.Type())
                    << ", queryLen = " << queryLen << ", targetLen = " << targetLen
                    << ", CIGAR: " << cigar.ToStdString();
                throw std::runtime_error(oss.str());
            }

            // Move down.
            queryPos += op.Length();

        } else if (op.Type() == Data::CigarOperationType::REFERENCE_SKIP) {
            // Sanity check.
            if ((targetPos + op.Length()) > targetLen) {
                std::ostringstream oss;
                oss << "Invalid CIGAR string (REFERENCE_SKIP): "
                    << "coordinates out of bounds! "
                    << "queryPos = " << queryPos << ", targetPos = " << targetPos
                    << ", offending CIGAR op: " << op.Length() << op.TypeToChar(op.Type())
                    << ", queryLen = " << queryLen << ", targetLen = " << targetLen
                    << ", CIGAR: " << cigar.ToStdString();
                throw std::runtime_error(oss.str());
            }

            // Move down.
            targetPos += op.Length();

        } else if (op.Type() == Data::CigarOperationType::HARD_CLIP) {
            // Do nothing.

        } else {
            std::ostringstream oss;
            oss << "CIGAR operation '" << op.TypeToChar(op.Type())
                << "' not supported by ExtractVariantString function.";
            throw std::runtime_error(oss.str());
        }
    }

    std::swap(retQueryVariants, varStrQuery);
    std::swap(retTargetVariants, varStrTarget);
    std::swap(retDiffsPerBase, diffsPerBase);
    std::swap(retDiffsPerEvent, diffsPerEvent);
}

DiffCounts ComputeMaskedDiffCounts(const Data::Cigar& cigar, const std::string_view queryVariants,
                                   const std::string_view targetVariants,
                                   const bool throwOnPartiallyMaskedIndels)
{

    int32_t aVarPos = 0;
    int32_t bVarPos = 0;
    const int32_t aVarLen = queryVariants.size();
    const int32_t bVarLen = targetVariants.size();

    DiffCounts diffs;

    for (const auto& op : cigar) {
        const int32_t opLen = op.Length();

        if (op.Type() == Data::CigarOperationType::SEQUENCE_MATCH) {
            // Move down.
            diffs.numEq += opLen;

        } else if (op.Type() == Data::CigarOperationType::SEQUENCE_MISMATCH) {
            if ((aVarPos + opLen) > aVarLen || (bVarPos + opLen) > bVarLen) {
                std::ostringstream oss;
                oss << "Variant position out of bounds. CIGAR op: " << op.Length()
                    << op.TypeToChar(op.Type()) << ", aVarPos = " << aVarPos
                    << ", bVarPos = " << bVarPos << ", aVarLen = " << aVarLen
                    << ", bVarLen = " << bVarLen;
                throw std::runtime_error(oss.str());
            }

            for (int32_t i = 0; i < opLen; ++i) {
                const int aLower = islower(queryVariants[aVarPos]);
                const int bLower = islower(targetVariants[bVarPos]);

                // Sanity check.
                if (aLower != bLower) {
                    std::ostringstream oss;
                    oss << "Incorrect variant masking, variant is uppercase in one instance and "
                           "lowercase in the other. "
                        << ", CIGAR op: " << op.Length() << op.TypeToChar(op.Type());
                    throw std::runtime_error(oss.str());
                }
                // Skip masked variants.
                if (aLower == 0 && bLower == 0) {
                    ++diffs.numX;
                }
                ++aVarPos;
                ++bVarPos;
            }

        } else if (op.Type() == Data::CigarOperationType::INSERTION) {
            if ((aVarPos + opLen) > aVarLen) {
                std::ostringstream oss;
                oss << "Variant position out of bounds. CIGAR op: " << op.Length()
                    << op.TypeToChar(op.Type()) << ", aVarPos = " << aVarPos
                    << ", bVarPos = " << bVarPos << ", aVarLen = " << aVarLen
                    << ", bVarLen = " << bVarLen;
                throw std::runtime_error(oss.str());
            }

            int32_t numMasked = 0;
            for (int32_t i = 0; i < opLen; ++i) {
                if (islower(queryVariants[aVarPos + i])) {
                    ++numMasked;
                }
            }
            if (throwOnPartiallyMaskedIndels && numMasked > 0 && numMasked < opLen) {
                std::ostringstream oss;
                oss << "Some positions in an insertion variant are masked, but not all. CIGAR op: "
                    << op.Length() << op.TypeToChar(op.Type()) << ", aVarPos = " << aVarPos
                    << ", bVarPos = " << bVarPos << ", aVarLen = " << aVarLen
                    << ", bVarLen = " << bVarLen << ", variant = '"
                    << queryVariants.substr(aVarPos, opLen) << "'";
                throw std::runtime_error(oss.str());
            }
            diffs.numI += (opLen - numMasked);
            aVarPos += opLen;

        } else if (op.Type() == Data::CigarOperationType::DELETION) {
            if ((bVarPos + opLen) > bVarLen) {
                std::ostringstream oss;
                oss << "Variant position out of bounds. CIGAR op: " << op.Length()
                    << op.TypeToChar(op.Type()) << ", aVarPos = " << aVarPos
                    << ", bVarPos = " << bVarPos << ", aVarLen = " << aVarLen
                    << ", bVarLen = " << bVarLen;
                throw std::runtime_error(oss.str());
            }

            int32_t numMasked = 0;
            for (int32_t i = 0; i < opLen; ++i) {
                if (islower(targetVariants[bVarPos + i])) {
                    ++numMasked;
                }
            }
            if (throwOnPartiallyMaskedIndels && numMasked > 0 && numMasked < opLen) {
                std::ostringstream oss;
                oss << "Some positions in an insertion variant are masked, but not all. CIGAR op: "
                    << op.Length() << op.TypeToChar(op.Type()) << ", aVarPos = " << aVarPos
                    << ", bVarPos = " << bVarPos << ", aVarLen = " << aVarLen
                    << ", bVarLen = " << bVarLen << ", variant = '"
                    << targetVariants.substr(bVarPos, opLen) << "'";
                throw std::runtime_error(oss.str());
            }
            diffs.numD += (opLen - numMasked);
            bVarPos += opLen;

        } else if (op.Type() == Data::CigarOperationType::SOFT_CLIP) {
            // Do nothing.
        } else if (op.Type() == Data::CigarOperationType::REFERENCE_SKIP) {
            // Do nothing.
        } else if (op.Type() == Data::CigarOperationType::HARD_CLIP) {
            // Do nothing.

        } else {
            std::ostringstream oss;
            oss << "CIGAR operation '" << op.TypeToChar(op.Type())
                << "' not supported by the ComputeDiffCounts function.";
            throw std::runtime_error(oss.str());
        }
    }

    return diffs;
}

int32_t FindTargetPosFromCigar(const Data::Cigar& cigar, const int32_t queryPos)
{
    if (cigar.empty()) {
        throw std::runtime_error("Empty CIGAR given to FindTargetPosFromCigar!");
    }
    if (queryPos < 0) {
        std::ostringstream oss;
        oss << "The queryPos should be >= 0, value " << queryPos
            << " was given to FindTargetPosFromCigar.";
        throw std::runtime_error(oss.str());
    }
    int32_t currQueryPos = 0;
    int32_t currTargetPos = 0;
    for (const auto& op : cigar) {
        int32_t opLen = op.Length();
        if (op.Type() == Data::CigarOperationType::SEQUENCE_MATCH ||
            op.Type() == Data::CigarOperationType::SEQUENCE_MISMATCH ||
            op.Type() == Data::CigarOperationType::ALIGNMENT_MATCH) {
            if (queryPos < (currQueryPos + opLen)) {
                const int32_t diff = queryPos - currQueryPos;
                return currTargetPos + diff;
            }
            currQueryPos += opLen;
            currTargetPos += opLen;
        } else if (op.Type() == Data::CigarOperationType::INSERTION) {
            if (queryPos < (currQueryPos + opLen)) {
                // By convention, insertions come after an actual base.
                return currTargetPos - 1;
            }
            currQueryPos += opLen;
        } else if (op.Type() == Data::CigarOperationType::SOFT_CLIP) {
            if (queryPos < (currQueryPos + opLen)) {
                std::ostringstream oss;
                oss << "Given query position is located in a soft clipped region! queryPos = "
                    << queryPos << ", CIGAR: " << cigar.ToStdString();
                throw std::runtime_error(oss.str());
            }
            currQueryPos += opLen;
        } else if (op.Type() == Data::CigarOperationType::DELETION ||
                   op.Type() == Data::CigarOperationType::REFERENCE_SKIP) {
            // if (queryPos < (currQueryPos + op.Length())) {
            //     return currTargetPos;
            // }
            // Don't report alignment position within a deletion - wait for an actual base.
            currTargetPos += opLen;
        }
    }

    std::ostringstream oss;
    oss << "Coordinate queryPos = " << queryPos
        << " is out of bounds of the supplied CIGAR alignment: " << cigar.ToStdString();
    throw std::runtime_error(oss.str());

    return -1;
}

void NormalizeM5AlignmentInPlace(std::string& queryAln, std::string& targetAln)
{
    /*
     * This function normalizes gaps by pushing them towards the ends of the
     * query and target sequences.
     * It also takes care of mismatches, shifting gaps through them.
     * Example of what this function does.
     * TTGACACT       TTGACACT
     * ||| X|||   ->  |||X |||
     * TTG-TACT       TTGT-ACT
    */

    if (queryAln.size() != targetAln.size()) {
        std::ostringstream oss;
        oss << "Invalid input alignment strings because size differs, queryAln.size() = "
            << queryAln.size() << ", targetAln.size() = " << targetAln.size();
        throw std::runtime_error(oss.str());
    }

    const int64_t len = queryAln.size();

    // Avoid getters for speed.
    char* query = queryAln.data();
    char* target = targetAln.data();

    for (int64_t i = 0; i < (len - 1); ++i) {
        if (query[i] == '-' && target[i] == '-') {
            continue;
        } else if (target[i] == '-') {
            for (int64_t j = (i + 1); j < len; ++j) {
                const char c = target[j];
                if (c == '-') {
                    continue;
                }
                if (c == query[i] || target[j] != query[j]) {
                    target[i] = c;
                    target[j] = '-';
                }
                break;
            }
        } else if (query[i] == '-') {
            for (int64_t j = (i + 1); j < len; ++j) {
                const char c = query[j];
                if (c == '-') {
                    continue;
                }
                if (c == target[i] || target[j] != query[j]) {
                    query[i] = c;
                    query[j] = '-';
                }
                break;
            }
        }
    }
}

void ConvertCigarToM5(const std::string_view query, const std::string_view target,
                      const Data::Cigar& cigar, std::string& retQueryAln, std::string& retTargetAln)
{
    // Clear the output.
    retQueryAln.clear();
    retTargetAln.clear();

    // Sanity check.
    if (cigar.empty()) {
        return;
    }

    const int64_t queryLen = query.size();
    const int64_t targetLen = target.size();

    // Compute diffs to know how many columns we need.
    const DiffCounts diffs = CigarDiffCounts(cigar);
    const int64_t querySpan = diffs.numEq + diffs.numX + diffs.numI;
    const int64_t targetSpan = diffs.numEq + diffs.numX + diffs.numD;

    // Sanity check.
    if (querySpan != queryLen || targetSpan != targetLen) {
        std::ostringstream oss;
        oss << "Invalid CIGAR string, query or target span do not match. CIGAR: "
            << cigar.ToStdString() << ", queryLen = " << queryLen << ", targetLen = " << targetLen
            << ", querySpan = " << querySpan << ", targetSpan = " << targetSpan;
        throw std::runtime_error(oss.str());
    }

    // Preallocate space.
    retQueryAln.resize(diffs.numEq + diffs.numX + diffs.numI + diffs.numD);
    retTargetAln.resize(diffs.numEq + diffs.numX + diffs.numI + diffs.numD);

    int64_t qPos = 0;
    int64_t tPos = 0;
    int64_t alnPos = 0;

    for (auto& cigarOp : cigar) {
        const auto op = cigarOp.Type();
        const int32_t count = cigarOp.Length();

        if (op == Data::CigarOperationType::ALIGNMENT_MATCH ||
            op == Data::CigarOperationType::SEQUENCE_MATCH ||
            op == Data::CigarOperationType::SEQUENCE_MISMATCH) {
            for (int32_t opPos = 0; opPos < count; ++opPos, ++alnPos) {
                retQueryAln[alnPos] = query[qPos];
                retTargetAln[alnPos] = target[tPos];
                ++qPos;
                ++tPos;
            }
        } else if (op == Data::CigarOperationType::INSERTION ||
                   op == Data::CigarOperationType::SOFT_CLIP) {
            for (int32_t opPos = 0; opPos < count; ++opPos, ++alnPos) {
                retQueryAln[alnPos] = query[qPos];
                retTargetAln[alnPos] = '-';
                ++qPos;
            }

        } else if (op == Data::CigarOperationType::DELETION ||
                   op == Data::CigarOperationType::REFERENCE_SKIP) {
            for (int32_t opPos = 0; opPos < count; ++opPos, ++alnPos) {
                retQueryAln[alnPos] = '-';
                retTargetAln[alnPos] = target[tPos];
                ++tPos;
            }
        } else {
            throw std::runtime_error{"ERROR: Unknown/unsupported CIGAR op: " +
                                     std::to_string(cigarOp.Char())};
        }
    }
}

Data::Cigar ConvertM5ToCigar(const std::string_view queryAln, const std::string_view targetAln)
{
    if (queryAln.size() != targetAln.size()) {
        std::ostringstream oss;
        oss << "Query and target M5 strings do not match in length! queryAln.size() = "
            << queryAln.size() << ", targetAln.size() = " << targetAln.size();
        throw std::runtime_error(oss.str());
    }

    Data::Cigar cigar;

    const char* queryAlnC = queryAln.data();
    const char* targetAlnC = targetAln.data();

    for (size_t alnPos = 0; alnPos < queryAln.size(); ++alnPos) {
        Data::CigarOperationType newOp;
        if (queryAlnC[alnPos] == targetAlnC[alnPos] && queryAlnC[alnPos] != '-') {
            newOp = Data::CigarOperationType::SEQUENCE_MATCH;
        } else if (queryAlnC[alnPos] != targetAlnC[alnPos] && queryAlnC[alnPos] != '-' &&
                   targetAlnC[alnPos] != '-') {
            newOp = Data::CigarOperationType::SEQUENCE_MISMATCH;
        } else if (queryAlnC[alnPos] == '-' && targetAlnC[alnPos] != '-') {
            newOp = Data::CigarOperationType::DELETION;
        } else if (queryAlnC[alnPos] != '-' && targetAlnC[alnPos] == '-') {
            newOp = Data::CigarOperationType::INSERTION;
        } else {
            // Both are '-'.
            continue;
        }
        AppendToCigar(cigar, newOp, 1);
    }

    return cigar;
}

Data::Cigar NormalizeCigar(const std::string_view query, const std::string_view target,
                           const Data::Cigar& cigar)
{
    std::string queryAln;
    std::string targetAln;

    PacBio::Pancake::ConvertCigarToM5(query, target, cigar, queryAln, targetAln);

    NormalizeM5AlignmentInPlace(queryAln, targetAln);

    return PacBio::Pancake::ConvertM5ToCigar(queryAln, targetAln);
}

bool TrimCigar(const Data::Cigar& cigar, const int32_t windowSize, const int32_t minMatches,
               const bool clipOnFirstMatch, Data::Cigar& retTrimmedCigar, TrimmingInfo& retTrimming)
{
    // Hardcode the max window size so that we can allocate on stack.
    static const int32_t MAX_WINDOW_SIZE = 512;

    // Sanity check.
    if (windowSize >= MAX_WINDOW_SIZE) {
        std::ostringstream oss;
        oss << "Too large window size. Requested: " << windowSize
            << ", max allowed: " << MAX_WINDOW_SIZE;
        throw std::runtime_error(oss.str());
    }

    // Reset the return values.
    retTrimming = TrimmingInfo();
    retTrimmedCigar.clear();

    // Temporary storage until the end, so that we don't return partial results.
    TrimmingInfo trimInfo;

    const auto ProcessCigarOp =
        [](const Data::Cigar& _cigar, const int32_t opId, const int32_t _windowSize,
           const int32_t _minMatches, const bool _clipOnFirstMatch,
           std::array<std::pair<int32_t, int32_t>, 512>& buff, int32_t& buffStart, int32_t& buffEnd,
           int32_t& matchCount, int32_t& foundOpId, int32_t& foundOpInternalId, int32_t& posQuery,
           int32_t& posTarget) -> bool {
        /*
         * This function processes a single CIGAR operation and adds it to the circular buffer.
         * The circular buffer represents the window.
         * Every base of the CIGAR event is processed and added separately to the window, and the
         * amount of matches in the window are maintained.
         * Once the window is filled to its size, we check if it is valid or it needs to be trimmed.
         * \returns true if a valid window was found. Also, the ID of the CIGAR operation and the
         *          internal ID of that CIGAR operation are returned via parameters.
        */

        const auto& op = _cigar[opId];
        int32_t opLen = op.Length();
        int32_t buffSize = buff.size();
        for (int32_t i = 0; i < opLen; ++i) {
            const int32_t currWindowSize =
                (buffEnd >= buffStart) ? (buffEnd - buffStart) : (buffSize - buffStart + buffEnd);

            // This happens only after the window has been filled.
            if (currWindowSize >= _windowSize) {
                // Get the start operation, which will be pushed outside of the window.
                const auto& windowOpPair = buff[buffStart];
                const int32_t startOpId = windowOpPair.first;
                const int32_t startOpInternalId = windowOpPair.second;
                const auto startOpType = _cigar[startOpId].Type();
                buffStart = (buffStart + 1) % buffSize;

                // Check if we found our target window.
                if (matchCount >= _minMatches &&
                    (_clipOnFirstMatch == false ||
                     (_clipOnFirstMatch &&
                      startOpType == Data::CigarOperationType::SEQUENCE_MATCH))) {
                    foundOpId = startOpId;
                    foundOpInternalId = startOpInternalId;
                    return true;
                }

                // Move window down and maintain the match count.
                if (startOpType == Data::CigarOperationType::SEQUENCE_MATCH) {
                    // If the start operation was a match, reduce the count as it leaves the window.
                    matchCount = std::max(matchCount - 1, 0);
                    ++posQuery;
                    ++posTarget;
                } else if (startOpType == Data::CigarOperationType::SEQUENCE_MISMATCH) {
                    ++posQuery;
                    ++posTarget;
                } else if (startOpType == Data::CigarOperationType::INSERTION) {
                    ++posQuery;
                } else if (startOpType == Data::CigarOperationType::DELETION) {
                    ++posTarget;
                } else {
                    throw std::runtime_error(
                        "Unsupported CIGAR operation when trimming the alignment: '" +
                        std::string(Data::CigarOperation::TypeToChar(startOpType), 1) + "'");
                }
            }

            // Add to the window.
            buff[buffEnd].first = opId;
            buff[buffEnd].second = i;
            buffEnd = (buffEnd + 1) % buffSize;
            if (op.Type() == Data::CigarOperationType::SEQUENCE_MATCH) {
                ++matchCount;
            }
        }
        return false;
    };

    // Circular buffer.
    std::array<std::pair<int32_t, int32_t>, MAX_WINDOW_SIZE> buff;

    // Clipping information.
    PacBio::Data::CigarOperation prefixOp;
    int32_t infixOpIdStart = 0;
    PacBio::Data::CigarOperation suffixOp;
    int32_t infixOpIdEnd = 0;

    int32_t prefixOpId = 0;
    int32_t prefixOpInternalId = 0;
    int32_t suffixOpId = 0;
    int32_t suffixOpInternalId = 0;

    // Find clipping of the front part.
    {
        int32_t buffStart = 0;
        int32_t buffEnd = 0;
        int32_t matchCount = 0;
        int32_t posQuery = 0;
        int32_t posTarget = 0;
        int32_t foundOpId = 0;
        int32_t foundOpInternalId = 0;
        bool foundGoodWindow = false;

        for (int32_t opId = 0; opId < static_cast<int32_t>(cigar.size()); ++opId) {
            foundGoodWindow = ProcessCigarOp(cigar, opId, windowSize, minMatches, clipOnFirstMatch,
                                             buff, buffStart, buffEnd, matchCount, foundOpId,
                                             foundOpInternalId, posQuery, posTarget);
            if (foundGoodWindow) {
                break;
            }
        }

        // If we cannot find a good window, just return.
        // This means that we looped through the entire CIGAR string, and it was bad in its entirety.
        if (foundGoodWindow == false) {
            return false;
        }

        // The window may begin within a CIGAR operation, so handle the first CIGAR operation
        // separately, as a "prefixOp". Other operations (except the last one) will be included
        // in their entirety. These are called "infix" operations.
        const auto& foundOp = cigar[foundOpId];
        prefixOpId = foundOpId;
        prefixOpInternalId = foundOpInternalId;
        prefixOp = PacBio::Data::CigarOperation(
            foundOp.Type(), static_cast<int32_t>(foundOp.Length()) - foundOpInternalId);
        infixOpIdStart = foundOpId + 1;
        trimInfo.queryFront = posQuery;
        trimInfo.targetFront = posTarget;
    }

    // Find clipping of the back part.
    {
        int32_t buffStart = 0;
        int32_t buffEnd = 0;
        int32_t matchCount = 0;
        int32_t posQuery = 0;
        int32_t posTarget = 0;
        int32_t foundOpId = 0;
        int32_t foundOpInternalId = 0;
        bool foundGoodWindow = false;

        for (int32_t opId = (static_cast<int32_t>(cigar.size()) - 1); opId >= 0; --opId) {
            foundGoodWindow = ProcessCigarOp(cigar, opId, windowSize, minMatches, clipOnFirstMatch,
                                             buff, buffStart, buffEnd, matchCount, foundOpId,
                                             foundOpInternalId, posQuery, posTarget);
            if (foundGoodWindow) {
                break;
            }
        }

        // If we cannot find a good window, just return.
        // This means that we looped through the entire CIGAR string, and it was bad in it's entirety.
        if (foundGoodWindow == false) {
            return false;
        }

        // The last window may end within a CIGAR operation, so handle the last CIGAR operation
        // separately, as a "suffixOp". Other operations (except the first one) will be included
        // in their entirety. These are called "infix" operations.
        const auto& foundOp = cigar[foundOpId];
        suffixOpId = foundOpId;
        suffixOpInternalId = static_cast<int32_t>(foundOp.Length()) - foundOpInternalId;
        suffixOp = PacBio::Data::CigarOperation(foundOp.Type(), suffixOpInternalId);
        infixOpIdEnd = foundOpId;
        trimInfo.queryBack = posQuery;
        trimInfo.targetBack = posTarget;
    }

    // Create the new trimmed CIGAR string.
    retTrimmedCigar.clear();
    if (prefixOpId == suffixOpId) {
        // There is no infix and start and end operation are the same.
        const auto& foundOp = cigar[prefixOpId];
        auto newOp =
            PacBio::Data::CigarOperation(foundOp.Type(), suffixOpInternalId - prefixOpInternalId);
        retTrimmedCigar.emplace_back(newOp);

    } else {
        if (prefixOp.Length() > 0) {
            retTrimmedCigar.emplace_back(prefixOp);
        }
        if (infixOpIdEnd > infixOpIdStart) {
            retTrimmedCigar.insert(retTrimmedCigar.end(), cigar.begin() + infixOpIdStart,
                                   cigar.begin() + infixOpIdEnd);
        }
        if (suffixOp.Length() > 0) {
            retTrimmedCigar.emplace_back(suffixOp);
        }
    }

    // Set the clipping results.
    retTrimming = trimInfo;

    return true;
}

int32_t ScoreCigarAlignment(const Data::Cigar& cigar, const int32_t match, const int32_t mismatch,
                            const int32_t gapOpen, const int32_t gapExt)
{
    int64_t score = 0;
    for (const auto& op : cigar) {
        const int32_t count = op.Length();
        switch (op.Type()) {
            case Data::CigarOperationType::SEQUENCE_MATCH:
                // Scores are positive.
                score += match * count;
                break;
            case Data::CigarOperationType::SEQUENCE_MISMATCH:
                // Penalties are positive.
                score -= mismatch * count;
                break;
            case Data::CigarOperationType::INSERTION:
                // Penalties are positive.
                score -= (gapOpen + gapExt * (count - 1));
                break;
            case Data::CigarOperationType::DELETION:
                // Penalties are positive.
                score -= (gapOpen + gapExt * (count - 1));
                break;
            default:
                break;
        }
    }
    return score;
}

std::pair<int32_t, PacBio::Pancake::DiffCounts> ScoreCigarAlignment(
    const Data::Cigar& cigar, const int32_t match, const int32_t mismatch, const int32_t gapOpen1,
    const int32_t gapExt1, const int32_t gapOpen2, const int32_t gapExt2)
{
    DiffCounts diffs;
    int64_t score = 0;

    const int32_t longThreshold =
        (gapExt1 != gapExt2) ? ((gapOpen2 - gapOpen1) / (gapExt1 - gapExt2) - 1) : 0;

    for (const auto& op : cigar) {
        const int32_t count = op.Length();
        switch (op.Type()) {
            case Data::CigarOperationType::SEQUENCE_MATCH:
                // Scores are positive.
                score += match * count;
                diffs.numEq += count;
                break;
            case Data::CigarOperationType::SEQUENCE_MISMATCH:
                // Penalties are positive.
                score -= mismatch * count;
                diffs.numX += count;
                break;
            case Data::CigarOperationType::INSERTION:
                // Penalties are positive.
                score -= (count < longThreshold) ? (gapOpen1 + gapExt1 * (count - 1))
                                                 : (gapOpen2 + gapExt2 * (count - 1));
                diffs.numI += count;
                break;
            case Data::CigarOperationType::DELETION:
                // Penalties are positive.
                score -= (count < longThreshold) ? (gapOpen1 + gapExt1 * (count - 1))
                                                 : (gapOpen2 + gapExt2 * (count - 1));
                ;
                diffs.numD += count;
                break;
            default:
                break;
        }
    }
    return {score, diffs};
}

void MergeCigars(PacBio::Data::Cigar& dest, const PacBio::Data::Cigar& src)
{
    if (src.empty()) {
        return;
    }
    if (dest.size() > 0 && src.front().Type() == dest.back().Type()) {
        dest.back() = PacBio::Data::CigarOperation(src.front().Type(),
                                                   dest.back().Length() + src.front().Length());
    } else {
        dest.emplace_back(src.front());
    }
    dest.insert(dest.end(), src.begin() + 1, src.end());
}

std::vector<uint8_t> ComputeSimpleRepeatMask(const std::string_view seq, int32_t maxWindowSize)
{
    const int32_t seqLen = Utility::Ssize(seq);

    if (maxWindowSize <= 0) {
        return std::vector<uint8_t>(seqLen, 0);
    }
    if (maxWindowSize > 7) {
        maxWindowSize = 7;
        assert(false && "maxWindowSize is > 7");
    }

    // Bitmasks of different lengths.
    constexpr std::array<uint64_t, 8> masks = {
        0,
        (1 << 2) - 1,
        (1 << 4) - 1,
        (1 << 6) - 1,
        (1 << 8) - 1,
        (1 << 10) - 1,
        (1 << 12) - 1,
        (1 << 14) - 1,
    };

    // Compute the repeat masks.
    std::vector<uint8_t> ret(seqLen, 0);
    uint64_t window = 0;
    int32_t beginI = 0;
    for (int32_t i = 0; i < seqLen; ++i) {
        const int32_t base = seq[i];
        const uint64_t baseTwobit = BASE_TO_TWO_BIT[base];

        // Check if non-ACTG base occurred.
        if (baseTwobit > 3) {
            window = 0;
            beginI = i + 1;
            continue;
        }

        window = (window << 2) | baseTwobit;

        const int32_t distToBegin = i - beginI;
        const int32_t maxSpan = std::min((distToBegin / 2) + (distToBegin % 2), maxWindowSize);

        for (int32_t span = 1; span <= maxSpan; ++span) {
            const uint64_t prev = window >> (span * 2) & masks[span];
            const int8_t isSame = (window & masks[span]) == prev;
            const int8_t flag = (isSame << (span - 1));

            for (int32_t k = (i - span * 2 + 1); k <= i; ++k) {
                ret[k] |= flag;
            }
        }
    }

    return ret;
}

bool CheckAlignmentOutOfBand(const PacBio::Data::Cigar& cigar, const int32_t bandwidth)
{
    const int32_t fuzz = 1;
    const int32_t upperDiag = std::max(0, bandwidth - fuzz);
    const int32_t lowerDiag = -upperDiag;
    int32_t qpos = 0;
    int32_t tpos = 0;
    for (const PacBio::Data::CigarOperation& op : cigar) {
        qpos += (PacBio::Data::ConsumesQuery(op.Type()) ? op.Length() : 0);
        tpos += (PacBio::Data::ConsumesReference(op.Type()) ? op.Length() : 0);
        const int32_t diag = tpos - qpos;
        if (diag != 0 && (diag >= upperDiag || diag <= lowerDiag)) {
            return true;
        }
    }
    return false;
}

}  // namespace Pancake
}  // namespace PacBio
