// Authors: Ivan Sovic

#include <pancake/util/RunLengthEncoding.hpp>

namespace PacBio {
namespace Pancake {

void RunLengthEncoding(const std::string& seq, std::string& encodedSeq,
                       std::vector<int32_t>& runLengths)
{
    RunLengthEncoding(seq.c_str(), seq.size(), encodedSeq, runLengths);
}

void RunLengthEncoding(const char* seq, int64_t seqLen, std::string& encodedSeq,
                       std::vector<int32_t>& runLengths)
{
    encodedSeq = std::string(seq, seqLen);
    runLengths.clear();
    if (seq == NULL || seqLen == 0) return;

    int64_t outSeqLen = RunLengthEncoding(encodedSeq.data(), seqLen, runLengths);
    encodedSeq[outSeqLen] = '\n';
    encodedSeq.resize(outSeqLen);
}

int64_t RunLengthEncoding(char* seq, int64_t seqLen, std::vector<int32_t>& runLengths)
{
    runLengths.clear();
    if (seq == NULL || seqLen == 0) return 0;

    runLengths.resize(seqLen);
    int64_t prev = 0;
    int64_t count = 0;

    for (int64_t i = 0; i < seqLen; ++i) {
        if (seq[prev] != seq[i]) {
            ++prev;
            count = 0;
            seq[prev] = seq[i];
        }
        ++count;
        runLengths[prev] = count;
    }

    ++prev;
    runLengths.resize(prev);

    return prev;
}

void RunLengthEncoding(const std::string& seq, std::string& encodedSeq,
                       std::vector<int32_t>& seqToHPCCoords, std::vector<int32_t>& hpcToSeqCoords)
{
    RunLengthEncoding(seq.c_str(), seq.size(), encodedSeq, seqToHPCCoords, hpcToSeqCoords);
}

void RunLengthEncoding(const char* seq, int64_t seqLen, std::string& encodedSeq,
                       std::vector<int32_t>& seqToHPCCoords, std::vector<int32_t>& hpcToSeqCoords)
{
    encodedSeq = std::string(seq, seqLen);
    seqToHPCCoords.clear();
    hpcToSeqCoords.clear();
    if (seq == NULL || seqLen == 0) return;

    int64_t outSeqLen =
        RunLengthEncoding(encodedSeq.data(), seqLen, seqToHPCCoords, hpcToSeqCoords);
    encodedSeq[outSeqLen] = '\n';
    encodedSeq.resize(outSeqLen);
}

int64_t RunLengthEncoding(char* seq, int64_t seqLen, std::vector<int32_t>& seqToHPCCoords,
                          std::vector<int32_t>& hpcToSeqCoords)
{
    seqToHPCCoords.clear();
    hpcToSeqCoords.clear();
    if (seq == NULL || seqLen == 0) return 0;

    seqToHPCCoords.resize(seqLen, 0);
    hpcToSeqCoords.resize(seqLen, 0);
    int64_t prev = 0;
    int64_t count = 0;

    for (int64_t i = 0; i < seqLen; ++i) {
        if (seq[prev] != seq[i]) {
            ++prev;
            count = 0;
            seq[prev] = seq[i];
        }
        ++count;
        seqToHPCCoords[i] = prev;
        hpcToSeqCoords[prev] = i;
    }

    ++prev;
    hpcToSeqCoords.resize(prev);

    return prev;
}

int64_t RunLengthEncoding(const char* seq, int64_t seqLen, std::vector<char>& destHPC,
                          int32_t& hpcLen, std::vector<int32_t>& seqToHPCCoords,
                          std::vector<int32_t>& hpcToSeqCoords)
{
    if (seq == NULL || seqLen == 0) return 0;

    // All buffers should have the same length for speed.
    // If not, clear them so that the below if can only compare one of them.
    if (seqToHPCCoords.size() != hpcToSeqCoords.size() || destHPC.size() != seqToHPCCoords.size() ||
        destHPC.size() != hpcToSeqCoords.size()) {
        destHPC.clear();
        seqToHPCCoords.clear();
        hpcToSeqCoords.clear();
    }
    if (seqLen > static_cast<int64_t>(destHPC.size())) {
        // Reduce the amount of resizing. Only increase the size of the buffer when required.
        destHPC.resize(seqLen);
        seqToHPCCoords.resize(seqLen);
        hpcToSeqCoords.resize(seqLen);
    }

    int64_t prev = 0;
    int64_t count = 0;
    char* destRaw = destHPC.data();

    for (int64_t i = 0; i < seqLen; ++i) {
        if (destRaw[prev] != seq[i]) {
            ++prev;
            count = 0;
            destRaw[prev] = seq[i];
        }
        ++count;
        seqToHPCCoords[i] = prev;
        hpcToSeqCoords[prev] = i;
    }

    hpcLen = prev;

    return prev;
}

}  // namespace Pancake
}  // namespace PacBio
