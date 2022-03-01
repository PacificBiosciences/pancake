// Authors: Ivan Sovic

#include <pancake/AlignerBatchGPUEdelweiss.hpp>

#include <pancake/AlignerKSW2.hpp>
#include <pancake/AlignmentTools.hpp>
#include <pancake/util/Util.hpp>

#include <pbcopper/utility/Ssize.h>
#include <pbcopper/utility/Stopwatch.h>

#include <cstdint>
#include <iostream>
#include <sstream>
#include <string_view>

namespace PacBio {
namespace Pancake {

AlignerBatchGPUEdelweiss::AlignerBatchGPUEdelweiss(const int32_t numThreads,
                                                   const AlignmentParameters& alnParams,
                                                   int32_t minBandwidth, int32_t maxBandwidth,
                                                   const int32_t deviceId,
                                                   const int64_t maxGPUMemoryCap)
    : AlignerBatchGPUEdelweiss(nullptr, alnParams, deviceId, minBandwidth, maxBandwidth,
                               maxGPUMemoryCap)
{
    fafFallback_ = std::make_unique<Parallel::FireAndForget>(numThreads);
    faf_ = fafFallback_.get();
}

AlignerBatchGPUEdelweiss::AlignerBatchGPUEdelweiss(Parallel::FireAndForget* faf,
                                                   const AlignmentParameters& alnParams,
                                                   int32_t minBandwidth, int32_t maxBandwidth,
                                                   const int32_t deviceId,
                                                   const int64_t maxGPUMemoryCap)
    : faf_{faf}
    , fafFallback_{nullptr}
    , aligner_{deviceId, maxGPUMemoryCap}
    , deviceId_(deviceId)
    , maxGPUMemoryCap_(maxGPUMemoryCap)
    , minBandwidth_(minBandwidth)
    , maxBandwidth_(maxBandwidth)
    , alnParams_(alnParams)
{}

AlignerBatchGPUEdelweiss::~AlignerBatchGPUEdelweiss()
{
    if (fafFallback_) {
        fafFallback_->Finalize();
    }
}

void AlignerBatchGPUEdelweiss::Clear()
{
    concatQueries_.clear();
    queryStarts_.clear();
    concatTargets_.clear();
    targetStarts_.clear();
    querySpans_.clear();
    targetSpans_.clear();
}

StatusAddSequencePair AlignerBatchGPUEdelweiss::AddSequencePairForExtensionAlignment(
    const std::string_view /*qseq*/, const std::string_view /*tseq*/)
{
    std::ostringstream oss;
    oss << "The Edelweiss aligner does not currently support extension alignment.";
    throw std::runtime_error(oss.str());
    return {};
}

StatusAddSequencePair AlignerBatchGPUEdelweiss::AddSequencePairForGlobalAlignment(
    const std::string_view qseq, const std::string_view tseq)
{
    return AddSequencePairForGlobalAlignment_(qseq, tseq);
}

StatusAddSequencePair AlignerBatchGPUEdelweiss::AddSequencePairForGlobalAlignment_(
    const std::string_view qseq, const std::string_view tseq)
{
    auto AppendSeq = [](std::vector<char>& destSeq, std::vector<int64_t>& destSeqStarts,
                        std::vector<int32_t>& destSeqSpans, const std::string_view seq) {
        const int64_t beginPos = destSeq.size();
        destSeqStarts.emplace_back(Utility::Ssize(destSeq));
        destSeqSpans.emplace_back(Utility::Ssize(seq));
        if ((destSeq.size() + seq.size()) > destSeq.capacity()) {
            const int64_t newCapacity = 2 * (Utility::Ssize(destSeq) + Utility::Ssize(seq));
            destSeq.reserve(newCapacity);
        }
        destSeq.resize(destSeq.size() + seq.size());
        std::memcpy(destSeq.data() + beginPos, seq.data(), seq.size());
    };

    AppendSeq(concatQueries_, queryStarts_, querySpans_, qseq);
    AppendSeq(concatTargets_, targetStarts_, targetSpans_, tseq);

    return StatusAddSequencePair::OK;
}

void AlignerBatchGPUEdelweiss::AlignAll()
{
    assert(queryStarts_.size() == targetStarts_.size());
    assert(querySpans_.size() == targetSpans_.size());
    assert(queryStarts_.size() == querySpans_.size());

    alnResults_.clear();

    const PacBio::Edelweiss::AlignerParameters edelweissAlnParams{
        alnParams_.matchScore, alnParams_.mismatchPenalty, alnParams_.gapOpen1,
        alnParams_.gapExtend1, alnParams_.gapOpen2,        alnParams_.gapExtend2};

    // Edelweiss string view objects need to be initialized here, because adding seqs will realloc.
    std::vector<Edelweiss::Aligner::SeqPair> edelweissSeqs;
    for (size_t i = 0; i < queryStarts_.size(); ++i) {
        const char* queryPtr = concatQueries_.data() + queryStarts_[i];
        const size_t querySpan = querySpans_[i];
        const char* targetPtr = concatTargets_.data() + targetStarts_[i];
        const size_t targetSpan = targetSpans_[i];
        edelweissSeqs.emplace_back(std::string_view(targetPtr, targetSpan),
                                   std::string_view(queryPtr, querySpan), -1);
    }

    // Edelweiss::Aligner aligner{deviceId_, maxGPUMemoryCap_};
    const std::vector<Edelweiss::Aligner::OutputAlignment> edelResults =
        aligner_.Align(edelweissSeqs, minBandwidth_, maxBandwidth_, edelweissAlnParams);

    PopulateResults_(querySpans_, targetSpans_, alnParams_, edelResults, alnResults_, faf_);
}

AlignmentResult AlignerBatchGPUEdelweiss::ConvertEdelweissResultsToPancake_(
    const Edelweiss::Aligner::OutputAlignment& edelweissResult, int32_t queryLen, int32_t targetLen,
    const AlignmentParameters& alnParams)
{
    AlignmentResult ret;

    const auto [score, diffs] = ScoreCigarAlignment(
        edelweissResult.cigar, alnParams.matchScore, alnParams.mismatchPenalty, alnParams.gapOpen1,
        alnParams.gapExtend1, alnParams.gapOpen2, alnParams.gapExtend2);
    const bool isSuboptimal =
        CheckAlignmentOutOfBand(edelweissResult.cigar, alnParams.alignBandwidth);
    ret.valid = (edelweissResult.status == Edelweiss::Aligner::OutputAlignment::Status::SUCCESS) &&
                (!isSuboptimal);
    ret.lastQueryPos = queryLen;
    ret.lastTargetPos = targetLen;
    ret.maxQueryPos = queryLen - 1;
    ret.maxTargetPos = targetLen - 1;
    ret.score = score;
    ret.diffs = diffs;
    ret.maxScore = ret.score;
    ret.zdropped = false;
    if (ret.valid) {
        ret.cigar = edelweissResult.cigar;
    }
    return ret;
}

void AlignerBatchGPUEdelweiss::PopulateResults_(
    const std::vector<int32_t>& querySpans, const std::vector<int32_t>& targetSpans,
    const AlignmentParameters& alnParams,
    const std::vector<Edelweiss::Aligner::OutputAlignment>& edelweissResults,
    std::vector<AlignmentResult>& alnResults, Parallel::FireAndForget* faf)
{
    alnResults.clear();
    alnResults.resize(edelweissResults.size());

    // Determine how many records should land in each thread, spread roughly evenly.
    const int32_t numRecords = edelweissResults.size();
    const std::vector<std::pair<int32_t, int32_t>> jobsPerThread =
        PacBio::Pancake::DistributeJobLoad<int32_t>(faf ? faf->NumThreads() : 1, numRecords);

    // Run the mapping in parallel.
    const auto Submit = [&querySpans, &targetSpans, &alnParams, &alnResults, &edelweissResults,
                         &jobsPerThread](int32_t jobId) {
        const int32_t jobStart = jobsPerThread[jobId].first;
        const int32_t jobEnd = jobsPerThread[jobId].second;

        for (int32_t i = jobStart; i < jobEnd; ++i) {
            alnResults[i] = ConvertEdelweissResultsToPancake_(edelweissResults[i], querySpans[i],
                                                              targetSpans[i], alnParams);
        }
    };
    Parallel::Dispatch(faf, jobsPerThread.size(), Submit);
}

}  // namespace Pancake
}  // namespace PacBio
