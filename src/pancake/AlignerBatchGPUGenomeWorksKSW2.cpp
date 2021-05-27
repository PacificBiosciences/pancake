// Authors: Ivan Sovic

#include <pacbio/alignment/AlignmentTools.h>
#include <pacbio/pancake/AlignerBatchGPUGenomeWorksKSW2.h>
#include <pacbio/pancake/AlignerKSW2.h>
#include <pacbio/util/Util.h>
#include <pbcopper/utility/Stopwatch.h>
#include <cstdint>
#include <iostream>
#include <sstream>

namespace PacBio {
namespace Pancake {

static void ksw_gen_simple_mat(int m, int8_t* mat, int8_t a, int8_t b, int8_t sc_ambi)
{
    int i, j;
    a = a < 0 ? -a : a;
    b = b > 0 ? -b : b;
    sc_ambi = sc_ambi > 0 ? -sc_ambi : sc_ambi;
    for (i = 0; i < m - 1; ++i) {
        for (j = 0; j < m - 1; ++j)
            mat[i * m + j] = i == j ? a : b;
        mat[i * m + m - 1] = sc_ambi;
    }
    for (j = 0; j < m; ++j)
        mat[(m - 1) * m + j] = sc_ambi;
}

AlignerBatchGPUGenomeWorksKSW2::AlignerBatchGPUGenomeWorksKSW2(const int32_t numThreads,
                                                               const AlignmentParameters& alnParams,
                                                               const uint32_t deviceId,
                                                               const int64_t maxGPUMemoryCap)
    : AlignerBatchGPUGenomeWorksKSW2(nullptr, alnParams, deviceId, maxGPUMemoryCap)
{
    fafFallback_ = std::make_unique<Parallel::FireAndForget>(numThreads);
    faf_ = fafFallback_.get();
}

AlignerBatchGPUGenomeWorksKSW2::AlignerBatchGPUGenomeWorksKSW2(Parallel::FireAndForget* faf,
                                                               const AlignmentParameters& alnParams,
                                                               const uint32_t /*deviceId*/,
                                                               const int64_t maxGPUMemoryCap)
    : faf_{faf}
    , fafFallback_{nullptr}
    , alnParams_(alnParams)
    , cudaStream_(claraparabricks::genomeworks::make_cuda_stream())
    , currentlyConsumedMem_(0)
    , maxGPUMemoryCap_(maxGPUMemoryCap)
{
    deviceAllocator_ = claraparabricks::genomeworks::create_default_device_allocator(
        claraparabricks::genomeworks::cudautils::find_largest_contiguous_device_memory_section(),
        cudaStream_.get());
}

AlignerBatchGPUGenomeWorksKSW2::~AlignerBatchGPUGenomeWorksKSW2()
{
    if (fafFallback_) {
        fafFallback_->Finalize();
    }
}

void AlignerBatchGPUGenomeWorksKSW2::Clear()
{
    currentlyConsumedMem_ = 0;

    concatQueries_.clear();
    queryStarts_.clear();
    concatTargets_.clear();
    targetStarts_.clear();
    gwResults_.clear();

    alnResults_.clear();
    querySpans_.clear();
    targetSpans_.clear();
}

StatusAddSequencePair AlignerBatchGPUGenomeWorksKSW2::AddSequencePairForExtensionAlignment(
    const char* /*query*/, int32_t /*queryLen*/, const char* /*target*/, int32_t /*targetLen*/)
{
    std::ostringstream oss;
    oss << "The GenomeWorks Cudaaligner does not support extension alignment at this "
           "point.";
    throw std::runtime_error(oss.str());
    return {};
}

StatusAddSequencePair AlignerBatchGPUGenomeWorksKSW2::AddSequencePairForGlobalAlignment(
    const char* query, int32_t queryLen, const char* target, int32_t targetLen)
{
    if (queryLen == 0 || targetLen == 0) {
        // Add a dummy sequence pair to keep the place for the alignment.
        StatusAddSequencePair rv = AddSequencePairForGlobalAlignment_("A", 1, "A", 1);
        // The spans have been added by the function call one line up.
        querySpans_.back() = queryLen;
        targetSpans_.back() = targetLen;
        return rv;
    }

    return AddSequencePairForGlobalAlignment_(query, queryLen, target, targetLen);
}

StatusAddSequencePair AlignerBatchGPUGenomeWorksKSW2::AddSequencePairForGlobalAlignment_(
    const char* query, int32_t queryLen, const char* target, int32_t targetLen)
{
    if (queryLen < 0 || targetLen < 0) {
        return StatusAddSequencePair::SEQUENCE_LEN_BELOW_ZERO;
    }

    int64_t requiredSize = 0;
    // Size required for the Z-matrix:
    // requiredSize = (2 * std::min(queryLen, targetLen) - 1 +
    //                 std::min(alnParams_.alignBandwidth,
    //                          std::max(queryLen, targetLen) - std::min(queryLen, targetLen)) *
    //                     (alnParams_.alignBandwidth + 1));
    // requiredSize = (2 * std::min(queryLen, targetLen) - 1 +
    //                 std::min(alnParams_.alignBandwidth,
    //                          std::max(queryLen, targetLen)) *
    //                     (alnParams_.alignBandwidth + 1)) * 6;
    requiredSize = 2 * std::min(queryLen, targetLen) - 1 +
                   (alnParams_.alignBandwidth + 1) * std::max(queryLen, targetLen) * 3;

    // Size required for the query, the target, plus extra.
    requiredSize += queryLen + targetLen + (queryLen + 1) * 12 + targetLen * 20 + 40;

    // Check if it fits.
    // std::cerr << "currentlyConsumedMem_ + requiredSize = " << (currentlyConsumedMem_ + requiredSize) << ", maxGPUMemoryCap_ = " << maxGPUMemoryCap_ << " (numPairs = " << querySpans_.size() << ", currentlyConsumedMem_ = " << currentlyConsumedMem_ << ", requiredSize = " << requiredSize << ")\n";
    if ((currentlyConsumedMem_ + requiredSize) > maxGPUMemoryCap_) {
        //        std::cerr << "Filled up."
        //                  << " currentlyConsumedMem_ + requiredSize = "
        //                  << (currentlyConsumedMem_ + requiredSize)
        //                  << ", maxGPUMemoryCap_ = " << maxGPUMemoryCap_
        //                  << " (numPairs = " << querySpans_.size()
        //                  << ", currentlyConsumedMem_ = " << currentlyConsumedMem_
        //                  << ", requiredSize = " << requiredSize << ")\n";
        // std::cerr << "------------------------------\n";
        return StatusAddSequencePair::EXCEEDED_MAX_ALIGNMENTS;
    }
    currentlyConsumedMem_ += requiredSize;

    const std::vector<uint8_t> qseqInt =
        AlignerKSW2::ConvertSeqAlphabet(query, queryLen, &BaseToTwobit[0]);
    const std::vector<uint8_t> tseqInt =
        AlignerKSW2::ConvertSeqAlphabet(target, targetLen, &BaseToTwobit[0]);

    queryStarts_.emplace_back(concatQueries_.size());
    querySpans_.emplace_back(queryLen);
    concatQueries_.insert(concatQueries_.end(), qseqInt.begin(), qseqInt.end());

    targetStarts_.emplace_back(concatTargets_.size());
    targetSpans_.emplace_back(targetLen);
    concatTargets_.insert(concatTargets_.end(), tseqInt.begin(), tseqInt.end());

    return StatusAddSequencePair::OK;
}

AlignmentResult ConvertKSW2ResultsToPancake(const gw_ksw_extz_t& ez, const uint8_t* qseqInt,
                                            int32_t qseqIntLen, const uint8_t* tseqInt,
                                            int32_t tseqIntLen, int32_t bandwidth)
{
    AlignmentResult ret;

    PacBio::Data::Cigar currCigar;
    int32_t qAlnLen = 0;
    int32_t tAlnLen = 0;
    AlignerKSW2::ConvertMinimap2CigarToPbbam(ez.cigar, ez.n_cigar, qseqInt, qseqIntLen, tseqInt,
                                             tseqIntLen, currCigar, qAlnLen, tAlnLen, ret.diffs);
    const bool isSuboptimal = CheckAlignmentOutOfBand(currCigar, bandwidth);
    ret.valid = !isSuboptimal && !ez.zdropped && (ez.n_cigar != 0) && (qAlnLen == qseqIntLen) &&
                (tAlnLen == tseqIntLen);
    ret.lastQueryPos = qseqIntLen;
    ret.lastTargetPos = tseqIntLen;
    ret.maxQueryPos = ez.max_q;
    ret.maxTargetPos = ez.max_t;
    ret.score = ez.score;
    ret.maxScore = ez.max;
    ret.zdropped = ez.zdropped;
    if (ret.valid) {
        ret.cigar = std::move(currCigar);
    }

    return ret;
}

std::pair<int64_t, int64_t> AlignerBatchGPUGenomeWorksKSW2::AlignAll()
{
    int64_t cpuTime = 0;
    int64_t gpuTime = 0;
    PacBio::Utility::Stopwatch timer;
    alnResults_.clear();
    cpuTime += timer.ElapsedNanoseconds();
    timer.Reset();

    if (querySpans_.empty()) {
        return {};
    }

    const int32_t scoringMatrixSize = 5;
    std::vector<int8_t> scoringMatrix(25);
    ksw_gen_simple_mat(5, &scoringMatrix[0], alnParams_.matchScore, alnParams_.mismatchPenalty, 1);

    gwResults_.clear();
    gwResults_.resize(queryStarts_.size(), gw_ksw_extz_t());

    // GW_NVTX_RANGE(profiler, "function_call");
    // Intentionally turn off zdrop and endBonus for global alignment.
    claraparabricks::genomeworks::cudaalignerdoublegap::generate_matrix(
        deviceAllocator_, concatQueries_, queryStarts_, concatTargets_, targetStarts_,
        scoringMatrixSize, scoringMatrix, alnParams_.gapOpen1, alnParams_.gapExtend1,
        alnParams_.gapOpen2, alnParams_.gapExtend2, alnParams_.alignBandwidth, -1, 0,
        KSW_EZ_APPROX_MAX, gwResults_, cudaStream_.get());

    gpuTime += timer.ElapsedNanoseconds();
    timer.Reset();

    assert(querySpans_.size() == gwResults_.size());

    PopulateResults_(concatQueries_, queryStarts_, querySpans_, concatTargets_, targetStarts_,
                     targetSpans_, alnParams_, gwResults_, alnResults_, faf_);

    cpuTime += timer.ElapsedNanoseconds();
    return std::make_pair(cpuTime, gpuTime);
}

void AlignerBatchGPUGenomeWorksKSW2::PopulateResults_(
    const std::vector<uint8_t>& concatQueries, const std::vector<int64_t>& queryStarts,
    const std::vector<int32_t>& querySpans, const std::vector<uint8_t>& concatTargets,
    const std::vector<int64_t>& targetStarts, const std::vector<int32_t>& targetSpans,
    const AlignmentParameters& alnParams, const std::vector<gw_ksw_extz_t>& gwResults,
    std::vector<AlignmentResult>& alnResults, Parallel::FireAndForget* faf)
{
    alnResults.clear();
    alnResults.resize(gwResults.size());

    // Determine how many records should land in each thread, spread roughly evenly.
    const int32_t numRecords = gwResults.size();
    const std::vector<std::pair<int32_t, int32_t>> jobsPerThread =
        PacBio::Pancake::DistributeJobLoad<int32_t>(faf ? faf->NumThreads() : 1, numRecords);

    // Run the mapping in parallel.
    const auto Submit = [&concatQueries, &queryStarts, &querySpans, &concatTargets, &targetStarts,
                         &targetSpans, &alnParams, &alnResults, &gwResults,
                         &jobsPerThread](int32_t jobId) {
        const int32_t jobStart = jobsPerThread[jobId].first;
        const int32_t jobEnd = jobsPerThread[jobId].second;

        for (int32_t i = jobStart; i < jobEnd; ++i) {
            // Handle the edge cases which the aligner does not handle properly.
            if (querySpans[i] == 0 || targetSpans[i] == 0) {
                alnResults[i] = EdgeCaseAlignmentResult(
                    querySpans[i], targetSpans[i], alnParams.matchScore, alnParams.mismatchPenalty,
                    alnParams.gapOpen1, alnParams.gapExtend1);
                continue;
            }

            const uint8_t* qseqInt = concatQueries.data() + queryStarts[i];
            const uint8_t* tseqInt = concatTargets.data() + targetStarts[i];
            alnResults[i] =
                ConvertKSW2ResultsToPancake(gwResults[i], qseqInt, querySpans[i], tseqInt,
                                            targetSpans[i], alnParams.alignBandwidth);
        }
    };
    Parallel::Dispatch(faf, jobsPerThread.size(), Submit);
}

}  // namespace Pancake
}  // namespace PacBio
