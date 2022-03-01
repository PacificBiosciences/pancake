// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNER_BATCH_GPU_EDELWEISS_HPP
#define PANCAKE_ALIGNER_BATCH_GPU_EDELWEISS_HPP

#include <pancake/AlignerBatchBase.hpp>
#include <pancake/AlignerFactory.hpp>

#include <pbbam/Cigar.h>
#include <pbcopper/parallel/FireAndForget.h>
#include <pbedelweiss/Align.hpp>

#include <cstdint>
#include <cstring>
#include <memory>
#include <vector>

namespace PacBio {
namespace Pancake {

class AlignerBatchGPUEdelweiss : public AlignerBatchBase
{
public:
    AlignerBatchGPUEdelweiss(int32_t numThreads, const AlignmentParameters& alnParams,
                             int32_t minBandwidth, int32_t maxBandwidth, int32_t deviceId,
                             int64_t maxGPUMemoryCap);

    AlignerBatchGPUEdelweiss(Parallel::FireAndForget* faf, const AlignmentParameters& alnParams,
                             int32_t minBandwidth, int32_t maxBandwidth, int32_t deviceId,
                             int64_t maxGPUMemoryCap);

    ~AlignerBatchGPUEdelweiss() override;

    /*
     * Clears the internal states (sequence data and results).
    */
    void Clear() override;

    /*
     * Adds a single sequence pair for global alignment to the internal state (modifies the state),
     * but does not align the sequences right away. Alignment will be performed only when the
     * AlignAll function is called.
     * The data will be copied internally, so the query and target pointers do not have to
     * remain valid after this function has been called.
     * The return value is either StatusAddSequencePair::OK if the sequences were added successfully,
     * or another value describing the reason for rejecting the sequence pair.
     * Calling this function will clear the internal alignment results.
    */
    StatusAddSequencePair AddSequencePairForGlobalAlignment(std::string_view qseq,
                                                            std::string_view tseq) override;

    /*
     * Adds a single sequence pair for extension alignment to the internal state (modifies the state),
     * but does not align the sequences right away. Alignment will be performed only when the
     * AlignAll function is called.
     * The data will be copied internally, so the query and target pointers do not have to
     * remain valid after this function has been called.
     * The return value is either StatusAddSequencePair::OK if the sequences were added successfully,
     * or another value describing the reason for rejecting the sequence pair.
     * Calling this function will clear the internal alignment results.
    */
    StatusAddSequencePair AddSequencePairForExtensionAlignment(std::string_view qseq,
                                                               std::string_view tseq) override;

    /*
     * Aligns all the sequence pairs added to the aligner, in parallel.
     * This function modifies the internal state, because the alignment results are stored internally.
    */
    void AlignAll() override;

    /*
     * Const getter for the alignment results. Returns an empty vector if the sequences have not
     * been aligned yet.
    */
    const std::vector<AlignmentResult>& GetAlnResults() const override { return alnResults_; }

    /*
     * Non-const getter for the alignment results. Returns an empty vector if the sequences have not
     * been aligned yet.
    */
    std::vector<AlignmentResult>& GetAlnResults() override { return alnResults_; }

    /**
     * \brief Current number of sequence pairs in the batch.
    */
    size_t BatchSize() const override { return querySpans_.size(); }

    /*
     * Reset the maximum bandwidth cap for alignment.
    */
    void ResetMaxBandwidth(int32_t maxBandwidth) { maxBandwidth_ = maxBandwidth; }

private:
    Parallel::FireAndForget* faf_;
    std::unique_ptr<Parallel::FireAndForget> fafFallback_;
    Edelweiss::Aligner aligner_;
    int32_t deviceId_ = 0;
    int64_t maxGPUMemoryCap_ = 0;
    int32_t minBandwidth_ = 64;
    int32_t maxBandwidth_ = 256;
    AlignmentParameters alnParams_;
    std::vector<char> concatQueries_;
    std::vector<int64_t> queryStarts_;
    std::vector<char> concatTargets_;
    std::vector<int64_t> targetStarts_;
    std::vector<int32_t> querySpans_;
    std::vector<int32_t> targetSpans_;
    std::vector<AlignmentResult> alnResults_;

    StatusAddSequencePair AddSequencePairForGlobalAlignment_(std::string_view qseq,
                                                             std::string_view tseq);

    static AlignmentResult ConvertEdelweissResultsToPancake_(
        const Edelweiss::Aligner::OutputAlignment& edelweissResult, int32_t queryLen,
        int32_t targetLen, const AlignmentParameters& alnParams);

    static void PopulateResults_(
        const std::vector<int32_t>& querySpans, const std::vector<int32_t>& targetSpans,
        const AlignmentParameters& alnParams,
        const std::vector<Edelweiss::Aligner::OutputAlignment>& edelweissResults,
        std::vector<AlignmentResult>& alnResults, Parallel::FireAndForget* faf);
};
}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ALIGNER_BATCH_GPU_EDELWEISS_HPP
