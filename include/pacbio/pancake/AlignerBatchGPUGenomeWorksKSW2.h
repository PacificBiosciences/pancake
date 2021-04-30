// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNER_BATCH_GPU_KSW2_H
#define PANCAKE_ALIGNER_BATCH_GPU_KSW2_H

#include <pacbio/pancake/AlignerBatchBase.h>
#include <pacbio/pancake/AlignerFactory.h>
#include <pbbam/Cigar.h>
#include <claraparabricks/genomeworks/cudaalignerdoublegap/alignerdoublegap.cuh>
#include <claraparabricks/genomeworks/utils/cudautils.hpp>
#include <cstdint>
#include <memory>
#include <vector>

namespace PacBio {
namespace Pancake {

class AlignerBatchGPUGenomeWorksKSW2 : public AlignerBatchBase
{
public:
    AlignerBatchGPUGenomeWorksKSW2(const AlignmentParameters& alnParams, uint32_t deviceId,
                                   int64_t maxGPUMemoryCap);

    ~AlignerBatchGPUGenomeWorksKSW2() override;

    /*
     * Clears the internal states (sequences for alignment and results).
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
    StatusAddSequencePair AddSequencePairForGlobalAlignment(const char* query, int32_t queryLen,
                                                            const char* target,
                                                            int32_t targetLen) override;

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
    StatusAddSequencePair AddSequencePairForExtensionAlignment(const char* query, int32_t queryLen,
                                                               const char* target,
                                                               int32_t targetLen) override;

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

    size_t BatchSize() const override { return querySpans_.size(); }

    void ResetMaxBandwidth(int32_t maxBandwidth) { alnParams_.alignBandwidth = maxBandwidth; }

private:
    AlignmentParameters alnParams_;
    claraparabricks::genomeworks::CudaStream cudaStream_;
    claraparabricks::genomeworks::DefaultDeviceAllocator deviceAllocator_;

    int64_t currentlyConsumedMem_;
    int64_t maxGPUMemoryCap_;

    std::vector<uint8_t> concatQueries_;
    std::vector<int64_t> queryStarts_;
    std::vector<uint8_t> concatTargets_;
    std::vector<int64_t> targetStarts_;
    std::vector<gw_ksw_extz_t> gwResults_;

    std::vector<int32_t> querySpans_;
    std::vector<int32_t> targetSpans_;
    std::vector<AlignmentResult> alnResults_;

    StatusAddSequencePair AddSequencePairForGlobalAlignment_(const char* query, int32_t queryLen,
                                                             const char* target, int32_t targetLen);
};
}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ALIGNER_BATCH_GPU_H
