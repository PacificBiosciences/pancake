// Copyright (c) 2019, Pacific Biosciences of California, Inc.
// All rights reserved.
// See LICENSE.txt.
//
// Contributions from NVIDIA are Copyright (c) 2021, NVIDIA Corporation.
// All rights reserved.
// SPDX-License-Identifier: BSD-3-Clause-Clear
//
// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNER_BATCH_GPU_HPP
#define PANCAKE_ALIGNER_BATCH_GPU_HPP

#include <pancake/AlignerBatchBase.hpp>
#include <pancake/AlignerFactory.hpp>

#include <pbcopper/data/Cigar.h>
#include <claraparabricks/genomeworks/cudaaligner/aligner.hpp>
#include <claraparabricks/genomeworks/cudaaligner/alignment.hpp>
#include <claraparabricks/genomeworks/cudaaligner/cudaaligner.hpp>
#include <claraparabricks/genomeworks/utils/cudautils.hpp>

#include <cstdint>
#include <memory>
#include <vector>

namespace PacBio {
namespace Pancake {

int64_t ComputeMaxGPUMemory(int64_t cudaalignerBatches, double maxGPUMemoryFraction);

class AlignerBatchGPU : public AlignerBatchBase
{
public:
    AlignerBatchGPU(const AlignmentParameters& alnParams, uint32_t maxBandwidth, uint32_t deviceId,
                    int64_t maxGPUMemoryCap);

    ~AlignerBatchGPU() override;

    /*
     * Clears the internal states (sequences for alignment and results).
    */
    void Clear() override;

    /*
     * Reset the cuda aligner to tthe provided maxiumum bandwidth.
    */
    void ResetMaxBandwidth(int32_t maxBandwidth);

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

    size_t BatchSize() const override { return querySpans_.size(); }

    struct AlignerBatchGPUHostBuffers;

private:
    AlignmentParameters alnParams_;

    // Make sure gpuStream is the first CUDA-related data element,
    // such that it gets destroyed after the other CUDA-related
    // data elements which depend on it.
    claraparabricks::genomeworks::CudaStream gpuStream_;
    std::unique_ptr<claraparabricks::genomeworks::cudaaligner::FixedBandAligner> aligner_;
    std::unique_ptr<AlignerBatchGPUHostBuffers> gpuHostBuffers_;

    std::vector<int32_t> querySpans_;
    std::vector<int32_t> targetSpans_;
    std::vector<AlignmentResult> alnResults_;
};
}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ALIGNER_BATCH_GPU_HPP
