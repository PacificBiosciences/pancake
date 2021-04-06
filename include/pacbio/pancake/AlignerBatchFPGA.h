// Copyright (c) 2019, Pacific Biosciences of California, Inc.
// All rights reserved.
// See LICENSE.txt.
//
// Contributions from NVIDIA are Copyright (c) 2021, NVIDIA Corporation.
// All rights reserved.
// SPDX-License-Identifier: BSD-3-Clause-Clear
//
// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNER_BATCH_FPGA_H
#define PANCAKE_ALIGNER_BATCH_FPGA_H

#include <pacbio/pancake/AlignerBase.h>
#include <pacbio/pancake/AlignerFactory.h>
#include <pacbio/pancake/AlignmentResult.h>
#include <pbbam/Cigar.h>
#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

namespace PacBio {
namespace Pancake {

class AlignerBatchFPGA
{
public:
    AlignerBatchFPGA(const AlignmentParameters& alnParams, uint32_t deviceId);

    ~AlignerBatchFPGA();

    /*
     * Clears the internal states (sequences for alignment and results).
    */
    void Clear();

    void AddSequencePair(const char* query, int32_t queryLen, const char* target,
                         int32_t targetLen);

    /*
     * Aligns all the sequence pairs added to the aligner, in parallel.
     * This function modifies the internal state, because the alignment results are stored internally.
    */
    std::pair<int64_t, int64_t> AlignAll(const bool extension);

    /*
     * Const getter for the alignment results. Returns an empty vector if the sequences have not
     * been aligned yet.
    */
    const std::vector<AlignmentResult>& GetAlnResults() const { return alnResults_; }

    /*
     * Non-const getter for the alignment results. Returns an empty vector if the sequences have not
     * been aligned yet.
    */
    std::vector<AlignmentResult>& GetAlnResults() { return alnResults_; }

    size_t BatchSize() const { return fpgaBatch_.size(); }

    struct AlignerBatchFPGAHostBuffers;

private:
    AlignmentParameters alnParams_;
    uint32_t deviceId_;

    std::vector<std::pair<std::vector<uint8_t>, std::vector<uint8_t>>> fpgaBatch_;
    std::vector<AlignmentResult> alnResults_;
};
}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ALIGNER_BATCH_FPGA_H
