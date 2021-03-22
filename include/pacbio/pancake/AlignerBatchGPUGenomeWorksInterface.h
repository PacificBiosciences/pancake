// Copyright (c) 2019, Pacific Biosciences of California, Inc.
// All rights reserved.
// See LICENSE.txt.
//
// PBCopper is Copyright (c) 2016-2018, Pacific Biosciences of California, Inc.
// see respective LICENSE.txt
//
// Contributions from NVIDIA are Copyright (c) 2021, NVIDIA Corporation.
// All rights reserved.
// SPDX-License-Identifier: BSD-3-Clause-Clear

#ifndef PANCAKE_ALIGNER_BATCH_GPU_GENOMEWORKS_INTERFACE_H
#define PANCAKE_ALIGNER_BATCH_GPU_GENOMEWORKS_INTERFACE_H

#include <cuda_runtime_api.h>
#include <pbcopper/data/Cigar.h>
#include <claraparabricks/genomeworks/cudaaligner/aligner.hpp>
#include <cstdint>

namespace PacBio {

namespace Pancake {

namespace GWInterface {

void RunConvertToPacBioCigarAndScoreGpuKernel(
    PacBio::Data::CigarOperation* cigar, int4* diffs, int64_t* scores,
    const claraparabricks::genomeworks::cudaaligner::DeviceAlignmentsPtrs& aln_ptrs,
    int32_t match_score, int32_t mismatch_penalty, int32_t gap_open_penalty,
    int32_t gap_ext_penalty, cudaStream_t stream, int32_t device_id);

}  // namespace GWInterface

}  // namespace Pancake

}  // namespace PacBio

#endif  // PANCAKE_ALIGNER_BATCH_GPU_GENOMEWORKS_INTERFACE_H
