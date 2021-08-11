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

#include <pacbio/pancake/AlignerBatchGPUGenomeWorksInterface.h>

#include <claraparabricks/genomeworks/cudaaligner/cudaaligner.hpp>
#include <claraparabricks/genomeworks/utils/cudautils.hpp>
#include <claraparabricks/genomeworks/utils/mathutils.hpp>
#include <limits>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <cuda/atomic>
#pragma GCC diagnostic pop

namespace PacBio {
namespace Pancake {

namespace GWInterface {

template <typename RandomAccessIterator, typename ValueType>
__device__ RandomAccessIterator upper_bound(RandomAccessIterator lower_bound,
                                            RandomAccessIterator upper_bound, ValueType query)
{
    assert(upper_bound >= lower_bound);
    while (upper_bound > lower_bound) {
        RandomAccessIterator mid = lower_bound + (upper_bound - lower_bound) / 2;
        const auto mid_value = *mid;
        if (mid_value <= query)
            lower_bound = mid + 1;
        else
            upper_bound = mid;
    }
    return lower_bound;
}

__device__ PacBio::Data::CigarOperationType gw_cigar_to_pb_cigar_op(int8_t action)
{
    using claraparabricks::genomeworks::cudaaligner::AlignmentState;
    using PacBio::Data::CigarOperationType;
    switch (action) {
        case AlignmentState::match:
            return CigarOperationType::SEQUENCE_MATCH;
        case AlignmentState::mismatch:
            return CigarOperationType::SEQUENCE_MISMATCH;
        case AlignmentState::insertion:
            return CigarOperationType::INSERTION;
        case AlignmentState::deletion:
            return CigarOperationType::DELETION;
        default:
            return CigarOperationType::UNKNOWN_OP;
    };
}

__device__ int64_t compute_score(int8_t action, int32_t runlength, int32_t match_score,
                                 int32_t mismatch_penalty, int32_t gap_open_penalty,
                                 int32_t gap_ext_penalty)
{
    using claraparabricks::genomeworks::cudaaligner::AlignmentState;
    switch (action) {
        case AlignmentState::match:
            // Scores are positive.
            return match_score * runlength;
        case AlignmentState::mismatch:
            // Penalties are positive.
            return -(mismatch_penalty * runlength);
        case AlignmentState::insertion:
            // Penalties are positive.
            return -(gap_open_penalty + gap_ext_penalty * (runlength - 1));
        case AlignmentState::deletion:
            // Penalties are positive.
            return -(gap_open_penalty + gap_ext_penalty * (runlength - 1));
        default:
            assert(false);  // action should always be one of the cases above.
            return 0;
    }
}

__global__ void convert_to_pacbio_kernel(PacBio::Data::CigarOperation* cigar, int4* diffs,
                                         cuda::atomic<int64_t, cuda::thread_scope_device>* scores,
                                         const int64_t* starts, const int32_t* lengths,
                                         const int8_t* actions, const int32_t* runlengths,
                                         int64_t total_length, int32_t n_alignments,
                                         int32_t match_score, int32_t mismatch_penalty,
                                         int32_t gap_open_penalty, int32_t gap_ext_penalty)
{
    // This kernel reverses the raw alignment results from cudaaligner
    // and translates them into sequences of PacBio::Data::CigarOperation,
    // i.e. PacBio Cigars.
    using claraparabricks::genomeworks::cudaaligner::AlignmentState;
    const int64_t i = static_cast<int64_t>(blockIdx.x) * blockDim.x + threadIdx.x;

    // Find alignment index:
    const int32_t alignment = upper_bound(starts, starts + n_alignments, i) - starts - 1;

    // Get start and length of the alignment data:
    const int64_t start = starts[alignment];
    const int32_t len = abs(lengths[alignment]);
    // For all data belonging to the alignment...
    if (i < start + len) {
        // Get cudaaligner RLE data
        const int8_t c = actions[2 * start + len - i - 1];
        const int32_t rl = runlengths[2 * start + len - i - 1];
        assert(c >= 0);
        assert(c < 4);
        assert(rl >= 0);
        // Translate
        cigar[i] = {gw_cigar_to_pb_cigar_op(c), static_cast<uint32_t>(rl)};

        // Compute score contribution
        int64_t score =
            compute_score(c, rl, match_score, mismatch_penalty, gap_open_penalty, gap_ext_penalty);

        // Add the information to the alignment's diff and score
        atomicAdd(&diffs[alignment].x, c == AlignmentState::match ? rl : 0);
        atomicAdd(&diffs[alignment].y, c == AlignmentState::mismatch ? rl : 0);
        atomicAdd(&diffs[alignment].z, c == AlignmentState::insertion ? rl : 0);
        atomicAdd(&diffs[alignment].w, c == AlignmentState::deletion ? rl : 0);
        scores[alignment].fetch_add(score, cuda::memory_order_relaxed);
    }
}

void RunConvertToPacBioCigarAndScoreGpuKernel(
    PacBio::Data::CigarOperation* cigar, int4* diffs, int64_t* scores,
    const claraparabricks::genomeworks::cudaaligner::DeviceAlignmentsPtrs& aln_ptrs,
    int32_t match_score, int32_t mismatch_score, int32_t gap_open_score, int32_t gap_ext_score,
    cudaStream_t stream, int32_t device_id)
{
    using claraparabricks::genomeworks::ceiling_divide;
    static_assert(sizeof(int64_t) == sizeof(cuda::atomic<int64_t, cuda::thread_scope_device>),
                  "the cast from int64_t -> cuda::atomic<int64> hack is definitely not supported.");
    claraparabricks::genomeworks::scoped_device_switch device(device_id);
    const int64_t total_length = aln_ptrs.total_length;
    const int32_t n_alignments = aln_ptrs.n_alignments;
    constexpr int32_t n_threads = 256;
    const int64_t n_blocks = ceiling_divide<int64_t>(total_length, n_threads);
    assert(n_alignments > 0);
    assert(n_blocks <= static_cast<int64_t>(std::numeric_limits<uint32_t>::max()));
    GW_CU_CHECK_ERR(cudaMemsetAsync(diffs, 0, n_alignments * sizeof(int4), stream));
    GW_CU_CHECK_ERR(cudaMemsetAsync(scores, 0, n_alignments * sizeof(int64_t), stream));
    convert_to_pacbio_kernel<<<n_blocks, n_threads, 0, stream>>>(
        cigar, diffs,
        reinterpret_cast<cuda::atomic<int64_t, cuda::thread_scope_device>*>(
            scores) /* please replace by atomic_ref once it is available */,
        aln_ptrs.starts, aln_ptrs.lengths, aln_ptrs.actions, aln_ptrs.runlengths, total_length,
        n_alignments, match_score, mismatch_score, gap_open_score, gap_ext_score);
    GW_CU_CHECK_ERR(cudaPeekAtLastError());
}

}  // namespace GWInterface
}  // namespace Pancake
}  // namespace PacBio
