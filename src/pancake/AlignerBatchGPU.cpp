// Copyright (c) 2019, Pacific Biosciences of California, Inc.
// All rights reserved.
// See LICENSE.txt.
//
// Contributions from NVIDIA are Copyright (c) 2021, NVIDIA Corporation.
// All rights reserved.
// SPDX-License-Identifier: BSD-3-Clause-Clear
//
// Authors: Ivan Sovic

#include <claraparabricks/genomeworks/utils/device_buffer.hpp>
#include <claraparabricks/genomeworks/utils/pinned_host_vector.hpp>
#include <cstdint>
#include <iostream>
#include <limits>
#include <pancake/AlignerBatchGPU.hpp>
#include <pancake/AlignerBatchGPUGenomeWorksInterface.hpp>
#include <pancake/AlignmentTools.hpp>
#include <pancake/util/Util.hpp>
#include <sstream>

namespace PacBio {
namespace Pancake {

namespace cpgw = claraparabricks::genomeworks;

struct AlignerBatchGPU::AlignerBatchGPUHostBuffers
{
    cpgw::pinned_host_vector<int32_t> cigarOffsets;
    cpgw::pinned_host_vector<uint32_t> metaData;
    cpgw::pinned_host_vector<int4> diffs;
    cpgw::pinned_host_vector<int64_t> scores;
    cpgw::pinned_host_vector<PacBio::Data::CigarOperation> cigarOpsBuffer;

    void ResizeForNumAlignments(int32_t numAlns)
    {
        assert(numAlns >= 0);
        if (metaData.size() < static_cast<size_t>(numAlns)) {
            try {
                cigarOffsets.clear();
                metaData.clear();
                diffs.clear();
                scores.clear();
                cigarOffsets.shrink_to_fit();
                metaData.shrink_to_fit();
                diffs.shrink_to_fit();
                scores.shrink_to_fit();
                const int32_t newSize = numAlns * 12 / 10;
                cigarOffsets.resize(newSize + 1);
                metaData.resize(newSize);
                diffs.resize(newSize);
                scores.resize(newSize);
            } catch (const std::bad_alloc&) {
                // Free all host buffers and try again with minimum size.
                // If this fails it fails.
                free();
                const int32_t newSize = numAlns;
                cigarOffsets.resize(newSize + 1);
                metaData.resize(newSize);
                diffs.resize(newSize);
                scores.resize(newSize);
            }
        }
    }

    void free()
    {
        cigarOpsBuffer.clear();
        cigarOpsBuffer.shrink_to_fit();
        cigarOffsets.clear();
        metaData.clear();
        diffs.clear();
        scores.clear();
        cigarOffsets.shrink_to_fit();
        metaData.shrink_to_fit();
        diffs.shrink_to_fit();
        scores.shrink_to_fit();
    }
};

static int64_t FindLongestCigarLength(cpgw::pinned_host_vector<int32_t>::const_iterator begin,
                                      cpgw::pinned_host_vector<int32_t>::const_iterator end)
{
    int64_t longestCigarLength = 0;
    if (begin != end) {
        for (auto it = begin + 1; it != end; ++it) {
            const int64_t cigarLength = *it - *(it - 1);
            assert(cigarLength >= 0);
            longestCigarLength = std::max(longestCigarLength, cigarLength);
        }
    }
    return longestCigarLength;
}

template <typename T>
static int64_t GetMaxSize(const cpgw::device_buffer<T>& buffer)
{
    return buffer.get_allocator().get_size_of_largest_free_memory_block();
}

template <typename T>
static int64_t GetMaxSize(const cpgw::pinned_host_vector<T>&)
{
    return std::numeric_limits<int64_t>::max();
}

template <typename T>
static void Reallocate(cpgw::device_buffer<T>& buffer, int64_t size)
{
    buffer.clear_and_resize(size);
}

template <typename T>
static void Reallocate(cpgw::pinned_host_vector<T>& buffer, int64_t size)
{
    buffer.clear();
    buffer.shrink_to_fit();
    buffer.resize(size);
}

template <typename Buffer>
static bool TryReallocate(Buffer& buffer, int64_t size)
{
    try {
        Reallocate(buffer, size);
    } catch (const std::bad_alloc&) {
        return false;
    } catch (const cpgw::device_memory_allocation_exception&) {
        return false;
    }
    return true;
}

template <typename Buffer, typename Iterator>
static void ResizeToEither(Buffer& buffer, Iterator begin, Iterator end)
{
    assert(std::is_sorted(begin, end, std::greater<int64_t>()));
    if (begin == end) {
        return;
    }

    Reallocate(buffer, 0);
    const int64_t maxSize = GetMaxSize(buffer) / sizeof(typename Buffer::value_type);
    Iterator last = end - 1;
    for (Iterator it = begin; it != last; ++it) {
        const int64_t size = *it;
        if (size > maxSize) {
            continue;
        }
        if (TryReallocate(buffer, size)) {
            return;
        }
    }
    // Try to allocate the last size (assuming it is the smallest).
    // If it cannot be allocated => let it throw an exception
    Reallocate(buffer, *last);
}

static std::vector<int64_t> computeSizesToTry(int64_t maxSize, int64_t minSize)
{
    assert(minSize >= 0);
    assert(maxSize >= minSize);
    std::vector<int64_t> sizesToTry;
    sizesToTry.reserve(64);
    for (int64_t size = maxSize; size > minSize; size /= 2) {
        sizesToTry.push_back(size);
    }
    sizesToTry.push_back(minSize);
    return sizesToTry;
}

void RetrieveResultsAsPacBioCigar(
    AlignerBatchGPU::AlignerBatchGPUHostBuffers& hostBuffers,
    const std::vector<int32_t>& querySpans, const std::vector<int32_t>& targetSpans,
    claraparabricks::genomeworks::cudaaligner::FixedBandAligner& aligner,
    std::vector<AlignmentResult>& results, int32_t matchScore, int32_t mismatchScore,
    int32_t gapOpenScore, int32_t gapExtScore)
{
    GW_NVTX_RANGE(profiler, "pancake retrieve results");
    cudaStream_t stream = aligner.get_stream();
    auto allocator = aligner.get_device_allocator();
    const int32_t numberOfAlignments = aligner.num_alignments();

    // Number of alignments should be the same as number of sequences added.
    if (querySpans.size() != static_cast<size_t>(numberOfAlignments)) {
        throw std::runtime_error(
            "Number of alignments doesn't match number of input sequences, in " +
            std::string(__func__) + ".");
    }

    results.clear();

    if (numberOfAlignments == 0) {
        return;
    }

    hostBuffers.ResizeForNumAlignments(numberOfAlignments);

    GW_CU_CHECK_ERR(cudaStreamSynchronize(stream));
    aligner.free_temporary_device_buffers();  // such that we can reuse the memory here.

    cpgw::cudaaligner::DeviceAlignmentsPtrs aln = aligner.get_alignments_device();
    // aln.cigar_offsets contains the start and end of the different cigars in aln.cigar_operations.
    // The offset entry at position i is at the same time the start of cigar i and the end of the
    // previous cigar (except for the first enty).
    // The last entry aln.cigar_offsets[numberOfAlignments] is the end of the last cigar.
    // Therefore, the array has numberOfAlignments + 1 entries in total.
    cpgw::cudautils::device_copy_n_async(aln.cigar_offsets, numberOfAlignments + 1,
                                         hostBuffers.cigarOffsets.data(), stream);

    cpgw::device_buffer<int4> diffsDevice(numberOfAlignments, allocator, stream);
    cpgw::device_buffer<int64_t> scoresDevice(numberOfAlignments, allocator, stream);
    cpgw::device_buffer<PacBio::Data::CigarOperation> pacbioCigarsDevice(0, allocator, stream);

    if (pacbioCigarsDevice.size() < aln.total_length) {
        std::array<int64_t, 3> sizesToTry = {aln.total_length + aln.total_length / 5,
                                             aln.total_length, 0};
        ResizeToEither(pacbioCigarsDevice, begin(sizesToTry), end(sizesToTry));
    }
    if (static_cast<int64_t>(hostBuffers.cigarOpsBuffer.size()) < aln.total_length) {
        std::array<int64_t, 3> sizesToTry = {aln.total_length + aln.total_length / 5,
                                             aln.total_length, 0};
        ResizeToEither(hostBuffers.cigarOpsBuffer, begin(sizesToTry), end(sizesToTry));
    }
    int64_t longestCigarLength = 0;
    if (pacbioCigarsDevice.size() < aln.total_length ||
        static_cast<int64_t>(hostBuffers.cigarOpsBuffer.size()) < aln.total_length) {
        GW_CU_CHECK_ERR(cudaStreamSynchronize(stream));  // wait for cigarOffsets
        longestCigarLength =
            FindLongestCigarLength(begin(hostBuffers.cigarOffsets), end(hostBuffers.cigarOffsets));
        cpgw::cudautils::device_copy_n_async(aln.metadata, numberOfAlignments,
                                             hostBuffers.metaData.data(), stream);
        const std::vector<int64_t> sizesToTry =
            computeSizesToTry(aln.total_length, longestCigarLength);
        if (pacbioCigarsDevice.size() < aln.total_length) {
            ResizeToEither(pacbioCigarsDevice, begin(sizesToTry), end(sizesToTry));
        }
        if (static_cast<int64_t>(hostBuffers.cigarOpsBuffer.size()) < pacbioCigarsDevice.size()) {
            ResizeToEither(hostBuffers.cigarOpsBuffer, begin(sizesToTry), end(sizesToTry));
        }
    } else {
        cpgw::cudautils::device_copy_n_async(aln.metadata, numberOfAlignments,
                                             hostBuffers.metaData.data(), stream);
    }

    results.resize(numberOfAlignments);
    GW_CU_CHECK_ERR(cudaStreamSynchronize(stream));
    const int64_t cigarBufSize = std::min({static_cast<int64_t>(hostBuffers.cigarOpsBuffer.size()),
                                           pacbioCigarsDevice.size(), aln.total_length});
    int32_t completedAlignments = 0;
    while (completedAlignments < numberOfAlignments) {
        GWInterface::RunConvertToPacBioCigarAndScoreGpuKernel(
            pacbioCigarsDevice.data(), diffsDevice.data(), scoresDevice.data(), aln, cigarBufSize,
            completedAlignments, matchScore, mismatchScore, gapOpenScore, gapExtScore, stream,
            aligner.get_device());

        cpgw::cudautils::device_copy_n_async(diffsDevice.data(), numberOfAlignments,
                                             hostBuffers.diffs.data(), stream);
        cpgw::cudautils::device_copy_n_async(scoresDevice.data(), numberOfAlignments,
                                             hostBuffers.scores.data(), stream);
        cpgw::cudautils::device_copy_n_async(pacbioCigarsDevice.data(), cigarBufSize,
                                             hostBuffers.cigarOpsBuffer.data(), stream);

        GW_CU_CHECK_ERR(cudaStreamSynchronize(stream));
        int32_t i = completedAlignments;
        for (; i < numberOfAlignments; ++i) {
            const int32_t idx =
                hostBuffers.metaData[i] & cpgw::cudaaligner::DeviceAlignmentsPtrs::index_mask;
            const bool isOptimal = hostBuffers.metaData[i] >> 31;
            const int64_t startPos =
                hostBuffers.cigarOffsets[i] - hostBuffers.cigarOffsets[completedAlignments];
            const int64_t len = hostBuffers.cigarOffsets[i + 1] - hostBuffers.cigarOffsets[i];
            if (startPos + len > cigarBufSize) {
                if (i == completedAlignments) {
                    throw std::runtime_error(
                        "Could not allocate enough device or pinned host memory!");
                }
                break;
            }
            auto& res = results[idx];
            res.cigar.clear();
            res.cigar.reserve(len);
            res.cigar.insert(end(res.cigar), begin(hostBuffers.cigarOpsBuffer) + startPos,
                             begin(hostBuffers.cigarOpsBuffer) + startPos + len);
            res.lastQueryPos = querySpans[idx];
            res.lastTargetPos = targetSpans[idx];
            res.maxQueryPos = querySpans[idx];
            res.maxTargetPos = targetSpans[idx];
            const int32_t numEq = hostBuffers.diffs[idx].x;
            const int32_t numX = hostBuffers.diffs[idx].y;
            const int32_t numI = hostBuffers.diffs[idx].z;
            const int32_t numD = hostBuffers.diffs[idx].w;
            const int32_t querySpan = numEq + numX + numI;
            const int32_t targetSpan = numEq + numX + numD;
            assert(!isOptimal || querySpan == 0 || querySpan == querySpans[idx]);
            assert(!isOptimal || targetSpan == 0 || targetSpan == targetSpans[idx]);
            res.valid =
                (querySpan == querySpans[idx] && targetSpan == targetSpans[idx] && isOptimal);
            res.score = hostBuffers.scores[idx];
            res.maxScore = hostBuffers.scores[idx];
            res.zdropped = false;
            res.diffs.numEq = numEq;
            res.diffs.numX = numX;
            res.diffs.numI = numI;
            res.diffs.numD = numD;
        }
        completedAlignments += i;
    }
}

int64_t ComputeMaxGPUMemory(int64_t cudaalignerBatches, double maxGPUMemoryFraction)
{
    size_t freeMem = 0, totalMem = 0;

    GW_CU_CHECK_ERR(cudaMemGetInfo(&freeMem, &totalMem));
    const size_t freeUsableMemory = static_cast<double>(freeMem) * maxGPUMemoryFraction;
    const int64_t usableMemoryPerAligner = freeUsableMemory / cudaalignerBatches;

    std::cerr << "cudaalignerBatches = " << cudaalignerBatches << "\n";
    std::cerr << "freeMem = " << freeMem / (1024.0 * 1024.0) << " MB\n";
    std::cerr << "totalMem = " << totalMem / (1024.0 * 1024.0) << " MB\n";
    std::cerr << "usableMemoryPerAligner = " << usableMemoryPerAligner / (1024.0 * 1024.0)
              << " MB \n";

    return usableMemoryPerAligner;
}

AlignerBatchGPU::AlignerBatchGPU(const AlignmentParameters& alnParams, uint32_t maxBandwidth,
                                 uint32_t deviceId, int64_t maxGPUMemoryCap)
    : alnParams_(alnParams)
    , gpuStream_(claraparabricks::genomeworks::make_cuda_stream())
    , gpuHostBuffers_(std::make_unique<AlignerBatchGPUHostBuffers>())
{
    aligner_ = claraparabricks::genomeworks::cudaaligner::create_aligner(
        claraparabricks::genomeworks::cudaaligner::AlignmentType::global_alignment, maxBandwidth,
        gpuStream_.get(), deviceId, maxGPUMemoryCap);
}

AlignerBatchGPU::~AlignerBatchGPU() = default;

void AlignerBatchGPU::Clear()
{
    aligner_->reset();
    alnResults_.clear();
    querySpans_.clear();
    targetSpans_.clear();
}

void AlignerBatchGPU::ResetMaxBandwidth(int32_t maxBandwidth)
{
    aligner_->reset_max_bandwidth(maxBandwidth);
}

StatusAddSequencePair AlignerBatchGPU::AddSequencePairForExtensionAlignment(const char* /*query*/,
                                                                            int32_t /*queryLen*/,
                                                                            const char* /*target*/,
                                                                            int32_t /*targetLen*/)
{
    std::ostringstream oss;
    oss << "The GenomeWorks Cudaaligner does not support extension alignment at this "
           "point.";
    throw std::runtime_error(oss.str());
    return {};
}

StatusAddSequencePair AlignerBatchGPU::AddSequencePairForGlobalAlignment(const char* query,
                                                                         int32_t queryLen,
                                                                         const char* target,
                                                                         int32_t targetLen)
{
    if (queryLen < 0 || targetLen < 0) {
        return StatusAddSequencePair::SEQUENCE_LEN_BELOW_ZERO;
    }

    const claraparabricks::genomeworks::cudaaligner::StatusType s =
        aligner_->add_alignment(target, targetLen, query, queryLen);

    if (s == claraparabricks::genomeworks::cudaaligner::StatusType::exceeded_max_alignments) {
        return StatusAddSequencePair::EXCEEDED_MAX_ALIGNMENTS;

    } else if (s == claraparabricks::genomeworks::cudaaligner::StatusType::
                        exceeded_max_alignment_difference ||
               s == claraparabricks::genomeworks::cudaaligner::StatusType::exceeded_max_length) {
        // Do nothing as this case will be handled by CPU aligner.
        throw std::runtime_error(
            "Not implemented yet: s == StatusType::exceeded_max_alignment_difference || s == "
            "StatusType::exceeded_max_length.\n");

    } else if (s != claraparabricks::genomeworks::cudaaligner::StatusType::success) {
        // fprintf(stderr, "Unknown error in cuda aligner!\n");
        throw std::runtime_error("Unknown error in cuda aligner!\n");

    } else {
        querySpans_.emplace_back(queryLen);
        targetSpans_.emplace_back(targetLen);
        alnResults_.clear();
    }

    return StatusAddSequencePair::OK;
}

void AlignerBatchGPU::AlignAll()
{
    alnResults_.clear();
    aligner_->align_all();
    RetrieveResultsAsPacBioCigar(*gpuHostBuffers_.get(), querySpans_, targetSpans_, *aligner_.get(),
                                 alnResults_, alnParams_.matchScore, alnParams_.mismatchPenalty,
                                 alnParams_.gapOpen1, alnParams_.gapExtend1);
}

}  // namespace Pancake
}  // namespace PacBio
