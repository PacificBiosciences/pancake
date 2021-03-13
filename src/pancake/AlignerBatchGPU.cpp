// Copyright (c) 2019, Pacific Biosciences of California, Inc.
// All rights reserved.
// See LICENSE.txt.
//
// Contributions from NVIDIA are Copyright (c) 2021, NVIDIA Corporation.
// All rights reserved.
// SPDX-License-Identifier: BSD-3-Clause-Clear
//
// Authors: Ivan Sovic

#include <pacbio/alignment/AlignmentTools.h>
#include <pacbio/pancake/AlignerBatchGPU.h>
#include <pacbio/util/Util.h>
#include <pbcopper/utility/Stopwatch.h>
#include <pacbio/pancake/AlignerBatchGPUGenomeWorksInterface.h>
#include <claraparabricks/genomeworks/utils/device_buffer.hpp>
#include <claraparabricks/genomeworks/utils/pinned_host_vector.hpp>
#include <cstdint>
#include <iostream>
#include <sstream>

namespace PacBio {
namespace Pancake {

struct AlignerBatchGPU::AlignerBatchGPUHostBuffers
{
    claraparabricks::genomeworks::pinned_host_vector<int64_t> starts;
    claraparabricks::genomeworks::pinned_host_vector<int32_t> lengths;
    claraparabricks::genomeworks::pinned_host_vector<int4> diffs;
    claraparabricks::genomeworks::pinned_host_vector<int64_t> scores;
    claraparabricks::genomeworks::pinned_host_vector<PacBio::Data::CigarOperation> cigarOpsBuffer;

    void clearAndResize(int64_t cigarOpsBufferLen, int32_t numberAlns)
    {
        if(starts.size() < static_cast<size_t>(numberAlns))
        {
            // add 20% to avoid later reallocations
            const int32_t newSize = numberAlns * 12 / 10;
            starts.clear();
            starts.resize(newSize);
            lengths.clear();
            lengths.resize(newSize);
            diffs.clear();
            diffs.resize(newSize);
            scores.clear();
            scores.resize(newSize);
        }
        if(cigarOpsBuffer.size() < static_cast<size_t>(cigarOpsBufferLen))
        {
            // add 20% to avoid later reallocations
            const int32_t newSize = cigarOpsBufferLen * 12 / 10;
            cigarOpsBuffer.clear();
            cigarOpsBuffer.resize(newSize);
        }
    }
};

void RetrieveResultsAsPacBioCigar(AlignerBatchGPU::AlignerBatchGPUHostBuffers* hostBuffers, const std::vector<int32_t>& querySpans, const std::vector<int32_t>& targetSpans, const claraparabricks::genomeworks::cudaaligner::FixedBandAligner* aligner, std::vector<AlignmentResult>& results, int32_t matchScore, int32_t mismatchScore, int32_t gapOpenScore, int32_t gapExtScore)
{
    GW_NVTX_RANGE(profiler, "pancake retrieve results");
    namespace gw = claraparabricks::genomeworks;
    if(hostBuffers == nullptr || aligner == nullptr)
    {
        throw std::runtime_error(
                "hostBuffers or aligner should not be nullptr " +
            std::string(__func__) + ".");
    }
    cudaStream_t stream = aligner->get_stream();
    auto allocator = aligner->get_device_allocator();
    gw::cudaaligner::DeviceAlignmentsPtrs aln = aligner->get_alignments_device();
    const int32_t numberOfAlignments = aln.n_alignments;

    // Number of alignments should be the same as number of sequences added.
    if (querySpans.size() != static_cast<size_t>(numberOfAlignments)){
        throw std::runtime_error(
            "Number of alignments doesn't match number of input sequences, in " +
            std::string(__func__) + ".");
    }

    gw::device_buffer<PacBio::Data::CigarOperation> pacbio_cigars_device(aln.total_length, allocator, stream);
    gw::device_buffer<int4> diffs_device(numberOfAlignments, allocator, stream);
    gw::device_buffer<int64_t> scores_device(numberOfAlignments, allocator, stream);
    GWInterface::RunConvertToPacBioCigarAndScoreGpuKernel(pacbio_cigars_device.data(), diffs_device.data(), scores_device.data(), aln, matchScore, mismatchScore, gapOpenScore, gapExtScore, stream, aligner->get_device());

    hostBuffers->clearAndResize(aln.total_length, numberOfAlignments);

    gw::cudautils::device_copy_n_async(aln.starts, numberOfAlignments, hostBuffers->starts.data(), stream);
    gw::cudautils::device_copy_n_async(aln.lengths, numberOfAlignments, hostBuffers->lengths.data(), stream);
    gw::cudautils::device_copy_n_async(diffs_device.data(), numberOfAlignments, hostBuffers->diffs.data(), stream);
    gw::cudautils::device_copy_n_async(scores_device.data(), numberOfAlignments, hostBuffers->scores.data(), stream);
    gw::cudautils::device_copy_n_async(pacbio_cigars_device.data(), aln.total_length, hostBuffers->cigarOpsBuffer.data(), stream);

    GW_CU_CHECK_ERR(cudaStreamSynchronize(stream));
    for(int32_t i=0; i < numberOfAlignments; ++i)
    {
        auto& res = results[i];
        const int32_t len = std::abs(hostBuffers->lengths[i]);
        const bool is_optimal = (hostBuffers->lengths[i] >= 0); // if a alignment is optimal is encoded in the sign of length.
        res.cigar.clear();
        res.cigar.insert(end(res.cigar), begin(hostBuffers->cigarOpsBuffer) + hostBuffers->starts[i], begin(hostBuffers->cigarOpsBuffer) + hostBuffers->starts[i] + len);
        res.lastQueryPos = querySpans[i];
        res.lastTargetPos = targetSpans[i];
        res.maxQueryPos = querySpans[i];
        res.maxTargetPos = targetSpans[i];
        const int32_t numEq = hostBuffers->diffs[i].x;
        const int32_t numX  = hostBuffers->diffs[i].y;
        const int32_t numI  = hostBuffers->diffs[i].z;
        const int32_t numD  = hostBuffers->diffs[i].w;
        const int32_t querySpan  = numEq + numX + numI;
        const int32_t targetSpan = numEq + numX + numD;
        res.valid = (querySpan == querySpans[i] && targetSpan == targetSpans[i] && is_optimal);
        res.score    = hostBuffers->scores[i];
        res.maxScore = hostBuffers->scores[i];
        res.zdropped = false;
        res.diffs.numEq = hostBuffers->diffs[i].x;
        res.diffs.numX  = hostBuffers->diffs[i].y;
        res.diffs.numI  = hostBuffers->diffs[i].z;
        res.diffs.numD  = hostBuffers->diffs[i].w;
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

AlignerBatchGPU::AlignerBatchGPU(const AlignmentParameters& alnParams,
                                 uint32_t maxBandwidth, uint32_t deviceId, int64_t maxGPUMemoryCap)
    : alnParams_(alnParams), gpuStream_(claraparabricks::genomeworks::make_cuda_stream()),
    gpuHostBuffers_(std::make_unique<AlignerBatchGPUHostBuffers>())
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

StatusAddSequencePair AlignerBatchGPU::AddSequencePair(const char* query, int32_t queryLen,
                                                        const char* target, int32_t targetLen)
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

std::pair<int64_t, int64_t> AlignerBatchGPU::AlignAll()
{
    int64_t cpuTime = 0;
    int64_t gpuTime = 0;
    PacBio::Utility::Stopwatch timer;
    alnResults_.clear();
    cpuTime += timer.ElapsedNanoseconds();
    timer.Reset();

    aligner_->align_all();
    alnResults_.resize(querySpans_.size());
    RetrieveResultsAsPacBioCigar(gpuHostBuffers_.get(), querySpans_, targetSpans_, aligner_.get(), alnResults_, alnParams_.matchScore, alnParams_.mismatchPenalty, alnParams_.gapOpen1, alnParams_.gapExtend1);
    gpuTime += timer.ElapsedNanoseconds();
    timer.Reset();
    cpuTime += timer.ElapsedNanoseconds();
    return std::make_pair(cpuTime, gpuTime);
}

}  // namespace Pancake
}  // namespace PacBio
