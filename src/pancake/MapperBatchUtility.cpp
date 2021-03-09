// Authors: Ivan Sovic

#include <pacbio/alignment/AlignmentTools.h>
#include <pacbio/pancake/MapperBatchUtility.h>
#include <pacbio/pancake/OverlapWriterBase.h>
#include <pbcopper/logging/Logging.h>

namespace PacBio {
namespace Pancake {

std::vector<PacBio::Pancake::MapperBatchChunk> ConstructBatchData(
    const std::vector<std::pair<std::vector<std::string>, std::vector<std::string>>>& inData)
{
    std::vector<PacBio::Pancake::MapperBatchChunk> batchData;

    for (const auto& dataPair : inData) {
        const auto& queries = dataPair.first;
        const auto& targets = dataPair.second;

        PacBio::Pancake::MapperBatchChunk chunk;

        // Add the target sequences to the chunk.
        for (size_t i = 0; i < targets.size(); ++i) {
            const auto& seq = targets[i];
            const int32_t seqId = i;
            auto seqCache = PacBio::Pancake::FastaSequenceCached(std::to_string(seqId), seq.c_str(),
                                                                 seq.size(), seqId);
            chunk.targetSeqs.emplace_back(std::move(seqCache));
        }

        // Add the query sequences to the chunk.
        for (size_t i = 0; i < queries.size(); ++i) {
            const auto& seq = queries[i];
            const int32_t seqId = i;
            auto seqCache = PacBio::Pancake::FastaSequenceCached(std::to_string(seqId), seq.c_str(),
                                                                 seq.size(), seqId);
            chunk.querySeqs.emplace_back(std::move(seqCache));
        }

        batchData.emplace_back(std::move(chunk));
    }

    return batchData;
}

void PrepareSequencesForBatchAlignment(
    const std::vector<MapperBatchChunk>& batchChunks,
    const std::vector<std::vector<std::string>>& querySeqsRev,
    const std::vector<std::vector<MapperBaseResult>>& mappingResults,
    const MapperSelfHitPolicy selfHitPolicy, std::vector<PairForBatchAlignment>& retPartsGlobal,
    std::vector<PairForBatchAlignment>& retPartsSemiglobal,
    std::vector<AlignmentStitchInfo>& retAlnStitchInfo, int32_t& retLongestSequence)
{
    retPartsGlobal.clear();
    retPartsSemiglobal.clear();
    retAlnStitchInfo.clear();
    retLongestSequence = 0;

    // Results are a vector for every chunk (one chunk is one ZMW).
    for (size_t resultId = 0; resultId < mappingResults.size(); ++resultId) {
        const auto& result = mappingResults[resultId];
        const auto& chunk = batchChunks[resultId];

        // One chunk can have multiple queries (subreads).
        for (size_t qId = 0; qId < result.size(); ++qId) {
            // Prepare the query data in fwd and rev.
            const char* qSeqFwd = chunk.querySeqs[qId].c_str();
            const std::string& qSeqRevString = querySeqsRev[resultId][qId];
            const char* qSeqRev = qSeqRevString.c_str();

            // Each query can have multiple mappings.
            for (size_t mapId = 0; mapId < result[qId].mappings.size(); ++mapId) {
                const auto& mapping = result[qId].mappings[mapId];
                const char* tSeq = chunk.targetSeqs[mapping->mapping->Bid].c_str();

                const auto& aln = mapping->mapping;

                // Skip self-hits unless the default policy is used, in which case align all.
                if (selfHitPolicy != MapperSelfHitPolicy::DEFAULT && aln->Aid == aln->Bid) {
                    continue;
                }

                // AlignmentStitchVector singleQueryStitches;
                AlignmentStitchInfo singleAlnStitches(resultId, qId, mapId);

                // Each mapping is split into regions in between seed hits for alignment.
                for (size_t regId = 0; regId < mapping->regionsForAln.size(); ++regId) {
                    const auto& region = mapping->regionsForAln[regId];

                    // Prepare the sequences for alignment.
                    const char* qSeqInStrand = region.queryRev ? qSeqRev : qSeqFwd;
                    const char* tSeqInStrand = tSeq;
                    int32_t qStart = region.qStart;
                    int32_t tStart = region.tStart;
                    const int32_t qSpan = region.qSpan;
                    const int32_t tSpan = region.tSpan;

                    PairForBatchAlignment part{qSeqInStrand + qStart, qSpan,
                                               tSeqInStrand + tStart, tSpan,
                                               region.type,           region.queryRev};

                    retLongestSequence = std::max(retLongestSequence, std::max(qSpan, tSpan));

                    if (region.type == RegionType::GLOBAL) {
                        singleAlnStitches.parts.emplace_back(AlignmentStitchPart{
                            region.type, static_cast<int64_t>(retPartsGlobal.size()),
                            static_cast<int64_t>(regId)});
                        retPartsGlobal.emplace_back(std::move(part));
                    } else {
                        singleAlnStitches.parts.emplace_back(AlignmentStitchPart{
                            region.type, static_cast<int64_t>(retPartsSemiglobal.size()),
                            static_cast<int64_t>(regId)});
                        retPartsSemiglobal.emplace_back(std::move(part));
                    }
                }
                retAlnStitchInfo.emplace_back(std::move(singleAlnStitches));
            }
        }
    }
}

OverlapPtr StitchSingleAlignment(const OverlapPtr& aln,
                                 const std::vector<AlignmentRegion>& regionsForAln,
                                 const std::vector<AlignmentResult>& internalAlns,
                                 const std::vector<AlignmentResult>& flankAlns,
                                 const std::vector<AlignmentStitchPart>& parts)
{
    if (parts.empty()) {
        return nullptr;
    }

    auto ret = createOverlap(aln);
    ret->Cigar.clear();

    int32_t newQueryStart = -1;
    int32_t newTargetStart = -1;
    int32_t newQueryEnd = 0;
    int32_t newTargetEnd = 0;

    Alignment::DiffCounts diffs;
    int32_t score = 0;

    for (const auto& part : parts) {
        const auto& region = regionsForAln[part.regionId];

        if (part.regionType == RegionType::FRONT) {
            const auto& partAln = flankAlns[part.partId];
            if (partAln.valid == false) {
                return nullptr;
            }
            PacBio::Data::Cigar cigar = partAln.cigar;
            std::reverse(cigar.begin(), cigar.end());
            MergeCigars(ret->Cigar, cigar);
            newQueryStart = region.qStart + region.qSpan - partAln.lastQueryPos;
            newTargetStart = region.tStart + region.tSpan - partAln.lastTargetPos;
            newQueryEnd = region.qStart + region.qSpan;
            newTargetEnd = region.tStart + region.tSpan;
            diffs += partAln.diffs;
            score += partAln.score;

        } else if (part.regionType == RegionType::BACK) {
            const auto& partAln = flankAlns[part.partId];
            if (partAln.valid == false) {
                return nullptr;
            }
            MergeCigars(ret->Cigar, partAln.cigar);
            if (newQueryStart < 0) {
                newQueryStart = region.qStart;
                newTargetStart = region.tStart;
            }
            newQueryEnd = region.qStart + partAln.lastQueryPos;
            newTargetEnd = region.tStart + partAln.lastTargetPos;
            diffs += partAln.diffs;
            score += partAln.score;

        } else {
            const auto& partAln = internalAlns[part.partId];
            if (partAln.valid == false) {
                return nullptr;
            }
            MergeCigars(ret->Cigar, partAln.cigar);
            if (newQueryStart < 0) {
                newQueryStart = region.qStart;
                newTargetStart = region.tStart;
            }
            newQueryEnd = region.qStart + partAln.lastQueryPos;
            newTargetEnd = region.tStart + partAln.lastTargetPos;
            diffs += partAln.diffs;
            score += partAln.score;
        }
    }

    // Skip if the alignment is not valid.
    if (ret == nullptr) {
        return ret;
    }

    // If there is no alignment, reset the overlap.
    if (ret->Cigar.empty()) {
        return nullptr;
    }

    ret->Astart = newQueryStart;
    ret->Aend = newQueryEnd;
    ret->Bstart = newTargetStart;
    ret->Bend = newTargetEnd;

    // Reverse the CIGAR and the coordinates if needed.
    if (ret->Brev) {
        // CIGAR reversal.
        std::reverse(ret->Cigar.begin(), ret->Cigar.end());

        // Reverse the query coordinates.
        std::swap(ret->Astart, ret->Aend);
        ret->Astart = ret->Alen - ret->Astart;
        ret->Aend = ret->Alen - ret->Aend;

        // Get the forward-oriented target coordinates.
        std::swap(ret->Bstart, ret->Bend);
        ret->Bstart = ret->Blen - ret->Bstart;
        ret->Bend = ret->Blen - ret->Bend;
    }

    // Set the alignment identity and edit distance.
    // Alignment::DiffCounts diffs = CigarDiffCounts(ret->Cigar);
    diffs.Identity(false, false, ret->Identity, ret->EditDistance);
    ret->Score = score;

    // std::cerr << "Testing: " << *ret
    //           << "\n";
    // std::cerr << "    - qSpan = " << (diffs.numEq + diffs.numX + diffs.numI) << "\n";
    // std::cerr << "    - tSpan = " << (diffs.numEq + diffs.numX + diffs.numD) << "\n";
    // std::cerr << "    - diffs = " << diffs << "\n";
    // std::cerr << "\n";

    return ret;
}

void StitchAlignments(std::vector<std::vector<MapperBaseResult>>& mappingResults,
                      const std::vector<MapperBatchChunk>& batchChunks,
                      const std::vector<std::vector<std::string>>& querySeqsRev,
                      const std::vector<AlignmentResult>& internalAlns,
                      const std::vector<AlignmentResult>& flankAlns,
                      const std::vector<AlignmentStitchInfo>& alnStitchInfo)
{
    for (const auto& singleAlnInfo : alnStitchInfo) {
        // Not initialized for some reason, skip it.
        if (singleAlnInfo.batchId < 0 || singleAlnInfo.queryId < 0 || singleAlnInfo.mapId < 0) {
            continue;
        }

        // Sanity check.
        if (mappingResults[singleAlnInfo.batchId][singleAlnInfo.queryId]
                .mappings[singleAlnInfo.mapId] == nullptr) {
            continue;
        }
        auto& mapping = mappingResults[singleAlnInfo.batchId][singleAlnInfo.queryId]
                            .mappings[singleAlnInfo.mapId];

        // Sanity check.
        if (mapping->mapping == nullptr) {
            continue;
        }
        auto& aln = mapping->mapping;

        // Do the stitching, and swap.
        OverlapPtr newAln = StitchSingleAlignment(aln, mapping->regionsForAln, internalAlns,
                                                  flankAlns, singleAlnInfo.parts);
        std::swap(aln, newAln);

        if (aln == nullptr) {
            continue;
        }

        {  // Validation of the final alignment.
            const char* querySeq =
                (aln->Brev)
                    ? querySeqsRev[singleAlnInfo.batchId][singleAlnInfo.queryId].c_str()
                    : batchChunks[singleAlnInfo.batchId].querySeqs[singleAlnInfo.queryId].c_str();
            const char* targetSeq = batchChunks[singleAlnInfo.batchId].targetSeqs[aln->Bid].c_str();

            const int32_t qStart = (aln->Brev) ? (aln->Alen - aln->Aend) : aln->Astart;
            const int32_t tStart = aln->BstartFwd();

            PacBio::BAM::Cigar revCigar;
            if (aln->Brev) {
                revCigar.insert(revCigar.end(), aln->Cigar.rbegin(), aln->Cigar.rend());
            }
            PacBio::BAM::Cigar& cigarInStrand = (aln->Brev) ? revCigar : aln->Cigar;

            try {
                ValidateCigar(querySeq + qStart, aln->ASpan(), targetSeq + tStart, aln->BSpan(),
                              cigarInStrand, "Full length validation, fwd.");

            } catch (std::exception& e) {
                PBLOG_WARN << "[Note: Exception caused by ValidateCigar in StitchAlignments] "
                           << e.what() << "\n";
                PBLOG_DEBUG << "singleAlnInfo.batchId = " << singleAlnInfo.batchId;
                PBLOG_DEBUG << "singleAlnInfo.queryId = " << singleAlnInfo.queryId;
                PBLOG_DEBUG << "Aligned: \"" << *newAln << "\"";
                PBLOG_DEBUG << mappingResults[singleAlnInfo.batchId][singleAlnInfo.queryId] << "\n";
                aln = nullptr;
                continue;
            }
        }
    }
}

void StitchAlignmentsInParallel(std::vector<std::vector<MapperBaseResult>>& mappingResults,
                                const std::vector<MapperBatchChunk>& batchChunks,
                                const std::vector<std::vector<std::string>>& querySeqsRev,
                                const std::vector<AlignmentResult>& internalAlns,
                                const std::vector<AlignmentResult>& flankAlns,
                                const std::vector<AlignmentStitchInfo>& alnStitchInfo,
                                Parallel::FireAndForget* faf)
{
    // Determine how many records should land in each thread, spread roughly evenly.
    const int32_t numThreads = faf ? faf->NumThreads() : 1;
    const int32_t numRecords = alnStitchInfo.size();
    const std::vector<std::pair<int32_t, int32_t>> jobsPerThread =
        PacBio::Pancake::DistributeJobLoad<int32_t>(numThreads, numRecords);

    const auto Submit = [&jobsPerThread, &batchChunks, &querySeqsRev, &internalAlns, &flankAlns,
                         &alnStitchInfo, &mappingResults](int32_t idx) {
        const int32_t jobStart = jobsPerThread[idx].first;
        const int32_t jobEnd = jobsPerThread[idx].second;
        for (int32_t jobId = jobStart; jobId < jobEnd; ++jobId) {
            const AlignmentStitchInfo& singleAlnInfo = alnStitchInfo[jobId];

            // Not initialized for some reason, skip it.
            if (singleAlnInfo.batchId < 0 || singleAlnInfo.queryId < 0 || singleAlnInfo.mapId < 0) {
                continue;
            }

            // Sanity check.
            if (mappingResults[singleAlnInfo.batchId][singleAlnInfo.queryId]
                    .mappings[singleAlnInfo.mapId] == nullptr) {
                continue;
            }
            auto& mapping = mappingResults[singleAlnInfo.batchId][singleAlnInfo.queryId]
                                .mappings[singleAlnInfo.mapId];

            // Sanity check.
            if (mapping->mapping == nullptr) {
                continue;
            }
            auto& aln = mapping->mapping;

            // Do the stitching, and swap.
            OverlapPtr newAln = StitchSingleAlignment(aln, mapping->regionsForAln, internalAlns,
                                                      flankAlns, singleAlnInfo.parts);
            std::swap(aln, newAln);

            if (aln == nullptr) {
                continue;
            }

            {  // Validation of the final alignment.
                const char* querySeq =
                    (aln->Brev) ? querySeqsRev[singleAlnInfo.batchId][singleAlnInfo.queryId].c_str()
                                : batchChunks[singleAlnInfo.batchId]
                                      .querySeqs[singleAlnInfo.queryId]
                                      .c_str();
                const char* targetSeq =
                    batchChunks[singleAlnInfo.batchId].targetSeqs[aln->Bid].c_str();

                const int32_t qStart = (aln->Brev) ? (aln->Alen - aln->Aend) : aln->Astart;
                const int32_t tStart = aln->BstartFwd();

                PacBio::BAM::Cigar revCigar;
                if (aln->Brev) {
                    revCigar.insert(revCigar.end(), aln->Cigar.rbegin(), aln->Cigar.rend());
                }
                PacBio::BAM::Cigar& cigarInStrand = (aln->Brev) ? revCigar : aln->Cigar;

                try {
                    ValidateCigar(querySeq + qStart, aln->ASpan(), targetSeq + tStart, aln->BSpan(),
                                  cigarInStrand, "Full length validation, fwd.");

                } catch (std::exception& e) {
                    PBLOG_WARN << "[Note: Exception caused by ValidateCigar in StitchAlignments] "
                               << e.what() << "\n";
                    PBLOG_DEBUG << "singleAlnInfo.batchId = " << singleAlnInfo.batchId;
                    PBLOG_DEBUG << "singleAlnInfo.queryId = " << singleAlnInfo.queryId;
                    PBLOG_DEBUG << "Aligned: \"" << *newAln << "\"";
                    PBLOG_DEBUG << mappingResults[singleAlnInfo.batchId][singleAlnInfo.queryId]
                                << "\n";
                    aln = nullptr;
                    continue;
                }
            }
        }
    };
    Parallel::Dispatch(faf, jobsPerThread.size(), Submit);
}

void SetUnalignedAndMockedMappings(std::vector<std::vector<MapperBaseResult>>& mappingResults,
                                   const bool mockPerfectAlignment,
                                   const int32_t matchScoreForMockAlignment)
{
    for (size_t chunkId = 0; chunkId < mappingResults.size(); ++chunkId) {
        auto& result = mappingResults[chunkId];

        // One chunk can have multiple queries (subreads).
        for (size_t qId = 0; qId < result.size(); ++qId) {
            // Each query can have multiple alignments.
            for (size_t mapId = 0; mapId < result[qId].mappings.size(); ++mapId) {
                if (result[qId].mappings[mapId] == nullptr ||
                    result[qId].mappings[mapId]->mapping == nullptr) {
                    continue;
                }
                OverlapPtr& aln = result[qId].mappings[mapId]->mapping;

                if (mockPerfectAlignment && aln->Aid == aln->Bid) {
                    aln = CreateMockedAlignment(aln, matchScoreForMockAlignment);
                }
                if (aln->Cigar.empty()) {
                    aln = nullptr;
                }
            }
        }
    }
}

std::vector<std::vector<std::string>> ComputeReverseComplements(
    const std::vector<MapperBatchChunk>& batchChunks,
    const std::vector<std::vector<MapperBaseResult>>& mappingResults, Parallel::FireAndForget* faf)
{
    std::vector<std::vector<uint8_t>> shouldReverse(batchChunks.size());
    for (size_t i = 0; i < batchChunks.size(); ++i) {
        shouldReverse[i].resize(batchChunks[i].querySeqs.size(), false);
    }

    // Results are a vector for every chunk (one chunk is one ZMW).
    for (size_t chunkId = 0; chunkId < mappingResults.size(); ++chunkId) {
        auto& result = mappingResults[chunkId];
        // One chunk can have multiple queries (subreads).
        for (size_t qId = 0; qId < result.size(); ++qId) {
            for (size_t mapId = 0; mapId < result[qId].mappings.size(); ++mapId) {
                if (result[qId].mappings[mapId] == nullptr ||
                    result[qId].mappings[mapId]->mapping == nullptr) {
                    continue;
                }
                const OverlapPtr& aln = result[qId].mappings[mapId]->mapping;
                shouldReverse[chunkId][qId] |= aln->Brev;
            }
        }
    }

    // Determine how many records should land in each thread, spread roughly evenly.
    const int32_t numThreads = faf ? faf->NumThreads() : 1;
    const int32_t numRecords = batchChunks.size();
    const std::vector<std::pair<int32_t, int32_t>> jobsPerThread =
        PacBio::Pancake::DistributeJobLoad<int32_t>(numThreads, numRecords);

    std::vector<std::vector<std::string>> querySeqsRev(batchChunks.size());

    const auto Submit = [&batchChunks, &jobsPerThread, &shouldReverse, &querySeqsRev](int32_t idx) {
        const int32_t jobStart = jobsPerThread[idx].first;
        const int32_t jobEnd = jobsPerThread[idx].second;
        for (int32_t chunkId = jobStart; chunkId < jobEnd; ++chunkId) {
            auto& revSeqs = querySeqsRev[chunkId];
            for (size_t qId = 0; qId < batchChunks[chunkId].querySeqs.size(); ++qId) {
                if (shouldReverse[chunkId][qId]) {
                    const auto& query = batchChunks[chunkId].querySeqs[qId];
                    revSeqs.emplace_back(
                        PacBio::Pancake::ReverseComplement(query.c_str(), 0, query.size()));
                } else {
                    revSeqs.emplace_back("");
                }
            }
        }
    };
    Parallel::Dispatch(faf, jobsPerThread.size(), Submit);

    return querySeqsRev;
}

}  // namespace Pancake
}  // namespace PacBio
