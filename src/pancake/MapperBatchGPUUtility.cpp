// Authors: Ivan Sovic

#include <pancake/MapperBatchGPUUtility.hpp>

#include <pbcopper/logging/Logging.h>

#include <cstdint>
#include <string>

namespace PacBio {
namespace Pancake {

int32_t AlignPartsOnGPU(const std::vector<PairForBatchAlignment>& parts, AlignerBatchBase& aligner,
                        int32_t maxAllowedGap, std::vector<AlignmentResult>& retAlns)
{
    auto PartBandwidthAboveLimit = [](const PairForBatchAlignment& part,
                                      const int32_t maxAllowedGap_) {
        if (maxAllowedGap_ < 0) {
            return false;
        }
        const int32_t lenDiff = std::abs(part.queryLen - part.targetLen);
        return (lenDiff > maxAllowedGap_) || (part.maxGap > maxAllowedGap_);
    };

    const std::string logPrefix("[" + std::string(__FUNCTION__) + "]");

    retAlns.resize(parts.size());

    int32_t totalNumNotValid = 0;
    size_t partId = 0;
    int32_t numSkipped = 0;
    while (partId < parts.size()) {
        aligner.Clear();

        std::vector<size_t> partIds;
        partIds.reserve(parts.size());

        PBLOG_TRACE << logPrefix << " Preparing sequences for GPU alignment.";
        for (; partId < parts.size(); ++partId) {
            const PairForBatchAlignment& part = parts[partId];

            if (retAlns[partId].valid) {
                continue;
            }

            // Skip parts with bandwidth above the bin size.
            if (PartBandwidthAboveLimit(part, maxAllowedGap)) {
                ++numSkipped;
                continue;
            }

            partIds.emplace_back(partId);

            StatusAddSequencePair rv;
            if (part.regionType == RegionType::FRONT) {
                // Reverse the sequences for front flank alignment. No need to complement.
                std::string query(part.query, part.queryLen);
                std::reverse(query.begin(), query.end());
                std::string target(part.target, part.targetLen);
                std::reverse(target.begin(), target.end());
                rv = aligner.AddSequencePairForExtensionAlignment(query, target);
            } else {
                if (part.regionType == RegionType::GLOBAL) {
                    rv = aligner.AddSequencePairForGlobalAlignment(
                        std::string_view(part.query, part.queryLen),
                        std::string_view(part.target, part.targetLen));
                } else {
                    rv = aligner.AddSequencePairForExtensionAlignment(
                        std::string_view(part.query, part.queryLen),
                        std::string_view(part.target, part.targetLen));
                }
            }

            if (rv == StatusAddSequencePair::EXCEEDED_MAX_ALIGNMENTS) {
                break;

            } else if (rv != StatusAddSequencePair::OK) {
                throw std::runtime_error(
                    logPrefix +
                    " Error occurred while trying to add sequences for batch alignment.");
            }
        }

        PBLOG_TRACE << logPrefix << " Aligning batch of " << aligner.BatchSize()
                    << " sequence pairs.";
        aligner.AlignAll();

        std::vector<AlignmentResult>& partInternalAlns = aligner.GetAlnResults();

        int32_t numNotValid = numSkipped;
        for (size_t i = 0; i < partInternalAlns.size(); ++i) {
            const auto& aln = partInternalAlns[i];
            if (aln.valid == false) {
                ++numNotValid;
            }
            retAlns[partIds[i]] = std::move(partInternalAlns[i]);
        }
        totalNumNotValid += numNotValid;
    }
    PBLOG_TRACE << logPrefix << " Total not valid: " << totalNumNotValid << " / " << retAlns.size();
    return totalNumNotValid;
}

}  // namespace Pancake
}  // namespace PacBio
