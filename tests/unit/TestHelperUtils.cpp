// Author: Ivan Sovic

#include "TestHelperUtils.hpp"
#include <pbbam/FastaReader.h>
#include <fstream>
#include <sstream>

namespace PacBio {
namespace PancakeTests {

std::vector<PacBio::BAM::FastaSequence> HelperLoadFasta(const std::string& inFasta)
{
    std::vector<PacBio::BAM::FastaSequence> ret;
    PacBio::BAM::FastaReader inReader{inFasta};
    PacBio::BAM::FastaSequence record;
    while (inReader.GetNext(record))
        ret.emplace_back(record);
    return ret;
}

std::vector<PacBio::Pancake::FastaSequenceId> HelperLoadFastaWithId(const std::string& inFasta,
                                                                    const int32_t seqIdOffset)
{
    std::vector<PacBio::Pancake::FastaSequenceId> ret;
    PacBio::BAM::FastaReader inReader{inFasta};
    PacBio::BAM::FastaSequence record;
    int32_t seqId = 0;
    while (inReader.GetNext(record)) {
        ret.emplace_back(PacBio::Pancake::FastaSequenceId(record, seqId + seqIdOffset));
        ++seqId;
    }
    return ret;
}

std::string HelperLoadFastaAsString(const std::string& inFasta)
{
    std::ostringstream oss;
    auto records = HelperLoadFasta(inFasta);
    for (const auto& record : records)
        oss << ">" << record.Name() << "\n" << record.Bases() << "\n";
    return oss.str();
}

std::vector<std::string> HelperLoadFastaAsStringVector(const std::string& inFasta)
{
    std::vector<std::string> ret;
    PacBio::BAM::FastaReader inReader{inFasta};
    PacBio::BAM::FastaSequence record;
    while (inReader.GetNext(record)) {
        ret.emplace_back(record.Bases());
    }
    return ret;
}

std::vector<std::string> HelperLoadFile(const std::string& inFile)
{
    std::vector<std::string> ret;
    std::string line;
    std::ifstream ifs(inFile);
    if (ifs.is_open() == false) {
        throw std::runtime_error("Cannot open file " + inFile + " for reading!");
    }
    while (std::getline(ifs, line)) {
        ret.emplace_back(std::move(line));
    }
    return ret;
}

void HelperLoadBatchData(
    const std::vector<std::pair<std::string, std::string>>& batchDataSequenceFiles,
    const int32_t seqIdOffset, const PacBio::Pancake::MapperCLRMapSettings& mapSettings,
    std::vector<PacBio::Pancake::MapperBatchChunk>& retBatchData,
    std::vector<PacBio::BAM::FastaSequence>& retAllSeqs)
{
    /*
     * This function takes a vector of pairs, where each pair is a pair of filename paths, one
     * for the target sequences and one for the query sequences; and then loads the sequences.
     * One pair is what a single MapperCLR::MapAndAlign would be run on.
     * Here, the MapperBatchCPU will take a vector of those chunks and align them at once.
     *
     * This function returns a flat vector of all sequences that were loaded (targets and queries),
     * and a vector of MapperBatchChunk.
     * The flat vector is needed because MapperBatchChunk uses FastaSequenceCached which uses pointers
     * to those sequences, so the life span of the original data needs to be ensured.
    */

    retBatchData.clear();
    retAllSeqs.clear();

    for (const auto& vals : batchDataSequenceFiles) {
        const auto& targetFile = vals.first;
        const auto& queryFile = vals.second;
        PacBio::Pancake::MapperBatchChunk bd;

        // Load target sequences and construct the FastaSequenceCached objects.
        std::vector<PacBio::BAM::FastaSequence> targetSeqs =
            PacBio::PancakeTests::HelperLoadFasta(targetFile);
        for (size_t seqId = 0; seqId < targetSeqs.size(); ++seqId) {
            const auto& seq = targetSeqs[seqId];
            // const int32_t seqId = retAllSeqs.size();
            retAllSeqs.emplace_back(std::move(seq));
            auto newFsc = PacBio::Pancake::FastaSequenceCached(
                std::to_string(seqId), retAllSeqs.back().Bases().c_str(),
                retAllSeqs.back().Bases().size(), seqId + seqIdOffset);
            bd.targetSeqs.AddRecord(std::move(newFsc));
        }

        // Load query sequences and construct the FastaSequenceCached objects.
        std::vector<PacBio::BAM::FastaSequence> querySeqs =
            PacBio::PancakeTests::HelperLoadFasta(queryFile);
        for (size_t seqId = 0; seqId < querySeqs.size(); ++seqId) {
            const auto& seq = querySeqs[seqId];
            // const int32_t seqId = retAllSeqs.size();
            retAllSeqs.emplace_back(std::move(seq));
            auto newFsc = PacBio::Pancake::FastaSequenceCached(
                std::to_string(seqId), retAllSeqs.back().Bases().c_str(),
                retAllSeqs.back().Bases().size(), seqId + seqIdOffset);
            bd.querySeqs.AddRecord(std::move(newFsc));
        }

        bd.mapSettings = mapSettings;

        retBatchData.emplace_back(std::move(bd));
    }
}

void HelperLoadBatchData(
    const std::vector<std::pair<std::string, std::string>>& batchDataSequenceFiles,
    const int32_t seqIdOffset, const double freqPercentile,
    const PacBio::Pancake::SeedDBParameters& seedParamsPrimary,
    const PacBio::Pancake::SeedDBParameters& seedParamsFallback,
    std::vector<PacBio::Pancake::MapperBatchChunk>& retBatchData,
    std::vector<PacBio::BAM::FastaSequence>& retAllSeqs)
{
    /*
     * This function takes a vector of pairs, where each pair is a pair of filename paths, one
     * for the target sequences and one for the query sequences; and then loads the sequences.
     * One pair is what a single MapperCLR::MapAndAlign would be run on.
     * Here, the MapperBatchCPU will take a vector of those chunks and align them at once.
     *
     * This function returns a flat vector of all sequences that were loaded (targets and queries),
     * and a vector of MapperBatchChunk.
     * The flat vector is needed because MapperBatchChunk uses FastaSequenceCached which uses pointers
     * to those sequences, so the life span of the original data needs to be ensured.
    */

    PacBio::Pancake::MapperCLRMapSettings mapSettings;
    mapSettings.freqPercentile = freqPercentile;
    mapSettings.seedParams = seedParamsPrimary;
    mapSettings.seedParamsFallback = seedParamsFallback;
    HelperLoadBatchData(batchDataSequenceFiles, seqIdOffset, mapSettings, retBatchData, retAllSeqs);
}

void HelperLoadBatchData(
    const std::vector<std::tuple<std::string, std::string, PacBio::Pancake::MapperCLRMapSettings>>&
        batchDataSequenceFiles,
    const int32_t seqIdOffset, std::vector<PacBio::Pancake::MapperBatchChunk>& retBatchData,
    std::vector<PacBio::BAM::FastaSequence>& retAllSeqs)
{
    /*
     * This function takes a vector of pairs, where each pair is a pair of filename paths, one
     * for the target sequences and one for the query sequences; and then loads the sequences.
     * One pair is what a single MapperCLR::MapAndAlign would be run on.
     * Here, the MapperBatchCPU will take a vector of those chunks and align them at once.
     *
     * This function returns a flat vector of all sequences that were loaded (targets and queries),
     * and a vector of MapperBatchChunk.
     * The flat vector is needed because MapperBatchChunk uses FastaSequenceCached which uses pointers
     * to those sequences, so the life span of the original data needs to be ensured.
    */

    retBatchData.clear();
    retAllSeqs.clear();

    for (const auto& vals : batchDataSequenceFiles) {
        const auto& targetFile = std::get<0>(vals);
        const auto& queryFile = std::get<1>(vals);
        const auto& mapSettings = std::get<2>(vals);
        PacBio::Pancake::MapperBatchChunk bd;

        // Load target sequences and construct the FastaSequenceCached objects.
        std::vector<PacBio::BAM::FastaSequence> targetSeqs =
            PacBio::PancakeTests::HelperLoadFasta(targetFile);
        for (size_t seqId = 0; seqId < targetSeqs.size(); ++seqId) {
            const auto& seq = targetSeqs[seqId];
            // const int32_t seqId = retAllSeqs.size();
            retAllSeqs.emplace_back(std::move(seq));
            auto newFsc = PacBio::Pancake::FastaSequenceCached(
                std::to_string(seqId), retAllSeqs.back().Bases().c_str(),
                retAllSeqs.back().Bases().size(), seqId + seqIdOffset);
            bd.targetSeqs.AddRecord(std::move(newFsc));
        }

        // Load query sequences and construct the FastaSequenceCached objects.
        std::vector<PacBio::BAM::FastaSequence> querySeqs =
            PacBio::PancakeTests::HelperLoadFasta(queryFile);
        for (size_t seqId = 0; seqId < querySeqs.size(); ++seqId) {
            const auto& seq = querySeqs[seqId];
            // const int32_t seqId = retAllSeqs.size();
            retAllSeqs.emplace_back(std::move(seq));
            auto newFsc = PacBio::Pancake::FastaSequenceCached(
                std::to_string(seqId), retAllSeqs.back().Bases().c_str(),
                retAllSeqs.back().Bases().size(), seqId + seqIdOffset);
            bd.querySeqs.AddRecord(std::move(newFsc));
        }

        bd.mapSettings = mapSettings;

        retBatchData.emplace_back(std::move(bd));
    }
}

std::vector<std::vector<std::string>> HelperFormatBatchMappingResults(
    const std::vector<std::vector<PacBio::Pancake::MapperBaseResult>>& results)
{
    std::vector<std::vector<std::string>> resultsStr;
    for (size_t i = 0; i < results.size(); ++i) {
        const auto& resultsForBatchElement = results[i];
        std::vector<std::string> chunkResults;
        for (size_t j = 0; j < resultsForBatchElement.size(); ++j) {
            const auto& queryMappings = resultsForBatchElement[j];
            for (const auto& mapping : queryMappings.mappings) {
                // std::cerr << PacBio::Pancake::OverlapWriterBase::PrintOverlapAsM4(
                //                  *mapping->mapping, "", "", true, false)
                //           << "\n";

                chunkResults.emplace_back(PacBio::Pancake::OverlapWriterBase::PrintOverlapAsM4(
                    *mapping->mapping, "", "", true, false));

                std::cerr << "\"" << chunkResults.back() << "\",\n";
            }
        }
        resultsStr.emplace_back(std::move(chunkResults));
    }
    return resultsStr;
}

std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<std::string>>>
HelperFormatBatchMappingResultsForMockingFlags(
    const std::vector<std::vector<PacBio::Pancake::MapperBaseResult>>& results,
    const std::vector<std::tuple<std::string, std::string, PacBio::Pancake::MapperCLRMapSettings>>&
        batchData,
    const PacBio::Pancake::MapperSelfHitPolicy alignSelfHitPolicy)
{
    std::vector<std::vector<std::string>> resultsStr;
    std::vector<std::vector<std::string>> expectedStr;

    // To be usesd as a label for the alignment policy.
    const std::string mockAlignments =
        (alignSelfHitPolicy == PacBio::Pancake::MapperSelfHitPolicy::PERFECT_ALIGNMENT)
            ? "mocked"
            : "computed";

    for (size_t i = 0; i < results.size(); ++i) {
        const auto& resultsForBatchElement = results[i];

        // To be usesd as a label for the mapping policy. Chunk specific.
        const std::string mockMappings = (std::get<2>(batchData[i]).selfHitPolicy ==
                                          PacBio::Pancake::MapperSelfHitPolicy::PERFECT_ALIGNMENT)
                                             ? "mocked"
                                             : "computed";

        // Results and expected alignments for a single batch.
        std::vector<std::string> chunkResults;
        std::vector<std::string> chunkExpected;

        // Collect all the results/expected alignments for this batch and append labels for mocking.
        for (size_t j = 0; j < resultsForBatchElement.size(); ++j) {
            const auto& queryMappings = resultsForBatchElement[j];
            for (const auto& mapping : queryMappings.mappings) {
                const std::string resultOvl = PacBio::Pancake::OverlapWriterBase::PrintOverlapAsM4(
                    *mapping->mapping, "", "", true, false);

                const bool isMockingCandidate = mapping->mapping->Aid == mapping->mapping->Bid &&
                                                mapping->mapping->Astart == 0 &&
                                                mapping->mapping->Aend == mapping->mapping->Alen &&
                                                mapping->mapping->Bstart == 0 &&
                                                mapping->mapping->Bend == mapping->mapping->Blen;

                const std::string expectedOvlWithMocking =
                    resultOvl + (isMockingCandidate ? (" " + mockMappings + " " + mockAlignments)
                                                    : " computed computed");

                const std::string resultOvlWithMocking =
                    resultOvl + (mapping->isMockedMapping ? " mocked" : " computed") +
                    (mapping->isMockedAlignment ? " mocked" : " computed");

                chunkExpected.emplace_back(expectedOvlWithMocking);
                chunkResults.emplace_back(resultOvlWithMocking);

                // // Debug output.
                // std::cerr << "[batch] expectedOvlWithMocking = " << expectedOvlWithMocking << "\n";
                // std::cerr << "[batch] resultOvlWithMocking   = " << resultOvlWithMocking << "\n";
                // std::cerr << "\n";
            }
        }

        resultsStr.emplace_back(std::move(chunkResults));
        expectedStr.emplace_back(std::move(chunkExpected));
    }
    return {expectedStr, resultsStr};
}

PacBio::Pancake::MapperCLRSettings HelperInitPancakeSettingsSubread()
{
    PacBio::Pancake::MapperCLRSettings settings;

    settings.map.seedParams.KmerSize = 15;
    settings.map.seedParams.MinimizerWindow = 5;
    settings.map.seedParams.Spacing = 0;
    settings.map.seedParams.UseHPCForSeedsOnly = true;

    settings.map.secondaryAllowedOverlapFractionQuery = 0.0;
    settings.map.secondaryAllowedOverlapFractionTarget = 0.05;

    settings.map.seedParamsFallback = settings.map.seedParams;
    settings.map.seedParamsFallback.KmerSize = 10;
    settings.map.seedParamsFallback.MinimizerWindow = 5;

    settings.map.freqPercentile = 0.000;
    settings.map.seedJoinDist = 10000;
    settings.map.longMergeBandwidth = 10000;
    settings.map.maxFlankExtensionDist = settings.map.seedJoinDist;
    settings.map.minAlignmentSpan = 200;

    settings.align.alnParamsGlobal.zdrop = 400;
    settings.align.alnParamsGlobal.zdrop2 = 200;
    settings.align.alnParamsGlobal.alignBandwidth = 500;
    settings.align.alnParamsGlobal.endBonus = 1000;
    settings.align.alnParamsGlobal.matchScore = 2;
    settings.align.alnParamsGlobal.mismatchPenalty = 4;
    settings.align.alnParamsGlobal.gapOpen1 = 4;
    settings.align.alnParamsGlobal.gapExtend1 = 2;
    settings.align.alnParamsGlobal.gapOpen2 = 24;
    settings.align.alnParamsGlobal.gapExtend2 = 1;

    settings.align.alignerTypeGlobal = PacBio::Pancake::AlignerType::KSW2;

    settings.align.alignerTypeExt = PacBio::Pancake::AlignerType::KSW2;
    settings.align.alnParamsExt = settings.align.alnParamsGlobal;

    return settings;
}

}  // namespace PancakeTests
}  // namespace PacBio
