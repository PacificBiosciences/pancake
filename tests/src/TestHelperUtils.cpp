// Author: Ivan Sovic

#include "TestHelperUtils.h"
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
    const int32_t seqIdOffset, const double freqPercentile,
    const PacBio::Pancake::SeedDB::SeedDBParameters& seedParamsPrimary,
    const PacBio::Pancake::SeedDB::SeedDBParameters& seedParamsFallback,
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

        // Set the seed parameter settings and create a mapper.
        bd.mapSettings.freqPercentile = freqPercentile;
        bd.mapSettings.seedParams = seedParamsPrimary;
        bd.mapSettings.seedParamsFallback = seedParamsFallback;

        retBatchData.emplace_back(std::move(bd));
    }
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
            }
        }
        resultsStr.emplace_back(std::move(chunkResults));
    }
    return resultsStr;
}

}  // namespace PancakeTests
}  // namespace PacBio
