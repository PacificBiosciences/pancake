// Author: Ivan Sovic

#ifndef PANCAKE_TEST_HELPER_UTILS_H
#define PANCAKE_TEST_HELPER_UTILS_H

#include <pbbam/FastaSequence.h>
#include <pancake/FastaSequenceId.hpp>
#include <pancake/MapperBatchUtility.hpp>
#include <pancake/OverlapWriterBase.hpp>
#include <string>
#include <vector>

namespace PacBio {
namespace PancakeTests {

std::vector<PacBio::BAM::FastaSequence> HelperLoadFasta(const std::string& inFasta);
std::vector<PacBio::Pancake::FastaSequenceId> HelperLoadFastaWithId(const std::string& inFasta,
                                                                    const int32_t seqIdOffset);
std::vector<std::string> HelperLoadFastaAsStringVector(const std::string& inFasta);
std::string HelperLoadFastaAsString(const std::string& inFasta);
std::vector<std::string> HelperLoadFile(const std::string& inFile);

void HelperLoadBatchData(
    const std::vector<std::pair<std::string, std::string>>& batchDataSequenceFiles,
    const int32_t seqIdOffset, const double freqPercentile,
    const PacBio::Pancake::SeedDBParameters& seedParamsPrimary,
    const PacBio::Pancake::SeedDBParameters& seedParamsFallback,
    std::vector<PacBio::Pancake::MapperBatchChunk>& retBatchData,
    std::vector<PacBio::BAM::FastaSequence>& retAllSeqs);

void HelperLoadBatchData(
    const std::vector<std::tuple<std::string, std::string, PacBio::Pancake::MapperCLRMapSettings>>&
        batchDataSequenceFiles,
    const int32_t seqIdOffset, std::vector<PacBio::Pancake::MapperBatchChunk>& retBatchData,
    std::vector<PacBio::BAM::FastaSequence>& retAllSeqs);

void HelperLoadBatchData(
    const std::vector<std::pair<std::string, std::string>>& batchDataSequenceFiles,
    const int32_t seqIdOffset, const PacBio::Pancake::MapperCLRMapSettings& mapSettings,
    std::vector<PacBio::Pancake::MapperBatchChunk>& retBatchData,
    std::vector<PacBio::BAM::FastaSequence>& retAllSeqs);

std::vector<std::vector<std::string>> HelperFormatBatchMappingResults(
    const std::vector<std::vector<PacBio::Pancake::MapperBaseResult>>& results);

std::tuple<std::vector<std::vector<std::string>>, std::vector<std::vector<std::string>>>
HelperFormatBatchMappingResultsForMockingFlags(
    const std::vector<std::vector<PacBio::Pancake::MapperBaseResult>>& results,
    const std::vector<std::tuple<std::string, std::string, PacBio::Pancake::MapperCLRMapSettings>>&
        batchData,
    const PacBio::Pancake::MapperSelfHitPolicy alignSelfHitPolicy);

PacBio::Pancake::MapperCLRSettings HelperInitPancakeSettingsSubread();

}  // namespace PancakeTests
}  // namespace PacBio

#endif  // PANCAKE_TEST_HELPER_UTILS_H