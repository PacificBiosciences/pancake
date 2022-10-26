// Authors: Ivan Sovic

#include <pancake/MapperBase.hpp>

#include <pbcopper/logging/Logging.h>

namespace PacBio {
namespace Pancake {

std::vector<MapperBaseResult> MapperBase::MapAndAlign(const std::vector<std::string>& targetSeqs,
                                                      const std::vector<std::string>& querySeqs)
{
    try {
        std::vector<FastaSequenceCached> targetSeqsCached;
        targetSeqsCached.reserve(targetSeqs.size());
        for (int32_t i = 0; i < Utility::Ssize(targetSeqs); ++i) {
            targetSeqsCached.emplace_back(FastaSequenceCached(
                std::to_string(i), targetSeqs[i].c_str(), targetSeqs[i].size(), i));
        }

        std::vector<FastaSequenceCached> querySeqsCached;
        querySeqsCached.reserve(querySeqs.size());
        for (int32_t i = 0; i < Utility::Ssize(querySeqs); ++i) {
            querySeqsCached.emplace_back(FastaSequenceCached(
                std::to_string(i), querySeqs[i].c_str(), querySeqs[i].size(), i));
        }

        return MapAndAlign(targetSeqsCached, querySeqsCached);
    } catch (const std::exception& e) {
        // Log, but do not fail. Important for clients of this class.
        // Return a vector of the size of the input, but with empty results for each query.
        PBLOG_DEBUG << "MapperBase caught an exception in " << std::string(__FUNCTION__)
                    << ". Message: " << e.what();
        return std::vector<MapperBaseResult>(querySeqs.size());
    }
}

std::vector<MapperBaseResult> MapperBase::MapAndAlign(
    const std::vector<FastaSequenceCached>& targetSeqs,
    const std::vector<FastaSequenceCached>& querySeqs)
{
    try {
        PacBio::Pancake::FastaSequenceCachedStore targetSeqsStore(targetSeqs);
        PacBio::Pancake::FastaSequenceCachedStore querySeqsStore(querySeqs);
        return MapAndAlign(targetSeqsStore, querySeqsStore);
    } catch (const std::exception& e) {
        // Log, but do not fail. Important for clients of this class.
        // Return a vector of the size of the input, but with empty results for each query.
        PBLOG_DEBUG << "MapperBase caught an exception in " << std::string(__FUNCTION__)
                    << ". Message: " << e.what();
        return std::vector<MapperBaseResult>(querySeqs.size());
    }
}

}  // namespace Pancake
}  // namespace PacBio