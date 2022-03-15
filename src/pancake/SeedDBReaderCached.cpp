// Authors: Ivan Sovic

#include <pancake/SeedDBReader.hpp>

#include <pancake/SeedDBReaderCached.hpp>

#include <sstream>

namespace PacBio {
namespace Pancake {

SeedDBReaderCached::SeedDBReaderCached(
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache>& seedDBCache, int32_t blockId)
    : seedDBIndexCache_(seedDBCache), blockId_(blockId)
{
    // Sanity check.
    ValidateSeedDBIndexCache(seedDBCache);

    // Create a SeedDB reader.
    PacBio::Pancake::SeedDBReader reader(seedDBCache);

    // Fetch the block of data.
    reader.GetBlock(records_, blockId);

    for (int32_t i = 0; i < static_cast<int32_t>(records_.size()); ++i) {
        headerToOrdinalId_[records_[i].Name()] = i;
        seqIdToOrdinalId_[records_[i].Id()] = i;
    }
}

SeedDBReaderCached::~SeedDBReaderCached() = default;

const SequenceSeeds& SeedDBReaderCached::GetSeedsForSequence(int32_t seqId) const
{
    const auto it = seqIdToOrdinalId_.find(seqId);
    if (it == seqIdToOrdinalId_.end()) {
        std::ostringstream oss;
        oss << "(SeedDBReaderCached) Invalid seqId, not found in block " << blockId_
            << ". seqId = " << seqId << ", records_.size() = " << records_.size();
        throw std::runtime_error(oss.str());
    }
    const int32_t ordinalId = it->second;
    return records_[ordinalId];
}

const SequenceSeeds& SeedDBReaderCached::GetSeedsForSequence(const std::string& seqName) const
{
    const auto it = headerToOrdinalId_.find(seqName);
    if (it == headerToOrdinalId_.end()) {
        std::ostringstream oss;
        oss << "(SeedDBReaderCached) Invalid seqName, not found in block " << blockId_
            << ". seqName = '" << seqName << ".";
        throw std::runtime_error(oss.str());
    }
    const int32_t ordinalId = it->second;
    return records_[ordinalId];
}

}  // namespace Pancake
}  // namespace PacBio
