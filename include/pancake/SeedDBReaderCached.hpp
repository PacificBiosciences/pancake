// Author: Ivan Sovic

#ifndef PANCAKE_SEEDDB_READER_CACHED_H
#define PANCAKE_SEEDDB_READER_CACHED_H

#include <pancake/SeedDBIndexCache.hpp>
#include <pancake/SequenceSeeds.hpp>
#include <pancake/util/Util.hpp>

#include <memory>
#include <ostream>
#include <string>

namespace PacBio {
namespace Pancake {

class SeedDBReaderCached
{
public:
    SeedDBReaderCached(std::shared_ptr<PacBio::Pancake::SeedDBIndexCache>& seedDBCache,
                       int32_t blockId);
    ~SeedDBReaderCached();

    const SequenceSeeds& GetSeedsForSequence(int32_t seqId) const;
    const SequenceSeeds& GetSeedsForSequence(const std::string& seqName) const;
    const std::vector<PacBio::Pancake::SequenceSeeds>& records() const { return records_; }

private:
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBIndexCache_;
    int32_t blockId_;
    std::vector<PacBio::Pancake::SequenceSeeds> records_;

    // Info to allow random access.
    std::unordered_map<std::string, int32_t> headerToOrdinalId_;
    std::unordered_map<int32_t, int32_t> seqIdToOrdinalId_;
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEQDB_READER_CACHED_H