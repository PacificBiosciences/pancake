// Author: Ivan Sovic

#ifndef PANCAKE_SEED_DB_READER_HPP
#define PANCAKE_SEED_DB_READER_HPP

#include <pancake/SeedDBIndexCache.hpp>
#include <pancake/SequenceSeeds.hpp>
#include <pancake/util/CommonTypes.hpp>
#include <pancake/util/Util.hpp>

#include <memory>
#include <ostream>
#include <string>

namespace PacBio {
namespace Pancake {

class SeedDBReader
{
public:
    SeedDBReader(std::shared_ptr<PacBio::Pancake::SeedDBIndexCache>& seedDBCache);
    ~SeedDBReader();

    bool GetSeedsForSequence(SequenceSeeds& record, int64_t seqId);
    bool GetSeedsForSequence(SequenceSeeds& record, const std::string& seqName);
    bool GetNext(SequenceSeeds& record);

    bool GetNextBatch(std::vector<SequenceSeeds>& records, int64_t batchSize);
    bool GetBlock(std::vector<SequenceSeeds>& records, int32_t blockId);

    bool JumpTo(int64_t seqId);
    bool JumpTo(const std::string& seqName);

private:
    using FilePtr = std::unique_ptr<FILE, FileDeleter>;
    class OpenFileHandler
    {
    public:
        FilePtr fp = nullptr;
        int32_t fileId = -1;
        int64_t pos = 0;
        int32_t nextOrdinalId = 0;
    };
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBIndexCache_;
    OpenFileHandler fileHandler_;
    HeaderLookupType headerToOrdinalId_;

    void AccessLocation_(OpenFileHandler& fileHandler,
                         const std::vector<PacBio::Pancake::SeedDBFileLine>& fileLines,
                         const std::string& indexParentFolder, int32_t fileId,
                         int32_t nextOrdinalId, int64_t offset) const;

    void LoadSeedsForSequence_(Pancake::SequenceSeeds& record, OpenFileHandler& fileHandler,
                               const std::vector<PacBio::Pancake::SeedDBFileLine>& fileLines,
                               const std::string& indexParentFolder, const SeedDBSeedsLine& sl,
                               int32_t ordinalId) const;
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEED_DB_READER_HPP
