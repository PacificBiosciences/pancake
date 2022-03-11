// Author: Ivan Sovic

#ifndef PANCAKE_SEEDDB_READER_RAW_BLOCK_H
#define PANCAKE_SEEDDB_READER_RAW_BLOCK_H

#include <memory>
#include <pancake/ContiguousFilePart.hpp>
#include <pancake/Seed.hpp>
#include <pancake/SeedDBIndexCache.hpp>
#include <pancake/util/Util.hpp>
#include <string>
#include <vector>

namespace PacBio {
namespace Pancake {

class SeedDBReaderRawBlock
{
public:
    SeedDBReaderRawBlock(const std::shared_ptr<PacBio::Pancake::SeedDBIndexCache>& seedDBCache);
    ~SeedDBReaderRawBlock() = default;

    std::vector<SeedDB::SeedRaw> GetBlock(int32_t blockId) const;

private:
    using FilePtr = std::unique_ptr<FILE, FileDeleter>;
    class OpenFileHandler
    {
    public:
        FilePtr fp = nullptr;
        int32_t fileId = -1;
        int64_t pos = -1;
    };

    const std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> seedDBIndexCache_;
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEEDDB_READER_RAW_BLOCK_H