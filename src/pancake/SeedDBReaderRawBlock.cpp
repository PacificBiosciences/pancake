// Authors: Ivan Sovic

#include <pancake/SeedDBReaderRawBlock.hpp>

#include <pancake/SeedDBReader.hpp>

#include <functional>
#include <iostream>
#include <sstream>

namespace PacBio {
namespace Pancake {

SeedDBReaderRawBlock::SeedDBReaderRawBlock(
    const std::shared_ptr<PacBio::Pancake::SeedDBIndexCache>& seedDBCache)
    : seedDBIndexCache_(seedDBCache)
{
    // Sanity check.
    ValidateSeedDBIndexCache(seedDBCache);
}

std::vector<SeedRaw> SeedDBReaderRawBlock::GetBlock(int32_t blockId) const
{
    // Sanity check for the sequence ID.
    if (blockId < 0 || blockId >= static_cast<int32_t>(seedDBIndexCache_->blockLines.size())) {
        std::ostringstream oss;
        oss << "Invalid blockId (SeedDBReader). blockId = " << blockId
            << ", blocks.size() = " << seedDBIndexCache_->blockLines.size();
        throw std::runtime_error(oss.str());
    }

    // Sequences in the block might not be stored contiguously in the file,
    // for example if a user has permuted or filtered the DB.
    // We will collect all contiguous stretches of bytes here, and then
    // fetch those parts later.
    const std::vector<PacBio::Pancake::ContiguousFilePart> contiguousParts =
        GetSeedDBContiguousParts(seedDBIndexCache_, blockId);

    int64_t totalBytes = 0;
    for (const auto& part : contiguousParts) {
        totalBytes += (part.endOffset - part.startOffset);
    }
    const int64_t totalItems = totalBytes / 16;

    // Preallocate the space for the loaded data.
    std::vector<SeedRaw> ret(totalItems);

    OpenFileHandler fh;
    int64_t readItems = 0;
    for (const auto& part : contiguousParts) {
        // Open a new file if required.
        if (part.fileId != fh.fileId) {
            if (part.fileId < 0 ||
                part.fileId >= static_cast<int32_t>(seedDBIndexCache_->fileLines.size())) {
                throw std::runtime_error("(SeedDBReaderRawBlock) Invalid fileId value: " +
                                         std::to_string(part.fileId));
            }
            const auto& fl = seedDBIndexCache_->fileLines[part.fileId];
            const std::string actualPath =
                JoinPath(seedDBIndexCache_->indexParentFolder, fl.filename);
            fh.fp = PacBio::Pancake::OpenFile(actualPath, "rb");
            fh.fileId = part.fileId;
            fh.pos = 0;
        }
        // Jump to a different offset in the file if required.
        if (part.startOffset != fh.pos) {
            const int32_t rv = std::fseek(fh.fp.get(), part.startOffset, SEEK_SET);
            if (rv) {
                throw std::runtime_error("(SeedDBReaderRawBlock) Could not fseek to position: " +
                                         std::to_string(part.startOffset));
            }
            fh.pos = part.startOffset;
        }

        // Load the bytes.
        const int64_t itemsToRead = (part.endOffset - part.startOffset) / 16;
        const int64_t n = std::fread(&ret[readItems], sizeof(SeedRaw), itemsToRead, fh.fp.get());
        readItems += n;

        // Sanity check.
        if (n != itemsToRead) {
            std::ostringstream oss;
            oss << "(SeedDBReaderRawBlock) Could not read seeds for the following part: "
                << "fileId = " << part.fileId << ", offsetStart = " << part.startOffset
                << ", offsetEnd = " << part.endOffset;
            throw std::runtime_error(oss.str());
        }
    }

    return ret;
}

}  // namespace Pancake
}  // namespace PacBio
