// Authors: Ivan Sovic

#include <pancake/SeedDBReaderCachedBlock.hpp>

#include <pancake/SeedDBReader.hpp>

#include <iostream>
#include <sstream>

namespace PacBio {
namespace Pancake {

SeedDBReaderCachedBlock::SeedDBReaderCachedBlock(
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache>& seedDBCache)
    : indexCache_(seedDBCache)
{
    // Sanity check.
    ValidateSeedDBIndexCache(seedDBCache);
}

SeedDBReaderCachedBlock::SeedDBReaderCachedBlock(
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache>& seedDBCache,
    const std::vector<int32_t>& blockIds)
    : indexCache_(seedDBCache), blockIds_(blockIds)
{
    // Sanity check.
    if (indexCache_->fileLines.empty())
        throw std::runtime_error("There are no file specifications in the input index file.");
    if (indexCache_->seedLines.empty())
        throw std::runtime_error("There are no sequences in the input index file.");
    if (indexCache_->blockLines.empty())
        throw std::runtime_error("There are no blocks in the input index file.");

    LoadBlock(blockIds);
}

SeedDBReaderCachedBlock::~SeedDBReaderCachedBlock() = default;

void SeedDBReaderCachedBlock::LoadBlock(const std::vector<int32_t>& blockIds)
{
    blockIds_ = blockIds;

    // Form the contiguous blocks for loading.
    // First, get all the spans for the first block. These will be stored directly in the
    // accumulator called 'parts'.
    std::vector<ContiguousFilePart> parts = GetSeedDBContiguousParts(indexCache_, blockIds_[0]);
    if (parts.empty()) {
        std::runtime_error(
            "Unknown error occurred, the GetSeqDBContiguousParts returned empty in "
            "LoadBlockUncompressed_.");
    }
    // Next, get all other parts for all other blocks. If possible, extend
    // the previously added parts, otherwise just append.
    for (size_t i = 1; i < blockIds_.size(); ++i) {
        const std::vector<ContiguousFilePart> newParts =
            GetSeedDBContiguousParts(indexCache_, blockIds_[i]);
        for (const auto& newPart : newParts) {
            if (newPart.CanAppendTo(parts.back())) {
                parts.back().ExtendWith(newPart);
            } else {
                parts.emplace_back(newPart);
            }
        }
    }

    // Preallocate the data size.
    int64_t totalSize = 0;
    for (const auto& blockId : blockIds_) {
        const auto& bl = indexCache_->GetBlockLine(blockId);
        totalSize += bl.numBytes;
    }
    totalSize /= 16;
    data_.resize(totalSize);

    // Preallocate the space for all the records.
    const auto& blFirst = indexCache_->GetBlockLine(blockIds_.front());
    const auto& blLast = indexCache_->GetBlockLine(blockIds_.back());
    records_.resize(blLast.endSeqId - blFirst.startSeqId);

    int64_t currDataPos = 0;
    int64_t currRecord = 0;
    for (const auto& part : parts) {
        // Open the file and position to the correct offset.
        const auto& fl = indexCache_->GetFileLine(part.fileId);
        const std::string actualPath = JoinPath(indexCache_->indexParentFolder, fl.filename);
        std::unique_ptr<FILE, FileDeleter> fp = PacBio::Pancake::OpenFile(actualPath.c_str(), "rb");
        const int32_t rv = fseek(fp.get(), part.startOffset, SEEK_SET);
        if (rv) {
            throw std::runtime_error("Could not fseek to position: " +
                                     std::to_string(part.startOffset));
        }

        // Load the bytes.
        const int64_t itemsToRead = (part.endOffset - part.startOffset) / 16;
        const int64_t numItemsRead =
            fread(&data_[currDataPos], sizeof(PacBio::Pancake::Int128t), itemsToRead, fp.get());
        if (itemsToRead != numItemsRead) {
            std::ostringstream oss;
            oss << "(SeedDBReaderCachedBlock) Could not read data for the following part: "
                << "fileId = " << part.fileId << ", offsetStart = " << part.startOffset
                << ", offsetEnd = " << part.endOffset << ", firstId = " << part.seqIds.front()
                << ", lastId = " << part.seqIds.back() << ", itemsToRead = " << itemsToRead
                << ", numItemsRead = " << numItemsRead;
            throw std::runtime_error(oss.str());
        }

        // Create the records, and add them to the lookups.
        int64_t seqStart = currDataPos;
        for (const auto& id : part.seqIds) {
            const auto& sl = indexCache_->GetSeedsLine(id);
            records_[currRecord] =
                SequenceSeedsCached(sl.header, &data_[seqStart], sl.numSeeds, sl.seqId);
            headerToOrdinalId_[sl.header] = currRecord;
            seqIdToOrdinalId_[sl.seqId] = currRecord;
            seqStart += sl.numSeeds;
            ++currRecord;
        }

        // Increment the storage location for the next part.
        currDataPos += numItemsRead;
    }
}

const SequenceSeedsCached& SeedDBReaderCachedBlock::GetSeedsForSequence(int32_t seqId) const
{
    const auto it = seqIdToOrdinalId_.find(seqId);
    if (it == seqIdToOrdinalId_.end()) {
        std::ostringstream oss;
        oss << "(SeedDBReaderCachedBlock) Invalid seqId, not found in blocks: {";
        for (const auto& blockId : blockIds_) {
            oss << blockId << ", ";
        }
        oss << "}. seqId = " << seqId << ", records_.size() = " << records_.size();
        throw std::runtime_error(oss.str());
    }
    const int32_t ordinalId = it->second;
    return records_[ordinalId];
}

const SequenceSeedsCached& SeedDBReaderCachedBlock::GetSeedsForSequence(
    const std::string& seqName) const
{
    const auto it = headerToOrdinalId_.find(seqName);
    if (it == headerToOrdinalId_.end()) {
        std::ostringstream oss;
        oss << "(SeedDBReaderCachedBlock) Invalid seqName, not found in blocks: {";
        for (const auto& blockId : blockIds_) {
            oss << blockId << ", ";
        }
        oss << "}. seqName = '" << seqName << ".";
        throw std::runtime_error(oss.str());
    }
    const int32_t ordinalId = it->second;
    return records_[ordinalId];
}

}  // namespace Pancake
}  // namespace PacBio
