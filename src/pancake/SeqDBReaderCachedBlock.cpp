// Authors: Ivan Sovic

#include <pancake/SeqDBReaderCachedBlock.hpp>

#include <pancake/SeqDBReader.hpp>
#include <pancake/Twobit.hpp>
#include <pancake/util/RunLengthEncoding.hpp>
#include <pancake/util/Util.hpp>

#include <iostream>
#include <sstream>

namespace PacBio {
namespace Pancake {

SeqDBReaderCachedBlock::SeqDBReaderCachedBlock(
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache>& seqDBCache, bool useHomopolymerCompression)
    : seqDBIndexCache_(seqDBCache), useHomopolymerCompression_(useHomopolymerCompression)
{
    ValidateSeqDBIndexCache(seqDBCache);
}

SeqDBReaderCachedBlock::~SeqDBReaderCachedBlock() = default;

void SeqDBReaderCachedBlock::LoadBlocks(const std::vector<int32_t>& blockIds)
{
    // Collect all sequence IDs for all blocks.
    std::vector<int32_t> seqIds;
    for (const auto& blockId : blockIds) {
        const auto& bl = seqDBIndexCache_->GetBlockLine(blockId);
        std::vector<int32_t> newSpan(bl.endSeqId - bl.startSeqId);
        std::iota(newSpan.begin(), newSpan.end(), bl.startSeqId);
        seqIds.insert(seqIds.end(), newSpan.begin(), newSpan.end());
    }

    // Get the contiguous file parts for loading.
    std::vector<ContiguousFilePart> parts = GetSeqDBContiguousParts(seqDBIndexCache_, seqIds);

    // Actually load the data.
    if (seqDBIndexCache_->compressionLevel == 0) {
        return LoadBlockUncompressed_(parts);
    }
    return LoadBlockCompressed_(parts);
}

void SeqDBReaderCachedBlock::LoadSequences(const std::vector<int32_t>& seqIds)
{
    // Get the contiguous file parts for loading.
    std::vector<ContiguousFilePart> parts = GetSeqDBContiguousParts(seqDBIndexCache_, seqIds);

    // Actually load the data.
    if (seqDBIndexCache_->compressionLevel == 0) {
        return LoadBlockUncompressed_(parts);
    }
    return LoadBlockCompressed_(parts);
}

void SeqDBReaderCachedBlock::LoadSequences(const std::vector<std::string>& seqNames)
{
    // Get the contiguous file parts for loading.
    std::vector<ContiguousFilePart> parts = GetSeqDBContiguousParts(seqDBIndexCache_, seqNames);

    // Actually load the data.
    if (seqDBIndexCache_->compressionLevel == 0) {
        return LoadBlockUncompressed_(parts);
    }
    return LoadBlockCompressed_(parts);
}

void SeqDBReaderCachedBlock::LoadBlockCompressed_(const std::vector<ContiguousFilePart>& parts)
{
    // Count the data size.
    int64_t totalBases = 0;
    for (const auto& part : parts) {
        for (const auto& sId : part.seqIds) {
            const auto& sl = seqDBIndexCache_->GetSeqLine(sId);
            totalBases += sl.numBases;
        }
    }

    // Preallocate the space for all the records.
    data_.resize(totalBases);
    recordStore_.Clear();

    // Position of a current record in the data_ vector.
    int64_t seqStart = 0;
    for (const auto& part : parts) {
        // Open the file and position to the correct offset.
        const auto& fl = seqDBIndexCache_->GetFileLine(part.fileId);
        const std::string actualPath = JoinPath(seqDBIndexCache_->indexParentFolder, fl.filename);
        std::unique_ptr<std::FILE, FileDeleter> fp = PacBio::Pancake::OpenFile(actualPath, "rb");
        const int32_t rv = std::fseek(fp.get(), part.startOffset, SEEK_SET);
        if (rv) {
            throw std::runtime_error("Could not fseek to position: " +
                                     std::to_string(part.startOffset));
        }

        // Load the bytes.
        const int64_t itemsToRead = (part.endOffset - part.startOffset);
        std::vector<uint8_t> tempData;
        tempData.resize(itemsToRead);
        const int64_t numItemsRead =
            std::fread(tempData.data(), sizeof(uint8_t), itemsToRead, fp.get());
        if (itemsToRead != numItemsRead) {
            std::ostringstream oss;
            oss << "(SeqDBReaderCachedBlock) Could not read data for the following part: "
                << "fileId = " << part.fileId << ", offsetStart = " << part.startOffset
                << ", offsetEnd = " << part.endOffset << ", frontId = " << part.seqIds.front()
                << ", backId = " << part.seqIds.back() << ", itemsToRead = " << itemsToRead
                << ", numItemsRead = " << numItemsRead;
            throw std::runtime_error(oss.str());
        }

        // File offset of the first record in the part. Needed to normalize
        // the access to the tempData vector.
        int64_t startByteOffset = part.startOffset;

        // Decompress and create the records, and add them to the lookups.
        for (const auto& id : part.seqIds) {
            const auto& sl = seqDBIndexCache_->GetSeqLine(id);
            const int64_t firstByte = sl.fileOffset - startByteOffset;
            DecompressSequenceCStyle({&tempData[firstByte], static_cast<size_t>(sl.numBytes)},
                                     sl.numBases, sl.ranges, &data_[seqStart]);
            recordStore_.AddRecord(
                FastaSequenceCached{sl.header, &data_[seqStart], sl.numBases, sl.seqId});
            seqStart += sl.numBases;
        }
    }

    if (useHomopolymerCompression_) {
        CompressHomopolymers_();
    }
}

void SeqDBReaderCachedBlock::LoadBlockUncompressed_(const std::vector<ContiguousFilePart>& parts)
{
    // Count the data size.
    int64_t totalBases = 0;
    for (const auto& part : parts) {
        for (const auto& sId : part.seqIds) {
            const auto& sl = seqDBIndexCache_->GetSeqLine(sId);
            totalBases += sl.numBases;
        }
    }

    // Preallocate the space for all the records.
    data_.resize(totalBases);
    recordStore_.Clear();

    int64_t currDataPos = 0;
    for (const auto& part : parts) {
        // Open the file and position to the correct offset.
        const auto& fl = seqDBIndexCache_->GetFileLine(part.fileId);
        const std::string actualPath = JoinPath(seqDBIndexCache_->indexParentFolder, fl.filename);
        std::unique_ptr<std::FILE, FileDeleter> fp = PacBio::Pancake::OpenFile(actualPath, "rb");
        const int32_t rv = std::fseek(fp.get(), part.startOffset, SEEK_SET);
        if (rv) {
            throw std::runtime_error("Could not fseek to position: " +
                                     std::to_string(part.startOffset));
        }

        // Load the bytes.
        const int64_t itemsToRead = (part.endOffset - part.startOffset);
        const int64_t numItemsRead =
            std::fread(&data_[currDataPos], sizeof(uint8_t), itemsToRead, fp.get());
        if (itemsToRead != numItemsRead) {
            std::ostringstream oss;
            oss << "(SeqDBReaderCachedBlock) Could not read data for the following part: "
                << "fileId = " << part.fileId << ", offsetStart = " << part.startOffset
                << ", offsetEnd = " << part.endOffset << ", firstId = " << part.seqIds.front()
                << ", lastId = " << part.seqIds.back() << ", itemsToRead = " << itemsToRead
                << ", numItemsRead = " << numItemsRead;
            throw std::runtime_error(oss.str());
        }

        // Create the records, and add them to the lookups.
        int64_t seqStart = currDataPos;
        for (const auto& id : part.seqIds) {
            const auto& sl = seqDBIndexCache_->GetSeqLine(id);
            recordStore_.AddRecord(
                FastaSequenceCached{sl.header, &data_[seqStart], sl.numBases, sl.seqId});
            seqStart += sl.numBases;
        }

        // Increment the storage location for the next part.
        currDataPos += numItemsRead;
    }

    if (useHomopolymerCompression_) {
        CompressHomopolymers_();
    }
}

void SeqDBReaderCachedBlock::CompressHomopolymers_()
{
    std::vector<int32_t> runLengths;
    for (size_t i = 0; i < recordStore_.records().size(); ++i) {
        auto& record = recordStore_.records()[i];
        int64_t comprLen = PacBio::Pancake::RunLengthEncoding(const_cast<char*>(record.c_str()),
                                                              record.size(), runLengths);
        record.Size(comprLen);
    }
}

}  // namespace Pancake
}  // namespace PacBio
