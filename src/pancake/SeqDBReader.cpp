// Authors: Ivan Sovic

#include <pancake/SeqDBReader.hpp>

#include <pancake/CompressedSequence.hpp>
#include <pancake/util/Util.hpp>

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>

namespace PacBio {
namespace Pancake {

SeqDBReader::SeqDBReader(std::shared_ptr<PacBio::Pancake::SeqDBIndexCache>& seqDBCache)
    : seqDBIndexCache_(seqDBCache)
{
    ValidateSeqDBIndexCache(seqDBIndexCache_);
    ComputeSeqDBIndexHeaderLookup(*seqDBIndexCache_, headerToOrdinalId_);
}

SeqDBReader::~SeqDBReader() = default;

const std::shared_ptr<PacBio::Pancake::SeqDBIndexCache>& SeqDBReader::GetSeqDBIndexCache() const
{
    return seqDBIndexCache_;
}

const SeqDBSequenceLine& SeqDBReader::GetSequenceInfo(int64_t seqId) const
{
    return seqDBIndexCache_->GetSeqLine(seqId);
}

const SeqDBSequenceLine& SeqDBReader::GetSequenceInfo(const std::string& seqName) const
{
    return GetSeqLine(*seqDBIndexCache_, headerToOrdinalId_, seqName);
}

bool SeqDBReader::GetSequence(Pancake::FastaSequenceId& record, int64_t seqId)
{
    // Clear the output.
    record.Bases("");
    record.Name("");

    // Access the SequenceLine object.
    const int32_t ordinalId = seqId;
    const auto& sl = seqDBIndexCache_->GetSeqLine(ordinalId);

    // Load the data and decompress if required, and update fileHandler_->nextOrdinalId.
    LoadAndUnpackSequence_(record, fileHandler_, seqDBIndexCache_->fileLines,
                           seqDBIndexCache_->indexParentFolder, sl, ordinalId,
                           seqDBIndexCache_->compressionLevel);

    return true;
}

bool SeqDBReader::GetSequence(Pancake::FastaSequenceId& record, const std::string& seqName)
{
    // Clear the output.
    record.Bases("");
    record.Name("");

    // Find the sequence.
    auto it = headerToOrdinalId_.find(seqName);
    if (it == headerToOrdinalId_.end()) return false;
    int32_t ordinalId = it->second;

    // Access the SequenceLine object.
    const auto& sl = seqDBIndexCache_->seqLines[ordinalId];

    // Load the data and decompress if required, and update fileHandler_->nextOrdinalId.
    LoadAndUnpackSequence_(record, fileHandler_, seqDBIndexCache_->fileLines,
                           seqDBIndexCache_->indexParentFolder, sl, ordinalId,
                           seqDBIndexCache_->compressionLevel);

    return true;
}

bool SeqDBReader::JumpTo(int64_t seqId)
{
    // Access the SequenceLine object.
    const int32_t ordinalId = seqId;
    const auto& sl = seqDBIndexCache_->GetSeqLine(ordinalId);

    // Jump to the correct file and offset, and update the fileHandler_->nextOrdinalId.
    AccessLocation_(fileHandler_, seqDBIndexCache_->fileLines, seqDBIndexCache_->indexParentFolder,
                    sl.fileId, ordinalId, sl.fileOffset);

    return true;
}

bool SeqDBReader::JumpTo(const std::string& seqName)
{
    // Find the sequence.
    auto it = headerToOrdinalId_.find(seqName);
    if (it == headerToOrdinalId_.end())
        throw std::runtime_error(
            "Invalid sequence name, it does not exist in the SeqDB. seqName = " + seqName);
    int32_t ordinalId = it->second;

    // Access the SequenceLine object.
    const auto& sl = seqDBIndexCache_->seqLines[ordinalId];

    // Jump to the correct file and offset, and update the fileHandler_->nextOrdinalId.
    AccessLocation_(fileHandler_, seqDBIndexCache_->fileLines, seqDBIndexCache_->indexParentFolder,
                    sl.fileId, ordinalId, sl.fileOffset);

    return true;
}

bool SeqDBReader::GetNext(Pancake::FastaSequenceId& record)
{
    // Clear the output.
    record.Bases("");
    record.Name("");

    // Sanity check for the sequence ID.
    if (fileHandler_.nextOrdinalId < 0) throw std::runtime_error("Invalid nextSeqId < 0.");

    // Can't go to the next sequence, we loaded all of them.
    if (fileHandler_.nextOrdinalId >= static_cast<int32_t>(seqDBIndexCache_->seqLines.size()))
        return false;

    // Access the SequenceLine object.
    const auto& sl = seqDBIndexCache_->seqLines[fileHandler_.nextOrdinalId];

    // Load the data and decompress if required, and update fileHandler_->nextOrdinalId.
    LoadAndUnpackSequence_(record, fileHandler_, seqDBIndexCache_->fileLines,
                           seqDBIndexCache_->indexParentFolder, sl, fileHandler_.nextOrdinalId,
                           seqDBIndexCache_->compressionLevel);

    return true;
}

bool SeqDBReader::GetNextBatch(std::vector<Pancake::FastaSequenceId>& records, int64_t batchSize)
{
    records.clear();

    // Sanity check for the sequence ID.
    if (fileHandler_.nextOrdinalId < 0) throw std::runtime_error("Invalid nextSeqId < 0.");

    // Batch size < 0 loads everything as one batch.
    batchSize = (batchSize < 0) ? std::numeric_limits<int64_t>::max() : batchSize;

    int32_t numSeqLines = seqDBIndexCache_->seqLines.size();
    int64_t loadedBases = 0;
    while (fileHandler_.nextOrdinalId < numSeqLines) {
        // Access the SequenceLine object.
        const auto& sl = seqDBIndexCache_->seqLines[fileHandler_.nextOrdinalId];

        // Load the data and decompress if required.
        Pancake::FastaSequenceId record;
        // Load the data and decompress if required, and update fileHandler_->nextOrdinalId.
        LoadAndUnpackSequence_(record, fileHandler_, seqDBIndexCache_->fileLines,
                               seqDBIndexCache_->indexParentFolder, sl, fileHandler_.nextOrdinalId,
                               seqDBIndexCache_->compressionLevel);

        // Store the record and increase the counters.
        records.emplace_back(record);
        loadedBases += static_cast<int64_t>(record.Bases().size());

        if (loadedBases >= batchSize) break;
    }

    if (records.empty()) return false;

    return true;
}

bool SeqDBReader::GetBlock(std::vector<Pancake::FastaSequenceId>& records, int32_t blockId)
{
    records.clear();

    // Sanity check for the sequence ID.
    if (blockId < 0 || blockId >= static_cast<int32_t>(seqDBIndexCache_->blockLines.size())) {
        std::ostringstream oss;
        oss << "Invalid blockId (b). blockId = " << blockId
            << ", blocks.size() = " << seqDBIndexCache_->blockLines.size();
        throw std::runtime_error(oss.str());
    }

    const auto& block = seqDBIndexCache_->blockLines[blockId];

    // Sanity check that the block's range is good.
    int32_t numSeqLines = seqDBIndexCache_->seqLines.size();
    if (block.startSeqId < 0 || block.endSeqId <= block.startSeqId ||
        block.startSeqId >= numSeqLines || block.endSeqId > numSeqLines) {
        std::ostringstream oss;
        oss << "Invalid block values: start = " << block.startSeqId << ", end = " << block.endSeqId
            << ", numSeqLines = " << numSeqLines << ".";
        throw std::runtime_error(oss.str());
    }

    fileHandler_.nextOrdinalId = block.startSeqId;

    while (fileHandler_.nextOrdinalId < block.endSeqId) {
        // Access the SequenceLine object.
        const auto& sl = seqDBIndexCache_->seqLines[fileHandler_.nextOrdinalId];

        // Load the data and decompress if required.
        Pancake::FastaSequenceId record;
        // Load the data and decompress if required, and update fileHandler_->nextOrdinalId.
        LoadAndUnpackSequence_(record, fileHandler_, seqDBIndexCache_->fileLines,
                               seqDBIndexCache_->indexParentFolder, sl, fileHandler_.nextOrdinalId,
                               seqDBIndexCache_->compressionLevel);

        // Store the record.
        records.emplace_back(record);
    }

    if (records.empty()) return false;

    return true;
}

void SeqDBReader::AccessLocation_(OpenFileHandler& fileHandler,
                                  const std::vector<PacBio::Pancake::SeqDBFileLine>& fileLines,
                                  const std::string& indexParentFolder, int32_t fileId,
                                  int32_t nextOrdinalId, int64_t offset) const
{
    if (fileId < 0 || fileId >= static_cast<int32_t>(fileLines.size()))
        throw std::runtime_error("Invalid fileId value: " + std::to_string(fileId));
    if (offset < 0)
        throw std::runtime_error("Invalid file offset value: " + std::to_string(offset));

    const auto& fl = fileLines[fileId];

    if (fileHandler.fileId != fileId) {
        std::string actualPath = JoinPath(indexParentFolder, fl.filename);
        fileHandler.fp = PacBio::Pancake::OpenFile(actualPath.c_str(), "rb");
        fileHandler.fileId = fileId;
        fileHandler.pos = 0;
    }
    if (offset != fileHandler.pos) {
        int32_t rv = std::fseek(fileHandler.fp.get(), offset, SEEK_SET);
        if (rv) throw std::runtime_error("Could not fseek to position: " + std::to_string(offset));
        fileHandler.pos = offset;
    }
    fileHandler.nextOrdinalId = nextOrdinalId;
}

void SeqDBReader::LoadAndUnpackSequence_(
    Pancake::FastaSequenceId& record, OpenFileHandler& fileHandler,
    const std::vector<PacBio::Pancake::SeqDBFileLine>& fileLines,
    const std::string& indexParentFolder, const SeqDBSequenceLine& sl, int32_t ordinalId,
    bool isCompressed) const
{

    // Jump to the correct file and offset.
    AccessLocation_(fileHandler, fileLines, indexParentFolder, sl.fileId, ordinalId + 1,
                    sl.fileOffset);

    if (isCompressed) {
        // Load and decompress the sequence.
        std::vector<uint8_t> data(sl.numBytes);
        int32_t n = std::fread(data.data(), sizeof(uint8_t), sl.numBytes, fileHandler.fp.get());
        fileHandler.pos += n;
        if (n != sl.numBytes) {
            std::ostringstream oss;
            oss << "Could not read sequence of '" << sl.header << "'. Num bytes read: " << n
                << ", expected: " << sl.numBytes;
            oss << ", file_id = " << sl.fileId << ", file_offset = " << sl.fileOffset
                << ", ftell = " << std::ftell(fileHandler.fp.get());
            throw std::runtime_error(oss.str());
        }
        std::string bases;
        DecompressSequence(data, sl.numBases, sl.ranges, bases);
        record.Name(sl.header);
        record.Bases(bases);
        record.Id(sl.seqId);

    } else {
        std::string bases(sl.numBytes, '\0');
        int32_t n = std::fread(bases.data(), sizeof(char), sl.numBytes, fileHandler.fp.get());
        fileHandler.pos += n;
        if (n != sl.numBytes) {
            std::ostringstream oss;
            oss << "Could not read sequence of '" << sl.header << "'. Num bytes read: " << n
                << ", expected: " << sl.numBytes;
            throw std::runtime_error(oss.str());
        }
        record.Name(sl.header);
        record.Bases(bases);
        record.Id(sl.seqId);
    }
}

}  // namespace Pancake
}  // namespace PacBio
