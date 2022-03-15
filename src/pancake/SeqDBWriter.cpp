// Authors: Ivan Sovic

#include <pancake/SeqDBWriter.hpp>

#include <pancake/CompressedSequence.hpp>
#include <pancake/util/Util.hpp>

#include <cmath>
#include <iostream>
#include <sstream>

namespace PacBio {
namespace Pancake {

std::unique_ptr<SeqDBWriter> CreateSeqDBWriter(const std::string& filenamePrefix,
                                               bool useCompression, int64_t flushSize,
                                               int64_t blockSize, bool splitBlocks)
{
    return std::make_unique<SeqDBWriter>(filenamePrefix, useCompression, flushSize, blockSize,
                                         splitBlocks);
}

SeqDBWriter::SeqDBWriter(const std::string& filenamePrefix, bool useCompression, int64_t flushSize,
                         int64_t blockSize, bool splitBlocks)
    : filenamePrefix_(filenamePrefix)
    , flushSizeBytes_(flushSize)
    , blockSize_(blockSize)
    , splitBlocks_(splitBlocks)
{
    cache_.version = version_;
    cache_.compressionLevel = useCompression;

    // Perform the actual throwable work here, so that the constructor doesn't throw.
    SplitPath(filenamePrefix_, parentFolder_, basenamePrefix_);
    seqBuffer_.reserve(flushSizeBytes_);
    OpenNewSequenceFile_();
    OpenNewIndexFile_();
}

SeqDBWriter::~SeqDBWriter()
{
    if (currentBlock_.Span() > 0) {
        cache_.blockLines.emplace_back(currentBlock_);
    }
    FlushSequenceBuffer();
    WriteIndex();
}

void SeqDBWriter::AddSequence(const std::string& header, const std::string& seq)
{
    // We need access to the open file to count the number of sequences, bases and bytes.
    if (cache_.fileLines.empty()) {
        throw std::runtime_error("There are no output sequence files open.");
    }

    // Only open a new file before writing to it. Otherwise, we'll always end up
    // with an extra empty file at the end of each block.
    if (openNewSequenceFileUponNextWrite_) {
        OpenNewSequenceFile_();
    }

    openNewSequenceFileUponNextWrite_ = false;

    int32_t numBytes = 0;
    std::vector<Range> ranges;
    int64_t numUncompressedBases = 0;
    int64_t numCompressedBases = 0;

    // Add the bases (either compressed or uncompressed), and initialize the
    // byte and length values properly.
    if (cache_.compressionLevel) {
        // Compress the sequence.
        const PacBio::Pancake::CompressedSequence compressed =
            PacBio::Pancake::CompressedSequence(seq);
        const auto& bytes = compressed.GetTwobit();
        // Add the compressed bytes to the buffer.
        seqBuffer_.insert(seqBuffer_.end(), bytes.begin(), bytes.end());
        // Set the numeric values used for the index.
        numBytes = static_cast<int32_t>(bytes.size());
        ranges = compressed.GetRanges();
        numUncompressedBases = compressed.GetNumUncompressedBases();
        numCompressedBases = compressed.GetNumCompressedBases();

    } else {
        // Convert the bytes to the internal type.
        const uint8_t* bytes = reinterpret_cast<const uint8_t*>(seq.data());
        // Add the bytes to the buffer.
        seqBuffer_.insert(seqBuffer_.end(), bytes, bytes + seq.size());
        // Set the numeric values used for the index.
        numBytes = static_cast<int32_t>(seq.size());
        ranges = {Range{0, static_cast<int32_t>(seq.size())}};
        numUncompressedBases = numBytes;
        numCompressedBases = numBytes;
    }

    // Create a new index registry object.
    SeqDBSequenceLine sl;
    sl.seqId = static_cast<int32_t>(cache_.seqLines.size());
    sl.header = header;
    sl.numBytes = numBytes;
    sl.numBases = static_cast<int32_t>(seq.size());
    sl.fileId = cache_.fileLines.back().fileId;
    sl.fileOffset = cache_.fileLines.back().numBytes;
    sl.ranges = ranges;
    cache_.seqLines.emplace_back(sl);

    // Extend the current block.
    currentBlock_.startSeqId = (currentBlock_.startSeqId < 0) ? sl.seqId : currentBlock_.startSeqId;
    currentBlock_.endSeqId = sl.seqId + 1;
    currentBlock_.numBytes += sl.numBytes;
    currentBlock_.numBases += sl.numBases;

    // Increase counts for the current file.
    cache_.fileLines.back().numBytes += numBytes;
    ++cache_.fileLines.back().numSequences;
    cache_.fileLines.back().numUncompressedBases += numUncompressedBases;
    cache_.fileLines.back().numCompressedBases += numCompressedBases;

    // Increase the global counts.
    ++totalOutSeqs_;
    totalOutBytes_ += numBytes;

    // Flush the sequences if we reached the buffer size.
    if (static_cast<int64_t>(seqBuffer_.size()) > flushSizeBytes_) {
        FlushSequenceBuffer();
    }

    if (currentBlock_.numBases >= blockSize_ && currentBlock_.Span() > 0) {
        FlushSequenceBuffer();
        cache_.blockLines.emplace_back(currentBlock_);
        currentBlock_ = SeqDBBlockLine();
        currentBlock_.blockId = static_cast<int32_t>(cache_.blockLines.size());
        openNewSequenceFileUponNextWrite_ = false;
        if (splitBlocks_) {
            openNewSequenceFileUponNextWrite_ = true;
        }
    }
}

void SeqDBWriter::OpenNewIndexFile_()
{
    outIndexFilename_ = filenamePrefix_ + ".seqdb";
    fpOutIndex_ = PacBio::Pancake::OpenFile(outIndexFilename_.c_str(), "w");
}

void SeqDBWriter::OpenNewSequenceFile_()
{
    // Register a new file object.
    SeqDBFileLine fileLine;
    fileLine.fileId = static_cast<int32_t>(cache_.fileLines.size());
    fileLine.filename =
        basenamePrefix_ + ".seqdb." + std::to_string(cache_.fileLines.size()) + ".seq";
    cache_.fileLines.emplace_back(fileLine);

    // Open the new file and return the pointer.
    fpOutSeqs_ =
        PacBio::Pancake::OpenFile(JoinPath(parentFolder_, fileLine.filename).c_str(), "wb");
}

void SeqDBWriter::FlushSequenceBuffer()
{
    WriteSequences();
    ClearSequenceBuffer();
}

bool SeqDBWriter::WriteSequences()
{
    // An output sequence file should be open at all times, starting from construction.
    if (fpOutSeqs_ == nullptr) {
        throw std::runtime_error("Cannot write sequences because a sequence file is not open.");
    }

    // Write the actual sequences.
    const size_t num =
        fwrite(seqBuffer_.data(), sizeof(int8_t), seqBuffer_.size(), fpOutSeqs_.get());
    return num == seqBuffer_.size();
}

void SeqDBWriter::WriteIndex() { WriteSeqDBIndexCache(fpOutIndex_.get(), cache_); }

void SeqDBWriter::ClearSequenceBuffer() { seqBuffer_.clear(); }

void SeqDBWriter::CloseFiles()
{
    fpOutIndex_ = nullptr;
    fpOutSeqs_ = nullptr;
}

}  // namespace Pancake
}  // namespace PacBio
