// Authors: Ivan Sovic

#include <pancake/FastaSequenceCached.hpp>
#include <pancake/SeedDBWriter.hpp>
#include <pancake/util/Util.hpp>

#include <cinttypes>
#include <cmath>
#include <iostream>
#include <sstream>

namespace PacBio {
namespace Pancake {

std::unique_ptr<SeedDBWriter> CreateSeedDBWriter(
    const std::string& filenamePrefix, bool splitBlocks,
    const PacBio::Pancake::SeedDB::SeedDBParameters& params)
{
    return std::make_unique<SeedDBWriter>(filenamePrefix, splitBlocks, params);
}

SeedDBWriter::SeedDBWriter(const std::string& filenamePrefix, bool splitBlocks,
                           const PacBio::Pancake::SeedDB::SeedDBParameters& params)
    : params_(params)
    , filenamePrefix_(filenamePrefix)
    , splitBlocks_(splitBlocks)
    , openNewSeedsFileUponNextWrite_(false)
{
    // Perform the actuall throwable work here, so that the constructor doesn't throw.
    SplitPath(filenamePrefix_, parentFolder_, basenamePrefix_);
    OpenNewSeedsFile_();
    OpenNewIndexFile_();
}

SeedDBWriter::~SeedDBWriter()
{
    MarkBlockEnd();
    WriteIndex();
}

void SeedDBWriter::OpenNewIndexFile_()
{
    outIndexFilename_ = filenamePrefix_ + ".seeddb";
    fpOutIndex_ = PacBio::Pancake::OpenFile(outIndexFilename_.c_str(), "w");
}

void SeedDBWriter::OpenNewSeedsFile_()
{
    // Register a new file object.
    SeedDBFileLine fileLine;
    fileLine.fileId = static_cast<int32_t>(fileLines_.size());
    fileLine.filename = basenamePrefix_ + ".seeddb." + std::to_string(fileLines_.size()) + ".seeds";
    fileLines_.emplace_back(fileLine);

    // Open the new file and return the pointer.
    fpOutSeeds_ =
        PacBio::Pancake::OpenFile(JoinPath(parentFolder_, fileLine.filename).c_str(), "wb");
}

void SeedDBWriter::WriteSeeds(const std::string& seqName, int32_t seqId, int32_t seqLen,
                              const std::vector<PacBio::Pancake::Int128t>& seeds)
{
    int64_t numBytes = static_cast<int64_t>(seeds.size() * 16);

    // Only open a new file before writing to it. Otherwise, we'll always end up
    // with an extra empty file at the end of each block.
    if (openNewSeedsFileUponNextWrite_) {
        OpenNewSeedsFile_();
    }

    openNewSeedsFileUponNextWrite_ = false;

    // Create a new index registry object.
    SeedDBSeedsLine sl;
    sl.seqId = seqId;
    sl.header = seqName;
    sl.fileId = fileLines_.back().fileId;
    sl.fileOffset = fileLines_.back().numBytes;
    sl.numBytes = numBytes;
    sl.numBases = seqLen;
    sl.numSeeds = static_cast<int32_t>(seeds.size());
    seedsLines_.emplace_back(sl);

    currentBlock_.numBytes += numBytes;
    currentBlock_.startSeqId = (currentBlock_.startSeqId < 0) ? seqId : currentBlock_.startSeqId;
    currentBlock_.endSeqId = seqId + 1;

    // Increase counts for the current file.
    fileLines_.back().numBytes += numBytes;
    ++fileLines_.back().numSequences;

    // Write the actual sequences.
    const int64_t numWritten = seeds.empty()
                                   ? 0
                                   : fwrite(reinterpret_cast<const uint8_t*>(seeds.data()),
                                            sizeof(uint8_t), seeds.size() * 16, fpOutSeeds_.get());

    // Sanity check.
    if (numWritten != numBytes) {
        std::ostringstream oss;
        oss << "Could not write all seed bytes. Written " << numWritten << " of " << numBytes
            << " bytes.";
        throw std::runtime_error(oss.str());
    }
}

void SeedDBWriter::WriteSeeds(const std::vector<PacBio::Pancake::FastaSequenceId>& seqs,
                              const std::vector<std::vector<PacBio::Pancake::Int128t>>& seeds)
{
    if (seqs.size() != seeds.size()) {
        std::ostringstream oss;
        oss << "The seqs and seeds vector are of different lengths (seqs.size() = " << seqs.size()
            << ", seeds.size() = " << seeds.size() << ").";
        throw std::runtime_error(oss.str());
    }

    for (size_t i = 0; i < seqs.size(); ++i) {
        WriteSeeds(seqs[i].Name(), static_cast<int32_t>(seqs[i].Id()),
                   static_cast<int32_t>(seqs[i].Bases().size()), seeds[i]);
    }
}

void SeedDBWriter::WriteSeeds(const std::vector<PacBio::Pancake::FastaSequenceCached>& seqs,
                              const std::vector<std::vector<PacBio::Pancake::Int128t>>& seeds)
{
    if (seqs.size() != seeds.size()) {
        std::ostringstream oss;
        oss << "The seqs and seeds vector are of different lengths (seqs.size() = " << seqs.size()
            << ", seeds.size() = " << seeds.size() << ").";
        throw std::runtime_error(oss.str());
    }

    for (size_t i = 0; i < seqs.size(); ++i) {
        WriteSeeds(seqs[i].Name(), seqs[i].Id(), seqs[i].size(), seeds[i]);
    }
}

void SeedDBWriter::MarkBlockEnd()
{
    if (currentBlock_.endSeqId > currentBlock_.startSeqId) {
        blockLines_.emplace_back(currentBlock_);
        openNewSeedsFileUponNextWrite_ = false;
        if (splitBlocks_) {
            openNewSeedsFileUponNextWrite_ = true;
        }
    }
    currentBlock_ = SeedDBBlockLine();
    currentBlock_.blockId = static_cast<int32_t>(blockLines_.size());
}

void SeedDBWriter::WriteIndex()
{
    // An output index file should be open at all times, starting from construction.
    if (fpOutIndex_ == nullptr) {
        throw std::runtime_error("Cannot write the index because an output file is not open.");
    }

    // Write the version and compression information.
    fprintf(fpOutIndex_.get(), "V\t%s\n", version_.c_str());

    // Write the parameters used to compute the seeds.
    fprintf(fpOutIndex_.get(), "P\tk=%d,w=%d,s=%d,hpc=%d,rc=%d\n", params_.KmerSize,
            params_.MinimizerWindow, params_.Spacing, params_.UseHPC, params_.UseRC);

    // Write all the files and their sizes.
    for (const auto& f : fileLines_) {
        fprintf(fpOutIndex_.get(), "F\t%d\t%s\t%d\t%" PRIi64 "\n", f.fileId, f.filename.c_str(),
                f.numSequences, f.numBytes);
    }

    // Write the indexes of all sequences.
    for (size_t i = 0; i < seedsLines_.size(); ++i) {
        fprintf(fpOutIndex_.get(), "S\t%d\t%s\t%d\t%" PRIi64 "\t%" PRIi64 "\t%d\t%d\n",
                seedsLines_[i].seqId, seedsLines_[i].header.c_str(), seedsLines_[i].fileId,
                seedsLines_[i].fileOffset, seedsLines_[i].numBytes, seedsLines_[i].numBases,
                seedsLines_[i].numSeeds);
    }

    // Write the blocks of all sequences.
    for (size_t i = 0; i < blockLines_.size(); ++i) {
        fprintf(fpOutIndex_.get(), "B\t%d\t%d\t%d\t%" PRIi64 "\n", blockLines_[i].blockId,
                blockLines_[i].startSeqId, blockLines_[i].endSeqId, blockLines_[i].numBytes);
    }
}

void SeedDBWriter::CloseFiles()
{
    fpOutIndex_ = nullptr;
    fpOutSeeds_ = nullptr;
}

}  // namespace Pancake
}  // namespace PacBio
