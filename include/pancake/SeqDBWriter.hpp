// Author: Ivan Sovic

#ifndef PANCAKE_SEQDB_WRITER_COMPRESSED_H
#define PANCAKE_SEQDB_WRITER_COMPRESSED_H

#include <cstdint>
#include <memory>
#include <pancake/CompressedSequence.hpp>
#include <pancake/SeqDBIndexCache.hpp>
#include <pancake/SeqDBWriterBase.hpp>
#include <pancake/util/Util.hpp>
#include <string>

/*
    # Sequence DB
    1. Metadata file: <prefix>.seqdb
    2. One or more sequence files: <prefix>.<file_id>.seq

    ## Metadata file:
    Text file containing the following fields:
    ```
    V <string:semantic_version>
    C <int32_t:compression_level>
    F <int32_t:file_id> <string:filename> <int32_t:num_seqs> <int64_t:file_size_in_bytes> <int64_t:num_compressed_bases>
    S <int32_t:seq_id> <string:header> <int32_t:file_id> <int64_t:file_offset> <int32_t:byte_size> <int32_t:num_bases> <int32_t:num_ranges> <int64_t:start_1> <int64_t:end_1> [<int64_t:start_2> <int64_t:end_2> ...]
    ```

    ## Sequence file:
    Binary file. Contains all bytes concatenated together, no headers, no new line chars.
    It can be either compressed (2-bit) or uncompressed (1-byte per base).
*/

namespace PacBio {
namespace Pancake {

class SeqDBWriter : SeqDBWriterBase
{
public:
    SeqDBWriter(const std::string& filenamePrefix, bool useCompression, int64_t flushSize,
                int64_t blockSize, bool splitBlocks);
    ~SeqDBWriter() override;

    void AddSequence(const std::string& header, const std::string& seq) override;
    bool WriteSequences() override;
    void WriteIndex() override;
    void ClearSequenceBuffer() override;
    void FlushSequenceBuffer() override;
    void CloseFiles() override;

private:
    void OpenNewSequenceFile_();
    void OpenNewIndexFile_();

    const std::string version_{"0.1.0"};
    std::string filenamePrefix_;
    std::string parentFolder_;
    std::string basenamePrefix_;
    int64_t flushSizeBytes_ = 1024 * 1024 * 1024;  // 1GB
    int64_t blockSize_ = 0;
    bool splitBlocks_ = false;
    std::vector<uint8_t> seqBuffer_;
    SeqDBIndexCache cache_;
    SeqDBBlockLine currentBlock_;
    bool openNewSequenceFileUponNextWrite_ = false;
    int64_t totalOutBytes_ = 0;
    int32_t totalOutSeqs_ = 0;
    // Output file handlers.
    std::unique_ptr<FILE, FileDeleter> fpOutIndex_{nullptr};
    std::string outIndexFilename_;
    std::unique_ptr<FILE, FileDeleter> fpOutSeqs_{nullptr};
};

std::unique_ptr<SeqDBWriter> CreateSeqDBWriter(const std::string& filenamePrefix,
                                               bool useCompression, int64_t flushSize,
                                               int64_t blockSize, bool splitBlocks);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEQDB_WRITER_H