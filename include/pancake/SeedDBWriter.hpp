// Author: Ivan Sovic

#ifndef PANCAKE_SEED_DB_WRITER_HPP
#define PANCAKE_SEED_DB_WRITER_HPP

#include <pancake/FastaSequenceCached.hpp>
#include <pancake/FastaSequenceId.hpp>
#include <pancake/SeedDBIndexCache.hpp>
#include <pancake/SeedDBParameters.hpp>
#include <pancake/util/CommonTypes.hpp>
#include <pancake/util/Util.hpp>

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

/*
    # Seed DB
    1. Metadata file: <prefix>.seeddb
    2. One or more files with seeds: <prefix>.<file_id>.seeds

    ## Metadata file:
    Text file containing the following fields:
    ```
    V <string:semantic_version>
    P <param1=val1,param2=val2,...>
    F <int32_t:file_id> <string:filename> <int32_t:num_seqs> <int64_t:file_size_in_bytes>
    S <int32_t:seq_id> <string:header> <int32_t:file_id> <int64_t:file_offset> <int64_t:byte_size> <int32_t:num_bases> <int32_t:num_seeds>
    B <int32_t:block_id> <int32_t:start_seq_id> <int32_t:end_seq_id> <int64_t:byte_size>
    ```

    ## Seed file:
    Binary file. Contains all bytes concatenated together, no headers, no new line chars.
*/

namespace PacBio {
namespace Pancake {

class SeedDBWriter
{
public:
    SeedDBWriter(const std::string& filenamePrefix, bool splitBlocks,
                 const PacBio::Pancake::SeedDB::SeedDBParameters& params);
    ~SeedDBWriter();

    void WriteSeeds(const std::string& seqName, int32_t seqId, int32_t seqLen,
                    const std::vector<PacBio::Pancake::Int128t>& seeds);
    void WriteSeeds(const std::vector<PacBio::Pancake::FastaSequenceId>& seqs,
                    const std::vector<std::vector<PacBio::Pancake::Int128t>>& seeds);
    void WriteSeeds(const std::vector<PacBio::Pancake::FastaSequenceCached>& seqs,
                    const std::vector<std::vector<PacBio::Pancake::Int128t>>& seeds);
    void MarkBlockEnd();
    void WriteIndex();
    void CloseFiles();

private:
    void OpenNewSeedsFile_();
    void OpenNewIndexFile_();

    const std::string version_{"0.1.0"};
    PacBio::Pancake::SeedDB::SeedDBParameters params_;
    std::string filenamePrefix_;
    std::string parentFolder_;
    std::string basenamePrefix_;
    bool splitBlocks_;
    std::vector<SeedDBFileLine> fileLines_;
    std::vector<SeedDBSeedsLine> seedsLines_;
    std::vector<SeedDBBlockLine> blockLines_;
    SeedDBBlockLine currentBlock_;
    bool openNewSeedsFileUponNextWrite_;
    // Output file handlers.
    std::unique_ptr<FILE, FileDeleter> fpOutIndex_{nullptr};
    std::string outIndexFilename_;
    std::unique_ptr<FILE, FileDeleter> fpOutSeeds_{nullptr};
};

std::unique_ptr<SeedDBWriter> CreateSeedDBWriter(
    const std::string& filenamePrefix, bool splitBlocks,
    const PacBio::Pancake::SeedDB::SeedDBParameters& params);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEED_DB_WRITER_HPP
