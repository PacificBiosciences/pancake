// Author: Ivan Sovic

#ifndef PANCAKE_SEQ_DB_INDEX_CACHE_HPP
#define PANCAKE_SEQ_DB_INDEX_CACHE_HPP

#include <pancake/ContiguousFilePart.hpp>
#include <pancake/Range.hpp>
#include <pancake/third-party/flat_hash_map/flat_hash_map.hpp>
#include <pancake/util/CommonTypes.hpp>
#include <pancake/util/Conversion.hpp>

#include <cstdint>
#include <cstdio>
#include <fstream>
#include <memory>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace PacBio {
namespace Pancake {

/// \brief      Container, describes a sequence file which accompanies the DB index.
///
class SeqDBFileLine
{
public:
    int32_t fileId = 0;
    std::string filename;
    int32_t numSequences = 0;
    int64_t numBytes = 0;
    int64_t numUncompressedBases = 0;
    int64_t numCompressedBases = 0;
};

/// \brief      Container, index information for a particular sequence.
///
class SeqDBSequenceLine
{
public:
    int32_t seqId = 0;
    std::string header;
    int32_t numBytes = 0;
    int32_t numBases = 0;
    int32_t fileId = 0;
    int64_t fileOffset = 0;
    std::vector<PacBio::Pancake::Range> ranges;
};

class SeqDBBlockLine
{
public:
    int32_t blockId = 0;
    int32_t startSeqId = -1;
    int32_t endSeqId = -1;
    int64_t numBytes = 0;
    int64_t numBases = 0;

    int32_t Span() const { return endSeqId - startSeqId; }
};

/// \brief      Container used to store all the indexing information from the DB.
///
class SeqDBIndexCache
{
public:
    std::string indexFilename;
    std::string indexParentFolder;
    std::string indexBasename;
    std::string version{"unknown"};
    int32_t compressionLevel = 0;
    std::vector<SeqDBFileLine> fileLines;
    std::vector<SeqDBSequenceLine> seqLines;
    std::vector<SeqDBBlockLine> blockLines;

    SeqDBIndexCache() = default;
    ~SeqDBIndexCache() = default;
    friend std::ostream& operator<<(std::ostream& os, const PacBio::Pancake::SeqDBIndexCache& r);

    const SeqDBSequenceLine& GetSeqLine(int32_t seqId) const;
    const SeqDBSequenceLine& GetSeqLine(const std::string& header) const;
    const SeqDBBlockLine& GetBlockLine(int32_t blockId) const;
    const SeqDBFileLine& GetFileLine(int32_t fileId) const;

    const HeaderLookupType& GetHeaderLookup() const;
    bool IsHeaderLookupConstructed() const;

    /// \brief Checks that the fileLines, seqLines and blockLines are not empty.
    ///        Throws if they are.
    void Validate() const;

    /// \brief Constructs a header lookup from this index object and stores it
    ///         to the private members.
    void ConstructHeaderLookup();

private:
    HeaderLookupType headerToOrdinalId_;
    bool headerToOrdinalIdConstructed_ = false;
};

/// \brief Utility function to fetch a sequence by header.
const SeqDBSequenceLine& GetSeqLine(const SeqDBIndexCache& indexCache,
                                    const HeaderLookupType& headerToOrdinalId,
                                    const std::string& header);

/// \brief Validates the pointer and then calls indexCache->Validate.
void ValidateSeqDBIndexCache(const std::shared_ptr<PacBio::Pancake::SeqDBIndexCache>& indexCache);

/// \brief Writes the .seqdb file to disk.
void WriteSeqDBIndexCache(std::FILE* fpOut, const SeqDBIndexCache& cache);

void ComputeSeqDBIndexHeaderLookup(const PacBio::Pancake::SeqDBIndexCache& dbCache,
                                   HeaderLookupType& headerToOrdinalId);

/// \brief Given a vector of SeqDBSequenceLine objects, this function groups them in
///        creates blocks of given block size.
std::vector<SeqDBBlockLine> CreateSeqDBBlocks(const std::vector<SeqDBSequenceLine>& seqLines,
                                              int64_t blockSize);

/// \brief When a SeqDBIndexCache is subsampled or manually constructed, the seqID of each sequence
///        might not correspond to the sequence's ordinal ID. This function updates the seqID to
///        match, and also constructs the blocks of given block size.
///        It modifies the cache object in place.
void NormalizeSeqDBIndexCache(SeqDBIndexCache& cache, int64_t blockSize);

/// \brief Loads the SeqDB index from file.
///
/// \note  This function also computes the ranges of blocks give a block size
///        in bytes. If the block size is < 0, then the entire index is in
///        a single block. If block size is 0, each sequence will be in one block.

/// \param[in]  indexFilename           Path to the .seqdb index file.
/// \param[out] blockSizeInBytes        Block size.
///
/// \returns Pointer to the parsed index.
///
std::unique_ptr<PacBio::Pancake::SeqDBIndexCache> LoadSeqDBIndexCache(
    const std::string& indexFilename);

std::unique_ptr<PacBio::Pancake::SeqDBIndexCache> LoadSeqDBIndexCache(
    std::istream& is, const std::string& indexFilename);

std::unique_ptr<PacBio::Pancake::SeqDBIndexCache> LoadSeqDBIndexCache(
    std::FILE* fpIn, const std::string& indexFilename);

std::vector<ContiguousFilePart> GetSeqDBContiguousParts(
    const std::shared_ptr<PacBio::Pancake::SeqDBIndexCache>& seqDBIndexCache, int32_t blockId);

std::vector<ContiguousFilePart> GetSeqDBContiguousParts(
    const std::shared_ptr<PacBio::Pancake::SeqDBIndexCache>& seqDBIndexCache,
    std::vector<int32_t> seqIdsToFetch);

std::vector<ContiguousFilePart> GetSeqDBContiguousParts(
    const std::shared_ptr<PacBio::Pancake::SeqDBIndexCache>& seqDBIndexCache,
    const std::vector<std::string>& seqNamesToFetch);

/// \brief Filters the list of sequence lines from a SeqDB index.
///
/// \param outSeqLines Resulting list of sequence lines.
/// \param inSeqLines Input list of sequence lines.
/// \param samplingType Optional sampling of the input DB: linear, random or none.
/// \param sampledBases If sampling is not none, then this many bases will be sampled.
/// \param randomSeed Seed for random sampling, if random sampling is applied.
/// \param filterList A set of sequence names to keep/filter, depending on the filterType.
/// \param filterType If specified, the filterList will be used. Can be: whitelist, blacklist or none.
void PerformSeqDBSequenceLineSampling(std::vector<SeqDBSequenceLine>& outSeqLines,
                                      const std::vector<SeqDBSequenceLine>& inSeqLines,
                                      SamplingType sampling, int64_t sampledBases,
                                      int64_t randomSeed,
                                      const std::unordered_set<std::string>& filterList,
                                      FilterListType filterType);

/// \brief Filters the SeqDB index and returns a filtered index.
///
/// \param inSeqDBCache Input SeqDB index for filtering.
/// \param samplingType Optional sampling of the input DB: linear, random or none.
/// \param sampledBases If sampling is not none, then this many bases will be sampled.
/// \param randomSeed Seed for random sampling, if random sampling is applied.
/// \param filterType If specified, the filterList will be used. Can be: whitelist, blacklist or none.
/// \param filterList A set of sequence names to keep/filter, depending on the filterType.
/// \param doNormalization If true, the sequence IDs and blocks will be reset. Otherwise, they
///                        will remain the same as in the input DB. Normalization should not be
///                        applied if the filtered DB will be used for random access by sequence ID of
///                        the original input DB.
/// \param normBlockSize Block size to be used during normalization.
/// \param outIndexFilename The full path (relative or absolute) where the filtered DB
///                         might be written. This function will not write the DB, but
///                         this infor is stored internally in the DB to facilitate finding
///                         of the .seq files.
/// \returns Pointer to the newly constructed filtered DB.
std::unique_ptr<PacBio::Pancake::SeqDBIndexCache> FilterSeqDBIndexCache(
    const SeqDBIndexCache& inSeqDBCache, SamplingType samplingType, int64_t sampledBases,
    int64_t randomSeed, FilterListType filterType,
    const std::unordered_set<std::string>& filterList, bool doNormalization, int32_t normBlockSize,
    const std::string& outIndexFilename = "");

/// \brief Takes a sequence header and returns the ID of the sequence. In case the header is already
///         a numeric ID stored as a string, this function will only parse the number. Otherwise, it will
///         look-up the name in the SeqDB.
///
/// \param header The sequence header to be converted to ID.
/// \param headerIsNumeric If true, then the function will consider the header string as a string
///                         that contains only a number, and parse this number as the ID.
/// \param seqDBCache The SeqDB index which contains the header to ID mapping. It should have been
///                     initialized using "seqDBCache->ConstructHeaderLookup()" before the call to this function.
///
/// \returns ID of the sequence with the specified name. Otherwise, it throws if the name cannot be found,
///             or if the ID cannot be parsed.
int32_t GetSequenceIdFromHeader(const std::string& header, bool headerIsNumeric,
                                const Pancake::SeqDBIndexCache& seqDBCache);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEQ_DB_INDEX_CACHE_HPP
