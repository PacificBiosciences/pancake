// Author: Ivan Sovic

#ifndef PANCAKE_SEEDDB_INDEX_CACHE_H
#define PANCAKE_SEEDDB_INDEX_CACHE_H

#include <pacbio/pancake/ContiguousFilePart.h>
#include <pacbio/pancake/Range.h>
#include <pacbio/pancake/SeedDBParameters.h>
#include <pacbio/util/CommonTypes.h>
#include <cstdint>
#include <fstream>
#include <lib/flat_hash_map/flat_hash_map.hpp>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace PacBio {
namespace Pancake {

/// \brief      Container, describes a seeds file which accompanies the SeedDB index.
///
class SeedDBFileLine
{
public:
    int32_t fileId = 0;
    std::string filename;
    int32_t numSequences = 0;
    int64_t numBytes = 0;
};

/// \brief      Container, index information for a particular sequence's set of seeds.
///
class SeedDBSeedsLine
{
public:
    int32_t seqId = 0;
    std::string header;
    int32_t fileId = 0;
    int64_t fileOffset = 0;
    int64_t numBytes = 0;
    int32_t numBases = 0;
    int32_t numSeeds = 0;
};

/// \brief      Container, index information for a particular block of sequences.
///
class SeedDBBlockLine
{
public:
    int32_t blockId = 0;
    int32_t startSeqId = -1;
    int32_t endSeqId = -1;
    int64_t numBytes = 0;

    int32_t Span() const { return endSeqId - startSeqId; }
};

/// \brief      Container used to store all the indexing information from the DB.
///
class SeedDBIndexCache
{
public:
    std::string indexFilename;
    std::string indexParentFolder;
    std::string indexBasename;
    std::string version{"unknown"};
    PacBio::Pancake::SeedDB::SeedDBParameters seedParams;
    std::vector<SeedDBFileLine> fileLines;
    std::vector<SeedDBSeedsLine> seedLines;
    std::vector<SeedDBBlockLine> blockLines;

    SeedDBIndexCache() = default;
    ~SeedDBIndexCache() = default;
    friend std::ostream& operator<<(std::ostream& os, const PacBio::Pancake::SeedDBIndexCache& r);

    const SeedDBFileLine& GetFileLine(int32_t fileId) const;
    const SeedDBBlockLine& GetBlockLine(int32_t blockId) const;
    const SeedDBSeedsLine& GetSeedsLine(int32_t seqId) const;

    /// \brief Checks that the fileLines, seqLines and blockLines are not empty.
    ///        Throws if they are.
    void Validate() const;
};

void ComputeSeedDBIndexHeaderLookup(const PacBio::Pancake::SeedDBIndexCache& dbCache,
                                    HeaderLookupType& headerToOrdinalId);

/// \brief Loads the SeedDB index from file.
///
/// \param[in]  indexFilename           Path to the .seqdb index file.
///
/// \returns Pointer to the parsed index.
///
std::unique_ptr<PacBio::Pancake::SeedDBIndexCache> LoadSeedDBIndexCache(
    const std::string& indexFilename);

std::unique_ptr<PacBio::Pancake::SeedDBIndexCache> LoadSeedDBIndexCache(
    std::istream& is, const std::string& indexFilename);

std::unique_ptr<PacBio::Pancake::SeedDBIndexCache> LoadSeedDBIndexCache(
    FILE* fpIn, const std::string& indexFilename);

std::vector<ContiguousFilePart> GetSeedDBContiguousParts(
    const std::shared_ptr<PacBio::Pancake::SeedDBIndexCache>& seedDBIndexCache, int32_t blockId);

/// \brief Validates the pointer and then calls indexCache->Validate.
void ValidateSeedDBIndexCache(const std::shared_ptr<PacBio::Pancake::SeedDBIndexCache>& indexCache);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEQDB_CACHE_H