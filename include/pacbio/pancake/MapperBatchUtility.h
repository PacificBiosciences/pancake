// Author: Ivan Sovic

#ifndef PANCAKE_MAPPER_BATCH_UTILITY_H
#define PANCAKE_MAPPER_BATCH_UTILITY_H

#include <pacbio/pancake/AlignmentSeeded.h>
#include <pacbio/pancake/FastaSequenceCached.h>
#include <pacbio/pancake/MapperBase.h>
#include <pacbio/pancake/MapperCLR.h>
#include <pbcopper/parallel/FireAndForget.h>
#include <cstdint>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

namespace PacBio {
namespace Pancake {

struct MapperBatchChunk
{
    FastaSequenceCachedStore targetSeqs;
    FastaSequenceCachedStore querySeqs;
    MapperCLRMapSettings mapSettings;
};

enum class BatchAlignerRegionType
{
    GLOBAL,
    SEMIGLOBAL,
    BOTH,
};

enum class StatusPrepareSequences
{
    OK,
    BATCH_FULL,
    FAILED,
};

struct PairForBatchAlignment
{
    const char* query = NULL;
    int32_t queryLen = 0;
    const char* target = NULL;
    int32_t targetLen = 0;
    RegionType regionType = RegionType::GLOBAL;
    bool regionIsRev = false;
};

/*
 * \brief The AlignmentStitchPart is used to relate the alignment regions
 * created during mapping to the actual alignment parts.
*/
struct AlignmentStitchPart
{
    // RegionType is defined in AlignmentSeeded. It's either FRONT, BACK or GLOBAL.
    // If FRONT, the sequence of this part will have to be reversed.
    RegionType regionType = RegionType::GLOBAL;

    // The partId corresponds to the ID in either the internal or the
    // flank PairForBatchAlignment vectors, depending on the regionType value.
    int64_t partId = 0;

    // ID of the corresponding region in ChainedRegion::regionsForAln for an alignment
    // (it's a std::vector<AlignmentRegion>).
    int64_t regionId = 0;
};

/*
 * \brief Keeps the information needed to stitch the batch alignment for a particular query.
 * Each object is enough to reconstruct one query alignment. The vector of parts holds the
 * information about which internal or flank alignment portion should be used to splice
 * the alignment together, linearly.
*/
class AlignmentStitchInfo
{
public:
    AlignmentStitchInfo(int32_t _ordinalBatchId, int32_t _ordinalQueryId, int32_t _ordinalMapId)
        : ordinalBatchId(_ordinalBatchId)
        , ordinalQueryId(_ordinalQueryId)
        , ordinalMapId(_ordinalMapId)
    {
    }

    std::vector<AlignmentStitchPart> parts;
    int32_t ordinalBatchId = -1;
    int32_t ordinalQueryId = -1;
    int32_t ordinalMapId = -1;
};

inline std::ostream& operator<<(std::ostream& os, const AlignmentStitchInfo& b)
{
    os << "ordinalBatchId = " << b.ordinalBatchId << ", ordinalQueryId = " << b.ordinalQueryId
       << ", ordinalMapId = " << b.ordinalMapId;
    return os;
}

// typedef std::vector<AlignmentStitchPart> AlignmentStitchVector;

/*
 * This utility function takes a vector of pairs of query+target sequences, and reshapes them into
 * the MapperBatchChunk data structure. One instance of a MapperBatchChunk object represents a single
 * mapping/alignment job, while a vector of these objects is the batch job which will be executed by
 * the batch mapper.
 * Note: each mapping/alignment job can align multiple queries to multiple targets, that's why the
 * pair has two vectors in it.
*/
std::vector<PacBio::Pancake::MapperBatchChunk> ConstructBatchData(
    const std::vector<std::pair<std::vector<std::string>, std::vector<std::string>>>& inData);

void PrepareSequencesForBatchAlignment(
    const std::vector<MapperBatchChunk>& batchChunks,
    const std::vector<FastaSequenceCachedStore>& querySeqsRev,
    const std::vector<std::vector<MapperBaseResult>>& mappingResults,
    const MapperSelfHitPolicy selfHitPolicy, std::vector<PairForBatchAlignment>& retPartsGlobal,
    std::vector<PairForBatchAlignment>& retPartsSemiglobal,
    std::vector<AlignmentStitchInfo>& retAlnStitchInfo, int32_t& retLongestSequence);

OverlapPtr StitchSingleAlignment(const OverlapPtr& aln,
                                 const std::vector<AlignmentRegion>& regionsForAln,
                                 const std::vector<AlignmentResult>& internalAlns,
                                 const std::vector<AlignmentResult>& flankAlns,
                                 const std::vector<AlignmentStitchPart>& parts);

void StitchAlignmentsInParallel(std::vector<std::vector<MapperBaseResult>>& mappingResults,
                                const std::vector<MapperBatchChunk>& batchChunks,
                                const std::vector<FastaSequenceCachedStore>& querySeqsRev,
                                const std::vector<AlignmentResult>& internalAlns,
                                const std::vector<AlignmentResult>& flankAlns,
                                const std::vector<AlignmentStitchInfo>& alnStitchInfo,
                                Parallel::FireAndForget* faf);

void SetUnalignedAndMockedMappings(std::vector<std::vector<MapperBaseResult>>& mappingResults,
                                   const bool mockPerfectAlignment,
                                   const int32_t matchScoreForMockAlignment);

/**
 * \brief This function computes the reverse complements of the query sequences.
 *
 * As an optimization, if the onlyWhenRequired == true then the reverse complement for
 * a query will be computed only if there is a mapping that maps the reverse strand
 * of a query.
 * Otherwise, an entry in the return vector will be generated, but the sequence will be
 * an empty string.
*/
std::vector<std::vector<FastaSequenceId>> ComputeQueryReverseComplements(
    const std::vector<MapperBatchChunk>& batchChunks,
    const std::vector<std::vector<MapperBaseResult>>& mappingResults, const bool onlyWhenRequired,
    Parallel::FireAndForget* faf);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_MAPPER_BATCH_UTILITY_H
