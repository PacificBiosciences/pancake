// Author: Ivan Sovic

#ifndef PANCAKE_OVERLAPHIFI_OVERLAPPER_H
#define PANCAKE_OVERLAPHIFI_OVERLAPPER_H

#include <pacbio/alignment/SesResults.h>
#include <pacbio/overlaphifi/Overlap.h>
#include <pacbio/overlaphifi/OverlapHifiSettings.h>
#include <pacbio/overlaphifi/SeedIndex.h>
#include <pacbio/seeddb/Seed.h>
#include <pacbio/seeddb/SeedDBIndexCache.h>
#include <pacbio/seeddb/SequenceSeedsCached.h>
#include <pacbio/seqdb/FastaSequenceCached.h>
#include <pacbio/seqdb/SeqDBReaderCachedBlock.h>
#include <pacbio/util/CommonTypes.h>
#include <cstdint>
#include <memory>
#include <unordered_map>
#include <vector>

namespace PacBio {
namespace Pancake {

class MapperResult
{
public:
    std::vector<PacBio::Pancake::OverlapPtr> overlaps;
};

class Mapper
{
public:
    Mapper(const OverlapHifiSettings& settings)
        : settings_{settings}, sesScratch_{std::make_shared<Pancake::Alignment::SESScratchSpace>()}
    {
    }
    ~Mapper() = default;

    /// \brief Maps a single query to a given set of targets. The targets
    ///         are provided with their sequences (targetSeqs) and the seeds (index).
    ///
    /// \param targetSeqs Cached target sequences (i.e. a single block of sequences).
    /// \param index The seed index for seed lookup.
    /// \param querySeq The query sequence which will be mapped.
    /// \param querySeeds Precomputed seeds for the query sequence. The seeds will not be
    ///                     computed internally like in other mappers.
    /// \param freqCutoff Maximum allowed frequency of any particular seed to retain it.
    /// \returns An object which contains a vector of all found overlaps.
    ///
    MapperResult Map(const PacBio::Pancake::SeqDBReaderCachedBlock& targetSeqs,
                     const PacBio::Pancake::SeedIndex& index,
                     const PacBio::Pancake::FastaSequenceCached& querySeq,
                     const PacBio::Pancake::SequenceSeedsCached& querySeeds,
                     int64_t freqCutoff) const;

private:
    OverlapHifiSettings settings_;
    std::shared_ptr<PacBio::Pancake::Alignment::SESScratchSpace> sesScratch_;

    /// \brief Writes the seed hits to a specified file, in a CSV format, useful for visualization.
    /// The header line contains:
    ///     <queryName> <queryStart> <querySpan> <targetName> <targetStart> <targetSpan> 0.0
    /// Every other line:
    ///     <queryPos> <targetPos> <targetID>
    ///
    /// \param outPath Output file where to write to. If the file cannot be open, the function will
    ///                 not throw but simply silently continue (intentionally).
    /// \param hits The vector of SeedHits which to write out.
    /// \param seedLen Length of any particular seed, will be used to write out the endpoints of seed hits.
    /// \param queryName Name of the query sequence.
    /// \param queryLen Length of the query sequence.
    /// \param targetName Name of the target sequence.
    /// \param targetLen Length of the target sequence.
    ///
    static void DebugWriteSeedHits_(const std::string& outPath, const std::vector<SeedHit>& hits,
                                    int32_t seedLen, const std::string& queryName, int64_t queryLen,
                                    const std::string& targetName, int64_t targetLen);

    /// \brief Forms anchors simply by binning seeds in narrow diagonals.
    /// \param sortedHits Hits should be sorted in the following order of priority:
    ///                   (targetID, targetReverse, diagonal, targetPos, queryPos)
    /// \param querySeq The query sequence.
    /// \param indexCache is needed to fetch the length of the target sequences.
    /// \param chainBandwidth Allowed diagonal bandwidth to bin hits together.
    /// \param minNumSeeds Minimum number of seeds per diagonal bin to retain it.
    /// \param minChainSpan Minimum span (in either query or target coordinates) of the
    ///                     hits that are binned in a diagonal window.
    /// \param skipSelfHits Ignore hits where Aid == Bid.
    /// \param skipSymmetricOverlaps Do not produce an anchor if Aid > Bid.
    ///
    static std::vector<OverlapPtr> FormDiagonalAnchors_(
        const std::vector<SeedHit>& sortedHits,
        const PacBio::Pancake::FastaSequenceCached& querySeq,
        const std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> indexCache, int32_t chainBandwidth,
        int32_t minNumSeeds, int32_t minChainSpan, bool skipSelfHits, bool skipSymmetricOverlaps);

    /// \brief  Helper function used by FormDiagonalAnchors_ which creates a new overlap object
    ///         based on the minimum and maximum hit IDs.
    /// \param sortedHits Hits should be sorted in the following order of priority:
    ///                   (targetID, targetReverse, diagonal, targetPos, queryPos)
    /// \param querySeq The query sequence.
    /// \param indexCache is needed to fetch the length of the target sequences.
    /// \param beginId The ID of the first element in sortedHits corresponding to the current
    ///                diagonal bin for which we're constructing the overlap.
    /// \param endId The ID of the last element in sortedHits corresponding to the current
    ///              diagonal bin for which we're constructing the overlap. Non-inclusive
    ///              value (points to one after the last element).
    /// \param minTargetPosID The ID of the hit within the diagonal bin with the smallest target pos.
    /// \param maxTargetPosID The ID of the hit within the diagonal bin with the largest target pos.
    ///
    static OverlapPtr MakeOverlap_(
        const std::vector<SeedHit>& sortedHits,
        const PacBio::Pancake::FastaSequenceCached& querySeq,
        const std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> indexCache, int32_t beginId,
        int32_t endId, int32_t minTargetPosId, int32_t maxTargetPosId);

    /// \brief Performs alignment and alignment extension of a given vector of overlaps.
    ///         A thin wrapper around AlignOverlap_, simply calls it for each overlap.
    /// \param targetSeqs A cached sequence reader, to allow random access to sequence data.
    /// \param querySeq The query sequence.
    /// \param overlaps A vector of all overlaps to align.
    /// \param alignBandwidth The maximum allowed bandwidth for alignment. Used for the banded
    ///                       O(nd) algorithm
    /// \param alignMaxDiff The maximum number of diffs allowed between the query and target pair.
    ///                     This is a parameter of the O(nd) algorithm.
    /// \param useTraceback Runs alignment with traceback, for more accurate
    ///                     identity computation (in terms of mismatches) and CIGAR construction.
    /// \returns A new vector of overlaps with alignment information and modified coordinates.
    ///
    static std::vector<OverlapPtr> AlignOverlaps_(
        const PacBio::Pancake::SeqDBReaderCachedBlock& targetSeqs,
        const PacBio::Pancake::FastaSequenceCached& querySeq,
        const std::vector<OverlapPtr>& overlaps, double alignBandwidth, double alignMaxDiff,
        bool useTraceback, bool noSNPs, bool noIndels, bool maskHomopolymers,
        bool maskSimpleRepeats,
        std::shared_ptr<PacBio::Pancake::Alignment::SESScratchSpace> sesScratch);

    /// \brief Performs alignment and alignment extension of a single overlap. Uses the
    ///        banded O(nd) algorithm to align the overlap. The edit distance is
    ///        (pessimistically) estimated using the number of diffs computed by the O(nd).
    /// \param targetSeq The target sequence (B-read) for alignment.
    /// \param querySeq The query sequence (A-read) for alignment.
    /// \param reverseQuerySeq The full reverse-complemented query sequence.
    /// \param ovl The overlap which to align.
    /// \param alignBandwidth The maximum allowed bandwidth for alignment. Used for the banded
    ///                       O(nd) algorithm
    /// \param alignMaxDiff The maximum number of diffs allowed between the query and target pair.
    ///                     This is a parameter of the O(nd) algorithm.
    /// \param useTraceback Runs alignment with traceback, for more accurate
    ///                     identity computation (in terms of mismatches) and CIGAR construction.
    /// \returns A new vector overlap with alignment information and modified coordinates.
    ///
    static OverlapPtr AlignOverlap_(
        const PacBio::Pancake::FastaSequenceCached& targetSeq,
        const PacBio::Pancake::FastaSequenceCached& querySeq, const std::string reverseQuerySeq,
        const OverlapPtr& ovl, double alignBandwidth, double alignMaxDiff, bool useTraceback,
        bool noSNPs, bool noIndels, bool maskHomopolymers, bool maskSimpleRepeats,
        std::shared_ptr<PacBio::Pancake::Alignment::SESScratchSpace> sesScratch);

    static void NormalizeAndExtractVariantsInPlace_(
        OverlapPtr& ovl, const PacBio::Pancake::FastaSequenceCached& targetSeq,
        const PacBio::Pancake::FastaSequenceCached& querySeq, const std::string reverseQuerySeq,
        bool noSNPs, bool noIndels, bool maskHomopolymers, bool maskSimpleRepeats);

    /// \brief Filters overlaps based on the number of seeds, identity, mapped span or length.
    ///
    /// \param overlaps A vector of overlaps to filter.
    /// \param minNumSeeds Minimum allowed number of seeds available to form the initial overlap anchor.
    /// \param minIdentity Minimum allowed estimated identity of the aligned overlap, in percentage.
    /// \param minMappedSpan Minimum allowed span of the overlap, in either query or target coordinats.
    /// \param minQueryLen Minimum allowed query length.
    /// \param minTargetLen Minimum allowed target length.
    /// \param allowedDovetailDist Used to determine the type of the overlap (5', 3', contained/contains, internal).
    ///                            For accurate data, this should be 0 (or a value close to 0).
    /// \param allowedExtendDist Heuristically extend the overlaps into a dovetail form by augmenting the
    ///                          coordinates, but only if the unaligned flank is < allowedExtendDist.
    /// \param bestN Keep only best N overlaps. If bestN <= 0, all overlaps are kept. This keeps overall best
    ///               overlaps and it's not side specific (i.e. this bestN does not care about 5' or 3' ends).
    /// \returns A new vector of remaining overlaps.
    ///
    static std::vector<OverlapPtr> FilterOverlaps_(const std::vector<OverlapPtr>& overlaps,
                                                   int32_t minNumSeeds, float minIdentity,
                                                   int32_t minMappedSpan, int32_t minQueryLen,
                                                   int32_t minTargetLen,
                                                   int32_t allowedDovetailDist,
                                                   int32_t allowedExtendDist, int32_t bestN);
    /// \brief  Filters multiple overlaps for the same query-target pair, for example tandem repeats,
    ///         and keeps only the longest spanning overlap. The maximum of (querySpan, targetSpan)
    ///         is taken for a particular query-target pair for comparison.
    /// \param overlaps A vector of overlaps to filter.
    /// \return A vector of filtered overlaps.
    ///
    static std::vector<OverlapPtr> FilterTandemOverlaps_(const std::vector<OverlapPtr>& overlaps);

    /// \brief  Helper function which converts a SeedHit into a packed 128-bit integer.
    ///         Difference between this one and the one native to SeedHit is that here we pack
    ///         The diagonal of the seed hit, which is important for proper sorting of hits
    ///         before constructing anchors.
    ///
    /// \param sh SeedHit which to pack into an 128-bit integer.
    /// \returns A packed 128-bit integer composed of:
    ///             targetId:31, targetRev:1, diag:32, targetPos:32, queryPos:32.
    ///
    static PacBio::Pancake::Int128t PackSeedHitWithDiagonalTo128_(const SeedHit& sh);

    /// \brief  Helper function which extracts a subsequence from a given sequence, and reverse
    ///         complements if needed.
    /// \param targetSeq The full sequence from which a subsequence should be extracted.
    /// \param seqStart Start position (0-based) to extract the subsequence.
    /// \param seqEnd End position (0-based, non-inclusive) to extract the subsequence.
    /// \param revCmp True if the sequence should be reverse complemented.
    /// \return Returns the requested bases as a string.
    ///
    static std::string FetchTargetSubsequence_(
        const PacBio::Pancake::FastaSequenceCached& targetSeq, int32_t seqStart, int32_t seqEnd,
        bool revCmp);
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_OVERLAPHIFI_OVERLAPPER_H
