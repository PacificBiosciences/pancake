// Author: Ivan Sovic

#include "MapCLRSettings.hpp"

#include <pancake/Version.hpp>

namespace PacBio {
namespace Pancake {
namespace MapCLROptionNames {

// clang-format off

const CLI_v2::PositionalArgument TargetDBPrefix {
R"({
    "name" : "target_prefix",
    "description" : "Prefix of the target SeqDB and SeedDB files. It should match."
})"};

const CLI_v2::PositionalArgument QueryDBPrefix {
R"({
    "name" : "query_prefix",
    "description" : "Prefix of the query SeqDB and SeedDB files. It should match."
})"};

/////////////////////////////////
///// Input/Output options. /////
/////////////////////////////////
const CLI_v2::Option OutFormat{
R"({
    "names" : ["out-fmt"],
    "choices" : ["m4", "ipa", "paf", "sam"],
    "type" : "string",
    "default" : "paf",
    "description" : "Select the output format."
})", std::string("paf")};

const CLI_v2::Option WriteIds{
R"({
    "names" : ["write-ids"],
    "description" : "Output overlaps will contain numeric IDs for the A and B reads (instead of names).",
    "type" : "bool"
})", MapCLRSettings::Defaults::WriteIds};

const CLI_v2::Option NoCigar{
R"({
    "names" : ["no-cigar"],
    "description" : "Do not output the CIGAR even if alignment was computed.",
    "type" : "bool"
})", !MapCLRSettings::Defaults::WriteCigar};

const CLI_v2::Option TargetBlockId{
R"({
    "names" : ["target-block", "tb"],
    "description" : "Block ID from the target DB. Queries will be mapped only onto this block.",
    "type" : "int"
})", MapCLRSettings::Defaults::TargetBlockId};

const CLI_v2::Option QueryBlockStartId{
R"({
    "names" : ["query-block-start", "qbs"],
    "description" : "Start block ID for a range of blocks to map. Zero based.",
    "type" : "int"
})", MapCLRSettings::Defaults::QueryBlockStartId};

const CLI_v2::Option QueryBlockEndId{
R"({
    "names" : ["query-block-end", "qbe"],
    "description" : "Start block ID for a range of blocks to map. Zero based, non-inclusive. Value == 0 runs until the end block.",
    "type" : "int"
})", MapCLRSettings::Defaults::QueryBlockEndId};

const CLI_v2::Option CombineBlocks{
R"({
    "names" : ["combine"],
    "description" : "Combines this many query blocks into one larger block for processing.",
    "type" : "int"
})", MapCLRSettings::Defaults::CombineBlocks};

/////////////////////
///// Indexing. /////
/////////////////////
const CLI_v2::Option PrimaryKmerSize{
R"({
    "names" : ["k"],
    "description" : "Kmer size for the primary seeding parameters.",
    "type" : "int"
})", MapCLRSettings::Defaults::PrimaryKmerSize};

const CLI_v2::Option PrimaryMinimizerWindow{
R"({
    "names" : ["w"],
    "description" : "Minimizer window size for the primary seeding parameters.",
    "type" : "int"
})", MapCLRSettings::Defaults::PrimaryMinimizerWindow};

const CLI_v2::Option PrimarySpacing{
R"({
    "names" : ["s"],
    "description" : "Spacing for spaced seeds, for the primary seeding parameters.",
    "type" : "int"
})", MapCLRSettings::Defaults::PrimarySpacing};

const CLI_v2::Option PrimaryUseHomopolymerCompression{
R"({
    "names" : ["hpc"],
    "description" : "Use homopolymer compression for seed formation, primary seeding parameters.",
    "type" : "bool"
})", MapCLRSettings::Defaults::PrimaryUseHomopolymerCompression};

const CLI_v2::Option PrimaryUseReverseComplement{
R"({
    "names" : ["revcmp"],
    "description" : "Use reverse complement target sequences for indexing, primary seeding parameters.",
    "type" : "bool"
})", MapCLRSettings::Defaults::PrimaryUseReverseComplement};

const CLI_v2::Option FallbackKmerSize{
R"({
    "names" : ["k-fallback"],
    "description" : "Kmer size, fallback seeding parameters.",
    "type" : "int"
})", MapCLRSettings::Defaults::FallbackKmerSize};

const CLI_v2::Option FallbackMinimizerWindow{
R"({
    "names" : ["w-fallback"],
    "description" : "Minimizer window, fallback seeding parameters.",
    "type" : "int"
})", MapCLRSettings::Defaults::FallbackMinimizerWindow};

const CLI_v2::Option FallbackSpacing{
R"({
    "names" : ["s-fallback"],
    "description" : "Spacing for spaced seeds, fallback seeding parameters.",
    "type" : "int"
})", MapCLRSettings::Defaults::FallbackSpacing};

const CLI_v2::Option FallbackUseHomopolymerCompression{
R"({
    "names" : ["hpc-fallback"],
    "description" : "Use homopolymer compression for seed formation, fallback seeding parameters.",
    "type" : "bool"
})", MapCLRSettings::Defaults::FallbackUseHomopolymerCompression};

const CLI_v2::Option FallbackUseReverseComplement{
R"({
    "names" : ["revcmp-fallback"],
    "description" : "Use reverse complement target sequences for indexing, fallback seeding parameters.",
    "type" : "bool"
})", MapCLRSettings::Defaults::FallbackUseReverseComplement};

////////////////////
///// Mapping. /////
////////////////////
const CLI_v2::Option FreqPercentile{
R"({
    "names" : ["freq-percentile"],
    "description" : "Filter frequent kmers.",
    "type" : "double"
})", MapCLRSettings::Defaults::FreqPercentile};

const CLI_v2::Option SeedOccurrenceMin{
R"({
    "names" : ["occ-min"],
    "description" : "Minimum value for the occurrence threshold to keep a seed for mapping. If the frequency percentile is smaller than this, the threshold is pinned to this value.",
    "type" : "int"
})", MapCLRSettings::Defaults::SeedOccurrenceMin};

const CLI_v2::Option SeedOccurrenceMax{
R"({
    "names" : ["occ-max"],
    "description" : "Maximum allowed occurrence of a seed to keep for mapping. Value <= 0 turns off this threshold.",
    "type" : "int"
})", MapCLRSettings::Defaults::SeedOccurrenceMax};

const CLI_v2::Option SeedOccurrenceMaxMemory{
R"({
    "names" : ["occ-max-memory"],
    "description" : "Maximum allowed memory to be consumed by collected seed hits. This is used to dynamically compute the maximum occurrence cutoff based on the seed hit histogram. Seeds are chosen in the sorted order by their occurrence until the memory threshold is reached. Value <= 0 turns off this threshold.",
    "type" : "int"
})", MapCLRSettings::Defaults::SeedOccurrenceMaxMemory};

const CLI_v2::Option ChainMaxSkip{
R"({
    "names" : ["chain-max-skip"],
    "description" : "Continue to next seed if more than this many predecessors have been skipped.",
    "type" : "int"
})", MapCLRSettings::Defaults::ChainMaxSkip};

const CLI_v2::Option ChainMaxPredecessors{
R"({
    "names" : ["chain-max-pred"],
    "description" : "Maximum number of predecessors to check when chaining.",
    "type" : "int"
})", MapCLRSettings::Defaults::ChainMaxPredecessors};

const CLI_v2::Option ChainBandwidth{
R"({
    "names" : ["chain-bw"],
    "description" : "Diagonal bandwidth to merge seeds into chains.",
    "type" : "int"
})", MapCLRSettings::Defaults::ChainBandwidth};

const CLI_v2::Option SeedJoinDist{
R"({
    "names" : ["seed-join-dist"],
    "description" : "Maximum distance in either query or target coordinates to chain neighboring seeds.",
    "type" : "int"
})", MapCLRSettings::Defaults::SeedJoinDist};

const CLI_v2::Option LongMergeBandwidth{
R"({
    "names" : ["long-merge-bw"],
    "description" : "Maximum gap distance for chaining (abs(query_dist - target_dist)).",
    "type" : "int"
})", MapCLRSettings::Defaults::LongMergeBandwidth};

const CLI_v2::Option MinNumSeeds{
R"({
    "names" : ["min-num-seeds"],
    "description" : "Minimum number of seeds in an anchor.",
    "type" : "int"
})", MapCLRSettings::Defaults::MinNumSeeds};

const CLI_v2::Option MinCoveredBases{
R"({
    "names" : ["min-cov-bases"],
    "description" : "Minimum number of bases covered by kmers in an anchor.",
    "type" : "int"
})", MapCLRSettings::Defaults::MinCoveredBases};

const CLI_v2::Option MinDPScore{
R"({
    "names" : ["min-dp-score"],
    "description" : "Minimum DP score to keep a chain.",
    "type" : "int"
})", MapCLRSettings::Defaults::MinDPScore};

const CLI_v2::Option SecondaryAllowedOverlapFractionQuery{
R"({
    "names" : ["sec-ovl-frac-query"],
    "description" : "Maximum allowed query overlap between two chains in query coordinates to allow one to be called a supplementary. Overlap longer than this fraction will be marked as secondary.",
    "type" : "double"
})", MapCLRSettings::Defaults::SecondaryAllowedOverlapFractionQuery};

const CLI_v2::Option SecondaryAllowedOverlapFractionTarget{
R"({
    "names" : ["sec-ovl-frac-target"],
    "description" : "Maximum allowed overlap target between two chains in target coordinates to allow one to be called a supplementary. Overlap longer than this fraction will be marked as secondary.",
    "type" : "double"
})", MapCLRSettings::Defaults::SecondaryAllowedOverlapFractionTarget};

const CLI_v2::Option SecondaryMinScoreFraction{
R"({
    "names" : ["sec-min-score-frac"],
    "description" : "Consider only chains with a score higher than this fraction of the maximum score as potential secondary mappings, and ignore shorter ones.",
    "type" : "double"
})", MapCLRSettings::Defaults::SecondaryMinScoreFraction};

const CLI_v2::Option NoLIS{
R"({
    "names" : ["no-lis"],
    "description" : "Do not use LIS for chaining.",
    "type" : "bool"
})", MapCLRSettings::Defaults::NoLIS};

const CLI_v2::Option SelfHitPolicyMapping{
R"({
    "names" : ["map-self-hit-policy"],
    "choices" : ["DEFAULT", "SKIP", "PERFECT_ALIGNMENT"],
    "type" : "string",
    "default" : "DEFAULT",
    "description" : "Special treatment of self hits (when Aid == Bid): DEFAULT maps them as any other read pair, SKIP skips and does not align such cases, and PERFECT_ALIGNMENT mocks the alignment and reports a perfect alignment result for self hits."
})", std::string("DEFAULT")};

const CLI_v2::Option NoRefineSeedHits{
R"({
    "names" : ["no-refine-hits"],
    "description" : "Do not refine seed hits during the mapping process.",
    "type" : "bool"
})", !MapCLRSettings::Defaults::RefineSeedHits};

const CLI_v2::Option RefineMinGap1{
R"({
    "names" : ["refine-min-gap-1"],
    "description" : "Minimum gap length to select potential breakpoints in the first seed hit refinement heuristic.",
    "type" : "int"
})", MapCLRSettings::Defaults::RefineMinGap1};

const CLI_v2::Option RefineDiffThreshold{
R"({
    "names" : ["refine-diff-threshold"],
    "description" : "Maxium indel diff to keep seed hits, in the first seed hit refinement heuristic.",
    "type" : "int"
})", MapCLRSettings::Defaults::RefineDiffThreshold};

const CLI_v2::Option RefineMinGap2{
R"({
    "names" : ["refine-min-gap-2"],
    "description" : "Minimum gap length to select potential breakpoints in the second seed hit refinement heuristic.",
    "type" : "int"
})", MapCLRSettings::Defaults::RefineMinGap2};

//////////////////////
///// Filtering. /////
//////////////////////

const CLI_v2::Option SkipSymmetricOverlaps{
R"({
    "names" : ["skip-sym"],
    "description" : "If Aid < Bid, only compute alignment Aid->Bid and skip computing alignment for Bid->Aid.",
    "type" : "bool"
})", MapCLRSettings::Defaults::SkipSymmetricOverlaps};

const CLI_v2::Option MinQueryLen{
R"({
    "names" : ["min-qlen"],
    "description" : "Ignore queries shorter than this.",
    "type" : "int"
})", MapCLRSettings::Defaults::MinQueryLen};

const CLI_v2::Option BestNSecondary{
R"({
    "names" : ["bestn"],
    "description" : "Output only best N alignments.",
    "type" : "int"
})", MapCLRSettings::Defaults::BestNSecondary};

//////////////////////
///// Alignment. /////
//////////////////////

const CLI_v2::Option Align{
R"({
    "names" : ["align"],
    "description" : "Align the mappings.",
    "type" : "bool"
})", MapCLRSettings::Defaults::Align};

const CLI_v2::Option SelfHitPolicyAlignment{
R"({
    "names" : ["align-self-hit-policy"],
    "choices" : ["DEFAULT", "SKIP", "PERFECT_ALIGNMENT"],
    "type" : "string",
    "default" : "DEFAULT",
    "description" : "Special treatment of self hits (when Aid == Bid): DEFAULT aligns them as any other read pair, SKIP skips and does not align such cases, and PERFECT_ALIGNMENT mocks the alignment and reports a perfect alignment result for self hits."
})", std::string("DEFAULT")};

const CLI_v2::Option MaxFlankExtensionDistance{
R"({
    "names" : ["max-flank-ext"],
    "description" : "Extend flank alignments up to at most this distance from the last mapped seed.",
    "type" : "int"
})", MapCLRSettings::Defaults::MaxFlankExtensionDistance};

const CLI_v2::Option FlankExtensionFactor{
R"({
    "names" : ["flank-ext-factor"],
    "description" : "Take this much more of the longer flanking sequence for alignment, to allow for indel errors.",
    "type" : "double"
})", MapCLRSettings::Defaults::FlankExtensionFactor};

const CLI_v2::Option MinAlignmentSpan{
R"({
    "names" : ["min-aln-span"],
    "description" : "If two seeds are closer than this, take the next seed, unless there are only 2 seeds left.",
    "type" : "int"
})", MapCLRSettings::Defaults::MinAlignmentSpan};

const CLI_v2::Option AlignerTypeGlobal{
R"({
    "names" : ["aligner-global"],
    "choices" : ["KSW2", "EDLIB"],
    "type" : "string",
    "default" : "KSW2",
    "description" : "Aligner used for internal (global) alignment."
})", std::string("KSW2")};

const CLI_v2::Option AlignerTypeExt{
R"({
    "names" : ["aligner-ext"],
    "choices" : ["KSW2", "EDLIB"],
    "type" : "string",
    "default" : "KSW2",
    "description" : "Aligner used for extension (semiglobal) alignment."
})", std::string("KSW2")};

///////////////////////////////////////
/// Global alignment parameters.    ///
///////////////////////////////////////

const CLI_v2::Option AlignerGlobalZdrop {
R"({
    "names" : ["global-zdrop"],
    "description" : "Global aligner, zdrop parameter.",
    "type" : "int"
})", MapCLRSettings::Defaults::AlnParams.zdrop};

const CLI_v2::Option AlignerGlobalZdrop2 {
R"({
    "names" : ["global-zdrop2"],
    "description" : "Global aligner, zdrop2 parameter.",
    "type" : "int"
})", MapCLRSettings::Defaults::AlnParams.zdrop2};

const CLI_v2::Option AlignerGlobalBandwidth {
R"({
    "names" : ["global-bw"],
    "description" : "Global aligner, alignment bandwidth parameter.",
    "type" : "int"
})", MapCLRSettings::Defaults::AlnParams.alignBandwidth};

//// This parameter doesn'isn;t used for global alignment.
// const CLI_v2::Option AlignerGlobalEndBonus {
// R"({
//     "names" : ["global-end-bonus"],
//     "description" : "Global aligner, end bonus parameter.",
//     "type" : "int"
// })", MapCLRSettings::Defaults::AlnParams.endBonus};

const CLI_v2::Option AlignerGlobalMatchScore {
R"({
    "names" : ["global-match"],
    "description" : "Global aligner, alignment match score parameter.",
    "type" : "int"
})", MapCLRSettings::Defaults::AlnParams.matchScore};

const CLI_v2::Option AlignerGlobalMismatchPenalty {
R"({
    "names" : ["global-mismatch"],
    "description" : "Global aligner, alignment mismatch penalty parameter, positive value.",
    "type" : "int"
})", MapCLRSettings::Defaults::AlnParams.mismatchPenalty};

const CLI_v2::Option AlignerGlobalGapOpen1 {
R"({
    "names" : ["global-gap-open-1"],
    "description" : "Global aligner, alignment gap open 1 parameter (first affine function).",
    "type" : "int"
})", MapCLRSettings::Defaults::AlnParams.gapOpen1};

const CLI_v2::Option AlignerGlobalGapExtend1 {
R"({
    "names" : ["global-gap-ext-1"],
    "description" : "Global aligner, alignment gap extend 1 parameter (first affine function).",
    "type" : "int"
})", MapCLRSettings::Defaults::AlnParams.gapExtend1};

const CLI_v2::Option AlignerGlobalGapOpen2 {
R"({
    "names" : ["global-gap-open-2"],
    "description" : "Global aligner, alignment gap open 1 parameter (second affine function).",
    "type" : "int"
})", MapCLRSettings::Defaults::AlnParams.gapOpen2};

const CLI_v2::Option AlignerGlobalGapExtend2 {
R"({
    "names" : ["global-gap-ext-2"],
    "description" : "Global aligner, alignment gap extend 1 parameter (second affine function).",
    "type" : "int"
})", MapCLRSettings::Defaults::AlnParams.gapExtend2};

const CLI_v2::Option AlignerGlobalDynamicBandwidth {
R"({
    "names" : ["global-dynamic-bw"],
    "description" : "Global aligner, dynamically increase the bandwdith to max instead of applyin max right away; attempts to speed up alignment.",
    "type" : "bool"
})", MapCLRSettings::Defaults::AlnParams.dynamicBandwidth};

///////////////////////////////////////
/// Extension alignment parameters. ///
///////////////////////////////////////
const CLI_v2::Option AlignerExtZdrop {
R"({
    "names" : ["ext-zdrop"],
    "description" : "Global aligner, zdrop parameter.",
    "type" : "int"
})", MapCLRSettings::Defaults::AlnParams.zdrop};

const CLI_v2::Option AlignerExtZdrop2 {
R"({
    "names" : ["ext-zdrop2"],
    "description" : "Global aligner, zdrop2 parameter.",
    "type" : "int"
})", MapCLRSettings::Defaults::AlnParams.zdrop2};

const CLI_v2::Option AlignerExtBandwidth {
R"({
    "names" : ["ext-bw"],
    "description" : "Global aligner, alignment bandwidth parameter.",
    "type" : "int"
})", MapCLRSettings::Defaults::AlnParams.alignBandwidth};

const CLI_v2::Option AlignerExtEndBonus {
R"({
    "names" : ["ext-end-bonus"],
    "description" : "Global aligner, end bonus parameter.",
    "type" : "int"
})", MapCLRSettings::Defaults::AlnParams.endBonus};

const CLI_v2::Option AlignerExtMatchScore {
R"({
    "names" : ["ext-match"],
    "description" : "Global aligner, alignment match score parameter.",
    "type" : "int"
})", MapCLRSettings::Defaults::AlnParams.matchScore};

const CLI_v2::Option AlignerExtMismatchPenalty {
R"({
    "names" : ["ext-mismatch"],
    "description" : "Global aligner, alignment mismatch penalty parameter, positive value.",
    "type" : "int"
})", MapCLRSettings::Defaults::AlnParams.mismatchPenalty};

const CLI_v2::Option AlignerExtGapOpen1 {
R"({
    "names" : ["ext-gap-open-1"],
    "description" : "Global aligner, alignment gap open 1 parameter (first affine function).",
    "type" : "int"
})", MapCLRSettings::Defaults::AlnParams.gapOpen1};

const CLI_v2::Option AlignerExtGapExtend1 {
R"({
    "names" : ["ext-gap-ext-1"],
    "description" : "Global aligner, alignment gap extend 1 parameter (first affine function).",
    "type" : "int"
})", MapCLRSettings::Defaults::AlnParams.gapExtend1};

const CLI_v2::Option AlignerExtGapOpen2 {
R"({
    "names" : ["ext-gap-open-2"],
    "description" : "Global aligner, alignment gap open 1 parameter (second affine function).",
    "type" : "int"
})", MapCLRSettings::Defaults::AlnParams.gapOpen2};

const CLI_v2::Option AlignerExtGapExtend2 {
R"({
    "names" : ["ext-gap-ext-2"],
    "description" : "Global aligner, alignment gap extend 1 parameter (second affine function).",
    "type" : "int"
})", MapCLRSettings::Defaults::AlnParams.gapExtend2};

const CLI_v2::Option AlignerExtDynamicBandwidth {
R"({
    "names" : ["ext-dynamic-bw"],
    "description" : "Global aligner, dynamically increase the bandwdith to max instead of applyin max right away; attempts to speed up alignment.",
    "type" : "bool"
})", MapCLRSettings::Defaults::AlnParams.dynamicBandwidth};

//////////////////
///// Other. /////
//////////////////
///// Options below are intentionally left here as a reminder to implement later.

// const CLI_v2::Option MinIdentity{
// R"({
//     "names" : ["min-idt"],
//     "description" : "Minimum percent alignment identity allowed to report the alignment. This is an overall threshold which takes into account both indels and SNPs.",
//     "type" : "double"
// })", MapCLRSettings::Defaults::MinIdentity};
// Defaults::NoIndelsInIdentity};

// const CLI_v2::Option MinMappedLength{
// R"({
//     "names" : ["min-map-len"],
//     "description" : "Output only alignments above this length.",
//     "type" : "int"
// })", MapCLRSettings::Defaults::MinMappedLength};

// const CLI_v2::Option SkipSymmetricOverlaps{
// R"({
//     "names" : ["skip-sym"],
//     "description" : "If Aid < Bid, only compute overlap Aid->Bid and skip computing overlap for Bid->Aid.",
//     "type" : "bool"
// })", MapCLRSettings::Defaults::SkipSymmetricOverlaps};

// const CLI_v2::Option AllowSelfHits{
// R"({
//     "names" : ["allow-self-hits"],
//     "description" : "If both the query and the target DBs are the same and Aid == Bid then this is a self-hit. This option enables the output of such overlaps.",
//     "type" : "bool"
// })", false};

// const CLI_v2::Option OneHitPerTarget{
// R"({
//     "names" : ["one-hit-per-target"],
//     "description" : "Allow only one alignment per query/target pair.",
//     "type" : "bool"
// })", MapCLRSettings::Defaults::OneHitPerTarget};

// const CLI_v2::Option WriteReverseOverlaps{
// R"({
//     "names" : ["write-rev"],
//     "description" : "For eveery overlap, write out its reverse complement too.",
//     "type" : "bool"
// })", MapCLRSettings::Defaults::WriteReverseOverlaps};

// const CLI_v2::Option MarkSecondary{
// R"({
//     "names" : ["mark-secondary"],
//     "description" : "Mask homopolymer errors when traceback is generated. This will impact identity calculation.",
//     "type" : "bool"
// })", MapCLRSettings::Defaults::MarkSecondary};

// clang-format on

}  // namespace MapCLROptionNames

MapCLRSettings::MapCLRSettings() = default;

MapCLRSettings::MapCLRSettings(const PacBio::CLI_v2::Results& options)
    : TargetDBPrefix{options[MapCLROptionNames::TargetDBPrefix]}
    , QueryDBPrefix{options[MapCLROptionNames::QueryDBPrefix]}
    , NumThreads{options.NumThreads()}
    , TargetBlockId{options[MapCLROptionNames::TargetBlockId]}
    , QueryBlockStartId{options[MapCLROptionNames::QueryBlockStartId]}
    , QueryBlockEndId{options[MapCLROptionNames::QueryBlockEndId]}

    , WriteIds{options[MapCLROptionNames::WriteIds]}
    , WriteCigar{options[MapCLROptionNames::Align] && (!options[MapCLROptionNames::NoCigar])}
    , CombineBlocks{options[MapCLROptionNames::CombineBlocks]}
{
    OutFormat = ParseOverlapWriterFormat(std::string(options[MapCLROptionNames::OutFormat]));
    if (OutFormat == OverlapWriterFormat::Unknown) {
        throw std::runtime_error("Unknown output format: '" +
                                 std::string(options[MapCLROptionNames::OutFormat]) + "'.");
    }

    // Indexing options.
    MapperSettings.map.seedParams.KmerSize = options[MapCLROptionNames::PrimaryKmerSize];
    MapperSettings.map.seedParams.MinimizerWindow =
        options[MapCLROptionNames::PrimaryMinimizerWindow];
    MapperSettings.map.seedParams.Spacing = options[MapCLROptionNames::PrimarySpacing];
    MapperSettings.map.seedParams.UseHPCForSeedsOnly =
        options[MapCLROptionNames::PrimaryUseHomopolymerCompression];
    MapperSettings.map.seedParams.UseRC = options[MapCLROptionNames::PrimaryUseReverseComplement];
    MapperSettings.map.seedParamsFallback.KmerSize = options[MapCLROptionNames::FallbackKmerSize];
    MapperSettings.map.seedParamsFallback.MinimizerWindow =
        options[MapCLROptionNames::FallbackMinimizerWindow];
    MapperSettings.map.seedParamsFallback.Spacing = options[MapCLROptionNames::FallbackSpacing];
    MapperSettings.map.seedParamsFallback.UseHPCForSeedsOnly =
        options[MapCLROptionNames::FallbackUseHomopolymerCompression];
    MapperSettings.map.seedParamsFallback.UseRC =
        options[MapCLROptionNames::FallbackUseReverseComplement];

    // Initialize the MapperCLRMapSettings.
    MapperSettings.map.freqPercentile = options[MapCLROptionNames::FreqPercentile];
    MapperSettings.map.seedOccurrenceMin = options[MapCLROptionNames::SeedOccurrenceMin];
    MapperSettings.map.seedOccurrenceMax = options[MapCLROptionNames::SeedOccurrenceMax];
    MapperSettings.map.seedOccurrenceMaxMemory =
        options[MapCLROptionNames::SeedOccurrenceMaxMemory];
    MapperSettings.map.chainMaxSkip = options[MapCLROptionNames::ChainMaxSkip];
    MapperSettings.map.chainMaxPredecessors = options[MapCLROptionNames::ChainMaxPredecessors];
    MapperSettings.map.chainBandwidth = options[MapCLROptionNames::ChainBandwidth];
    MapperSettings.map.seedJoinDist = options[MapCLROptionNames::SeedJoinDist];
    MapperSettings.map.longMergeBandwidth = options[MapCLROptionNames::LongMergeBandwidth];
    MapperSettings.map.minNumSeeds = options[MapCLROptionNames::MinNumSeeds];
    MapperSettings.map.minCoveredBases = options[MapCLROptionNames::MinCoveredBases];
    MapperSettings.map.minDPScore = options[MapCLROptionNames::MinDPScore];
    MapperSettings.map.secondaryAllowedOverlapFractionQuery =
        options[MapCLROptionNames::SecondaryAllowedOverlapFractionQuery];
    MapperSettings.map.secondaryAllowedOverlapFractionTarget =
        options[MapCLROptionNames::SecondaryAllowedOverlapFractionTarget];
    MapperSettings.map.secondaryMinScoreFraction =
        options[MapCLROptionNames::SecondaryMinScoreFraction];
    MapperSettings.map.useLIS = (!options[MapCLROptionNames::NoLIS]);
    MapperSettings.map.skipSymmetricOverlaps = options[MapCLROptionNames::SkipSymmetricOverlaps];
    MapperSettings.map.selfHitPolicy = MapperSelfHitPolicyFromString(
        std::string(options[MapCLROptionNames::SelfHitPolicyMapping]));
    MapperSettings.map.minQueryLen = options[MapCLROptionNames::MinQueryLen];
    MapperSettings.map.bestNSecondary = options[MapCLROptionNames::BestNSecondary];
    MapperSettings.map.maxFlankExtensionDist =
        options[MapCLROptionNames::MaxFlankExtensionDistance];
    MapperSettings.map.flankExtensionFactor = options[MapCLROptionNames::FlankExtensionFactor];
    MapperSettings.map.minAlignmentSpan = options[MapCLROptionNames::MinAlignmentSpan];

    // Refining seed hits.
    MapperSettings.map.refineSeedHits = !options[MapCLROptionNames::NoRefineSeedHits];
    MapperSettings.map.refineMinGap1 = options[MapCLROptionNames::RefineMinGap1];
    MapperSettings.map.refineDiffThreshold = options[MapCLROptionNames::RefineDiffThreshold];
    MapperSettings.map.refineMinGap2 = options[MapCLROptionNames::RefineMinGap2];

    // Initialize the MapperCLRAlignSettings.
    MapperSettings.align.align = options[MapCLROptionNames::Align];
    MapperSettings.align.selfHitPolicy = MapperSelfHitPolicyFromString(
        std::string(options[MapCLROptionNames::SelfHitPolicyAlignment]));
    MapperSettings.align.alignerTypeGlobal =
        AlignerTypeFromString(options[MapCLROptionNames::AlignerTypeGlobal]);
    MapperSettings.align.alignerTypeExt =
        AlignerTypeFromString(options[MapCLROptionNames::AlignerTypeExt]);

    // Alignment parameters, global.
    MapperSettings.align.alnParamsGlobal.zdrop = options[MapCLROptionNames::AlignerGlobalZdrop];
    MapperSettings.align.alnParamsGlobal.zdrop2 = options[MapCLROptionNames::AlignerGlobalZdrop2];
    MapperSettings.align.alnParamsGlobal.alignBandwidth =
        options[MapCLROptionNames::AlignerGlobalBandwidth];
    MapperSettings.align.alnParamsGlobal.matchScore =
        options[MapCLROptionNames::AlignerGlobalMatchScore];
    MapperSettings.align.alnParamsGlobal.mismatchPenalty =
        options[MapCLROptionNames::AlignerGlobalMismatchPenalty];
    MapperSettings.align.alnParamsGlobal.gapOpen1 =
        options[MapCLROptionNames::AlignerGlobalGapOpen1];
    MapperSettings.align.alnParamsGlobal.gapExtend1 =
        options[MapCLROptionNames::AlignerGlobalGapExtend1];
    MapperSettings.align.alnParamsGlobal.gapOpen2 =
        options[MapCLROptionNames::AlignerGlobalGapOpen2];
    MapperSettings.align.alnParamsGlobal.gapExtend2 =
        options[MapCLROptionNames::AlignerGlobalGapExtend2];
    MapperSettings.align.alnParamsGlobal.dynamicBandwidth =
        options[MapCLROptionNames::AlignerGlobalDynamicBandwidth];

    // Alignment parameters, extension.
    MapperSettings.align.alnParamsExt.zdrop = options[MapCLROptionNames::AlignerExtZdrop];
    MapperSettings.align.alnParamsExt.zdrop2 = options[MapCLROptionNames::AlignerExtZdrop2];
    MapperSettings.align.alnParamsExt.alignBandwidth =
        options[MapCLROptionNames::AlignerExtBandwidth];
    MapperSettings.align.alnParamsExt.endBonus = options[MapCLROptionNames::AlignerExtEndBonus];
    MapperSettings.align.alnParamsExt.matchScore = options[MapCLROptionNames::AlignerExtMatchScore];
    MapperSettings.align.alnParamsExt.mismatchPenalty =
        options[MapCLROptionNames::AlignerExtMismatchPenalty];
    MapperSettings.align.alnParamsExt.gapOpen1 = options[MapCLROptionNames::AlignerExtGapOpen1];
    MapperSettings.align.alnParamsExt.gapExtend1 = options[MapCLROptionNames::AlignerExtGapExtend1];
    MapperSettings.align.alnParamsExt.gapOpen2 = options[MapCLROptionNames::AlignerExtGapOpen2];
    MapperSettings.align.alnParamsExt.gapExtend2 = options[MapCLROptionNames::AlignerExtGapExtend2];
    MapperSettings.align.alnParamsExt.dynamicBandwidth =
        options[MapCLROptionNames::AlignerExtDynamicBandwidth];
}

PacBio::CLI_v2::Interface MapCLRSettings::CreateCLI()
{
    PacBio::CLI_v2::Interface i{"pancake", "Mapping and alignment of PacBio reads.",
                                PacBio::Pancake::PancakeFormattedVersion()};

    // clang-format off
    i.AddOptionGroup("Input/Output Options", {
        MapCLROptionNames::OutFormat,
        MapCLROptionNames::WriteIds,
        MapCLROptionNames::NoCigar,
        MapCLROptionNames::CombineBlocks,
        MapCLROptionNames::TargetBlockId,
        MapCLROptionNames::QueryBlockStartId,
        MapCLROptionNames::QueryBlockEndId
    });

    i.AddOptionGroup("Indexing options (ignored for SeedDB input)", {
        MapCLROptionNames::PrimaryKmerSize,
        MapCLROptionNames::PrimaryMinimizerWindow,
        MapCLROptionNames::PrimarySpacing,
        MapCLROptionNames::PrimaryUseHomopolymerCompression,
        MapCLROptionNames::PrimaryUseReverseComplement,
        MapCLROptionNames::FallbackKmerSize,
        MapCLROptionNames::FallbackMinimizerWindow,
        MapCLROptionNames::FallbackSpacing,
        MapCLROptionNames::FallbackUseHomopolymerCompression,
        MapCLROptionNames::FallbackUseReverseComplement,
    });

    i.AddOptionGroup("Mapping Options", {
        MapCLROptionNames::ChainMaxSkip,
        MapCLROptionNames::ChainMaxPredecessors,
        MapCLROptionNames::ChainBandwidth,
        MapCLROptionNames::SeedJoinDist,
        MapCLROptionNames::LongMergeBandwidth,
        MapCLROptionNames::MinNumSeeds,
        MapCLROptionNames::MinCoveredBases,
        MapCLROptionNames::MinDPScore,
        MapCLROptionNames::SecondaryAllowedOverlapFractionQuery,
        MapCLROptionNames::SecondaryAllowedOverlapFractionTarget,
        MapCLROptionNames::SecondaryMinScoreFraction,
        MapCLROptionNames::NoLIS,
        MapCLROptionNames::FreqPercentile,
        MapCLROptionNames::SeedOccurrenceMin,
        MapCLROptionNames::SeedOccurrenceMax,
        MapCLROptionNames::SeedOccurrenceMaxMemory,
        MapCLROptionNames::SelfHitPolicyMapping,
        MapCLROptionNames::NoRefineSeedHits,
        MapCLROptionNames::RefineMinGap1,
        MapCLROptionNames::RefineDiffThreshold,
        MapCLROptionNames::RefineMinGap2,
    });
    i.AddOptionGroup("Alignment Options", {
        MapCLROptionNames::Align,
        MapCLROptionNames::SelfHitPolicyAlignment,
        MapCLROptionNames::MaxFlankExtensionDistance,
        MapCLROptionNames::FlankExtensionFactor,
        MapCLROptionNames::MinAlignmentSpan,
        MapCLROptionNames::AlignerTypeGlobal,
        MapCLROptionNames::AlignerTypeExt,

        MapCLROptionNames::AlignerGlobalZdrop,
        MapCLROptionNames::AlignerGlobalZdrop2,
        MapCLROptionNames::AlignerGlobalBandwidth,
        MapCLROptionNames::AlignerGlobalMatchScore,
        MapCLROptionNames::AlignerGlobalMismatchPenalty,
        MapCLROptionNames::AlignerGlobalGapOpen1,
        MapCLROptionNames::AlignerGlobalGapExtend1,
        MapCLROptionNames::AlignerGlobalGapOpen2,
        MapCLROptionNames::AlignerGlobalGapExtend2,
        MapCLROptionNames::AlignerGlobalDynamicBandwidth,

        MapCLROptionNames::AlignerExtZdrop,
        MapCLROptionNames::AlignerExtZdrop2,
        MapCLROptionNames::AlignerExtBandwidth,
        MapCLROptionNames::AlignerExtEndBonus,
        MapCLROptionNames::AlignerExtMatchScore,
        MapCLROptionNames::AlignerExtMismatchPenalty,
        MapCLROptionNames::AlignerExtGapOpen1,
        MapCLROptionNames::AlignerExtGapExtend1,
        MapCLROptionNames::AlignerExtGapOpen2,
        MapCLROptionNames::AlignerExtGapExtend2,
        MapCLROptionNames::AlignerExtDynamicBandwidth,

    });
    i.AddOptionGroup("Filtering Options", {
        MapCLROptionNames::MinQueryLen,
        MapCLROptionNames::BestNSecondary,
        MapCLROptionNames::SkipSymmetricOverlaps,
    });
    i.AddPositionalArguments({
        MapCLROptionNames::TargetDBPrefix,
        MapCLROptionNames::QueryDBPrefix,
    });

    // clang-format on
    return i;
}
}  // namespace Pancake
}  // namespace PacBio
