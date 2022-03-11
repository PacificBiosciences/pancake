// Author: Ivan Sovic

#include <pancake/OverlapHifiSettings.hpp>
#include <pancake/Version.hpp>

namespace PacBio {
namespace Pancake {
namespace OverlapHiFiOptionNames {

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

const CLI_v2::PositionalArgument TargetBlockId {
R"({
    "name" : "target_block",
    "type" : "int",
    "description" : "Block ID from the target DB. Queries will be mapped only onto this block."
})"};

const CLI_v2::PositionalArgument QueryBlockStartId {
R"({
    "name" : "query_block_start",
    "type" : "int",
    "description" : "Start block ID for a range of blocks to map. Zero based."
})"};

const CLI_v2::PositionalArgument QueryBlockEndId {
R"({
    "name" : "query_block_end",
    "type" : "int",
    "description" : "Start block ID for a range of blocks to map. Zero based, non-inclusive. Value == 0 runs until the end block."
})"};



const CLI_v2::Option OutFormat{
R"({
    "names" : ["out-fmt"],
    "choices" : ["m4", "ipa", "paf", "sam"],
    "type" : "string",
    "default" : "m4",
    "description" : "Select the output format."
})", std::string("m4")};

const CLI_v2::Option FreqPercentile{
R"({
    "names" : ["freq-percentile"],
    "description" : "Filter frequent kmers.",
    "type" : "double"
})", OverlapHifiSettings::Defaults::FreqPercentile};

const CLI_v2::Option MinQueryLen{
R"({
    "names" : ["min-qlen"],
    "description" : "Ignore queries shorter than this.",
    "type" : "int"
})", OverlapHifiSettings::Defaults::MinQueryLen};

const CLI_v2::Option MinTargetLen{
R"({
    "names" : ["min-tlen"],
    "description" : "Ignore targets shorter than this.",
    "type" : "int"
})", OverlapHifiSettings::Defaults::MinTargetLen};

const CLI_v2::Option MaxSeedDistance{
R"({
    "names" : ["max-seed-dist"],
    "description" : "Maximum distance between two seeds to join into an anchor.",
    "type" : "int"
})", OverlapHifiSettings::Defaults::MaxSeedDistance};

const CLI_v2::Option MinNumSeeds{
R"({
    "names" : ["min-num-seeds"],
    "description" : "Minimum number of seeds in an anchor.",
    "type" : "int"
})", OverlapHifiSettings::Defaults::MinNumSeeds};

const CLI_v2::Option MinCoveredBases{
R"({
    "names" : ["min-cov-bases"],
    "description" : "Minimum number of bases covered by kmers in an anchor.",
    "type" : "int"
})", OverlapHifiSettings::Defaults::MinCoveredBases};

const CLI_v2::Option MinChainSpan{
R"({
    "names" : ["min-anchor-span"],
    "description" : "Minimum chain span to retain it.",
    "type" : "int"
})", OverlapHifiSettings::Defaults::MinChainSpan};

const CLI_v2::Option ChainBandwidth{
R"({
    "names" : ["chain-bw"],
    "description" : "Diagonal bandwidth to merge seeds into chains.",
    "type" : "int"
})", OverlapHifiSettings::Defaults::ChainBandwidth};

const CLI_v2::Option AlignmentBandwidth{
R"({
    "names" : ["aln-bw"],
    "description" : "Bandwidth for alignment, fraction of the query span.",
    "type" : "double"
})", OverlapHifiSettings::Defaults::AlignmentBandwidth};

const CLI_v2::Option AlignmentMaxD{
R"({
    "names" : ["aln-diff-rate"],
    "description" : "Expected maximum diff rate between sequences.",
    "type" : "double"
})", OverlapHifiSettings::Defaults::AlignmentMaxD};

const CLI_v2::Option MinIdentity{
R"({
    "names" : ["min-idt"],
    "description" : "Minimum percent alignment identity allowed to report the alignment. This is an overall threshold which takes into account both indels and SNPs.",
    "type" : "double"
})", OverlapHifiSettings::Defaults::MinIdentity};

const CLI_v2::Option NoSNPsInIdentity{
R"({
    "names" : ["no-snps"],
    "description" : "Ignore SNPs when computing the identity for an overlap. This only works in the traceback mode.",
    "type" : "bool"
})", OverlapHifiSettings::Defaults::NoSNPsInIdentity};

const CLI_v2::Option NoIndelsInIdentity{
R"({
    "names" : ["no-indels"],
    "description" : "Ignore indels when computing the identity for an overlap. This only works in the traceback mode.",
    "type" : "bool"
})", OverlapHifiSettings::Defaults::NoIndelsInIdentity};

const CLI_v2::Option MinMappedLength{
R"({
    "names" : ["min-map-len"],
    "description" : "Output only alignments above this length.",
    "type" : "int"
})", OverlapHifiSettings::Defaults::MinMappedLength};

const CLI_v2::Option SkipSymmetricOverlaps{
R"({
    "names" : ["skip-sym"],
    "description" : "If Aid < Bid, only compute overlap Aid->Bid and skip computing overlap for Bid->Aid.",
    "type" : "bool"
})", OverlapHifiSettings::Defaults::SkipSymmetricOverlaps};

const CLI_v2::Option AllowSelfHits{
R"({
    "names" : ["allow-self-hits"],
    "description" : "If both the query and the target DBs are the same and Aid == Bid then this is a self-hit. This option enables the output of such overlaps.",
    "type" : "bool"
})", false};

const CLI_v2::Option OneHitPerTarget{
R"({
    "names" : ["one-hit-per-target"],
    "description" : "Allow only one alignment per query/target pair.",
    "type" : "bool"
})", OverlapHifiSettings::Defaults::OneHitPerTarget};

const CLI_v2::Option SmartHitPerTarget{
R"({
    "names" : ["smart-hit-per-target"],
    "description" : "Allow supplementary alignments on the same target (e.g. circular overlaps). Secondary alignments on the same target are not allowed.",
    "type" : "bool"
})", OverlapHifiSettings::Defaults::SmartHitPerTarget};

const CLI_v2::Option WriteReverseOverlaps{
R"({
    "names" : ["write-rev"],
    "description" : "For eveery overlap, write out its reverse complement too.",
    "type" : "bool"
})", OverlapHifiSettings::Defaults::WriteReverseOverlaps};
const CLI_v2::Option WriteIds{
R"({
    "names" : ["write-ids"],
    "description" : "Output overlaps will contain numeric IDs for the A and B reads (instead of names).",
    "type" : "bool"
})", OverlapHifiSettings::Defaults::WriteIds};
const CLI_v2::Option WriteCigar{
R"({
    "names" : ["write-cigar"],
    "description" : "Write the CIGAR string if the sensitive alignment mode is applied.",
    "type" : "bool"
})", OverlapHifiSettings::Defaults::WriteCigar};

const CLI_v2::Option AllowedDovetailDist{
R"({
    "names" : ["dt-dist"],
    "description" : "Allowed distance of an overlap from the beginning of the sequences to call the overlap a dovetail.",
    "type" : "int"
})", OverlapHifiSettings::Defaults::AllowedDovetailDist};

const CLI_v2::Option AllowedHeuristicExtendDist{
R"({
    "names" : ["ext-dist"],
    "description" : "Heuristically modify the coordinats of an overlap into a dovetail overlap if are within this distance from the edges of the reads.",
    "type" : "int"
})", OverlapHifiSettings::Defaults::AllowedHeuristicExtendDist};

const CLI_v2::Option CombineBlocks{
R"({
    "names" : ["combine"],
    "description" : "Combines this many query blocks into one larger block for processing.",
    "type" : "int"
})", OverlapHifiSettings::Defaults::CombineBlocks};

const CLI_v2::Option BestN{
R"({
    "names" : ["bestn"],
    "description" : "Output only best N alignments.",
    "type" : "int"
})", OverlapHifiSettings::Defaults::BestN};

const CLI_v2::Option UseHPC{
R"({
    "names" : ["use-hpc"],
    "description" : "Enable homopolymer compression."
})", OverlapHifiSettings::Defaults::UseHPC};

const CLI_v2::Option UseTraceback{
R"({
    "names" : ["traceback"],
    "description" : "Run alignment traceback and compute mismatches.",
    "type" : "bool"
})", OverlapHifiSettings::Defaults::UseTraceback};

const CLI_v2::Option MaskHomopolymers{
R"({
    "names" : ["mask-hp"],
    "description" : "Mask homopolymer errors when traceback is generated. This will impact identity calculation.",
    "type" : "bool"
})", OverlapHifiSettings::Defaults::MaskHomopolymers};

const CLI_v2::Option MaskSimpleRepeats{
R"({
    "names" : ["mask-repeats"],
    "description" : "Mask indels in simple exact repeats when traceback is generated. This will impact identity calculation.",
    "type" : "bool"
})", OverlapHifiSettings::Defaults::MaskSimpleRepeats};

const CLI_v2::Option MaskHomopolymerSNPs{
R"({
    "names" : ["mask-hp-snps"],
    "description" : "Mask mismatches which occur in the homopolymer sequences. Applied only when traceback is generated. This will impact identity calculation.",
    "type" : "bool"
})", OverlapHifiSettings::Defaults::MaskHomopolymerSNPs};

const CLI_v2::Option MaskHomopolymersArbitrary{
R"({
    "names" : ["mask-hp-arbitrary"],
    "description" : "Allows arbitrary bases to be inserted into the HP stretches (the bases don't have to match the HP). Only used in combination with '--mask-hp'.",
    "type" : "bool"
})", OverlapHifiSettings::Defaults::MaskHomopolymersArbitrary};

const CLI_v2::Option MarkSecondary{
R"({
    "names" : ["mark-secondary"],
    "description" : "Mask homopolymer errors when traceback is generated. This will impact identity calculation.",
    "type" : "bool"
})", OverlapHifiSettings::Defaults::MarkSecondary};

const CLI_v2::Option SecondaryAllowedOverlapFraction{
R"({
    "names" : ["secondary-min-ovl-frac"],
    "description" : "Minimum amount overlap with primary alignment to call an alignment secondary.",
    "type" : "double"
})", OverlapHifiSettings::Defaults::SecondaryAllowedOverlapFraction};

const CLI_v2::Option SecondaryMinScoreFraction{
R"({
    "names" : ["secondary-min-score-frac"],
    "description" : "Minimum secondary-to-primary score ratio.",
    "type" : "double"
})", OverlapHifiSettings::Defaults::SecondaryMinScoreFraction};

const CLI_v2::Option TrimAlignment{
R"({
    "names" : ["trim"],
    "description" : "Applies window-based trimming of the front and end of the alignment. Can be used only in combination with '--traceback'.",
    "type" : "bool"
})", OverlapHifiSettings::Defaults::TrimAlignment};

const CLI_v2::Option TrimWindowSize{
R"({
    "names" : ["trim-window-size"],
    "description" : "Window size for trimming.",
    "type" : "int"
})", OverlapHifiSettings::Defaults::TrimWindowSize};

const CLI_v2::Option TrimWindowMatchFraction{
R"({
    "names" : ["trim-match-frac"],
    "description" : "Minimum fraction in a trimming window of match bases to stop trimming.",
    "type" : "double"
})", OverlapHifiSettings::Defaults::TrimWindowMatchFraction};

const CLI_v2::Option TrimToFirstMatch{
R"({
    "names" : ["trim-to-first-match"],
    "description" : "When trimming is applied, this option ensures that the first non-trimmed base will be a match operation. Can be used only in combination with '--trim'.",
    "type" : "bool"
})", OverlapHifiSettings::Defaults::TrimToFirstMatch};

// clang-format on

}  // namespace OverlapHiFiOptionNames

OverlapHifiSettings::OverlapHifiSettings() = default;

OverlapWriterFormat ParseOutFormat(const std::string& val)
{
    if (val == "m4") {
        return OverlapWriterFormat::M4;
    } else if (val == "ipa") {
        return OverlapWriterFormat::IPAOvl;
    } else if (val == "paf") {
        return OverlapWriterFormat::PAF;
    } else if (val == "sam") {
        return OverlapWriterFormat::SAM;
    }
    return OverlapWriterFormat::Unknown;
}

OverlapHifiSettings::OverlapHifiSettings(const PacBio::CLI_v2::Results& options)
    : TargetDBPrefix{options[OverlapHiFiOptionNames::TargetDBPrefix]}
    , QueryDBPrefix{options[OverlapHiFiOptionNames::QueryDBPrefix]}
    , NumThreads{options.NumThreads()}
    , TargetBlockId{std::stoi(options[OverlapHiFiOptionNames::TargetBlockId])}
    , QueryBlockStartId{std::stoi(options[OverlapHiFiOptionNames::QueryBlockStartId])}
    , QueryBlockEndId{std::stoi(options[OverlapHiFiOptionNames::QueryBlockEndId])}

    , FreqPercentile{options[OverlapHiFiOptionNames::FreqPercentile]}
    , MinQueryLen{options[OverlapHiFiOptionNames::MinQueryLen]}
    , MinTargetLen{options[OverlapHiFiOptionNames::MinTargetLen]}
    , MaxSeedDistance{options[OverlapHiFiOptionNames::MaxSeedDistance]}
    , MinNumSeeds{options[OverlapHiFiOptionNames::MinNumSeeds]}
    , MinCoveredBases{options[OverlapHiFiOptionNames::MinCoveredBases]}
    , MinChainSpan{options[OverlapHiFiOptionNames::MinChainSpan]}
    , ChainBandwidth{options[OverlapHiFiOptionNames::ChainBandwidth]}
    , AlignmentBandwidth{options[OverlapHiFiOptionNames::AlignmentBandwidth]}
    , AlignmentMaxD{options[OverlapHiFiOptionNames::AlignmentMaxD]}
    , MinIdentity{options[OverlapHiFiOptionNames::MinIdentity]}
    , NoSNPsInIdentity{options[OverlapHiFiOptionNames::NoSNPsInIdentity]}
    , NoIndelsInIdentity{options[OverlapHiFiOptionNames::NoIndelsInIdentity]}
    , MinMappedLength{options[OverlapHiFiOptionNames::MinMappedLength]}
    , OneHitPerTarget{options[OverlapHiFiOptionNames::OneHitPerTarget]}
    , SmartHitPerTarget{options[OverlapHiFiOptionNames::SmartHitPerTarget]}
    , WriteReverseOverlaps{options[OverlapHiFiOptionNames::WriteReverseOverlaps]}
    , WriteIds{options[OverlapHiFiOptionNames::WriteIds]}
    , WriteCigar{options[OverlapHiFiOptionNames::WriteCigar]}
    , AllowedDovetailDist{options[OverlapHiFiOptionNames::AllowedDovetailDist]}
    , AllowedHeuristicExtendDist{options[OverlapHiFiOptionNames::AllowedHeuristicExtendDist]}
    , CombineBlocks{options[OverlapHiFiOptionNames::CombineBlocks]}
    , BestN{options[OverlapHiFiOptionNames::BestN]}
    , UseHPC{options[OverlapHiFiOptionNames::UseHPC]}
    , UseTraceback{options[OverlapHiFiOptionNames::UseTraceback]}
    , MaskHomopolymers{options[OverlapHiFiOptionNames::MaskHomopolymers]}
    , MaskSimpleRepeats{options[OverlapHiFiOptionNames::MaskSimpleRepeats]}
    , MaskHomopolymerSNPs{options[OverlapHiFiOptionNames::MaskHomopolymerSNPs]}
    , MaskHomopolymersArbitrary{options[OverlapHiFiOptionNames::MaskHomopolymersArbitrary]}
    , MarkSecondary{options[OverlapHiFiOptionNames::MarkSecondary]}
    , SecondaryAllowedOverlapFraction{options
                                          [OverlapHiFiOptionNames::SecondaryAllowedOverlapFraction]}
    , SecondaryMinScoreFraction{options[OverlapHiFiOptionNames::SecondaryMinScoreFraction]}
    , TrimAlignment{options[OverlapHiFiOptionNames::TrimAlignment]}
    , TrimWindowSize{options[OverlapHiFiOptionNames::TrimWindowSize]}
    , TrimWindowMatchFraction{options[OverlapHiFiOptionNames::TrimWindowMatchFraction]}
    , TrimToFirstMatch{options[OverlapHiFiOptionNames::TrimToFirstMatch]}
{
    if ((NoSNPsInIdentity || NoIndelsInIdentity || MaskHomopolymers || MaskSimpleRepeats ||
         MaskHomopolymerSNPs || MaskHomopolymersArbitrary) &&
        (UseTraceback == false)) {
        throw std::runtime_error(
            "The '--no-snps', '--no-indels', '--mask-hp' and '--mask-rep' can only be used "
            "The '--no-snps', '--no-indels', '--mask-hp', '--mask-hp-snps', '--mask-hp-arbitrary' "
            "and '--mask-rep' can "
            "only be used "
            "together with the '--traceback' "
            "option.");
    }
    if (NoSNPsInIdentity && NoIndelsInIdentity) {
        PBLOG_WARN << "Both --no-snps and --no-indels options are specified, which means that all "
                      "identity values will be 100%.";
    }
    if (MaskHomopolymersArbitrary == true && MaskHomopolymers == false) {
        throw std::runtime_error(
            "Option '--mask-hp-arbitrary' can only be used together with '--mask-hp'.");
    }

    OutFormat = ParseOutFormat(options[OverlapHiFiOptionNames::OutFormat]);
    if (OutFormat == OverlapWriterFormat::Unknown) {
        throw std::runtime_error("Unknown output format: '" +
                                 std::string(options[OverlapHiFiOptionNames::OutFormat]) + "'.");
    }

    SkipSelfHits = false;
    if (static_cast<bool>(options[OverlapHiFiOptionNames::AllowSelfHits]) == false &&
        QueryDBPrefix == TargetDBPrefix) {
        SkipSelfHits = true;
    }

    SkipSymmetricOverlaps = false;
    if (static_cast<bool>(options[OverlapHiFiOptionNames::SkipSymmetricOverlaps]) &&
        QueryDBPrefix == TargetDBPrefix) {
        SkipSymmetricOverlaps = true;
    }

    if (TrimToFirstMatch == true && TrimAlignment == false) {
        throw std::runtime_error(
            "The '--trim-to-first-match' option can only be used when '--trim' is specified.");
    }
}

PacBio::CLI_v2::Interface OverlapHifiSettings::CreateCLI()
{
    PacBio::CLI_v2::Interface i{"pancake ovl-hifi", "HiFi overlapping.",
                                PacBio::Pancake::PancakeFormattedVersion()};

    // clang-format off
    i.AddOptionGroup("Input/Output Options", {
        OverlapHiFiOptionNames::OutFormat,
    });
    i.AddOptionGroup("Algorithm Options", {
        OverlapHiFiOptionNames::FreqPercentile,
        OverlapHiFiOptionNames::MinQueryLen,
        OverlapHiFiOptionNames::MinTargetLen,
        OverlapHiFiOptionNames::MaxSeedDistance,
        OverlapHiFiOptionNames::MinNumSeeds,
        OverlapHiFiOptionNames::MinCoveredBases,
        OverlapHiFiOptionNames::MinChainSpan,
        OverlapHiFiOptionNames::ChainBandwidth,
        OverlapHiFiOptionNames::AlignmentBandwidth,
        OverlapHiFiOptionNames::AlignmentMaxD,
        OverlapHiFiOptionNames::MinIdentity,
        OverlapHiFiOptionNames::NoSNPsInIdentity,
        OverlapHiFiOptionNames::NoIndelsInIdentity,
        OverlapHiFiOptionNames::MinMappedLength,
        OverlapHiFiOptionNames::SkipSymmetricOverlaps,
        OverlapHiFiOptionNames::AllowSelfHits,
        OverlapHiFiOptionNames::OneHitPerTarget,
        OverlapHiFiOptionNames::SmartHitPerTarget,
        OverlapHiFiOptionNames::WriteReverseOverlaps,
        OverlapHiFiOptionNames::WriteIds,
        OverlapHiFiOptionNames::WriteCigar,
        OverlapHiFiOptionNames::AllowedDovetailDist,
        OverlapHiFiOptionNames::AllowedHeuristicExtendDist,
        OverlapHiFiOptionNames::CombineBlocks,
        OverlapHiFiOptionNames::BestN,
        OverlapHiFiOptionNames::UseHPC,
        OverlapHiFiOptionNames::UseTraceback,
        OverlapHiFiOptionNames::MaskHomopolymers,
        OverlapHiFiOptionNames::MaskSimpleRepeats,
        OverlapHiFiOptionNames::MaskHomopolymerSNPs,
        OverlapHiFiOptionNames::MaskHomopolymersArbitrary,
        OverlapHiFiOptionNames::MarkSecondary,
        OverlapHiFiOptionNames::SecondaryAllowedOverlapFraction,
        OverlapHiFiOptionNames::SecondaryMinScoreFraction,
        OverlapHiFiOptionNames::TrimAlignment,
        OverlapHiFiOptionNames::TrimWindowSize,
        OverlapHiFiOptionNames::TrimWindowMatchFraction,
        OverlapHiFiOptionNames::TrimToFirstMatch,
    });
    i.AddPositionalArguments({
        OverlapHiFiOptionNames::TargetDBPrefix,
        OverlapHiFiOptionNames::QueryDBPrefix,
        OverlapHiFiOptionNames::TargetBlockId,
        OverlapHiFiOptionNames::QueryBlockStartId,
        OverlapHiFiOptionNames::QueryBlockEndId
    });

    // clang-format on
    return i;
}
}  // namespace Pancake
}  // namespace PacBio
