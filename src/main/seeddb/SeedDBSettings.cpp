// Author: Ivan Sovic

#include "SeedDBSettings.h"
#include <pacbio/Version.h>
#include <limits>

namespace PacBio {
namespace Pancake {
namespace SeedDB {
namespace OptionNames {

// clang-format off

const CLI_v2::PositionalArgument InputFile {
R"({
    "name" : "input.seqdb",
    "description" : "Path to the SeqDB to process."
})"};

const CLI_v2::PositionalArgument OutputPrefix {
R"({
    "name" : "prefix",
    "description" : "The prefix of the output SeedDB files."
})"};

const CLI_v2::Option SplitBlocks{
R"({
    "names" : ["split-blocks"],
    "description" : "Write seeds for each block into a separate file."
})", SeedDBSettings::Defaults::SplitBlocks};

const CLI_v2::Option KmerSize{
R"({
    "names" : ["k", "kmer-size"],
    "type" : "int",
    "default" : 28,
    "description" : "Kmer size for indexing. Maximum size is 28."
})", SeedDBSettings::Defaults::KmerSize};

const CLI_v2::Option MinimizerWindow{
R"({
    "names" : ["w", "window"],
    "type" : "int",
    "default" : 80,
    "description" : "Minimizer window size for indexing. Maximum size is 512."
})", SeedDBSettings::Defaults::MinimizerWindow};

const CLI_v2::Option Spacing{
R"({
    "names" : ["s", "space"],
    "type" : "int",
    "default" : 0,
    "description" : "Seed spacing. Maximum spacing is 32."
})", SeedDBSettings::Defaults::Spacing};

const CLI_v2::Option UseHPC{
R"({
    "names" : ["use-hpc"],
    "description" : "Enable homopolymer compression. This option compresses the input sequences from the SeqDB, and the resulting seeds are in the compressed space."
})", SeedDBSettings::Defaults::UseHPC};

const CLI_v2::Option UseHPCForSeedsOnly{
R"({
    "names" : ["use-hpc-seeds-only"],
    "description" : "Enable homopolymer compression of seeds. Unlike '--use-hpc', here the input sequences are not physically compressed. Instead, the seed coordinates correspond to the original uncompressed sequences, but the kmers skip HP bases with their span."
})", SeedDBSettings::Defaults::UseHPCForSeedsOnly};

const CLI_v2::Option NoRevCmp{
R"({
    "names" : ["no-rc"],
    "description" : "Do not produce seeds from the reverse complement strand."
})", SeedDBSettings::Defaults::NoRevCmp};

// clang-format on

}  // namespace OptionNames

SeedDBSettings::SeedDBSettings() = default;

SeedDBSettings::SeedDBSettings(const PacBio::CLI_v2::Results& options)
    : InputFile{options[OptionNames::InputFile]}
    , OutputPrefix{options[OptionNames::OutputPrefix]}
    , NumThreads{options.NumThreads()}
    , SplitBlocks{options[OptionNames::SplitBlocks]}
    , SeedParameters{options[OptionNames::KmerSize],
                     options[OptionNames::MinimizerWindow],
                     options[OptionNames::Spacing],
                     options[OptionNames::UseHPC],
                     options[OptionNames::UseHPCForSeedsOnly],
                     !options[OptionNames::NoRevCmp]}
{
}

PacBio::CLI_v2::Interface SeedDBSettings::CreateCLI()
{
    PacBio::CLI_v2::Interface i{"pancake", "Compute seeds from a SeqDB.",
                                PacBio::Pancake::PancakeFormattedVersion()};

    // clang-format off
    i.AddOptionGroup("Algorithm Options", {
        OptionNames::SplitBlocks,
        OptionNames::KmerSize,
        OptionNames::MinimizerWindow,
        OptionNames::Spacing,
        OptionNames::UseHPC,
        OptionNames::UseHPCForSeedsOnly,
        OptionNames::NoRevCmp,
    });
    i.AddPositionalArguments({
        OptionNames::InputFile,
        OptionNames::OutputPrefix,
    });

    // clang-format on
    return i;
}

}  // namespace SeedDB
}  // namespace Pancake
}  // namespace PacBio