// Author: Ivan Sovic

#include "DBFilterSettings.hpp"

#include <pancake/LibraryInfo.hpp>

#include <limits>

namespace PacBio {
namespace Pancake {
namespace DBFilterOptionNames {

// clang-format off

const CLI_v2::PositionalArgument InputPrefix {
R"({
    "name" : "in_prefix",
    "description" : "The prefix of the input DB files."
})"};

const CLI_v2::PositionalArgument OutputPrefix {
R"({
    "name" : "out_prefix",
    "description" : "The prefix of the output DB files."
})"};

const CLI_v2::Option Sampling{
R"({
    "names" : ["sampling"],
    "choices" : ["none", "linear", "random"],
    "type" : "string",
    "default" : "none",
    "description" : "Select sampling type: none, linear, random."
})", std::string("none")};

const CLI_v2::Option SampleBases{
R"({
    "names" : ["sample-bases"],
    "description" : "Number of bases to sample.",
    "type" : "int"
})", DBFilterSettings::Defaults::SampleBases};

const CLI_v2::Option BlockSize{
R"({
    "names" : ["block-size"],
    "description" : "Block size in megabases. Value 0 means one sequnece per block, value < 0 all sequences in one block.",
    "type" : "float"
})", DBFilterSettings::Defaults::BlockSize};

const CLI_v2::Option RandomSeed{
R"({
    "names" : ["random-seed"],
    "description" : "Random seed for sampling.",
    "type" : "int"
})", DBFilterSettings::Defaults::RandomSeed};

const CLI_v2::Option FilterListPath{
R"({
    "names" : ["filter-list"],
    "description" : "A text file containing headers of blacklisted sequences, one per line.",
    "type" : "string",
    "default" : ""
})"};

const CLI_v2::Option FilterType{
R"({
    "names" : ["filter-type"],
    "choices" : ["none", "whitelist", "blacklist"],
    "type" : "string",
    "default" : "none",
    "description" : "Select from: whitelist, blacklist, none."
})", std::string("none")};

const CLI_v2::Option Consolidate{
R"({
    "names" : ["consolidate"],
    "description" : "Create the new data files in addition to filtering the DB.",
    "type" : "bool"
})"};

const CLI_v2::Option CompressionLevel {
R"({
    "names" : ["c", "compression"],
    "description" : "Compression level for output sequences.",
    "type" : "int"
})", DBFilterSettings::Defaults::CompressionLevel};

const CLI_v2::Option BufferSize{
R"({
    "names" : ["buffer-size"],
    "description" : "Sequence buffer size in megabytes. Has to be >= 0.0.",
    "type" : "float"
})", DBFilterSettings::Defaults::BufferSize};

const CLI_v2::Option SplitBlocks{
R"({
    "names" : ["split-blocks"],
    "description" : "Write seeds for each block into a separate file."
})", DBFilterSettings::Defaults::SplitBlocks};

// clang-format on

}  // namespace DBFilterOptionNames

SamplingType ParseSamplingType(const std::string& val)
{
    if (val == "none") {
        return SamplingType::None;
    } else if (val == "linear") {
        return SamplingType::Linear;
    } else if (val == "random") {
        return SamplingType::Random;
    }
    return SamplingType::Unknown;
}

FilterListType ParseFilterListType(const std::string& val)
{
    if (val == "whitelist") {
        return FilterListType::Whitelist;
    } else if (val == "blacklist") {
        return FilterListType::Blacklist;
    } else if (val == "none") {
        return FilterListType::None;
    }
    return FilterListType::Unknown;
}

DBFilterSettings::DBFilterSettings() = default;

DBFilterSettings::DBFilterSettings(const PacBio::CLI_v2::Results& options)
    : InputPrefix{options[DBFilterOptionNames::InputPrefix]}
    , OutputPrefix{options[DBFilterOptionNames::OutputPrefix]}
    , SampleBases{options[DBFilterOptionNames::SampleBases]}
    , BlockSize{options[DBFilterOptionNames::BlockSize]}
    , RandomSeed{options[DBFilterOptionNames::RandomSeed]}
    , FilterListPath(options[DBFilterOptionNames::FilterListPath])
    , Consolidate(options[DBFilterOptionNames::Consolidate])
    , CompressionLevel{options[DBFilterOptionNames::CompressionLevel]}
    , BufferSize{options[DBFilterOptionNames::BufferSize]}
    , SplitBlocks{options[DBFilterOptionNames::SplitBlocks]}
{
    Sampling = ParseSamplingType(options[DBFilterOptionNames::Sampling]);
    if (Sampling == SamplingType::Unknown) {
        throw std::runtime_error("Unknown sampling type: '" +
                                 std::string(options[DBFilterOptionNames::Sampling]) + "'.");
    }

    // Parse the filter type into an enum class.
    FilterType = ParseFilterListType(options[DBFilterOptionNames::FilterType]);
    if (FilterType == FilterListType::Unknown) {
        throw std::runtime_error("Unknown filter type: '" +
                                 std::string(options[DBFilterOptionNames::FilterType]) + "'.");
    }

    // Convert block and buffer sizes from MB to bytes.
    BlockSize *= (1000 * 1000);
    BufferSize *= (1024 * 1024);

    // Negative block size indicates that everything should be in one block.
    if (BlockSize < 0.0f) {
        BlockSize = static_cast<float>(std::numeric_limits<int64_t>::max()) / (1024.0f * 1024.0f);
    }

    // Buffer size can be zero, but not negative.
    if (BufferSize < 0.0f) {
        throw std::runtime_error("Buffer size cannot be a negative value.");
    }
}

PacBio::CLI_v2::Interface DBFilterSettings::CreateCLI()
{
    PacBio::CLI_v2::Interface i{"pancake dbfilter", "Filters the input SeqDB.",
                                PacBio::Pancake::PancakeFormattedVersion()};

    // clang-format off
    i.AddOptionGroup("Algorithm Options", {
        DBFilterOptionNames::Sampling,
        DBFilterOptionNames::SampleBases,
        DBFilterOptionNames::BlockSize,
        DBFilterOptionNames::RandomSeed,
        DBFilterOptionNames::FilterListPath,
        DBFilterOptionNames::FilterType,
        DBFilterOptionNames::Consolidate,
    });
    i.AddOptionGroup("Consolidation Options", {
        DBFilterOptionNames::CompressionLevel,
        DBFilterOptionNames::BufferSize,
        DBFilterOptionNames::SplitBlocks,
    });
    i.AddPositionalArguments({
        DBFilterOptionNames::InputPrefix,
        DBFilterOptionNames::OutputPrefix,
    });

    // clang-format on
    return i;
}
}  // namespace Pancake
}  // namespace PacBio
