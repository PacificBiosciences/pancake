// Author: Ivan Sovic

#include "SeqDBSettings.hpp"

#include <pancake/Version.hpp>

#include <limits>

namespace PacBio {
namespace Pancake {
namespace SeqDBOptionNames {

// clang-format off

const CLI_v2::PositionalArgument OutputPrefix {
R"({
    "name" : "prefix",
    "description" : "The prefix of the DB files."
})"};

const CLI_v2::PositionalArgument Input {
R"({
    "name" : "<input.fasta> [...]",
    "description" : "One or more input sequence files, in FASTA, FASTQ or FOFN formats."
})"};

const CLI_v2::Option CompressionLevel {
R"({
    "names" : ["c", "compression"],
    "description" : "Compression level for output sequences.",
    "type" : "int"
})", SeqDBSettings::Defaults::CompressionLevel};

const CLI_v2::Option BufferSize{
R"({
    "names" : ["buffer-size"],
    "description" : "Sequence buffer size in megabytes. Has to be >= 0.0.",
    "type" : "float"
})", SeqDBSettings::Defaults::BufferSize};

const CLI_v2::Option BlockSize{
R"({
    "names" : ["block-size"],
    "description" : "Block size in megabases. Value 0 means one sequnece per block, value < 0 all sequences in one block.",
    "type" : "float"
})", SeqDBSettings::Defaults::BlockSize};

const CLI_v2::Option SplitBlocks{
R"({
    "names" : ["split-blocks"],
    "description" : "Write seeds for each block into a separate file.",
    "type" : "bool"
})", SeqDBSettings::Defaults::SplitBlocks};

const CLI_v2::Option KeepOriginalCase{
R"({
    "names" : ["keep-case"],
    "description" : "Prevent conversion from lowercase to uppercase bases. Only valid for uncompressed DBs.",
    "type" : "bool"
})", SeqDBSettings::Defaults::KeepOriginalCase};

// clang-format on

}  // namespace SeqDBOptionNames

SeqDBSettings::SeqDBSettings() = default;

SeqDBSettings::SeqDBSettings(const PacBio::CLI_v2::Results& options)
    : OutputPrefix{options[SeqDBOptionNames::OutputPrefix]}
    , InputFiles{options[SeqDBOptionNames::Input]}
    , NumThreads{options.NumThreads()}
    , CompressionLevel{options[SeqDBOptionNames::CompressionLevel]}
    , BufferSize{options[SeqDBOptionNames::BufferSize]}
    , BlockSize{options[SeqDBOptionNames::BlockSize]}
    , SplitBlocks{options[SeqDBOptionNames::SplitBlocks]}
    , KeepOriginalCase{options[SeqDBOptionNames::KeepOriginalCase]}
{
    if (KeepOriginalCase && CompressionLevel > 0) {
        throw std::runtime_error("Original case can only be used with uncompressed DBs.");
    }

    // Allow multiple positional input arguments.
    const auto& files = options.PositionalArguments();
    if (files.size() < 2) {
        throw std::runtime_error{"Not enough input files specified, at least one required."};
    }

    OutputPrefix = files[0];
    InputFiles.clear();
    for (size_t i = 1; i < files.size(); ++i)
        InputFiles.push_back(files[i]);

    // Convert block and buffer sizes from MB to bytes.
    BlockSize *= (1000 * 1000);
    BufferSize *= (1024 * 1024);

    // Negative block size indicates that everything should be in one block.
    if (BlockSize < 0.0f)
        BlockSize = static_cast<float>(std::numeric_limits<int64_t>::max()) / (1024.0f * 1024.0f);

    // Buffer size can be zero, but not negative.
    if (BufferSize < 0.0f) throw std::runtime_error("Buffer size cannot be a negative value.");
}

PacBio::CLI_v2::Interface SeqDBSettings::CreateCLI()
{
    PacBio::CLI_v2::Interface i{"pancake seqdb", "Convert FASTA/FASTQ sequences to SeqDB.",
                                PacBio::Pancake::PancakeFormattedVersion()};

    // clang-format off
    i.AddOptionGroup("Algorithm Options", {
        SeqDBOptionNames::CompressionLevel,
        SeqDBOptionNames::BufferSize,
        SeqDBOptionNames::BlockSize,
        SeqDBOptionNames::SplitBlocks,
        SeqDBOptionNames::KeepOriginalCase,
    });
    i.AddPositionalArguments({
        SeqDBOptionNames::OutputPrefix,
        SeqDBOptionNames::Input,
    });

    // clang-format on
    return i;
}
}  // namespace Pancake
}  // namespace PacBio
