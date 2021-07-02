// Author: Ivan Sovic

#include "SeqDBSettings.h"
#include <pacbio/Version.h>
#include <limits>

namespace PacBio {
namespace Pancake {
namespace OptionNames {

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
    "description" : "Write seeds for each block into a separate file."
})", SeqDBSettings::Defaults::SplitBlocks};

// clang-format on

}  // namespace OptionNames

SeqDBSettings::SeqDBSettings() = default;

SeqDBSettings::SeqDBSettings(const PacBio::CLI_v2::Results& options)
    : OutputPrefix{options[OptionNames::OutputPrefix]}
    , InputFiles{options[OptionNames::Input]}
    , NumThreads{options.NumThreads()}
    , CompressionLevel{options[OptionNames::CompressionLevel]}
    , BufferSize{options[OptionNames::BufferSize]}
    , BlockSize{options[OptionNames::BlockSize]}
    , SplitBlocks{options[OptionNames::SplitBlocks]}
{
    // Allow multiple positional input arguments.
    const auto& files = options.PositionalArguments();
    if (files.size() < 2)
        throw std::runtime_error{"Not enough input files specified, at least one required."};
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
    PacBio::CLI_v2::Interface i{"pancake", "Convert FASTA/FASTQ sequences to SeqDB.",
                                PacBio::Pancake::PancakeFormattedVersion()};

    // clang-format off
    i.AddOptionGroup("Algorithm Options", {
        OptionNames::CompressionLevel,
        OptionNames::BufferSize,
        OptionNames::BlockSize,
        OptionNames::SplitBlocks,
    });
    i.AddPositionalArguments({
        OptionNames::OutputPrefix,
        OptionNames::Input,
    });

    // clang-format on
    return i;
}
}  // namespace Pancake
}  // namespace PacBio