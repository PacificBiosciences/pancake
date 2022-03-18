// Author: Ivan Sovic

#include "SeqFetchSettings.hpp"

#include <pancake/Version.hpp>

#include <limits>

namespace PacBio {
namespace Pancake {
namespace SeqFetchOptionNames {

// clang-format off

const CLI_v2::PositionalArgument OutputFile {
R"({
    "name" : "out_fn",
    "description" : "Output file for the fetched sequences."
})"};

const CLI_v2::PositionalArgument InputFetchListFile {
R"({
    "name" : "fetch_list",
    "description" : "List of sequences to fetch, one per line."
})"};

const CLI_v2::PositionalArgument InputFiles {
R"({
    "name" : "<input.fasta/fastq/bam/fofn> [...]",
    "description" : "One or more input sequence files, in FASTA, FASTQ, BAM, SeqDB or FOFN formats."
})"};

const CLI_v2::Option OutputFormat{
R"({
    "names" : ["out-fmt"],
    "choices" : ["fasta", "fastq"],
    "type" : "string",
    "default" : "fasta",
    "description" : "Output format. If an input file is FASTA and out format is FASTQ, dummy QVs will be added."
})", std::string("fasta")};

const CLI_v2::Option DummyQV{
R"({
    "names" : ["dummy-qv"],
    "type" : "string",
    "default" : "!",
    "description" : "Dummy QV to be added to sequences when input format is FASTA and output FASTQ."
})", std::string("!")};

const CLI_v2::Option AliasSeqDBFile{
R"({
    "names" : ["alias"],
    "type" : "string",
    "default" : "",
    "description" : "SeqDB file path. If provided, the SeqDB will be used to look-up the provided sequences by their IDs."
})", std::string("")};

const CLI_v2::Option FailOnMissingQueries {
R"({
    "names" : ["fail"],
    "description" : "Exit non-zero if not all seqs are found.",
    "type" : "bool"
})", SeqFetchSettings::Defaults::FailOnMissingQueries};

const CLI_v2::Option WriteIds {
R"({
    "names" : ["write-ids"],
    "description" : "The output sequence names will be replaced by their IDs in the SeqDB, if the SeqDB was provided as input.",
    "type" : "bool"
})", SeqFetchSettings::Defaults::WriteIds};

const CLI_v2::Option UseHPC{
R"({
    "names" : ["use-hpc"],
    "description" : "Fetch homopolymer compressed sequences."
})", SeqFetchSettings::Defaults::UseHPC};

const CLI_v2::Option UseRLE{
R"({
    "names" : ["use-rle"],
    "description" : "Write a run-length-encoded file alongside to the output file. The RLE file contains conversion coordinates from the HPC space to the original space instead of the run-length-encoding. This option does not write the HPC sequence, for that please specify '--user-hpc'."
})", SeqFetchSettings::Defaults::UseRLE};

// clang-format on

}  // namespace SeqFetchOptionNames

SeqFetchOutFormat ParseSeqFetchOutFormat(const std::string& val)
{
    if (val == "fasta") {
        return SeqFetchOutFormat::Fasta;
    } else if (val == "fastq") {
        return SeqFetchOutFormat::Fastq;
    } else {
        throw std::runtime_error("Unknown output format '" + val + "'.");
    }
    return SeqFetchOutFormat::Fasta;
}

SeqFetchSettings::SeqFetchSettings() = default;

SeqFetchSettings::SeqFetchSettings(const PacBio::CLI_v2::Results& options)
    : OutputFile{options[SeqFetchOptionNames::OutputFile]}
    , InputFetchListFile{options[SeqFetchOptionNames::InputFetchListFile]}
    , InputFiles{options[SeqFetchOptionNames::InputFiles]}
    , AliasSeqDBFile(options[SeqFetchOptionNames::AliasSeqDBFile])
    , FailOnMissingQueries(options[SeqFetchOptionNames::FailOnMissingQueries])
    , WriteIds(options[SeqFetchOptionNames::WriteIds])
    , UseHPC(options[SeqFetchOptionNames::UseHPC])
    , UseRLE(options[SeqFetchOptionNames::UseRLE])
{
    // Allow multiple positional input arguments.
    const auto& files = options.PositionalArguments();
    if (files.size() < 3) throw std::runtime_error{"Not enough positional arguments specified."};
    OutputFile = files[0];
    InputFiles.clear();
    for (size_t i = 2; i < files.size(); ++i)
        InputFiles.push_back(files[i]);

    OutputFormat = ParseSeqFetchOutFormat(options[SeqFetchOptionNames::OutputFormat]);

    std::string tempDummyQV = options[SeqFetchOptionNames::DummyQV];
    if (tempDummyQV.size() != 1) {
        throw std::runtime_error("The dummyQV needs to be exactly 1 character in size.");
    }
    DummyQV = tempDummyQV[0];

    if (UseHPC && OutputFormat == SeqFetchOutFormat::Fastq) {
        throw std::runtime_error(
            "Fastq output format is not supported with the homopolymer compression option.");
    }
}

PacBio::CLI_v2::Interface SeqFetchSettings::CreateCLI()
{
    PacBio::CLI_v2::Interface i{"pancake seqfetch",
                                "Fetches a set of sequences in random access from a list of "
                                "specified indexed sequence files.",
                                PacBio::Pancake::PancakeFormattedVersion()};

    i.DisableNumThreadsOption();

    // clang-format off
    i.AddOptionGroup("Algorithm Options", {
        SeqFetchOptionNames::OutputFormat,
        SeqFetchOptionNames::DummyQV,
        SeqFetchOptionNames::AliasSeqDBFile,
        SeqFetchOptionNames::FailOnMissingQueries,
        SeqFetchOptionNames::WriteIds,
        SeqFetchOptionNames::UseHPC,
        SeqFetchOptionNames::UseRLE,
    });
    i.AddPositionalArguments({
        SeqFetchOptionNames::OutputFile,
        SeqFetchOptionNames::InputFetchListFile,
        SeqFetchOptionNames::InputFiles,
    });

    // clang-format on
    return i;
}
}  // namespace Pancake
}  // namespace PacBio
