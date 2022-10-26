// Author: Ivan Sovic

#include "SeqDBDumpSettings.hpp"

#include <pancake/LibraryInfo.hpp>

#include <limits>

namespace PacBio {
namespace Pancake {
namespace SeqDBDumpOptionNames {

// clang-format off
const CLI_v2::PositionalArgument InputSeqDB {
R"({
    "name" : "seqdb",
    "description" : "Input SeqDB."
})"};

const CLI_v2::PositionalArgument OutputFile {
R"({
    "name" : "out_fn",
    "description" : "Output file to write the sequences to, or '-' to write to stdout."
})"};

const CLI_v2::Option BlockId {
R"({
    "names" : ["block-id"],
    "description" : "Writes only the specified block to the output, or the entire DB if '-1' is given to this parameter.",
    "type" : "int"
})", SeqDBDumpSettings::Defaults::BlockId};

const CLI_v2::Option WriteIds {
R"({
    "names" : ["write-ids"],
    "description" : "The output sequence names will be replaced by their IDs in the SeqDB, if the SeqDB was provided as input.",
    "type" : "bool"
})", SeqDBDumpSettings::Defaults::WriteIds};

const CLI_v2::Option UseHPC{
R"({
    "names" : ["use-hpc"],
    "description" : "Fetch homopolymer compressed sequences."
})", SeqDBDumpSettings::Defaults::UseHPC};

// clang-format on

}  // namespace SeqDBDumpOptionNames

SeqDBDumpSettings::SeqDBDumpSettings() = default;

SeqDBDumpSettings::SeqDBDumpSettings(const PacBio::CLI_v2::Results& options)
    : InputSeqDB{options[SeqDBDumpOptionNames::InputSeqDB]}
    , OutputFile{options[SeqDBDumpOptionNames::OutputFile]}
    , BlockId(options[SeqDBDumpOptionNames::BlockId])
    , WriteIds(options[SeqDBDumpOptionNames::WriteIds])
    , UseHPC(options[SeqDBDumpOptionNames::UseHPC])
{}

PacBio::CLI_v2::Interface SeqDBDumpSettings::CreateCLI()
{
    PacBio::CLI_v2::Interface i{"pancake seqdb-dump",
                                "Writes the entire SeqDB as FASTA to file or stdout.",
                                PacBio::Pancake::PancakeFormattedVersion()};

    i.DisableNumThreadsOption();

    // clang-format off
    i.AddOptionGroup("Algorithm Options", {
        SeqDBDumpOptionNames::BlockId,
        SeqDBDumpOptionNames::WriteIds,
        SeqDBDumpOptionNames::UseHPC,
    });
    i.AddPositionalArguments({
        SeqDBDumpOptionNames::InputSeqDB,
        SeqDBDumpOptionNames::OutputFile,
    });

    // clang-format on
    return i;
}
}  // namespace Pancake
}  // namespace PacBio
