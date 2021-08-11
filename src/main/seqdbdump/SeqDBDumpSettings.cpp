// Author: Ivan Sovic

#include "SeqDBDumpSettings.h"

#include <limits>

#include <pacbio/Version.h>

namespace PacBio {
namespace Pancake {
namespace OptionNames {

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

}  // namespace OptionNames

SeqDBDumpSettings::SeqDBDumpSettings() = default;

SeqDBDumpSettings::SeqDBDumpSettings(const PacBio::CLI_v2::Results& options)
    : InputSeqDB{options[OptionNames::InputSeqDB]}
    , OutputFile{options[OptionNames::OutputFile]}
    , BlockId(options[OptionNames::BlockId])
    , WriteIds(options[OptionNames::WriteIds])
    , UseHPC(options[OptionNames::UseHPC])
{
}

PacBio::CLI_v2::Interface SeqDBDumpSettings::CreateCLI()
{
    PacBio::CLI_v2::Interface i{"pancake seqfetch",
                                "Fetches a set of sequences in random access from a list of "
                                "specified indexed sequence files.",
                                PacBio::Pancake::PancakeFormattedVersion()};

    i.DisableNumThreadsOption();

    // clang-format off
    i.AddOptionGroup("Algorithm Options", {
        OptionNames::BlockId,
        OptionNames::WriteIds,
        OptionNames::UseHPC,
    });
    i.AddPositionalArguments({
        OptionNames::InputSeqDB,
        OptionNames::OutputFile,
    });

    // clang-format on
    return i;
}
}  // namespace Pancake
}  // namespace PacBio
