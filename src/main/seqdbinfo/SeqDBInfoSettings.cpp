// Author: Ivan Sovic

#include "SeqDBInfoSettings.h"

#include <limits>

#include <pacbio/Version.h>

namespace PacBio {
namespace Pancake {
namespace SeqDBInfoOptionNames {

// clang-format off
const CLI_v2::PositionalArgument InputSeqDB {
R"({
    "name" : "seqdb",
    "description" : "Input SeqDB."
})"};

const CLI_v2::Option HumanReadableOutput {
R"({
    "names" : ["human"],
    "description" : "Write human-readable output instead of Json.",
    "type" : "bool"
})", SeqDBInfoSettings::Defaults::HumanReadableOutput};

const CLI_v2::Option GenomicUnit {
R"({
    "names" : ["unit"],
    "choices" : ["bp", "kbp", "Mbp", "Gbp"],
    "description" : "Unit for genomic length for the output stats. Valid choices: bp, kbp, Mbp, Gbp.",
    "type" : "string"
})", std::string("bp")};

// clang-format on

}  // namespace SeqDBInfoOptionNames

SeqDBInfoSettings::SeqDBInfoSettings() = default;

SeqDBInfoSettings::SeqDBInfoSettings(const PacBio::CLI_v2::Results& options)
    : InputSeqDB{options[SeqDBInfoOptionNames::InputSeqDB]}
    , HumanReadableOutput(options[SeqDBInfoOptionNames::HumanReadableOutput])
    , Unit(GenomicUnitFromString(options[SeqDBInfoOptionNames::GenomicUnit]))
{
    if (Unit == GenomicUnit::Unknown) {
        const std::string unitStr = options[SeqDBInfoOptionNames::GenomicUnit];
        throw std::runtime_error("Unknown genomic unit: '" + unitStr + "'.");
    }
}

PacBio::CLI_v2::Interface SeqDBInfoSettings::CreateCLI()
{
    PacBio::CLI_v2::Interface i{"pancake seqdb-info",
                                "Computes the sequence statistics in the SeqDB.",
                                PacBio::Pancake::PancakeFormattedVersion()};

    i.DisableNumThreadsOption();

    // clang-format off
    i.AddOptionGroup("Options", {
        SeqDBInfoOptionNames::HumanReadableOutput,
        SeqDBInfoOptionNames::GenomicUnit,
    });
    i.AddPositionalArguments({
        SeqDBInfoOptionNames::InputSeqDB,
    });

    // clang-format on
    return i;
}
}  // namespace Pancake
}  // namespace PacBio
