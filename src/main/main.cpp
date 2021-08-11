// Author: Ivan Sovic

#include <cstdlib>
#include <iostream>
#include <stdexcept>

#include <pbcopper/cli2/CLI.h>

#include <pacbio/Version.h>

#include <pacbio/overlaphifi/OverlapHifiSettings.h>
#include "dbfilter/DBFilterSettings.h"
#include "dbfilter/DBFilterWorkflow.h"
#include "overlaphifi/OverlapHifiWorkflow.h"
#include "seeddb/SeedDBSettings.h"
#include "seeddb/SeedDBWorkflow.h"
#include "seqdb/SeqDBSettings.h"
#include "seqdb/SeqDBWorkflow.h"
#include "seqdbdump/SeqDBDumpSettings.h"
#include "seqdbdump/SeqDBDumpWorkflow.h"
#include "seqdbinfo/SeqDBInfoSettings.h"
#include "seqdbinfo/SeqDBInfoWorkflow.h"
#include "seqfetch/SeqFetchSettings.h"
#include "seqfetch/SeqFetchWorkflow.h"

PacBio::CLI_v2::MultiToolInterface CreateMultiInterface()
{
    PacBio::CLI_v2::MultiToolInterface mi{"pancake", "PacBio HiFi overlapper.",
                                          PacBio::Pancake::PancakeFormattedVersion()};

    // clang-format off
    mi.AddTools(
    {
        {"seqdb",
            PacBio::Pancake::SeqDBSettings::CreateCLI(),
           &PacBio::Pancake::SeqDBWorkflow::Runner},
        {"seeddb",
            PacBio::Pancake::SeedDB::SeedDBSettings::CreateCLI(),
           &PacBio::Pancake::SeedDB::SeedDBWorkflow::Runner},
        {"ovl-hifi",
            PacBio::Pancake::OverlapHifiSettings::CreateCLI(),
           &PacBio::Pancake::OverlapHifiWorkflow::Runner},
        {"dbfilter",
            PacBio::Pancake::DBFilterSettings::CreateCLI(),
           &PacBio::Pancake::DBFilterWorkflow::Runner},
        {"seqfetch",
            PacBio::Pancake::SeqFetchSettings::CreateCLI(),
           &PacBio::Pancake::SeqFetchWorkflow::Runner},
        {"seqdb-dump",
            PacBio::Pancake::SeqDBDumpSettings::CreateCLI(),
           &PacBio::Pancake::SeqDBDumpWorkflow::Runner},
        {"seqdb-info",
            PacBio::Pancake::SeqDBInfoSettings::CreateCLI(),
           &PacBio::Pancake::SeqDBInfoWorkflow::Runner},
    });

    // clang-format on
    return mi;
}

int main(int argc, char* argv[])
{
    try {
        return PacBio::CLI_v2::Run(argc, argv, CreateMultiInterface());
    } catch (const std::exception& e) {
        std::cerr << "Pancake ERROR: " << e.what() << '\n';
        return EXIT_FAILURE;
    }
}