// Author: Ivan Sovic

#include <pancake/OverlapHifiSettings.hpp>
#include <pancake/Version.hpp>
#include "dbfilter/DBFilterSettings.hpp"
#include "dbfilter/DBFilterWorkflow.hpp"
#include "mapclr/MapCLRSettings.hpp"
#include "mapclr/MapCLRWorkflow.hpp"
#include "overlaphifi/OverlapHifiWorkflow.hpp"
#include "seeddb/SeedDBSettings.hpp"
#include "seeddb/SeedDBWorkflow.hpp"
#include "seqdb/SeqDBSettings.hpp"
#include "seqdb/SeqDBWorkflow.hpp"
#include "seqdbdump/SeqDBDumpSettings.hpp"
#include "seqdbdump/SeqDBDumpWorkflow.hpp"
#include "seqdbinfo/SeqDBInfoSettings.hpp"
#include "seqdbinfo/SeqDBInfoWorkflow.hpp"
#include "seqfetch/SeqFetchSettings.hpp"
#include "seqfetch/SeqFetchWorkflow.hpp"

#include <pbcopper/cli2/CLI.h>

#include <cstdlib>
#include <iostream>
#include <stdexcept>

PacBio::CLI_v2::MultiToolInterface CreateMultiInterface()
{
    PacBio::CLI_v2::MultiToolInterface mi{"pancake", "PacBio HiFi overlapper.",
                                          PacBio::Pancake::PancakeFormattedVersion()};

    // clang-format off
    mi.AddTools(
    {
        {"map",
            PacBio::Pancake::MapCLRSettings::CreateCLI(),
           &PacBio::Pancake::MapCLRWorkflow::Runner},
        {"seqdb",
            PacBio::Pancake::SeqDBSettings::CreateCLI(),
           &PacBio::Pancake::SeqDBWorkflow::Runner},
        {"seeddb",
            PacBio::Pancake::SeedDBSettings::CreateCLI(),
           &PacBio::Pancake::SeedDBWorkflow::Runner},
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