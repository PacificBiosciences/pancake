// Author: Ivan Sovic

#ifndef PANCAKE_SEEDDB_WORKFLOW_H
#define PANCAKE_SEEDDB_WORKFLOW_H

#include <pbcopper/cli2/CLI.h>

namespace PacBio {
namespace Pancake {
namespace SeedDB {

struct SeedDBWorkflow
{
    static int Runner(const PacBio::CLI_v2::Results& options);
};

}  // namespace SeedDB
}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEEDDB_WORKFLOW_H