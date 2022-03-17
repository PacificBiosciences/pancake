// Author: Ivan Sovic

#ifndef PANCAKE_DBFILTER_WORKFLOW_H
#define PANCAKE_DBFILTER_WORKFLOW_H

#include <pbcopper/cli2/CLI.h>

namespace PacBio {
namespace Pancake {

struct DBFilterWorkflow
{
    static int Runner(const PacBio::CLI_v2::Results& options);
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_DBFILTER_WORKFLOW_H