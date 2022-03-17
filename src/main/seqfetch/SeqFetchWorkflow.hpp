// Author: Ivan Sovic

#ifndef PANCAKE_SEQFETCH_WORKFLOW_H
#define PANCAKE_SEQFETCH_WORKFLOW_H

#include <pbcopper/cli2/CLI.h>

namespace PacBio {
namespace Pancake {

struct SeqFetchWorkflow
{
    static int Runner(const PacBio::CLI_v2::Results& options);
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEQFETCH_WORKFLOW_H