// Author: Ivan Sovic

#ifndef PANCAKE_SEQDB_WORKFLOW_H
#define PANCAKE_SEQDB_WORKFLOW_H

#include <pbcopper/cli2/CLI.h>

namespace PacBio {
namespace Pancake {

struct SeqDBWorkflow
{
    static int Runner(const PacBio::CLI_v2::Results& options);
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEQDB_WORKFLOW_H