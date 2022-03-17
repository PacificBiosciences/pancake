// Author: Ivan Sovic

#ifndef PANCAKE_OVERLAP_HIFI_WORKFLOW_H
#define PANCAKE_OVERLAP_HIFI_WORKFLOW_H

#include <pbcopper/cli2/CLI.h>

namespace PacBio {
namespace Pancake {

struct OverlapHifiWorkflow
{
    static int Runner(const PacBio::CLI_v2::Results& options);
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_OVERLAP_HIFI_WORKFLOW_H