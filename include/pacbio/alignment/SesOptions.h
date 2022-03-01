// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNMENT_SES_OPTIONS_H
#define PANCAKE_ALIGNMENT_SES_OPTIONS_H

namespace PacBio {
namespace Pancake {
namespace Alignment {

enum class SESAlignMode
{
    Global,      // Sequences are aligned end to end.
    Semiglobal,  // No penalty at the end of the query or target.
};

enum class SESTracebackMode
{
    Disabled,
    Enabled,
};

enum class SESTrimmingMode
{
    Disabled,
    Enabled,
};
}  // namespace Alignment
}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ALIGNMENT_SES_OPTIONS_H
