// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNMENT_PARAMETERS_H
#define PANCAKE_ALIGNMENT_PARAMETERS_H

#include <cstdint>
#include <memory>
#include <sstream>
#include <vector>

namespace PacBio {
namespace Pancake {

// clang-format off
class AlignmentParameters
{
public:
    int32_t zdrop = 100;                    // Zdrop for alignment extension.
    int32_t zdrop2 = 500;                   // In case the small Zdrop fails, this one is applied. ('zdrop_inv' in Minimap2)
    int32_t alignBandwidth = 500;           // Bandwidth used for alignment of non-long-gap regions.
    int32_t endBonus = 50;                  // Score used at the beginning of alignment extension.
    int32_t matchScore = 2;                 // 'a' in Minimap2.
    int32_t mismatchPenalty = 4;            // 'b' in Minimap2.
    int32_t gapOpen1 = 4;                   // 'q' in Minimap2.
    int32_t gapExtend1 = 2;                 // 'e' in Minimap2.
    int32_t gapOpen2 = 24;                  // 'q2' in Minimap2.
    int32_t gapExtend2 = 1;                 // 'e2' in Minimap2.
    bool dynamicBandwidth = false;
};
// clang-format on

inline std::ostream& operator<<(std::ostream& out, const AlignmentParameters& a)
{
    out << "zdrop = " << a.zdrop << "\n"
        << "zdrop2 = " << a.zdrop2 << "\n"
        << "alignBandwidth = " << a.alignBandwidth << "\n"
        << "endBonus = " << a.endBonus << "\n"
        << "matchScore = " << a.matchScore << "\n"
        << "mismatchPenalty = " << a.mismatchPenalty << "\n"
        << "gapOpen1 = " << a.gapOpen1 << "\n"
        << "gapExtend1 = " << a.gapExtend1 << "\n"
        << "gapOpen2 = " << a.gapOpen2 << "\n"
        << "gapExtend2 = " << a.gapExtend2 << "\n"
        << "dynamicBandwidth = " << std::boolalpha << a.dynamicBandwidth << "\n";
    return out;
}

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ALIGNMENT_PARAMETERS_H
