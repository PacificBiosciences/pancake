// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNER_BASE_HPP
#define PANCAKE_ALIGNER_BASE_HPP

#include <pancake/AlignmentParameters.hpp>
#include <pancake/AlignmentResult.hpp>

#include <cstdint>
#include <memory>
#include <string_view>
#include <vector>

namespace PacBio {
namespace Pancake {

class AlignerBase
{
public:
    virtual ~AlignerBase() = default;
    virtual AlignmentResult Global(std::string_view qseq, std::string_view tseq) = 0;
    virtual AlignmentResult Extend(std::string_view qseq, std::string_view tseq) = 0;
};

typedef std::shared_ptr<AlignerBase> AlignerBasePtr;

AlignmentResult EdgeCaseAlignmentResult(int32_t qlen, int32_t tlen, int32_t matchScore,
                                        int32_t mismatchPenalty, int32_t gapOpen,
                                        int32_t gapExtend);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ALIGNER_BASE_HPP
