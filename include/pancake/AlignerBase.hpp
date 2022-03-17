// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNER_BASE_HPP
#define PANCAKE_ALIGNER_BASE_HPP

#include <pancake/AlignmentParameters.hpp>
#include <pancake/AlignmentResult.hpp>

#include <cstdint>
#include <memory>
#include <vector>

namespace PacBio {
namespace Pancake {

class AlignerBase
{
public:
    virtual ~AlignerBase() = default;
    virtual AlignmentResult Global(const char* qseq, int64_t qlen, const char* tseq,
                                   int64_t tlen) = 0;
    virtual AlignmentResult Extend(const char* qseq, int64_t qlen, const char* tseq,
                                   int64_t tlen) = 0;
};

typedef std::shared_ptr<AlignerBase> AlignerBasePtr;

AlignmentResult EdgeCaseAlignmentResult(int32_t qlen, int32_t tlen, int32_t matchScore,
                                        int32_t mismatchPenalty, int32_t gapOpen,
                                        int32_t gapExtend);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ALIGNER_BASE_HPP
