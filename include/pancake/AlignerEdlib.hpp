// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNER_EDLIB_H
#define PANCAKE_ALIGNER_EDLIB_H

#include <pbbam/Cigar.h>
#include <array>
#include <cstdint>
#include <memory>
#include <pancake/AlignerBase.hpp>
#include <vector>

namespace PacBio {
namespace Pancake {

class AlignerEdlib;
std::shared_ptr<AlignerBase> CreateAlignerEdlib(const AlignmentParameters& opt);

class AlignerEdlib : public AlignerBase
{
public:
    AlignerEdlib(const AlignmentParameters& opt);
    ~AlignerEdlib() override;

    AlignmentResult Global(const char* qseq, int64_t qlen, const char* tseq, int64_t tlen) override;
    AlignmentResult Extend(const char* qseq, int64_t qlen, const char* tseq, int64_t tlen) override;

private:
    AlignmentParameters opt_;
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ALIGNER_EDLIB_H
