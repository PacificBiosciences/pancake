// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNER_SES2_HPP
#define PANCAKE_ALIGNER_SES2_HPP

#include <pancake/AlignerBase.hpp>
#include <pancake/AlignmentParameters.hpp>
#include <pancake/SesResults.hpp>

#include <pbbam/Cigar.h>

#include <array>
#include <cstdint>
#include <memory>
#include <vector>

namespace PacBio {
namespace Pancake {

class AlignerSES2;
std::shared_ptr<AlignerBase> CreateAlignerSES2(const AlignmentParameters& opt);

class AlignerSES2 : public AlignerBase
{
public:
    AlignerSES2(const AlignmentParameters& opt);
    ~AlignerSES2() override;

    AlignmentResult Global(const char* qseq, int64_t qlen, const char* tseq, int64_t tlen) override;
    AlignmentResult Extend(const char* qseq, int64_t qlen, const char* tseq, int64_t tlen) override;

private:
    AlignmentParameters opt_;
    std::shared_ptr<PacBio::Pancake::Alignment::SESScratchSpace> sesScratch_;
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ALIGNER_SES2_HPP
