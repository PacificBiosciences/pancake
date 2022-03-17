// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNER_SES1_HPP
#define PANCAKE_ALIGNER_SES1_HPP

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

class AlignerSES1;
std::shared_ptr<AlignerBase> CreateAlignerSES1(const AlignmentParameters& opt);

class AlignerSES1 : public AlignerBase
{
public:
    AlignerSES1(const AlignmentParameters& opt);
    ~AlignerSES1() override;

    AlignmentResult Global(const char* qseq, int64_t qlen, const char* tseq, int64_t tlen) override;
    AlignmentResult Extend(const char* qseq, int64_t qlen, const char* tseq, int64_t tlen) override;

private:
    AlignmentParameters opt_;
    std::shared_ptr<PacBio::Pancake::Alignment::SESScratchSpace> sesScratch_;
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ALIGNER_SES1_HPP
