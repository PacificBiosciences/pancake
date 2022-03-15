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

    AlignmentResult Global(std::string_view qseq, std::string_view tseq) override;
    AlignmentResult Extend(std::string_view qseq, std::string_view tseq) override;

private:
    AlignmentParameters opt_;
    std::shared_ptr<PacBio::Pancake::SESScratchSpace> sesScratch_;
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ALIGNER_SES1_HPP
