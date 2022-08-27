// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNER_SES2_HPP
#define PANCAKE_ALIGNER_SES2_HPP

#include <pancake/AlignerBase.hpp>
#include <pancake/AlignmentParameters.hpp>
#include <pancake/SesResults.hpp>

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

    AlignmentResult Global(std::string_view qseq, std::string_view tseq) override;
    AlignmentResult Extend(std::string_view qseq, std::string_view tseq) override;

private:
    AlignmentParameters opt_;
    std::shared_ptr<PacBio::Pancake::SESScratchSpace> sesScratch_;
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ALIGNER_SES2_HPP
