// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNER_EDLIB_HPP
#define PANCAKE_ALIGNER_EDLIB_HPP

#include <pancake/AlignerBase.hpp>
#include <pancake/AlignmentParameters.hpp>

#include <pbbam/Cigar.h>

#include <array>
#include <cstdint>
#include <memory>
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

    AlignmentResult Global(std::string_view qseq, std::string_view tseq) override;
    AlignmentResult Extend(std::string_view qseq, std::string_view tseq) override;

private:
    AlignmentParameters opt_;
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ALIGNER_EDLIB_HPP
