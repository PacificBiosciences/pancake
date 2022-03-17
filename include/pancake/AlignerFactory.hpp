// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNER_FACTORY_HPP
#define PANCAKE_ALIGNER_FACTORY_HPP

#include <pancake/AlignerBase.hpp>
#include <pancake/AlignerEdlib.hpp>
#include <pancake/AlignerKSW2.hpp>
#include <pancake/AlignerSES1.hpp>
#include <pancake/AlignerSES2.hpp>
#include <pancake/AlignmentParameters.hpp>

namespace PacBio {
namespace Pancake {

enum class AlignerType
{
    KSW2,
    EDLIB,
    SES1,
    SES2,
};

std::string AlignerTypeToString(const AlignerType& alignerType);
AlignerType AlignerTypeFromString(const std::string& alignerType);
std::shared_ptr<AlignerBase> AlignerFactory(const AlignerType& alignerType,
                                            const AlignmentParameters& alnParams);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ALIGNER_FACTORY_HPP
