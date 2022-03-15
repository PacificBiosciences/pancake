// Authors: Ivan Sovic

#include <pancake/AlignerFactory.hpp>

namespace PacBio {
namespace Pancake {

std::string AlignerTypeToString(const AlignerType alignerType)
{
    if (alignerType == AlignerType::SES2) {
        return "SES2";
    } else if (alignerType == AlignerType::KSW2) {
        return "KSW2";
    } else if (alignerType == AlignerType::EDLIB) {
        return "EDLIB";
    } else if (alignerType == AlignerType::SES1) {
        return "SES1";
    }
    return "Unknown";
}

AlignerType AlignerTypeFromString(const std::string& alignerType)
{
    if (alignerType == "KSW2") {
        return AlignerType::KSW2;
    } else if (alignerType == "EDLIB") {
        return AlignerType::EDLIB;
    } else if (alignerType == "SES1") {
        return AlignerType::SES1;
    } else if (alignerType == "SES2") {
        return AlignerType::SES2;
    }
    throw std::runtime_error("Unknown aligner type: '" + alignerType +
                             "' in AlignerTypeFromString.");
}

std::shared_ptr<AlignerBase> AlignerFactory(const AlignerType alignerType,
                                            const AlignmentParameters& alnParams)
{
    if (alignerType == AlignerType::KSW2) {
        return CreateAlignerKSW2(alnParams);

    } else if (alignerType == AlignerType::EDLIB) {
        return CreateAlignerEdlib(alnParams);

    } else if (alignerType == AlignerType::SES1) {
        return CreateAlignerSES1(alnParams);

    } else if (alignerType == AlignerType::SES2) {
        return CreateAlignerSES2(alnParams);

    } else {
        throw std::runtime_error("AlignerType " + AlignerTypeToString(alignerType) +
                                 " not supported yet!");
    }
    return nullptr;
}

}  // namespace Pancake
}  // namespace PacBio
