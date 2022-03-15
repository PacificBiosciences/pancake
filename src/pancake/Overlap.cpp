// Authors: Ivan Sovic

#include <pancake/Overlap.hpp>

namespace PacBio {
namespace Pancake {

OverlapType DetermineOverlapType(bool Arev, int32_t AstartFwd, int32_t AendFwd, int32_t Alen,
                                 bool Brev, int32_t BstartFwd, int32_t BendFwd, int32_t Blen,
                                 int32_t allowedDovetailDist)
{
    int32_t leftHangA = AstartFwd;
    int32_t rightHangA = Alen - AendFwd;
    int32_t leftHangB = BstartFwd;
    int32_t rightHangB = Blen - BendFwd;

    if (Arev) {
        throw std::runtime_error(
            "The A-read should always be forward oriented. (In DetermineOverlapType.)");
    }

    if (Brev) {
        std::swap(leftHangB, rightHangB);
    }

    OverlapType ovlType = OverlapType::Internal;

    if (leftHangA <= allowedDovetailDist && rightHangA <= allowedDovetailDist) {
        /*
            This is a valid containment which would get picked up by the 5' rule
            if grace > 0. For this reason, the containment check must come first.
            Dovetail 5'.
                        left_a            right_a
                  >|/////////////|<          >|<
            A:     o-------------=============>
            B:                  o=============>
                               >|<         >||<
                               left_b     right_b

                               left_a     right_a
                               >|<         >||<
            A:                  o=============>
            B:     o-------------=============>
                  >|/////////////|<          >|<
                        left_b            right_b
        */
        ovlType = OverlapType::Contained;

    } else if (leftHangB <= allowedDovetailDist && rightHangB <= allowedDovetailDist) {
        ovlType = OverlapType::Contains;

    } else if (leftHangA <= allowedDovetailDist && rightHangB <= allowedDovetailDist) {
        /*
        Dovetail 5'.
                         left_a            right_a
                         >|//|<          >|///////|<
        A:                o--=============-------->
        B:     o-------------============->
              >|/////////////|<        >|/|<
                    left_b            right_b
        */
        ovlType = OverlapType::FivePrime;

    } else if (rightHangA <= allowedDovetailDist && leftHangB <= allowedDovetailDist) {
        /*
            Dovetail 5'.
                        left_a            right_a
                  >|/////////////|<        >|/|<
            A:     o-------------============->
            B:                o--=============-------->
                              >|//|<          >|///////|<
                              left_b            right_b
        */
        ovlType = OverlapType::ThreePrime;
    }

    return ovlType;
}

OverlapType DetermineOverlapType(const Overlap& ovl, int32_t allowedDovetailDist)
{
    return DetermineOverlapType(ovl.Arev, ovl.AstartFwd(), ovl.AendFwd(), ovl.Alen, ovl.Brev,
                                ovl.BstartFwd(), ovl.BendFwd(), ovl.Blen, allowedDovetailDist);
}

void HeuristicExtendOverlapFlanks(OverlapPtr& ovl, int32_t allowedDist)
{
    /// Note: The Overlap coordinates are internally represented in the strand
    /// of the overlap.

    int32_t leftHangA = ovl->Astart;
    int32_t rightHangA = ovl->Alen - ovl->Aend;
    int32_t leftHangB = ovl->Bstart;
    int32_t rightHangB = ovl->Blen - ovl->Bend;

    int32_t minLeft = std::min(leftHangA, leftHangB);
    int32_t minRight = std::min(rightHangA, rightHangB);

    int32_t leftFix = (minLeft > allowedDist) ? 0 : minLeft;
    int32_t rightFix = (minRight > allowedDist) ? 0 : minRight;

    ovl->Astart -= leftFix;
    ovl->Bstart -= leftFix;
    ovl->Aend += rightFix;
    ovl->Bend += rightFix;
}

OverlapPtr HeuristicExtendOverlapFlanks(const OverlapPtr& ovl, int32_t allowedDist)
{
    auto newOvl = createOverlap(ovl);
    HeuristicExtendOverlapFlanks(newOvl, allowedDist);
    return newOvl;
}

OverlapPtr ParseM4OverlapFromString(const std::string& line)
{
    auto ovl = createOverlap();
    char type[500];
    int32_t Arev = 0;
    int32_t Brev = 0;
    int32_t n =
        sscanf(line.c_str(),
               "%d %d %f %f "
               "%d %d %d %d "
               "%d %d %d %d "
               "%s",
               &(ovl->Aid), &(ovl->Bid), &(ovl->Score), &(ovl->Identity), &Arev, &(ovl->Astart),
               &(ovl->Aend), &(ovl->Alen), &Brev, &(ovl->Bstart), &(ovl->Bend), &(ovl->Blen), type);

    ovl->Arev = Arev;
    ovl->Brev = Brev;
    ovl->Identity /= 100.0f;

    // Internally we represent the overlap in the strand of the overlap.
    // (The format specifies it always in the FWD strand.)
    if (ovl->Arev) {
        std::swap(ovl->Astart, ovl->Aend);
        ovl->Astart = ovl->Alen - ovl->Astart;
        ovl->Aend = ovl->Alen - ovl->Aend;
    }
    if (ovl->Brev) {
        std::swap(ovl->Bstart, ovl->Bend);
        ovl->Bstart = ovl->Blen - ovl->Bstart;
        ovl->Bend = ovl->Blen - ovl->Bend;
    }

    ovl->Atype = OverlapType::Unknown;
    ovl->Btype = OverlapType::Unknown;

    // If the type is specified in the line, parse it.
    if (n >= 13) {
        ovl->Atype = OverlapTypeFromString(type);
    }
    return ovl;
}

std::string OverlapTypeToString(const OverlapType& type)
{
    std::string ret = "*";
    switch (type) {
        case OverlapType::Unknown:
            ret = "*";
            break;
        case OverlapType::Contained:
            ret = "contained";
            break;
        case OverlapType::Contains:
            ret = "contains";
            break;
        case OverlapType::FivePrime:
            ret = "5";
            break;
        case OverlapType::ThreePrime:
            ret = "3";
            break;
        case OverlapType::Internal:
            ret = "u";
            break;
        default:
            ret = "*";
    }
    return ret;
}

std::string OverlapTypeToStringSingleChar(const OverlapType& type)
{
    std::string ret = "*";
    switch (type) {
        case OverlapType::Unknown:
            ret = "*";
            break;
        case OverlapType::Contained:
            ret = "c";
            break;
        case OverlapType::Contains:
            ret = "C";
            break;
        case OverlapType::FivePrime:
            ret = "5";
            break;
        case OverlapType::ThreePrime:
            ret = "3";
            break;
        case OverlapType::Internal:
            ret = "u";
            break;
        default:
            ret = "*";
    }
    return ret;
}

OverlapType OverlapTypeFromString(const std::string& typeStr)
{
    OverlapType type = OverlapType::Unknown;
    if (typeStr == "5") {
        type = OverlapType::FivePrime;
    } else if (typeStr == "3") {
        type = OverlapType::ThreePrime;
    } else if (typeStr == "contained" || typeStr == "c") {
        type = OverlapType::Contained;
    } else if (typeStr == "contains" || typeStr == "C") {
        type = OverlapType::Contains;
    } else if (typeStr == "u") {
        type = OverlapType::Internal;
    }
    return type;
}

}  // namespace Pancake
}  // namespace PacBio
