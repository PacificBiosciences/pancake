// Authors: Ivan Sovic

#include <pancake/Overlap.hpp>

#include <tuple>

namespace PacBio {
namespace Pancake {

Overlap::Overlap(const int32_t _Aid, const int32_t _Bid, const float _Score, const float _Identity,
                 const bool _Arev, const int32_t _Astart, const int32_t _Aend, const int32_t _Alen,
                 const bool _Brev, const int32_t _Bstart, const int32_t _Bend, const int32_t _Blen,
                 const int32_t _EditDistance, const int32_t _NumSeeds, const OverlapType _Atype,
                 const OverlapType _Btype, const PacBio::BAM::Cigar& _Cigar,
                 const std::string_view _Avars, const std::string_view _Bvars, bool _IsFlipped,
                 const bool _IsSupplementary, const bool _IsSecondary)
    : Aid(_Aid)
    , Arev(_Arev)
    , Astart(_Astart)
    , Aend(_Aend)
    , Alen(_Alen)
    , Bid(_Bid)
    , Brev(_Brev)
    , Bstart(_Bstart)
    , Bend(_Bend)
    , Blen(_Blen)
    , Score(_Score)
    , Identity(_Identity)
    , EditDistance(_EditDistance)
    , NumSeeds(_NumSeeds)
    , Atype(_Atype)
    , Btype(_Btype)
    , Cigar(_Cigar)
    , Avars(_Avars)
    , Bvars(_Bvars)
    , IsFlipped(_IsFlipped)
    , IsSupplementary(_IsSupplementary)
    , IsSecondary(_IsSecondary)
{}

Overlap::Overlap(const int32_t _Aid, const int32_t _Bid, const float _Score, const float _Identity,
                 const bool _Arev, const int32_t _Astart, const int32_t _Aend, const int32_t _Alen,
                 const bool _Brev, const int32_t _Bstart, const int32_t _Bend, const int32_t _Blen)
    : Aid(_Aid)
    , Arev(_Arev)
    , Astart(_Astart)
    , Aend(_Aend)
    , Alen(_Alen)
    , Bid(_Bid)
    , Brev(_Brev)
    , Bstart(_Bstart)
    , Bend(_Bend)
    , Blen(_Blen)
    , Score(_Score)
    , Identity(_Identity)
    , EditDistance(0)
    , NumSeeds(0)
    , Atype(OverlapType::Unknown)
    , Btype(OverlapType::Unknown)
    , Cigar()
    , Avars()
    , Bvars()
    , IsFlipped(false)
    , IsSupplementary(false)
    , IsSecondary(false)
{}

void Overlap::Flip()
{
    IsFlipped = !IsFlipped;

    std::swap(Aid, Bid);
    std::swap(Arev, Brev);
    std::swap(Alen, Blen);
    std::swap(Avars, Bvars);
    std::swap(Atype, Btype);

    // If the query/target context changed, then I/D operations need to
    // be updated.
    if (Cigar.size() > 0) {
        for (auto& op : Cigar) {
            if (op.Type() == PacBio::BAM::CigarOperationType::INSERTION) {
                op.Type(PacBio::BAM::CigarOperationType::DELETION);
            } else if (op.Type() == PacBio::BAM::CigarOperationType::DELETION) {
                op.Type(PacBio::BAM::CigarOperationType::INSERTION);
            }
        }
    }

    // Keep the A sequence in the fwd direction at all times.
    if (Arev) {
        // Reorient the coordinates. Internally, all coordinates are IN-STRAND
        // of the sequence, so if A orientation is being flipped, so should the
        // coordinates be.
        std::swap(Astart, Bend);
        std::swap(Aend, Bstart);
        Astart = Alen - Astart;
        Aend = Alen - Aend;
        Arev = !Arev;
        Bstart = Blen - Bstart;
        Bend = Blen - Bend;
        Brev = !Brev;

        // Reverse the CIGAR string.
        if (Cigar.size() > 0) {
            std::reverse(Cigar.begin(), Cigar.end());
        }

        // Reverse the variant positions.
        Avars = Pancake::ReverseComplement(Avars, 0, Avars.size());
        Bvars = Pancake::ReverseComplement(Bvars, 0, Bvars.size());
    } else {
        std::swap(Astart, Bstart);
        std::swap(Aend, Bend);
    }
}

void Overlap::NormalizeStrand()
{
    /*
     * If the A-read is reversed, this function normalizes the overlap so that
     * the A-read is always facing in the forward direction.
     */

    if (Arev == false) {
        return;
    }

    std::swap(Astart, Aend);
    Astart = Alen - Astart;
    Aend = Alen - Aend;
    Arev = !Arev;

    std::swap(Bstart, Bend);
    Bstart = Blen - Bstart;
    Bend = Blen - Bend;
    Brev = !Brev;

    std::reverse(Cigar.begin(), Cigar.end());
    std::reverse(Avars.begin(), Avars.end());
    std::reverse(Bvars.begin(), Bvars.end());
}

bool Overlap::operator==(const Overlap& rhs) const
{
    return std::tie(Aid, Bid, Score, Identity, EditDistance, NumSeeds, Atype, Btype, Arev, Astart,
                    Aend, Alen, Brev, Bstart, Bend, Blen, Cigar, Avars, Bvars, IsFlipped,
                    IsSupplementary, IsSecondary) ==
           std::tie(rhs.Aid, rhs.Bid, rhs.Score, rhs.Identity, rhs.EditDistance, rhs.NumSeeds,
                    rhs.Atype, rhs.Btype, rhs.Arev, rhs.Astart, rhs.Aend, rhs.Alen, rhs.Brev,
                    rhs.Bstart, rhs.Bend, rhs.Blen, rhs.Cigar, rhs.Avars, rhs.Bvars, rhs.IsFlipped,
                    rhs.IsSupplementary, rhs.IsSecondary);
}

std::unique_ptr<Overlap> CreateOverlap() { return std::unique_ptr<Overlap>(new Overlap()); }

std::unique_ptr<Overlap> CreateOverlap(
    const int32_t Aid, const int32_t Bid, const float score, const float identity, const bool Arev,
    const int32_t Astart, const int32_t Aend, const int32_t Alen, const bool Brev,
    const int32_t Bstart, const int32_t Bend, const int32_t Blen, const int32_t EditDistance,
    const int32_t NumSeeds, const OverlapType Atype, const OverlapType Btype,
    const PacBio::BAM::Cigar& Cigar, const std::string_view Avars, const std::string_view Bvars,
    const bool IsFlipped, const bool IsSupplementary, const bool IsSecondary)
{
    return std::unique_ptr<Overlap>(new Overlap(
        Aid, Bid, score, identity, Arev, Astart, Aend, Alen, Brev, Bstart, Bend, Blen, EditDistance,
        NumSeeds, Atype, Btype, Cigar, Avars, Bvars, IsFlipped, IsSupplementary, IsSecondary));
}

std::unique_ptr<Overlap> CreateOverlap(const int32_t Aid, const int32_t Bid, const float score,
                                       const float identity, const bool Arev, const int32_t Astart,
                                       const int32_t Aend, const int32_t Alen, bool Brev,
                                       int32_t Bstart, int32_t Bend, int32_t Blen,
                                       const int32_t EditDistance, const int32_t NumSeeds,
                                       const OverlapType Atype, const OverlapType Btype)
{
    return std::unique_ptr<Overlap>(new Overlap(Aid, Bid, score, identity, Arev, Astart, Aend, Alen,
                                                Brev, Bstart, Bend, Blen, EditDistance, NumSeeds,
                                                Atype, Btype, {}, {}, {}, false, false, false));
}

std::unique_ptr<Overlap> CreateOverlap(const Overlap& ovl)
{
    return std::unique_ptr<Overlap>(
        new Overlap(ovl.Aid, ovl.Bid, ovl.Score, ovl.Identity, ovl.Arev, ovl.Astart, ovl.Aend,
                    ovl.Alen, ovl.Brev, ovl.Bstart, ovl.Bend, ovl.Blen, ovl.EditDistance,
                    ovl.NumSeeds, ovl.Atype, ovl.Btype, ovl.Cigar, ovl.Avars, ovl.Bvars,
                    ovl.IsFlipped, ovl.IsSupplementary, ovl.IsSecondary));
}

std::unique_ptr<Overlap> CreateOverlap(const std::unique_ptr<Overlap>& ovl)
{
    return CreateOverlap(*ovl);
}

OverlapPtr CreateFlippedOverlap(const OverlapPtr& ovl)
{
    auto newOvl = CreateOverlap(ovl);
    newOvl->Flip();
    return newOvl;
}

OverlapType DetermineOverlapType(const bool Arev, const int32_t AstartFwd, const int32_t AendFwd,
                                 const int32_t Alen, const bool Brev, const int32_t BstartFwd,
                                 const int32_t BendFwd, const int32_t Blen,
                                 const int32_t allowedDovetailDist)
{
    const int32_t leftHangA = AstartFwd;
    const int32_t rightHangA = Alen - AendFwd;
    int32_t leftHangB = BstartFwd;
    int32_t rightHangB = Blen - BendFwd;

    if (Arev) {
        throw std::runtime_error(
            "The A-read should always be forward oriented. (In DetermineOverlapType.)");
    }

    if (Brev) {
        std::swap(leftHangB, rightHangB);
    }

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
        return OverlapType::Contained;

    } else if (leftHangB <= allowedDovetailDist && rightHangB <= allowedDovetailDist) {
        return OverlapType::Contains;

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
        return OverlapType::FivePrime;

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
        return OverlapType::ThreePrime;
    }

    return OverlapType::Internal;
}

OverlapType DetermineOverlapType(const Overlap& ovl, const int32_t allowedDovetailDist)
{
    return DetermineOverlapType(ovl.Arev, ovl.AstartFwd(), ovl.AendFwd(), ovl.Alen, ovl.Brev,
                                ovl.BstartFwd(), ovl.BendFwd(), ovl.Blen, allowedDovetailDist);
}

void HeuristicExtendOverlapFlanks(OverlapPtr& ovl, const int32_t allowedDist)
{
    /// Note: The Overlap coordinates are internally represented in the strand
    /// of the overlap.

    const int32_t leftHangA = ovl->Astart;
    const int32_t rightHangA = ovl->Alen - ovl->Aend;
    const int32_t leftHangB = ovl->Bstart;
    const int32_t rightHangB = ovl->Blen - ovl->Bend;

    const int32_t minLeft = std::min(leftHangA, leftHangB);
    const int32_t minRight = std::min(rightHangA, rightHangB);

    const int32_t leftFix = (minLeft > allowedDist) ? 0 : minLeft;
    const int32_t rightFix = (minRight > allowedDist) ? 0 : minRight;

    ovl->Astart -= leftFix;
    ovl->Bstart -= leftFix;
    ovl->Aend += rightFix;
    ovl->Bend += rightFix;
}

OverlapPtr HeuristicExtendOverlapFlanks(const OverlapPtr& ovl, const int32_t allowedDist)
{
    auto newOvl = CreateOverlap(ovl);
    HeuristicExtendOverlapFlanks(newOvl, allowedDist);
    return newOvl;
}

OverlapPtr ParseM4OverlapFromString(const std::string_view line)
{
    auto ovl = CreateOverlap();
    char type[500];
    int32_t Arev = 0;
    int32_t Brev = 0;
    int32_t n =
        sscanf(line.data(),
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

OverlapType OverlapTypeFromString(const std::string_view typeStr)
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
