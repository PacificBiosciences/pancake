// Author: Ivan Sovic

#ifndef PANCAKE_OVERLAP_HPP
#define PANCAKE_OVERLAP_HPP

#include <pancake/util/Util.hpp>

#include <pbbam/Cigar.h>

#include <algorithm>
#include <cstdint>
#include <memory>
#include <string>

namespace PacBio {
namespace Pancake {

class Overlap;

using OverlapPtr = std::unique_ptr<Overlap>;

enum class OverlapType
{
    Unknown,
    Internal,
    Contained,
    Contains,
    FivePrime,
    ThreePrime
};

OverlapType DetermineOverlapType(bool Arev, int32_t AstartFwd, int32_t AendFwd, int32_t Alen,
                                 bool Brev, int32_t BstartFwd, int32_t BendFwd, int32_t Blen,
                                 int32_t allowedDovetailDist);
OverlapType DetermineOverlapType(const Overlap& ovl, int32_t allowedDovetailDist);
void HeuristicExtendOverlapFlanks(OverlapPtr& ovl, int32_t allowedDist);
OverlapPtr HeuristicExtendOverlapFlanks(const OverlapPtr& ovl, int32_t allowedDist);
OverlapPtr ParseM4OverlapFromString(const std::string& line);
std::string OverlapTypeToString(const OverlapType& type);
std::string OverlapTypeToStringSingleChar(const OverlapType& type);
OverlapType OverlapTypeFromString(const std::string& typeStr);

class Overlap
{
public:
    int32_t Aid = -1;
    bool Arev = false;
    int32_t Astart = 0;
    int32_t Aend = 0;
    int32_t Alen = 0;

    int32_t Bid = -1;
    bool Brev = false;
    int32_t Bstart = 0;
    int32_t Bend = 0;
    int32_t Blen = 0;

    float Score = 0.0f;
    float Identity = 0.0f;
    int32_t EditDistance = -1;
    int32_t NumSeeds = -1;

    OverlapType Atype = OverlapType::Unknown;
    OverlapType Btype = OverlapType::Unknown;

    PacBio::BAM::Cigar Cigar;
    std::string Avars;
    std::string Bvars;

    // Important to mark whether an overlap was flipped, because the Aid and Bid contexts change.
    bool IsFlipped = false;

    bool IsSupplementary = false;
    bool IsSecondary = false;

public:
    Overlap() = default;
    ~Overlap() = default;

    Overlap(int32_t _Aid, int32_t _Bid, float _Score, float _Identity, bool _Arev, int32_t _Astart,
            int32_t _Aend, int32_t _Alen, bool _Brev, int32_t _Bstart, int32_t _Bend, int32_t _Blen,
            int32_t _EditDistance, int32_t _NumSeeds, OverlapType _Atype, OverlapType _Btype,
            const PacBio::BAM::Cigar& _Cigar, const std::string& _Avars, const std::string& _Bvars,
            bool _IsFlipped, bool _IsSupplementary, bool _IsSecondary)
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

    Overlap(int32_t _Aid, int32_t _Bid, float _Score, float _Identity, bool _Arev, int32_t _Astart,
            int32_t _Aend, int32_t _Alen, bool _Brev, int32_t _Bstart, int32_t _Bend, int32_t _Blen)
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

public:
    int32_t ASpan() const { return (Aend - Astart); }
    int32_t BSpan() const { return (Bend - Bstart); }
    int32_t AstartFwd() const { return (Arev ? (Alen - Aend) : Astart); }
    int32_t AendFwd() const { return (Arev ? (Alen - Astart) : Aend); }
    int32_t BstartFwd() const { return (Brev ? (Blen - Bend) : Bstart); }
    int32_t BendFwd() const { return (Brev ? (Blen - Bstart) : Bend); }

    void Flip()
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

    void NormalizeStrand()
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

public:
    bool operator==(const Overlap& rhs) const
    {
        return Aid == rhs.Aid && Bid == rhs.Bid && Score == rhs.Score && Identity == rhs.Identity &&
               EditDistance == rhs.EditDistance && NumSeeds == rhs.NumSeeds && Atype == rhs.Atype &&
               Btype == rhs.Btype && Arev == rhs.Arev && Astart == rhs.Astart && Aend == rhs.Aend &&
               Alen == rhs.Aend && Brev == rhs.Brev && Bstart == rhs.Bstart && Bend == rhs.Bend &&
               Blen == rhs.Bend && Cigar == rhs.Cigar && Avars == rhs.Avars && Bvars == rhs.Bvars &&
               IsFlipped == rhs.IsFlipped && IsSupplementary == rhs.IsSupplementary &&
               IsSecondary == rhs.IsSecondary;
    }
};

inline std::unique_ptr<Overlap> createOverlap() { return std::unique_ptr<Overlap>(new Overlap()); }

inline std::unique_ptr<Overlap> createOverlap(
    int32_t Aid, int32_t Bid, float score, float identity, bool Arev, int32_t Astart, int32_t Aend,
    int32_t Alen, bool Brev, int32_t Bstart, int32_t Bend, int32_t Blen, int32_t EditDistance,
    int32_t NumSeeds, OverlapType Atype, OverlapType Btype, const PacBio::BAM::Cigar& Cigar,
    const std::string& Avars, const std::string& Bvars, bool IsFlipped, bool IsSupplementary,
    bool IsSecondary)
{
    return std::unique_ptr<Overlap>(new Overlap(
        Aid, Bid, score, identity, Arev, Astart, Aend, Alen, Brev, Bstart, Bend, Blen, EditDistance,
        NumSeeds, Atype, Btype, Cigar, Avars, Bvars, IsFlipped, IsSupplementary, IsSecondary));
}

inline std::unique_ptr<Overlap> createOverlap(int32_t Aid, int32_t Bid, float score, float identity,
                                              bool Arev, int32_t Astart, int32_t Aend, int32_t Alen,
                                              bool Brev, int32_t Bstart, int32_t Bend, int32_t Blen,
                                              int32_t EditDistance, int32_t NumSeeds,
                                              OverlapType Atype, OverlapType Btype)
{
    return std::unique_ptr<Overlap>(new Overlap(Aid, Bid, score, identity, Arev, Astart, Aend, Alen,
                                                Brev, Bstart, Bend, Blen, EditDistance, NumSeeds,
                                                Atype, Btype, {}, {}, {}, false, false, false));
}

inline std::unique_ptr<Overlap> createOverlap(const Overlap& ovl)
{
    return std::unique_ptr<Overlap>(
        new Overlap(ovl.Aid, ovl.Bid, ovl.Score, ovl.Identity, ovl.Arev, ovl.Astart, ovl.Aend,
                    ovl.Alen, ovl.Brev, ovl.Bstart, ovl.Bend, ovl.Blen, ovl.EditDistance,
                    ovl.NumSeeds, ovl.Atype, ovl.Btype, ovl.Cigar, ovl.Avars, ovl.Bvars,
                    ovl.IsFlipped, ovl.IsSupplementary, ovl.IsSecondary));
}

inline std::unique_ptr<Overlap> createOverlap(const std::unique_ptr<Overlap>& ovl)
{
    return createOverlap(*ovl);
}

inline OverlapPtr CreateFlippedOverlap(const OverlapPtr& ovl)
{
    auto newOvl = createOverlap(ovl);
    newOvl->Flip();
    return newOvl;
}

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_OVERLAP_HPP
