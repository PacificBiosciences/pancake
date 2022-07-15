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
    float Identity = 0.0f;  // Range [0.0, 1.0].
    int32_t EditDistance = -1;
    int32_t NumSeeds = -1;

    OverlapType Atype = OverlapType::Unknown;
    OverlapType Btype = OverlapType::Unknown;

    PacBio::BAM::Cigar Cigar;
    std::string Avars;
    std::string Bvars;

    // This is important to mark whether an overlap was flipped, because the Aid and Bid contexts change.
    bool IsFlipped = false;

    bool IsSupplementary = false;
    bool IsSecondary = false;

public:
    Overlap() = default;

    ~Overlap() = default;

    Overlap(int32_t _Aid, int32_t _Bid, float _Score, float _Identity, bool _Arev, int32_t _Astart,
            int32_t _Aend, int32_t _Alen, bool _Brev, int32_t _Bstart, int32_t _Bend, int32_t _Blen,
            int32_t _EditDistance, int32_t _NumSeeds, OverlapType _Atype, OverlapType _Btype,
            const PacBio::BAM::Cigar& _Cigar, std::string_view _Avars, std::string_view _Bvars,
            bool _IsFlipped, bool _IsSupplementary, bool _IsSecondary);

    Overlap(int32_t _Aid, int32_t _Bid, float _Score, float _Identity, bool _Arev, int32_t _Astart,
            int32_t _Aend, int32_t _Alen, bool _Brev, int32_t _Bstart, int32_t _Bend,
            int32_t _Blen);

public:
    int32_t ASpan() const { return (Aend - Astart); }

    int32_t BSpan() const { return (Bend - Bstart); }

    int32_t AstartFwd() const { return (Arev ? (Alen - Aend) : Astart); }

    int32_t AendFwd() const { return (Arev ? (Alen - Astart) : Aend); }

    int32_t BstartFwd() const { return (Brev ? (Blen - Bend) : Bstart); }

    int32_t BendFwd() const { return (Brev ? (Blen - Bstart) : Bend); }

    void Flip();

    void NormalizeStrand();

public:
    bool operator==(const Overlap& rhs) const;
};

std::unique_ptr<Overlap> CreateOverlap();

std::unique_ptr<Overlap> CreateOverlap(int32_t Aid, int32_t Bid, float score, float identity,
                                       bool Arev, int32_t Astart, int32_t Aend, int32_t Alen,
                                       bool Brev, int32_t Bstart, int32_t Bend, int32_t Blen,
                                       int32_t EditDistance, int32_t NumSeeds, OverlapType Atype,
                                       OverlapType Btype, const PacBio::BAM::Cigar& Cigar,
                                       std::string_view Avars, std::string_view Bvars,
                                       bool IsFlipped, bool IsSupplementary, bool IsSecondary);

std::unique_ptr<Overlap> CreateOverlap(int32_t Aid, int32_t Bid, float score, float identity,
                                       bool Arev, int32_t Astart, int32_t Aend, int32_t Alen,
                                       bool Brev, int32_t Bstart, int32_t Bend, int32_t Blen,
                                       int32_t EditDistance, int32_t NumSeeds, OverlapType Atype,
                                       OverlapType Btype);

std::unique_ptr<Overlap> CreateOverlap(const Overlap& ovl);

std::unique_ptr<Overlap> CreateOverlap(const std::unique_ptr<Overlap>& ovl);

OverlapPtr CreateFlippedOverlap(const OverlapPtr& ovl);

OverlapType DetermineOverlapType(bool Arev, int32_t AstartFwd, int32_t AendFwd, int32_t Alen,
                                 bool Brev, int32_t BstartFwd, int32_t BendFwd, int32_t Blen,
                                 int32_t allowedDovetailDist);

OverlapType DetermineOverlapType(const Overlap& ovl, int32_t allowedDovetailDist);

void HeuristicExtendOverlapFlanks(OverlapPtr& ovl, int32_t allowedDist);

OverlapPtr HeuristicExtendOverlapFlanks(const OverlapPtr& ovl, int32_t allowedDist);

OverlapPtr ParseM4OverlapFromString(std::string_view line);

std::string OverlapTypeToString(const OverlapType& type);

std::string OverlapTypeToStringSingleChar(const OverlapType& type);

OverlapType OverlapTypeFromString(std::string_view typeStr);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_OVERLAP_HPP
