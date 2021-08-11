// Author: Ivan Sovic

#ifndef PANCAKE_OVERLAPHIFI_SEEDHIT_H
#define PANCAKE_OVERLAPHIFI_SEEDHIT_H

#include <pacbio/util/CommonTypes.h>
#include <cstdint>
#include <ostream>
#include <vector>

namespace PacBio {
namespace Pancake {

static const PacBio::Pancake::Int128t MASK128_LOW32bit = 0x000000000FFFFFFFF;
static const PacBio::Pancake::Int128t MASK128_LOW16bit = 0x0000000000000FFFF;
static const PacBio::Pancake::Int128t MASK128_LOW8bit = 0x000000000000000FF;
static const PacBio::Pancake::Int128t MASK128_LOW1bit = 0x00000000000000001;
static const int32_t SEED_HIT_FLAG_IGNORE_BIT_SET = 1 << 0;
static const int32_t SEED_HIT_FLAG_LONG_JOIN_BIT_SET = 1 << 1;

static const int32_t SEED_HIT_FLAG_IGNORE_BIT_UNSET = ~SEED_HIT_FLAG_IGNORE_BIT_SET;
static const int32_t SEED_HIT_FLAG_LONG_JOIN_BIT_UNSET = ~SEED_HIT_FLAG_LONG_JOIN_BIT_SET;

class SeedHit
{
public:
    SeedHit() = default;
    ~SeedHit() = default;

    SeedHit(const int32_t _targetId, const bool _targetRev, const int32_t _targetPos,
            const int32_t _queryPos, int32_t _targetSpan, const int32_t _querySpan, int32_t _flags)
        : targetId(_targetId)
        , targetRev(_targetRev)
        , targetPos(_targetPos)
        , queryPos(_queryPos)
        , targetSpan(_targetSpan)
        , querySpan(_querySpan)
        , flags(_flags)
    {
    }

    PacBio::Pancake::Int128t PackTo128() const
    {
        PacBio::Pancake::Int128t ret = 0;
        ret = ((static_cast<PacBio::Pancake::Int128t>(targetId) & MASK128_LOW32bit) << 97) |
              ((static_cast<PacBio::Pancake::Int128t>(targetRev) & MASK128_LOW1bit) << 96) |
              ((static_cast<PacBio::Pancake::Int128t>(targetPos) & MASK128_LOW32bit) << 64) |
              ((static_cast<PacBio::Pancake::Int128t>(queryPos) & MASK128_LOW32bit) << 32) |
              ((static_cast<PacBio::Pancake::Int128t>(targetSpan) & MASK128_LOW8bit) << 24) |
              ((static_cast<PacBio::Pancake::Int128t>(querySpan) & MASK128_LOW8bit) << 16) |
              ((static_cast<PacBio::Pancake::Int128t>(flags) & MASK128_LOW16bit) << 0);
        return ret;
    }
    bool operator<(const SeedHit& b) const { return this->PackTo128() < b.PackTo128(); }
    bool operator==(const SeedHit& b) const
    {
        return targetId == b.targetId && targetRev == b.targetRev && targetPos == b.targetPos &&
               targetSpan == b.targetSpan && querySpan == b.querySpan && flags == b.flags &&
               queryPos == b.queryPos;
    }
    int32_t Diagonal() const { return targetPos - queryPos; }

    /*
     * Flags.
    */
    void SetFlagIgnore(bool val)
    {
        if (val) {
            flags |= SEED_HIT_FLAG_IGNORE_BIT_SET;
        } else {
            flags &= SEED_HIT_FLAG_IGNORE_BIT_UNSET;
        }
    }
    void SetFlagIgnore() { flags |= SEED_HIT_FLAG_IGNORE_BIT_SET; }
    void UnsetFlagIgnore() { flags &= SEED_HIT_FLAG_IGNORE_BIT_UNSET; }
    bool CheckFlagIgnore() const { return (flags & SEED_HIT_FLAG_IGNORE_BIT_SET); }

    void SetFlagLongJoin(bool val)
    {
        if (val) {
            flags |= SEED_HIT_FLAG_IGNORE_BIT_SET;
        } else {
            flags &= SEED_HIT_FLAG_IGNORE_BIT_UNSET;
        }
    }
    void SetFlagLongJoin() { flags |= SEED_HIT_FLAG_IGNORE_BIT_SET; }
    void UnsetFlagLongJoin() { flags &= SEED_HIT_FLAG_IGNORE_BIT_UNSET; }
    bool CheckFlagLongJoin() const { return (flags & SEED_HIT_FLAG_IGNORE_BIT_SET); }

    /*
     * Data.
    */
    uint32_t targetId : 31;
    bool targetRev : 1;
    int32_t targetPos;
    int32_t queryPos;
    uint8_t targetSpan;
    uint8_t querySpan;
    uint16_t flags;
};

inline std::ostream& operator<<(std::ostream& os, const SeedHit& a)
{
    os << "tid = " << a.targetId << ", trev = " << a.targetRev << ", tpos = " << a.targetPos
       << ", qpos = " << a.queryPos << ", tspan = " << a.targetSpan << ", qspan = " << a.querySpan
       << ", flags = " << a.flags << ", diag = " << (a.targetPos - a.queryPos);
    return os;
}

inline std::tuple<int32_t, int32_t, int32_t, int32_t> PackSeedHitWithDiagonalToTuple(
    const SeedHit& sh)
{
    return std::make_tuple(((sh.targetId << 1) | sh.targetRev), (sh.targetPos - sh.queryPos),
                           sh.targetPos, sh.queryPos);
}

void CalcHitCoverage(const std::vector<SeedHit>& hits, int32_t hitsBegin, int32_t hitsEnd,
                     int32_t& coveredBasesQuery, int32_t& coveredBasesTarget);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_OVERLAPHIFI_SEEDHIT_H