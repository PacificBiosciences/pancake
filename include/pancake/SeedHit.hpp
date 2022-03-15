// Author: Ivan Sovic

#ifndef PANCAKE_SEED_HIT_HPP
#define PANCAKE_SEED_HIT_HPP

#include <pancake/util/CommonTypes.hpp>

#include <cstdint>
#include <ostream>
#include <vector>

namespace PacBio {
namespace Pancake {

constexpr PacBio::Pancake::Int128t MASK128_LOW32bit = 0x000000000FFFFFFFF;
constexpr PacBio::Pancake::Int128t MASK128_LOW16bit = 0x0000000000000FFFF;
constexpr PacBio::Pancake::Int128t MASK128_LOW8bit = 0x000000000000000FF;
constexpr PacBio::Pancake::Int128t MASK128_LOW1bit = 0x00000000000000001;

constexpr int32_t SEED_HIT_FLAG_IGNORE_BIT_SET = 1 << 0;
constexpr int32_t SEED_HIT_FLAG_LONG_JOIN_BIT_SET = 1 << 1;
constexpr int32_t SEED_HIT_FLAG_IGNORE_BIT_UNSET = ~SEED_HIT_FLAG_IGNORE_BIT_SET;
constexpr int32_t SEED_HIT_FLAG_LONG_JOIN_BIT_UNSET = ~SEED_HIT_FLAG_LONG_JOIN_BIT_SET;

class alignas(sizeof(PacBio::Pancake::Int128t)) SeedHit
{
public:
    uint32_t targetId : 31;
    bool targetRev : 1;
    int32_t targetPos;
    int32_t queryPos;
    uint8_t targetSpan;
    uint8_t querySpan;
    uint16_t flags;

public:
    SeedHit() = default;

    ~SeedHit() = default;

    SeedHit(int32_t _targetId, bool _targetRev, int32_t _targetPos, int32_t _queryPos,
            int32_t _targetSpan, int32_t _querySpan, int32_t _flags)
        : targetId(_targetId)
        , targetRev(_targetRev)
        , targetPos(_targetPos)
        , queryPos(_queryPos)
        , targetSpan(_targetSpan)
        , querySpan(_querySpan)
        , flags(_flags)
    {}

public:
    bool operator<(const SeedHit& b) const { return this->PackTo128() < b.PackTo128(); }

    bool operator==(const SeedHit& b) const
    {
        return targetId == b.targetId && targetRev == b.targetRev && targetPos == b.targetPos &&
               targetSpan == b.targetSpan && querySpan == b.querySpan && flags == b.flags &&
               queryPos == b.queryPos;
    }

public:
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

    void ParseFrom128(const PacBio::Pancake::Int128t vals)
    {
        targetId = (vals >> 97) & MASK128_LOW32bit;
        targetRev = (vals >> 96) & MASK128_LOW1bit;
        targetPos = (vals >> 64) & MASK128_LOW32bit;
        queryPos = (vals >> 32) & MASK128_LOW32bit;
        targetSpan = (vals >> 24) & MASK128_LOW8bit;
        querySpan = (vals >> 16) & MASK128_LOW8bit;
        flags = vals & MASK128_LOW16bit;
    }

    int32_t Diagonal() const { return targetPos - queryPos; }

public:
    /*
     * Flags.
    */
    void SetFlagIgnore(const bool val)
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
};

inline std::ostream& operator<<(std::ostream& os, const SeedHit& a)
{
    os << "tid = " << a.targetId << ", trev = " << a.targetRev << ", tpos = " << a.targetPos
       << ", qpos = " << a.queryPos << ", tspan = " << static_cast<int32_t>(a.targetSpan)
       << ", qspan = " << static_cast<int32_t>(a.querySpan) << ", flags = " << a.flags
       << ", diag = " << (a.targetPos - a.queryPos);
    return os;
}

/**
 * @brief Constructs a tuple composed of (targetID+rev, seed diagonal, targetPos, queryPos).
 *
 * @param sh Input seed hit.
 * @return std::tuple<int32_t, int32_t, int32_t, int32_t> Elements: (targetID+rev, diagonal, targetPos, queryPos).
 */
inline std::tuple<int32_t, int32_t, int32_t, int32_t> PackSeedHitWithDiagonalToTuple(
    const SeedHit& sh)
{
    return std::make_tuple(((sh.targetId << 1) | sh.targetRev), (sh.targetPos - sh.queryPos),
                           sh.targetPos, sh.queryPos);
}

/**
 * @brief Computes the number of bases covered by seed hits, in both query and target sequences.
 *
 * @param hits Input sorted seed hits.
 * @param hitsBegin ID of the first seed hit to begin iterating from.
 * @param hitsEnd ID of the end seed hit to finish iteration.
 * @return std::pair<int32_t, int32_t> First element: query coverage. Second element: target coverage.
 */
std::pair<int32_t, int32_t> CalcHitCoverage(const std::vector<SeedHit>& hits, int32_t hitsBegin,
                                            int32_t hitsEnd);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEED_HIT_HPP
