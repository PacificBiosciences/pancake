/*
 * seed.hpp
 *
 * This source file was obtained and adjusted from the Raptor
 * graph-based mapping tool codebase.
 *
 *  Created on: Oct 04, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_MINIMIZER_INDEX2_MINIMIZER_H_
#define SRC_MINIMIZER_INDEX2_MINIMIZER_H_

#include <pacbio/util/CommonTypes.h>
#include <cstdint>
#include <sstream>

namespace PacBio {
namespace Pancake {
namespace SeedDB {

static const int8_t MINIMIZER_FLAG_DEFAULT_FWD = 0x00;
static const int8_t MINIMIZER_FLAG_IS_REV = 0x01;
static const PacBio::Pancake::Int128t MINIMIZER_CODED_REV_BIT =
    (static_cast<PacBio::Pancake::Int128t>(1) << 32);

static const PacBio::Pancake::Int128t MINIMIZER_64bit_MASK = 0x0FFFFFFFFFFFFFFFF;
static const PacBio::Pancake::Int128t MINIMIZER_32bit_MASK = 0x0000000007FFFFFFF;
static const PacBio::Pancake::Int128t MINIMIZER_32bit_MASK_FULL = 0x000000000FFFFFFFF;
static const PacBio::Pancake::Int128t MINIMIZER_8bit_MASK = 0x000000000000000FF;
static const PacBio::Pancake::Int128t MINIMIZER_1bit_MASK = 0x00000000000000001;

using SeedRaw = PacBio::Pancake::Int128t;

class Seed
{
public:
    Seed() : key(0), span(0), seqID(0), seqRev(0), pos(0) {}

    Seed(uint64_t _key, int32_t _span, int32_t _seqID, int32_t _pos, bool _isRev)
        : key(_key), span(_span), seqID(_seqID), seqRev(_isRev), pos(_pos)
    {
        if (_span >= 256 || _span < 0) {
            span = 0;
        }
    }
    Seed(const PacBio::Pancake::Int128t& codedKeypos)
        : key(DecodeKey(codedKeypos))
        , span(DecodeSpan(codedKeypos))
        , seqID(DecodeSeqId(codedKeypos))
        , seqRev(DecodeIsRev(codedKeypos))
        , pos(DecodePos(codedKeypos))
    {
    }

    bool IsRev() const { return seqRev; }

    bool Valid() const { return span != 0; }

    void SetInvalid() { span = 0; }

    inline PacBio::Pancake::Int128t To128t()
    {
        PacBio::Pancake::Int128t ret = static_cast<PacBio::Pancake::Int128t>(key) << (64 + 8);
        ret |= (static_cast<PacBio::Pancake::Int128t>(span) & MINIMIZER_8bit_MASK) << 64;
        ret |= ((static_cast<PacBio::Pancake::Int128t>(seqID) << 1) & MINIMIZER_32bit_MASK) << 32;
        ret |= (static_cast<PacBio::Pancake::Int128t>(seqRev) & MINIMIZER_1bit_MASK) << 32;
        ret |= static_cast<PacBio::Pancake::Int128t>(pos) & MINIMIZER_32bit_MASK;
        return ret;
    }

    static inline PacBio::Pancake::Int128t Encode(uint64_t _key, int32_t _span, int32_t _seqID,
                                                  int32_t _pos, bool _isRev)
    {
        // If the specified span is not valid, set it to zero. This indicates that the seed
        // is not valid (zero-length seeds are not allowed).
        if (_span >= 256 || _span < 0) {
            _span = 0;
        }
        PacBio::Pancake::Int128t ret = static_cast<PacBio::Pancake::Int128t>(_key) << 72;
        ret |= (static_cast<PacBio::Pancake::Int128t>(_span) & MINIMIZER_8bit_MASK) << 64;
        ret |= ((static_cast<PacBio::Pancake::Int128t>(_seqID) << 1) & MINIMIZER_32bit_MASK) << 32;
        ret |= (static_cast<PacBio::Pancake::Int128t>(_isRev) & MINIMIZER_1bit_MASK) << 32;
        ret |= static_cast<PacBio::Pancake::Int128t>(_pos) & MINIMIZER_32bit_MASK;
        return ret;
    }

    static inline uint64_t DecodeKey(const PacBio::Pancake::Int128t& seed)
    {
        return ((seed >> (64 + 8)) & MINIMIZER_64bit_MASK);
    }

    static inline uint64_t DecodeSpan(const PacBio::Pancake::Int128t& seed)
    {
        return ((seed >> 64) & MINIMIZER_8bit_MASK);
    }

    static inline int32_t DecodePos(const PacBio::Pancake::Int128t& seed)
    {
        return (seed & MINIMIZER_32bit_MASK);
    }

    static inline int32_t DecodeSeqId(const PacBio::Pancake::Int128t& seed)
    {
        return ((seed >> 33) & MINIMIZER_32bit_MASK);
    }

    /*
     * Unlike DecodeSeqId, this method keeps the info about
     * of the strand still encoded in the return value.
     * The LSB is 0 for fwd, and 1 for rev.
     */
    static inline int32_t DecodeSeqIdWithRev(const PacBio::Pancake::Int128t& seed)
    {
        return ((seed >> 32) & MINIMIZER_32bit_MASK);
    }

    static inline bool DecodeIsRev(const PacBio::Pancake::Int128t& seed)
    {
        return (seed & MINIMIZER_CODED_REV_BIT);
    }

    std::string Verbose() const
    {
        std::stringstream ss;
        ss << "pos = " << pos << ", span = " << span << ", seqID = " << seqID
           << ", seqRev = " << seqRev << ", key = " << key;
        return ss.str();
    }

    uint64_t key : 56;
    uint64_t span : 8;
    uint32_t seqID : 31;
    bool seqRev : 1;
    uint32_t pos;
};

inline std::ostream& operator<<(std::ostream& os, const Seed& b)
{
    os << "pos = " << b.pos << ", span = " << b.span << ", seqID = " << b.seqID
       << ", seqRev = " << b.seqRev << ", key = " << b.key;
    return os;
}

}  // namespace SeedDB
}  // namespace Pancake
}  // namespace PacBio

#endif
