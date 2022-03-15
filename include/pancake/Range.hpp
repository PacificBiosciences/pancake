// Author: Ivan Sovic

#ifndef PANCAKE_RANGE_HPP
#define PANCAKE_RANGE_HPP

#include <cstdint>
#include <ostream>

namespace PacBio {
namespace Pancake {

class Range
{
public:
    int32_t start = 0;
    int32_t end = 0;

public:
    int32_t Span() const { return end - start; }

    bool operator==(const PacBio::Pancake::Range& b) const
    {
        return start == b.start && end == b.end;
    }

    friend std::ostream& operator<<(std::ostream& os, const PacBio::Pancake::Range& mr);
};

inline std::ostream& operator<<(std::ostream& os, const PacBio::Pancake::Range& r)
{
    os << "start = " << r.start << ", end = " << r.end;
    return os;
}

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_RANGE_HPP
