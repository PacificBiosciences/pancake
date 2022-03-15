// Author: Ivan Sovic

#include <pancake/SequenceSeedsCached.hpp>

namespace PacBio {
namespace Pancake {

SequenceSeedsCached::SequenceSeedsCached(std::string name, const PacBio::Pancake::Int128t* seeds,
                                         int64_t seedsSize, int64_t id)
    : name_{std::move(name)}, seeds_{seeds}, size_{seedsSize}, id_{id}
{}

const std::string& SequenceSeedsCached::Name() const { return name_; }

const PacBio::Pancake::Int128t* SequenceSeedsCached::Seeds() const { return seeds_; }

int64_t SequenceSeedsCached::Size() const { return size_; }

int64_t SequenceSeedsCached::Id() const { return id_; }

SequenceSeedsCached& SequenceSeedsCached::Name(std::string name)
{
    name_ = std::move(name);
    return *this;
}

SequenceSeedsCached& SequenceSeedsCached::Seeds(const PacBio::Pancake::Int128t* seeds)
{
    seeds_ = seeds;
    return *this;
}

SequenceSeedsCached& SequenceSeedsCached::Size(int64_t size)
{
    size_ = size;
    return *this;
}

SequenceSeedsCached& SequenceSeedsCached::Id(int64_t id)
{
    id_ = id;
    return *this;
}

}  // namespace Pancake
}  // namespace PacBio
