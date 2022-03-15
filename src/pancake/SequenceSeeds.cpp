// Author: Ivan Sovic

#include <cassert>
#include <cstdio>
#include <exception>
#include <pancake/SequenceSeeds.hpp>
#include <type_traits>

namespace PacBio {
namespace Pancake {

SequenceSeeds::SequenceSeeds(std::string name, std::vector<PacBio::Pancake::Int128t> seeds,
                             int64_t id)
    : name_{std::move(name)}, seeds_{std::move(seeds)}, id_{id}
{}

const std::string& SequenceSeeds::Name() const { return name_; }

SequenceSeeds& SequenceSeeds::Name(std::string name)
{
    name_ = std::move(name);
    return *this;
}

const std::vector<PacBio::Pancake::Int128t>& SequenceSeeds::Seeds() const { return seeds_; }

SequenceSeeds& SequenceSeeds::Seeds(std::vector<PacBio::Pancake::Int128t> seeds)
{
    seeds_ = std::move(seeds);
    return *this;
}

int64_t SequenceSeeds::Id() const { return id_; }

SequenceSeeds& SequenceSeeds::Id(int64_t id)
{
    id_ = id;
    return *this;
}

bool SequenceSeeds::operator==(const SequenceSeeds& other) const
{
    return (name_ == other.name_ && id_ == other.id_ && seeds_ == other.seeds_);
}

bool SequenceSeeds::operator!=(const SequenceSeeds& other) const { return !(*this == other); }

}  // namespace Pancake
}  // namespace PacBio
