// Author: Ivan Sovic

#include <pancake/FastaSequenceId.hpp>

#include <boost/algorithm/string.hpp>

#include <cassert>
#include <cstdio>
#include <exception>
#include <tuple>
#include <type_traits>

namespace PacBio {
namespace Pancake {

static_assert(std::is_copy_constructible<FastaSequenceId>::value,
              "Sequence(const Sequence&) is not = default");
static_assert(std::is_copy_assignable<FastaSequenceId>::value,
              "Sequence& operator=(const Sequence&) is not = default");

static_assert(std::is_nothrow_move_constructible<FastaSequenceId>::value,
              "Sequence(Sequence&&) is not = noexcept");
static_assert(std::is_nothrow_move_assignable<FastaSequenceId>::value ==
                  std::is_nothrow_move_assignable<std::string>::value,
              "");

FastaSequenceId::FastaSequenceId(std::string name, std::string bases, int64_t id)
    : BAM::FastaSequence(name, bases), id_{id}
{}

FastaSequenceId::FastaSequenceId(const BAM::FastaSequence& fastaSequence, int64_t id)
    : BAM::FastaSequence(fastaSequence), id_{id}
{}

int64_t FastaSequenceId::Id() const { return id_; }

FastaSequenceId& FastaSequenceId::Id(int64_t id)
{
    id_ = id;
    return *this;
}

bool FastaSequenceId::operator==(const FastaSequenceId& other) const
{
    return std::tie(Name(), id_, Bases()) == std::tie(other.Name(), id_, other.Bases());
}

bool FastaSequenceId::operator!=(const FastaSequenceId& other) const { return !(*this == other); }

}  // namespace Pancake
}  // namespace PacBio
