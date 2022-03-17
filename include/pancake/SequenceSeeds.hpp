// File Description
/// \file Sequence.h
/// \brief Defines the SequenceSeeds class. It contains all seeds for a given sequence,
//         tied with the sequence ID and length.
//
// Author: Ivan Sovic

#ifndef PANCAKE_SEQUENCE_SEEDS_HPP
#define PANCAKE_SEQUENCE_SEEDS_HPP

#include <pancake/util/CommonTypes.hpp>

#include <cstdint>
#include <string>
#include <vector>

namespace PacBio {
namespace Pancake {

///
/// \brief The SequenceSeeds class hold seeds for a given sequence.
///
class SequenceSeeds
{
public:
    /// \name Constructors & Related Methods
    /// \{

    ///
    /// \brief Sequence
    /// \param name
    /// \param seeds
    /// \param id
    ///
    explicit SequenceSeeds(std::string name, std::vector<PacBio::Pancake::Int128t> seeds,
                           int64_t id);

    SequenceSeeds() = default;

    /// \}

public:
    /// \name Attributes
    /// \{

    ///
    /// \brief Name
    /// \return
    ///
    const std::string& Name() const;

    ///
    /// \brief
    ///
    /// \param name
    /// \return SequenceSeeds&
    ///
    SequenceSeeds& Name(std::string name);

    ///
    /// \brief Seeds
    /// \return
    ///
    const std::vector<PacBio::Pancake::Int128t>& Seeds() const;

    ///
    /// \brief
    ///
    /// \param seeds
    /// \return SequenceSeeds&
    ///
    SequenceSeeds& Seeds(std::vector<PacBio::Pancake::Int128t> seeds);

    ///
    /// \brief
    ///
    /// \param id
    /// \return int64_t
    ///
    int64_t Id() const;

    ///
    /// \brief
    ///
    /// \param id
    /// \return int64_t
    ///
    SequenceSeeds& Id(int64_t id);

    /// \}

    bool operator==(const SequenceSeeds& other) const;
    bool operator!=(const SequenceSeeds& other) const;

private:
    std::string name_;
    std::vector<PacBio::Pancake::Int128t> seeds_;
    int64_t id_;
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEQUENCE_SEEDS_HPP
