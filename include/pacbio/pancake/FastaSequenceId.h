// File Description
/// \file Sequence.h
/// \brief Defines the Sequence class. Inherits from BAM::FastaSequence, but
//         differs from BAM::Sequence in that it contains the absolute ID
///        of the sequence as well.
//
// Author: Derek Barnett and Ivan Sovic.

#ifndef PANCAKE_FASTASEQUENCE_H
#define PANCAKE_FASTASEQUENCE_H

#include <pbbam/FastaSequence.h>
#include <string>

namespace PacBio {
namespace Pancake {

///
/// \brief The FastaSequenceId class represents a FASTA record (name & bases)
///
class FastaSequenceId : public BAM::FastaSequence
{
public:
    /// \name Constructors & Related Methods
    /// \{

    ///
    /// \brief Sequence
    /// \param name
    /// \param id
    /// \param bases
    ///
    explicit FastaSequenceId(std::string name, std::string bases, int64_t id);

    ///
    /// \brief Sequence
    /// \param fastaSequence
    /// \param id
    ///
    explicit FastaSequenceId(const BAM::FastaSequence& fastaSequence, int64_t id);

    FastaSequenceId() = default;

    /// \}

public:
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
    FastaSequenceId& Id(int64_t id);

    /// \}

    bool operator==(const FastaSequenceId& other) const;
    bool operator!=(const FastaSequenceId& other) const;

private:
    int64_t id_;
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // FASTASEQUENCE_H
