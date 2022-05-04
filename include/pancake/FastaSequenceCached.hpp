// Author: Ivan Sovic

#ifndef PANCAKE_FASTA_SEQUENCE_CACHED_HPP
#define PANCAKE_FASTA_SEQUENCE_CACHED_HPP

#include <pbcopper/data/Frames.h>

#include <cstdint>
#include <string>
#include <tuple>

namespace PacBio {
namespace Pancake {

/*
 * This is a container which does not store the actual data, it only points to
 * data existing somewhere else.
 * The FastaSequenceCached contained does not deallocate the memory!
*/
class FastaSequenceCached
{
public:
    FastaSequenceCached() : name_(), bases_(NULL), size_(0), id_(-1) {}

    FastaSequenceCached(std::string _name, const char* _bases, int64_t _size, int32_t _id,
                        const Data::Frames* _ipd = nullptr, const Data::Frames* _pw = nullptr)
        : name_(std::move(_name)), bases_(_bases), size_(_size), id_(_id), ipd_(_ipd), pw_(_pw)
    {}

    // Getters.
    const std::string& Name() const { return name_; }
    const char* Bases() const { return bases_; }
    int64_t Size() const { return size_; }
    int32_t Id() const { return id_; }

    // Setters.
    void Name(const std::string& val) { name_ = val; }
    void Bases(const char* val) { bases_ = val; }
    void Size(int64_t val) { size_ = val; }
    void Id(int32_t val) { id_ = val; }
    void IPD(const Data::Frames* ipd) { ipd_ = ipd; }
    void PW(const Data::Frames* pw) { pw_ = pw; }

    // Getters to make it compatible with the std::string.
    const char* data() const { return bases_; }
    const char* c_str() const { return bases_; }
    int64_t size() const { return size_; };

    const Data::Frames* IPD() const { return ipd_; }
    const Data::Frames* PW() const { return pw_; }

private:
    std::string name_;
    const char* bases_;
    int64_t size_;
    int32_t id_;

    const Data::Frames* ipd_;
    const Data::Frames* pw_;
};

/*
 * Note: the operator== here does not compare the actual sequence content, just the pointers.
 * Since this is a "view" class, it doesn't care about the actual data, just that it points to the
 * correct location.
*/
inline bool operator==(const FastaSequenceCached& lhs, const FastaSequenceCached& rhs) noexcept
{
    return std::tuple(lhs.Name(), lhs.Bases(), lhs.Size(), lhs.Id(), lhs.IPD(), lhs.PW()) ==
           std::tuple(rhs.Name(), rhs.Bases(), rhs.Size(), rhs.Id(), rhs.IPD(), rhs.PW());
}

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_FASTA_SEQUENCE_CACHED_HPP
