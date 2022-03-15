// Author: Ivan Sovic

#ifndef PANCAKE_SEED_DB_PARAMETERS_HPP
#define PANCAKE_SEED_DB_PARAMETERS_HPP

#include <cstdint>

namespace PacBio {
namespace Pancake {

// clang-format off
class SeedDBParameters
{
public:
    int32_t KmerSize = 28;
    int32_t MinimizerWindow = 80;
    int32_t Spacing = 0;
    bool UseHPC = false;                // This causes the input sequences from the DB to be HP-compressed.
    bool UseHPCForSeedsOnly = false;    // This takes the uncompressed sequences, and just skips HP bases when computing seeds.
    bool UseRC = true;

    SeedDBParameters() = default;
    SeedDBParameters(const int32_t kmerSize, const int32_t minimizerWindow, const int32_t spacing, const bool useHPC, const bool useHPCForSeedsOnly, const bool useRC)
        : KmerSize(kmerSize)
        , MinimizerWindow(minimizerWindow)
        , Spacing(spacing)
        , UseHPC(useHPC)
        , UseHPCForSeedsOnly(useHPCForSeedsOnly)
        , UseRC(useRC)
    {
    }
    ~SeedDBParameters() = default;

    bool operator==(const SeedDBParameters& rhs) const
    {
        return KmerSize == rhs.KmerSize && MinimizerWindow == rhs.MinimizerWindow &&
            Spacing == rhs.Spacing && UseHPC == rhs.UseHPC &&
            UseHPCForSeedsOnly == rhs.UseHPCForSeedsOnly &&
            UseRC == rhs.UseRC;
    }
    bool operator!=(const SeedDBParameters& rhs) const
    {
        return !((*this) == rhs);
    }
};
// clang-format on

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_SEED_DB_PARAMETERS_HPP
