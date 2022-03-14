// Author: Ivan Sovic

#ifndef PANCAKE_SEED_HIT_WRITER_H
#define PANCAKE_SEED_HIT_WRITER_H

#include <pancake/SeedHit.hpp>

#include <string>
#include <vector>

namespace PacBio {
namespace Pancake {

void WriteSeedHits(const std::string& outPath, const std::vector<SeedHit>& hits, size_t hitsStart,
                   size_t hitsEnd, int32_t hitsId, const std::string& queryName,
                   int64_t queryLength, const std::string& targetName, int64_t targetLength,
                   bool append);
}
}  // namespace PacBio

#endif  // PANCAKE_SEED_HIT_WRITER_H