// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNMENT_SEEDED_HPP
#define PANCAKE_ALIGNMENT_SEEDED_HPP

#include <pancake/AlignerFactory.hpp>
#include <pancake/AlignmentRegion.hpp>
#include <pancake/DPChain.hpp>
#include <pancake/Overlap.hpp>
#include <pancake/Seed.hpp>

#include <pbcopper/data/Cigar.h>

#include <cstdint>
#include <memory>
#include <string_view>
#include <vector>

namespace PacBio {
namespace Pancake {

class AlignRegionsGenericResult
{
public:
    Data::Cigar cigar;
    int32_t offsetFrontQuery = 0;
    int32_t offsetBackQuery = 0;
    int32_t offsetFrontTarget = 0;
    int32_t offsetBackTarget = 0;
    int32_t score = 0;
};

AlignmentResult AlignSingleRegion(std::string_view targetSeq, std::string_view querySeqFwd,
                                  std::string_view querySeqRev, AlignerBasePtr& alignerGlobal,
                                  AlignerBasePtr& alignerExt, const AlignmentRegion& region);

AlignRegionsGenericResult AlignRegionsGeneric(std::string_view targetSeq,
                                              std::string_view querySeqFwd,
                                              std::string_view querySeqRev,
                                              const std::vector<AlignmentRegion>& regions,
                                              AlignerBasePtr& alignerGlobal,
                                              AlignerBasePtr& alignerExt);

OverlapPtr AlignmentSeeded(const OverlapPtr& ovl, const std::vector<AlignmentRegion>& alnRegions,
                           std::string_view targetSeq, std::string_view querySeqFwd,
                           std::string_view querySeqRev, AlignerBasePtr& alignerGlobal,
                           AlignerBasePtr& alignerExt);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ALIGNMENT_SEEDED_HPP
