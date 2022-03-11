// Author: Ivan Sovic

#ifndef PANCAKE_CONTIGUOUS_FILE_PART_H
#define PANCAKE_CONTIGUOUS_FILE_PART_H

#include <cstdint>
#include <stdexcept>
#include <vector>

namespace PacBio {
namespace Pancake {

class ContiguousFilePart
{
public:
    int32_t fileId = 0;
    int64_t startOffset = 0;
    int64_t endOffset = 0;
    std::vector<int32_t> seqIds;

    bool operator==(const ContiguousFilePart& b) const
    {
        return fileId == b.fileId && startOffset == b.startOffset && endOffset == b.endOffset &&
               seqIds == b.seqIds;
    }

    bool CanAppendTo(const ContiguousFilePart& b) const
    {
        if (seqIds.empty() || b.seqIds.empty()) {
            throw std::runtime_error(
                "Cannot run CanAppendTo because one or both of the ContiguousFilePart objects has "
                "an empty seqIds vector");
        }
        if ((fileId == b.fileId) && (startOffset == b.endOffset) &&
            (seqIds.front() == b.seqIds.back()))
            return true;
        return false;
    }

    void ExtendWith(const ContiguousFilePart& b)
    {
        if (b.CanAppendTo(*this) == false) {
            throw std::runtime_error(
                "Attempted to extend a ContiguousFilePart with another part, but they are not "
                "sequential.");
        }
        endOffset = b.endOffset;
        seqIds.insert(seqIds.end(), b.seqIds.begin(), b.seqIds.end());
    }
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_CONTIGUOUS_FILE_PART_H