// Author: Ivan Sovic

#include <fstream>
#include <pancake/SeedHitWriter.hpp>
#include <sstream>

namespace PacBio {
namespace Pancake {

void WriteSeedHits(const std::string& outPath, const std::vector<SeedHit>& hits, size_t hitsStart,
                   size_t hitsEnd, int32_t hitsId, const std::string& queryName,
                   int64_t queryLength, const std::string& targetName, int64_t targetLength,
                   bool append)
{
    if (hitsStart >= hits.size() || hitsEnd > hits.size() || hitsStart > hitsEnd) {
        std::ostringstream oss;
        oss << "Invalid hitsStart and/or hitsEnd! hitsStart = " << hitsStart
            << ", hitsEnd = " << hitsEnd << ", hits.size() = " << hits.size();
        throw std::runtime_error(oss.str());
    }

    std::ofstream ofs;
    if (append) {
        ofs = std::ofstream(outPath, std::ios::app);
    } else {
        ofs = std::ofstream(outPath);
    }
    if (ofs.is_open() == false) {
        // Don't throw. This is a hidden feature which will only work if a user knows which folder to create.
        return;
        // throw std::runtime_error("Could not open file '" + outPath +
        //                          "' for writing! In MapperCLR::WriteSeedHits.");
    }
    if (append == false) {
        ofs << queryName.c_str() << "\t" << queryLength << "\t" << targetName.c_str() << "\t"
            << targetLength << "\n";
    }
    for (size_t j = hitsStart; j < hitsEnd; ++j) {
        ofs << hits[j].queryPos << "\t" << hits[j].targetPos << "\t" << hitsId << "\n";
    }
}
}  // namespace Pancake
}  // namespace PacBio