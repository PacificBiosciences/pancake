// Author: Ivan Sovic

#ifndef PANCAKE_OVERLAPHIFI_OVERLAP_WRITER_BASE_H
#define PANCAKE_OVERLAPHIFI_OVERLAP_WRITER_BASE_H

#include <pancake/FastaSequenceId.hpp>
#include <pancake/Overlap.hpp>
#include <pancake/SeqDBReaderCached.hpp>
#include <pancake/SeqDBReaderCachedBlock.hpp>
#include <pancake/util/Util.hpp>

#include <cstdint>
#include <memory>
#include <unordered_map>
#include <vector>

namespace PacBio {
namespace Pancake {

class OverlapWriterBase
{
public:
    virtual ~OverlapWriterBase() {}

    virtual void Write(const Overlap& ovl, const PacBio::Pancake::SeqDBReaderCached& targetSeqs,
                       const PacBio::Pancake::FastaSequenceId& querySeq) = 0;

    virtual void Write(const Overlap& ovl,
                       const PacBio::Pancake::SeqDBReaderCachedBlock& targetSeqs,
                       const PacBio::Pancake::FastaSequenceCached& querySeq) = 0;

    virtual void WriteHeader(const PacBio::Pancake::SeqDBReaderCached& targetSeqs) = 0;

    virtual void WriteHeader(const PacBio::Pancake::SeqDBReaderCachedBlock& targetSeqs) = 0;

    static void PrintOverlapAsIPAOvl(FILE* fpOut, const Overlap& ovl, const std::string& Aname,
                                     const std::string& Bname, bool writeIds, bool writeCigar);
    static void PrintOverlapAsM4(FILE* fpOut, const Overlap& ovl, const std::string& Aname,
                                 const std::string& Bname, bool writeIds, bool writeCigar);
    static void PrintOverlapAsPAF(FILE* fpOut, const Overlap& ovl, const std::string& Aname,
                                  const std::string& Bname, bool writeIds, bool writeCigar);
    static void PrintOverlapAsSAM(FILE* fpOut, const Overlap& ovl, const char* seq, int64_t seqLen,
                                  const std::string& Aname, const std::string& Bname, bool writeIds,
                                  bool writeCigar);
    static std::string PrintOverlapAsM4(const Overlap& ovl, const std::string& Aname,
                                        const std::string& Bname, bool writeIds, bool writeCigar);
    static std::string PrintOverlapAsM4(const Overlap& ovl, bool writeCigar = false);
};

constexpr char ConstexprTypeToChar(PacBio::BAM::CigarOperationType type)
{
    constexpr char lookup[11] = "MIDNSHP=XB";
    const int32_t x = static_cast<int>(type);
    return lookup[x];
}

std::ostream& operator<<(std::ostream& os, const Overlap& b);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_OVERLAPHIFI_OVERLAP_WRITER_BASE_H
