// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNER_KSW2_H
#define PANCAKE_ALIGNER_KSW2_H

#include <pacbio/pancake/AlignerBase.h>
#include <cstdint>
#include <memory>
#include <vector>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include <pbcopper/third-party/ksw2/ksw2.h>
#pragma GCC diagnostic pop

namespace PacBio {
namespace Pancake {

struct mm_tbuf_s
{
    void* km;
    int rep_len, frag_gap;
};
// memory buffer for thread-local storage during mapping
typedef struct mm_tbuf_s mm_tbuf_t;
typedef std::unique_ptr<mm_tbuf_t, void (*)(mm_tbuf_t*)> Minimap2ThreadBufferPtr;

class AlignerKSW2;
std::shared_ptr<AlignerBase> CreateAlignerKSW2(const AlignmentParameters& opt);

class AlignerKSW2 : public AlignerBase
{
public:
    AlignerKSW2(const AlignmentParameters& opt);
    ~AlignerKSW2() override;

    AlignmentResult Global(const char* qseq, int64_t qlen, const char* tseq, int64_t tlen) override;
    AlignmentResult Extend(const char* qseq, int64_t qlen, const char* tseq, int64_t tlen) override;

private:
    AlignmentParameters opt_;
    Minimap2ThreadBufferPtr buffer_;
    int8_t mat_[25];

    void ConvertMinimap2CigarToPbbam_(uint32_t* mm2Cigar, int32_t cigarLen,
                                      const std::vector<uint8_t>& qseq,
                                      const std::vector<uint8_t>& tseq,
                                      PacBio::Data::Cigar& retCigar, int32_t& retQueryAlignmentLen,
                                      int32_t& retTargetAlignmentLen,
                                      Alignment::DiffCounts& retDiffs);

    static std::vector<uint8_t> ConvertSeqAlphabet_(const char* seq, size_t seqlen,
                                                    const int8_t* conv_table);

    static void GenerateSimpleMatrix_(int m, int8_t* mat, int8_t a, int8_t b, int8_t scAmbi);

    static void AlignPair_(void* km, int qlen, const uint8_t* qseq, int tlen, const uint8_t* tseq,
                           const int8_t* mat, int w, int endBonus, int zdrop, int flag,
                           ksw_extz_t* ez, int q, int e, int q2, int e2);
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ALIGNER_KSW2_H
