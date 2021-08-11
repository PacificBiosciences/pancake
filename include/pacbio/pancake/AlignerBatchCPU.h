// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNER_BATCH_CPU_H
#define PANCAKE_ALIGNER_BATCH_CPU_H

#include <pacbio/pancake/AlignerBase.h>
#include <pacbio/pancake/AlignerBatchBase.h>
#include <pacbio/pancake/AlignerFactory.h>
#include <pacbio/pancake/Range.h>
#include <pbbam/Cigar.h>
#include <pbcopper/parallel/FireAndForget.h>
#include <array>
#include <cstdint>
#include <memory>
#include <vector>

namespace PacBio {
namespace Pancake {

class AlignerBatchCPU : public AlignerBatchBase
{
public:
    AlignerBatchCPU(int32_t numThreads, const AlignerType alnTypeGlobal,
                    const AlignmentParameters& alnParamsGlobal, const AlignerType alnTypeExt,
                    const AlignmentParameters& alnParamsExt);
    AlignerBatchCPU(Parallel::FireAndForget* faf, const AlignerType alnTypeGlobal,
                    const AlignmentParameters& alnParamsGlobal, const AlignerType alnTypeExt,
                    const AlignmentParameters& alnParamsExt);
    ~AlignerBatchCPU() override;

    void Clear() override;

    StatusAddSequencePair AddSequencePairForGlobalAlignment(const char* query, int32_t queryLen,
                                                            const char* target,
                                                            int32_t targetLen) override;

    StatusAddSequencePair AddSequencePairForExtensionAlignment(const char* query, int32_t queryLen,
                                                               const char* target,
                                                               int32_t targetLen) override;

    void AlignAll() override;

    const std::vector<AlignmentResult>& GetAlnResults() const override { return alnResults_; }

    std::vector<AlignmentResult>& GetAlnResults() override { return alnResults_; }

    size_t BatchSize() const override { return queries_.size(); }

private:
    Parallel::FireAndForget* faf_;
    std::unique_ptr<Parallel::FireAndForget> fafFallback_;

    AlignerType alnTypeGlobal_;
    AlignmentParameters alnParamsGlobal_;
    AlignerType alnTypeExt_;
    AlignmentParameters alnParamsExt_;

    std::vector<std::string> queries_;
    std::vector<std::string> targets_;
    std::vector<bool> isGlobalAlignment_;
    std::vector<AlignmentResult> alnResults_;

    StatusAddSequencePair AddSequencePair_(const char* query, int32_t queryLen, const char* target,
                                           int32_t targetLen, bool isGlobalAlignment);

    static void Worker_(const std::vector<std::string>& queries,
                        const std::vector<std::string>& targets,
                        const std::vector<bool>& isGlobalAlignment, int32_t alnStartId,
                        int32_t alnEndId, AlignerBasePtr& alignerGlobal, AlignerBasePtr& alignerExt,
                        std::vector<AlignmentResult>& alnResults);
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ALIGNER_BATCH_CPU_H
