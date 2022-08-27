// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNER_BATCH_CPU_HPP
#define PANCAKE_ALIGNER_BATCH_CPU_HPP

#include <pancake/AlignerBase.hpp>
#include <pancake/AlignerBatchBase.hpp>
#include <pancake/AlignerFactory.hpp>
#include <pancake/Range.hpp>

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
    AlignerBatchCPU(int32_t numThreads, AlignerType alnTypeGlobal,
                    const AlignmentParameters& alnParamsGlobal, AlignerType alnTypeExt,
                    const AlignmentParameters& alnParamsExt);
    AlignerBatchCPU(Parallel::FireAndForget* faf, AlignerType alnTypeGlobal,
                    const AlignmentParameters& alnParamsGlobal, AlignerType alnTypeExt,
                    const AlignmentParameters& alnParamsExt);
    ~AlignerBatchCPU() override;

    void Clear() override;

    StatusAddSequencePair AddSequencePairForGlobalAlignment(std::string_view qseq,
                                                            std::string_view tseq) override;

    StatusAddSequencePair AddSequencePairForExtensionAlignment(std::string_view qseq,
                                                               std::string_view tseq) override;

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

    StatusAddSequencePair AddSequencePair_(std::string_view qseq, std::string_view tseq,
                                           bool isGlobalAlignment);

    static void Worker_(const std::vector<std::string>& queries,
                        const std::vector<std::string>& targets,
                        const std::vector<bool>& isGlobalAlignment, int32_t alnStartId,
                        int32_t alnEndId, AlignerBasePtr& alignerGlobal, AlignerBasePtr& alignerExt,
                        std::vector<AlignmentResult>& alnResults);
};

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ALIGNER_BATCH_CPU_HPP
