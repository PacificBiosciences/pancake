// Author: Ivan Sovic

#ifndef PANCAKE_ALIGNER_BATCH_GPU_BASE_H
#define PANCAKE_ALIGNER_BATCH_GPU_BASE_H

#include <pacbio/pancake/AlignmentResult.h>
#include <cstdint>
#include <vector>

namespace PacBio {
namespace Pancake {

enum class StatusAddSequencePair
{
    OK,
    SEQUENCE_LEN_BELOW_ZERO,
    EXCEEDED_MAX_ALIGNMENTS,
    SEQUENCE_TOO_LONG,
};

class AlignerBatchBase
{
public:
    virtual ~AlignerBatchBase() = default;

    /*
     * Clears the internal states (sequences for alignment and results).
    */
    virtual void Clear() = 0;

    /*
     * Adds a single sequence pair for global alignment to the internal state (modifies the state),
     * but does not align the sequences right away. Alignment will be performed only when the
     * AlignAll function is called.
     * The data will be copied internally, so the query and target pointers do not have to
     * remain valid after this function has been called.
     * The return value is either StatusAddSequencePair::OK if the sequences were added successfully,
     * or another value describing the reason for rejecting the sequence pair.
     * Calling this function will clear the internal alignment results.
    */
    virtual StatusAddSequencePair AddSequencePairForGlobalAlignment(const char* query,
                                                                    int32_t queryLen,
                                                                    const char* target,
                                                                    int32_t targetLen) = 0;

    /*
     * Adds a single sequence pair for extension alignment to the internal state (modifies the state),
     * but does not align the sequences right away. Alignment will be performed only when the
     * AlignAll function is called.
     * The data will be copied internally, so the query and target pointers do not have to
     * remain valid after this function has been called.
     * The return value is either StatusAddSequencePair::OK if the sequences were added successfully,
     * or another value describing the reason for rejecting the sequence pair.
     * Calling this function will clear the internal alignment results.
    */
    virtual StatusAddSequencePair AddSequencePairForExtensionAlignment(const char* query,
                                                                       int32_t queryLen,
                                                                       const char* target,
                                                                       int32_t targetLen) = 0;

    /*
     * Aligns all the sequence pairs added to the aligner, in parallel.
     * This function modifies the internal state, because the alignment results are stored internally.
    */
    virtual void AlignAll() = 0;

    /*
     * Const getter for the alignment results. Returns an empty vector if the sequences have not
     * been aligned yet.
    */
    virtual const std::vector<AlignmentResult>& GetAlnResults() const = 0;

    /*
     * Non-const getter for the alignment results. Returns an empty vector if the sequences have not
     * been aligned yet.
    */
    virtual std::vector<AlignmentResult>& GetAlnResults() = 0;

    virtual size_t BatchSize() const = 0;
};
}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_ALIGNER_BATCH_GPU_BASE_H
