// Authors: Ivan Sovic

#include "SeedDBWorkflow.h"
#include <pacbio/pancake/FastaSequenceId.h>
#include <pacbio/pancake/Minimizers.h>
#include <pacbio/pancake/Seed.h>
#include <pacbio/pancake/SeedDBWriter.h>
#include <pacbio/pancake/SeqDBIndexCache.h>
#include <pacbio/pancake/SeqDBReaderCachedBlock.h>
#include "SeedDBSettings.h"

#include <pbcopper/parallel/FireAndForget.h>
#include <pbcopper/parallel/WorkQueue.h>

#include <iostream>

namespace PacBio {
namespace Pancake {
namespace SeedDB {

void Worker(const std::vector<FastaSequenceCached>& records, const SeedDBSettings& settings,
            int32_t start, int32_t end, int32_t startAbs,
            std::vector<std::vector<PacBio::Pancake::Int128t>>& seeds)
{
    const auto& sp = settings.SeedParameters;

    for (int32_t i = start; i < end; ++i) {
        const auto& record = records[i];
        const uint8_t* seq = reinterpret_cast<const uint8_t*>(record.c_str());
        int32_t seqLen = record.size();
        int rv =
            GenerateMinimizers(seeds[i], seq, seqLen, 0, record.Id(), sp.KmerSize,
                               sp.MinimizerWindow, sp.Spacing, sp.UseRC, sp.UseHPCForSeedsOnly);
        if (rv)
            throw std::runtime_error("Generating minimizers failed, startAbs = " +
                                     std::to_string(startAbs) + ", return code = " +
                                     std::to_string(rv));
    }
}

int SeedDBWorkflow::Runner(const PacBio::CLI_v2::Results& options)
{
    SeedDBSettings settings{options};

    // Load the DB.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(settings.InputFile);

    // Create a reader.
    PacBio::Pancake::SeqDBReaderCachedBlock reader(seqDBCache, settings.SeedParameters.UseHPC);

    int32_t numBlocks = seqDBCache->blockLines.size();
    int32_t absOffset = 0;

    auto writer = PacBio::Pancake::CreateSeedDBWriter(settings.OutputPrefix, settings.SplitBlocks,
                                                      settings.SeedParameters);

    for (int32_t blockId = 0; blockId < numBlocks; ++blockId) {
        // Load a block of records.
        reader.LoadBlocks({blockId});
        int32_t numRecords = reader.records().size();

        // Generate seeds in parallel.
        std::vector<std::vector<PacBio::Pancake::Int128t>> results(numRecords);
        PacBio::Parallel::FireAndForget faf(settings.NumThreads);
        for (int32_t i = 0; i < numRecords; ++i) {
            faf.ProduceWith(Worker, std::cref(reader.records()), std::cref(settings), i, i + 1,
                            i + absOffset, std::ref(results));
        }
        faf.Finalize();

        writer->WriteSeeds(reader.records(), results);
        writer->MarkBlockEnd();

        // Increase the abs offset counter.
        absOffset += numRecords;
    }

    return EXIT_SUCCESS;
}

}  // namespace SeedDB
}  // namespace Pancake
}  // namespace PacBio
