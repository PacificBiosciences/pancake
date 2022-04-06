// Authors: Ivan Sovic

#include "SeedDBWorkflow.hpp"
#include "SeedDBSettings.hpp"

#include <pancake/FastaSequenceId.hpp>
#include <pancake/Minimizers.hpp>
#include <pancake/Seed.hpp>
#include <pancake/SeedDBWriter.hpp>
#include <pancake/SeqDBIndexCache.hpp>
#include <pancake/SeqDBReaderCachedBlock.hpp>

#include <pbcopper/parallel/FireAndForget.h>
#include <pbcopper/parallel/WorkQueue.h>

#include <iostream>

namespace PacBio {
namespace Pancake {
namespace SeedDBCLI {

void Worker(const std::vector<FastaSequenceCached>& records, const SeedDBSettings& settings,
            int32_t start, int32_t end, int32_t startAbs,
            std::vector<std::vector<PacBio::Pancake::Int128t>>& seeds)
{
    const auto& sp = settings.SeedParameters;

    for (int32_t i = start; i < end; ++i) {
        const auto& record = records[i];
        const std::string_view seq(record.c_str(), record.size());
        int rv = GenerateMinimizers(seeds[i], seq, 0, record.Id(), sp.KmerSize, sp.MinimizerWindow,
                                    sp.Spacing, sp.UseRC, sp.UseHPCForSeedsOnly);
        if (rv)
            throw std::runtime_error(
                "Generating minimizers failed, startAbs = " + std::to_string(startAbs) +
                ", return code = " + std::to_string(rv));
    }
}
}  // namespace SeedDBCLI

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
            faf.ProduceWith(SeedDBCLI::Worker, std::cref(reader.records()), std::cref(settings), i,
                            i + 1, i + absOffset, std::ref(results));
        }
        faf.Finalize();

        writer->WriteSeeds(reader.records(), results);
        writer->MarkBlockEnd();

        // Increase the abs offset counter.
        absOffset += numRecords;
    }

    return EXIT_SUCCESS;
}

}  // namespace Pancake
}  // namespace PacBio
