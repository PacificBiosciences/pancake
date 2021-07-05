// Authors: Ivan Sovic

#include "SeqDBDumpWorkflow.h"
#include "SeqDBDumpSettings.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include <pacbio/pancake/SeqDBIndexCache.h>
#include <pacbio/pancake/SeqDBReader.h>
#include <pacbio/pancake/SeqDBReaderCachedBlock.h>
#include <pacbio/util/FileIO.h>
#include <pbcopper/logging/LogLevel.h>
#include <pbcopper/logging/Logging.h>

namespace PacBio {
namespace Pancake {

void WriteSeq(FILE* fp, const char* name, const size_t nameLen, const char* seq,
              const size_t seqLen)
{
    fprintf(fp, ">");
    fwrite(name, sizeof(char), nameLen, fp);
    fprintf(fp, "\n");
    fwrite(seq, sizeof(char), seqLen, fp);
    fprintf(fp, "\n");
}

int SeqDBDumpWorkflow::Runner(const PacBio::CLI_v2::Results& options)
{
    SeqDBDumpSettings settings{options};

    // Output file/stdout.
    FILE* fpOut = stdout;
    if (settings.OutputFile.size() > 0 && settings.OutputFile != "-") {
        fpOut = fopen(settings.OutputFile.c_str(), "w");
        if (fpOut == NULL) {
            throw std::runtime_error("Could not open file '" + settings.OutputFile +
                                     "' for writing!");
        }
        PBLOG_INFO << "Output is to file: " << settings.OutputFile;
    } else {
        PBLOG_INFO << "Output is to stdout.";
    }

    PBLOG_INFO << "Loading the SeqDB.";
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(settings.InputSeqDB);

    // Sanity check.
    const int32_t numBlocks = seqDBCache->blockLines.size();
    if (settings.BlockId >= numBlocks) {
        throw std::runtime_error("Specified block ID is too large, numBlocks = " +
                                 std::to_string(numBlocks) + ".");
    }

    // Open up a reader.
    PacBio::Pancake::SeqDBReaderCachedBlock reader(seqDBCache, settings.UseHPC);

    // Write the sequences.
    PBLOG_INFO << "Fetching the data.";
    const int32_t startBlockId = std::max(settings.BlockId, 0);
    const int32_t endBlockId = (settings.BlockId >= 0) ? (settings.BlockId + 1) : numBlocks;
    char nameIdBuffer[50];
    for (int32_t blockId = startBlockId; blockId < endBlockId; ++blockId) {
        reader.LoadBlocks({blockId});
        const FastaSequenceCachedStore& recordStore = reader.recordStore();
        const std::vector<FastaSequenceCached>& records = recordStore.records();

        for (const auto& record : records) {
            if (settings.WriteIds) {
                sprintf(nameIdBuffer, "%09d", record.Id());
                WriteSeq(fpOut, nameIdBuffer, strlen(nameIdBuffer), record.c_str(), record.size());
            } else {
                WriteSeq(fpOut, record.Name().c_str(), record.Name().size(), record.c_str(),
                         record.size());
            }
        }
    }

    if (fpOut && fpOut != stdout) {
        fclose(fpOut);
    }

    PBLOG_INFO << "Done!";

    return EXIT_SUCCESS;
}

}  // namespace Pancake
}  // namespace PacBio
