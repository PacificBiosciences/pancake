// Authors: Ivan Sovic

#include "DBFilterWorkflow.h"

#include <algorithm>
#include <iostream>
#include <random>
#include <string>
#include <unordered_set>
#include <vector>

#include <pbcopper/logging/LogLevel.h>
#include <pbcopper/logging/Logging.h>

#include <boost/algorithm/string/predicate.hpp>

#include <pancake/SeedDBIndexCache.hpp>
#include <pancake/SeqDBIndexCache.hpp>
#include <pancake/SeqDBReader.hpp>
#include <pancake/SeqDBWriter.hpp>
#include <pancake/util/TicToc.hpp>
#include <pancake/util/Util.hpp>

#include "DBFilterSettings.h"

namespace PacBio {
namespace Pancake {

std::unordered_set<std::string> ParseFilterList(const std::string& filterListPath)
{
    std::unordered_set<std::string> filterList;
    std::ifstream ifs(filterListPath);
    if (ifs.is_open() == false) {
        throw std::runtime_error("Could not open blacklist file '" + filterListPath + "'!");
    }
    std::string line;
    while (std::getline(ifs, line)) {
        filterList.emplace(line);
    }
    return filterList;
}

int DBFilterWorkflow::Runner(const PacBio::CLI_v2::Results& options)
{
    DBFilterSettings settings{options};

    std::string inSeqDBFile = settings.InputPrefix + ".seqdb";
    // std::string inSeedDBFile = settings.InputPrefix + ".seeddb";
    std::string outSeqDBFile = settings.OutputPrefix + ".seqdb";
    // std::string outSeedDBFile = settings.OutputPrefix + ".seeddb";

    PBLOG_INFO << "Loading the input DBs.";
    TicToc ttLoad;

    // Load the target SeqDB caches.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> inSeqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inSeqDBFile);
    PBLOG_INFO << "After loading target seq cache: " << ttLoad.VerboseSecs(true);

    // // Load the target SeedDB caches.
    // std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> inSeedDBCache =
    //     PacBio::Pancake::LoadSeedDBIndexCache(inSeedDBFile);
    // PBLOG_INFO << "After loading target seed cache: " << ttLoad.VerboseSecs(true);
    ttLoad.Stop();

    std::unordered_set<std::string> filterList;
    if (settings.FilterListPath.empty() == false) {
        filterList = ParseFilterList(settings.FilterListPath);
    }
    PBLOG_INFO << "Filter list size: " << filterList.size();

    // Construct a filtered SeqDB cache.
    // Do not normalize yet, we might need the original sequence IDs for fetching.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> filteredSeqDBCache = FilterSeqDBIndexCache(
        *inSeqDBCache, settings.Sampling, settings.SampleBases, settings.RandomSeed,
        settings.FilterType, filterList, false, 0, outSeqDBFile);

    PBLOG_INFO << "Original SeqDB sequences: " << inSeqDBCache->seqLines.size();
    PBLOG_INFO << "Filtered SeqDB sequences: " << filteredSeqDBCache->seqLines.size();

    // If the Consolidate option is not specified, then we only need to modify the
    // index cache. Otherwise, we need to create a copy of the data using the Writer.
    if (settings.Consolidate) {
        // Create a reader so we can fetch sequences.
        PacBio::Pancake::SeqDBReader reader(inSeqDBCache);

        // Create a new writer, which uses the actual output prefix.
        auto writer = PacBio::Pancake::CreateSeqDBWriter(
            settings.OutputPrefix, settings.CompressionLevel, settings.BufferSize,
            settings.BlockSize, settings.SplitBlocks);

        // Fetch and copy all the sequences into the new DB.
        Pancake::FastaSequenceId record;
        for (const auto& sl : filteredSeqDBCache->seqLines) {
            reader.GetSequence(record, sl.seqId);
            writer->AddSequence(record.Name(), record.Bases());
        }

    } else {
        NormalizeSeqDBIndexCache(*filteredSeqDBCache, settings.BlockSize);
        std::unique_ptr<FILE, FileDeleter> fpOutSeqDBCache =
            PacBio::Pancake::OpenFile(outSeqDBFile.c_str(), "w");
        WriteSeqDBIndexCache(fpOutSeqDBCache.get(), *filteredSeqDBCache);
    }

    PBLOG_INFO << "Done writing of the filtered DB.";

    return EXIT_SUCCESS;
}

}  // namespace Pancake
}  // namespace PacBio
