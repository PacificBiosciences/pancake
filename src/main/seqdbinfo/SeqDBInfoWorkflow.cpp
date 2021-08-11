// Authors: Ivan Sovic

#include "SeqDBInfoWorkflow.h"
#include "SeqDBInfoSettings.h"

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

#include <pacbio/pancake/SeqDBIndexCache.h>
#include <pacbio/util/SeqLengthStats.h>
#include <pbcopper/logging/LogLevel.h>
#include <pbcopper/logging/Logging.h>

namespace PacBio {
namespace Pancake {

int SeqDBInfoWorkflow::Runner(const PacBio::CLI_v2::Results& options)
{
    SeqDBInfoSettings settings{options};

    // Load the SeqDB index cache.
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(settings.InputSeqDB);

    // Collect and sort the length info.
    std::vector<int32_t> lengths(seqDBCache->seqLines.size(), 0);
    int64_t totalBytes = 0;
    for (size_t i = 0; i < seqDBCache->seqLines.size(); ++i) {
        const SeqDBSequenceLine& sl = seqDBCache->seqLines[i];
        lengths[i] = sl.numBases;
        totalBytes += sl.numBytes;
    }
    std::sort(lengths.begin(), lengths.end(), [](const auto& a, const auto& b) { return a > b; });

    // Compute stats.
    SeqLengthStats stats;
    ComputeSeqLengthStats(lengths, 0, stats);

    // Convert to the desired unit.
    stats.ChangeUnitOfLengths(settings.Unit);

    // Write the output.
    if (settings.HumanReadableOutput) {
        // clang-format off
        fprintf(stdout, "input\tunit\ttotal\tnum\tmin\tmax\tavg\tmedian\tAUC\tN10\tN10_n\tN25\tN25_n\tN50\tN50_n\tN75\tN75_n\tN90\tN90_n\tN100\tN100_n\n");

        fprintf(stdout, "%s\t%s\t%.2lf\t%d\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf"
                        "\t%.2lf\t%d\t%.2lf\t%d\t%.2lf\t%d\t%.2lf\t%d\t%.2lf\t%d\t%.2lf\t%d\n",
                        settings.InputSeqDB.c_str(),
                        GenomicUnitToString(stats.unit).c_str(),
                        stats.totalLength, stats.numSeqs, stats.lenMin,
                        stats.lenMax, stats.lenAvg, stats.lenMedian, stats.nxAUC,
                        std::get<1>(stats.Nx[10]), std::get<2>(stats.Nx[10]),
                        std::get<1>(stats.Nx[25]), std::get<2>(stats.Nx[25]),
                        std::get<1>(stats.Nx[50]), std::get<2>(stats.Nx[50]),
                        std::get<1>(stats.Nx[75]), std::get<2>(stats.Nx[75]),
                        std::get<1>(stats.Nx[90]), std::get<2>(stats.Nx[90]),
                        std::get<1>(stats.Nx[100]), std::get<2>(stats.Nx[100])
        );
        // clang-format on

    } else {
        PacBio::JSON::Json statsJson = SeqLengthStatsToJson(stats);
        statsJson["total_bytes"] = totalBytes;
        // statsJson["input"] = settings.InputSeqDB;
        std::cout << statsJson.dump(4);
    }

    return EXIT_SUCCESS;
}

}  // namespace Pancake
}  // namespace PacBio
