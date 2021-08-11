// Author: Ivan Sovic

#ifndef PANCAKE_UTIL_SEQ_LENGTH_STATS_H
#define PANCAKE_UTIL_SEQ_LENGTH_STATS_H

#include <pacbio/util/GenomicUnit.h>
#include <pbcopper/json/JSON.h>
#include <memory>
#include <ostream>
#include <vector>

namespace PacBio {
namespace Pancake {

class SeqLengthStats
{
public:
    double totalLength = 0;
    int32_t numSeqs = 0;
    double lenMin = 0;
    double lenMax = 0;
    double lenAvg = 0.0;
    double lenMedian = 0.0;
    double nxAUC = 0.0;
    // Tuple: Nx percent, Nx value, number of sequences.
    std::array<std::tuple<int32_t, double, int32_t>, 101> Nx;
    GenomicUnit unit = GenomicUnit::bp;

    SeqLengthStats();

    void Clear();
    void ChangeUnitOfLengths(GenomicUnit unit);

    bool operator==(const SeqLengthStats& b) const
    {
        return (totalLength == b.totalLength && numSeqs == b.numSeqs && lenMin == b.lenMin &&
                lenMax == b.lenMax && lenAvg == b.lenAvg && lenMedian == b.lenMedian &&
                nxAUC == b.nxAUC && Nx == b.Nx && unit == b.unit);
    }

private:
    void ScaleLengthsByFactor_(double factor);
};

inline std::ostream& operator<<(std::ostream& os, const SeqLengthStats& a)
{
    os << "unit = " << GenomicUnitToString(a.unit) << ", length = " << a.totalLength
       << ", num = " << a.numSeqs << ", min = " << a.lenMin << ", max = " << a.lenMax
       << ", avg = " << a.lenAvg << ", median = " << a.lenMedian << ", AUC = " << a.nxAUC << "\n";
    for (size_t i = 0; i < a.Nx.size(); ++i) {
        os << "  [x = " << std::get<0>(a.Nx[i]) << "] Nx = " << std::get<1>(a.Nx[i])
           << ", Nx_n = " << std::get<2>(a.Nx[i]) << "\n";
    }
    return os;
}

void ComputeSeqLengthStats(const std::vector<int32_t>& reverseSortedLengths, int64_t genomeSize,
                           SeqLengthStats& ret);

PacBio::JSON::Json SeqLengthStatsToJson(const SeqLengthStats& stats);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_UTIL_SEQ_LENGTH_STATS_H