// Authors: Ivan Sovic

#include <pacbio/util/SeqLengthStats.h>

namespace PacBio {
namespace Pancake {

SeqLengthStats::SeqLengthStats()
{
    // Initialize the Nx array to all zeros.
    Clear();
}

void SeqLengthStats::Clear()
{
    totalLength = numSeqs = lenMin = lenMax = 0;
    lenAvg = lenMedian = nxAUC = 0.0;
    unit = GenomicUnit::bp;
    for (size_t i = 0; i < Nx.size(); ++i) {
        Nx[i] = std::make_tuple(i, 0, 0);
    }
}

void SeqLengthStats::ChangeUnitOfLengths(GenomicUnit unitNew)
{
    const auto converter = GenomicUnitFromTo(unit, unitNew);
    ScaleLengthsByFactor_(converter.conversionFactor());
    unit = unitNew;
}

void SeqLengthStats::ScaleLengthsByFactor_(double factor)
{
    totalLength *= factor;
    lenMin *= factor;
    lenMax *= factor;
    lenAvg *= factor;
    lenMedian *= factor;
    nxAUC *= factor;
    for (size_t i = 0; i < Nx.size(); ++i) {
        std::get<1>(Nx[i]) *= factor;
    }
}

void ComputeSeqLengthStats(const std::vector<int32_t>& reverseSortedLengths, int64_t genomeSize,
                           SeqLengthStats& ret)
{
    ret.Clear();
    if (reverseSortedLengths.empty()) {
        return;
    }

    // Compute the total length and Nx AUC.
    ret.totalLength = 0;
    ret.nxAUC = 0.0;
    int32_t prevVal = reverseSortedLengths.front();
    for (const int32_t& val : reverseSortedLengths) {
        if (val > prevVal) {
            throw std::runtime_error(
                "Input vector 'reverseSortedLengths' is not reverse sorted!\n");
        }
        ret.totalLength += val;
        ret.nxAUC += val * val;
    }
    ret.nxAUC = (ret.totalLength == 0) ? 0.0 : ret.nxAUC / ret.totalLength;

    // Compute other basic stats.
    ret.unit = GenomicUnit::bp;
    ret.numSeqs = reverseSortedLengths.size();
    ret.lenMin = reverseSortedLengths.back();
    ret.lenMax = reverseSortedLengths.front();
    ret.lenAvg = (ret.numSeqs == 0) ? 0 : ret.totalLength / static_cast<double>(ret.numSeqs);
    ret.lenMedian = (ret.numSeqs == 0) ? 0
                                       : (reverseSortedLengths[ret.numSeqs / 2] +
                                          reverseSortedLengths[(ret.numSeqs - 1) / 2]) /
                                             2;

    // Compute the Nx graph for 1% increase in intervals.
    const int64_t nxGenomeSize = (genomeSize > 0) ? genomeSize : ret.totalLength;
    const double percStep = 0.01 * nxGenomeSize;
    double nextThreshold = 0;
    double subtotal = 0;
    int32_t x = 0;
    int32_t lastLen = reverseSortedLengths[0];
    int32_t lastLenId = 0;
    for (int32_t lenId = 0; lenId < static_cast<int32_t>(reverseSortedLengths.size()); ++lenId) {
        subtotal += reverseSortedLengths[lenId];
        while (subtotal >= nextThreshold && x < 101) {
            lastLen = reverseSortedLengths[lenId];
            lastLenId = lenId;
            ret.Nx[x] = std::make_tuple(x, lastLen, lastLenId + 1);
            nextThreshold += percStep;
            ++x;
        }
    }
    // Fill in any remaining bins.
    for (; x < static_cast<int32_t>(ret.Nx.size()); ++x) {
        ret.Nx[x] = std::make_tuple(x, lastLen, lastLenId + 1);
    }
}

PacBio::JSON::Json SeqLengthStatsToJson(const SeqLengthStats& stats)
{
    auto result = PacBio::JSON::Json::object();
    result["total_length"] = stats.totalLength;
    result["num_seqs"] = stats.numSeqs;
    result["len_min"] = stats.lenMin;
    result["len_max"] = stats.lenMax;
    result["len_avg"] = stats.lenAvg;
    result["len_median"] = stats.lenMedian;
    result["AUC"] = stats.nxAUC;
    result["Nx"] = stats.Nx;
    result["unit"] = GenomicUnitToString(stats.unit);
    return result;
}

}  // namespace Pancake
}  // namespace PacBio
