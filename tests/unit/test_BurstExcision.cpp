#include <pancake/BurstExcision.hpp>

#include <PancakeTestData.h>

#include <gtest/gtest.h>
#include <pbbam/BamReader.h>
#include <pbbam/BamRecord.h>
#include <pbbam/BamRecordImpl.h>
#include <pbbam/BamRecordTag.h>
#include <pbcopper/data/Frames.h>
#include <pbcopper/utility/Ssize.h>
#include <pancake/FastaSequenceCachedStore.hpp>
#include <pancake/FastaSequenceId.hpp>
#include <pancake/MapperCLR.hpp>

#include <memory>
#include <string>
#include <utility>

namespace {

struct Subread
{
    std::string Name;
    std::string Bases;
    PacBio::Data::Frames IPD;
    PacBio::Data::Frames PW;

    Subread(std::string name, std::string bases, PacBio::Data::Frames ipd, PacBio::Data::Frames pw)
        : Name(std::move(name)), Bases(std::move(bases)), IPD(std::move(ipd)), PW(std::move(pw))
    {}
};

std::vector<Subread> ParseSubreads(const std::string& path)
{
    std::vector<Subread> dst;

    PacBio::BAM::BamReader reader{path};
    for (const PacBio::BAM::BamRecord& record : reader) {
        PacBio::BAM::BamRecordImpl recordImpl{record.Impl()};
        dst.emplace_back(
            recordImpl.Name(), recordImpl.Sequence(),
            PacBio::Data::Frames::Decode(
                recordImpl.TagValue(PacBio::BAM::BamRecordTag::IPD).ToUInt8Array()),
            PacBio::Data::Frames::Decode(
                recordImpl.TagValue(PacBio::BAM::BamRecordTag::PULSE_WIDTH).ToUInt8Array()));
    }

    return dst;
}

PacBio::Pancake::FastaSequenceCachedStore CreateSubreadViews(const std::vector<Subread>& zmw)
{
    PacBio::Pancake::FastaSequenceCachedStore dst;

    int32_t id = 0;
    for (const auto& it : zmw) {
        dst.AddRecord(PacBio::Pancake::FastaSequenceCached{
            it.Name, it.Bases.c_str(), PacBio::Utility::Ssize(it.Bases), id++, &it.IPD, &it.PW});
    }

    return dst;
}

}  // namespace

TEST(BurstExcision, EmptyZmw)
{
    EXPECT_TRUE(PacBio::Pancake::BurstExcision(PacBio::Pancake::FastaSequenceCachedStore{},
                                               PacBio::Pancake::MapperCLRSettings{})
                    .empty());
}

TEST(BurstExcision, ZmwWithMissingKinetics)
{
    const std::string testData = PacBio::PancakeTestsConfig::Data_Dir +
                                 "/bursts/test-2-zmw-with-one-subread-having-two-bursts.bam";

    const std::vector<Subread> zmw = ParseSubreads(testData);

    PacBio::Pancake::FastaSequenceCachedStore zmwView = CreateSubreadViews(zmw);
    zmwView.records()[3].IPD(nullptr);

    const PacBio::Pancake::MapperCLRSettings settings{};

    const auto burstRegionsPerSubread = PacBio::Pancake::BurstExcision(zmwView, settings);

    EXPECT_EQ(std::size(zmw), std::size(burstRegionsPerSubread));

    for (const auto& burstRegions : burstRegionsPerSubread) {
        EXPECT_TRUE(burstRegions.empty());
    }
}

TEST(BurstExcision, ZmwWithoutBursts)
{
    const std::string testData =
        PacBio::PancakeTestsConfig::Data_Dir + "/bursts/test-1-zmw-without-bursts.bam";

    const std::vector<Subread> zmw = ParseSubreads(testData);

    const PacBio::Pancake::FastaSequenceCachedStore zmwView = CreateSubreadViews(zmw);

    const PacBio::Pancake::MapperCLRSettings settings{};

    const auto burstRegionsPerSubread =
        PacBio::Pancake::BurstExcision(zmwView, settings, 32, 0.666, 0.666, 0.666, 0.666, 0.666);

    EXPECT_EQ(std::size(zmw), std::size(burstRegionsPerSubread));

    for (const auto& burstRegions : burstRegionsPerSubread) {
        EXPECT_TRUE(burstRegions.empty());
    }
}

TEST(BurstExcision, ZmwWithTwoBursts)
{
    const std::string testData = PacBio::PancakeTestsConfig::Data_Dir +
                                 "/bursts/test-2-zmw-with-one-subread-having-two-bursts.bam";

    std::vector<Subread> zmw = ParseSubreads(testData);
    std::swap(zmw[0], zmw[2]);  // mimic CCS draft stage

    const PacBio::Pancake::FastaSequenceCachedStore zmwView = CreateSubreadViews(zmw);

    PacBio::Pancake::MapperCLRSettings settings{};
    {  // mimic CCS draft stage
        settings.map.seedParams.KmerSize = 15;
        settings.map.seedParams.MinimizerWindow = 5;
        settings.map.seedParams.Spacing = 0;
        settings.map.seedParams.UseHPCForSeedsOnly = true;
        settings.map.freqPercentile = 0.000;
        settings.map.seedOccurrenceMin = 10;
        settings.map.seedOccurrenceMax = 1000;
        settings.map.seedOccurrenceMaxMemory = 100'000'000;
        settings.map.secondaryAllowedOverlapFractionTarget = 0.05;
    }

    const auto burstRegionsPerSubread =
        PacBio::Pancake::BurstExcision(zmwView, settings, 32, 0.666, 0.666, 0.666, 0.666, 0.666);

    EXPECT_EQ(std::size(zmw), std::size(burstRegionsPerSubread));

    for (int32_t i = 0; i < 7; ++i) {
        EXPECT_TRUE(burstRegionsPerSubread[i].empty());
    }

    EXPECT_EQ(2ULL, std::size(burstRegionsPerSubread[7]));
    EXPECT_EQ(std::make_pair(666, 944), burstRegionsPerSubread[7][0]);
    EXPECT_EQ(std::make_pair(3818, 3902), burstRegionsPerSubread[7][1]);
}
