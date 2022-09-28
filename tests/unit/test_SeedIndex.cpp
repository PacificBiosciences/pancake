// Authors: Ivan Sovic

#include <PancakeTestData.h>
#include <gtest/gtest.h>
#include <pancake/Minimizers.hpp>
#include <pancake/Seed.hpp>
#include <pancake/SeedIndex.hpp>
#include <pancake/util/CommonTypes.hpp>
#include <sstream>
#include <tuple>

using namespace PacBio::Pancake;

TEST(SeedIndex, GetSeeds1)
{
    /*
     * Fetch the byte span of a single small block of 1 sequence.
     * It doesn't span more than 1 file or sequences.
    */

    const int32_t seqId = 123;
    const int32_t span = 15;
    std::vector<PacBio::Pancake::Int128t> inSeeds = {
        PacBio::Pancake::Seed::Encode(0, span, seqId, 0, false),
        PacBio::Pancake::Seed::Encode(0, span, seqId, 1, false),
        PacBio::Pancake::Seed::Encode(0, span, seqId, 2, false),
        PacBio::Pancake::Seed::Encode(0, span, seqId, 3, false),
        PacBio::Pancake::Seed::Encode(0, span, seqId, 4, false),
        PacBio::Pancake::Seed::Encode(0, span, seqId, 5, false),
    };
    uint64_t inKey = 0;

    // Expected results.
    std::vector<PacBio::Pancake::Int128t> expected = inSeeds;

    // Run the unit under test.
    PacBio::Pancake::SeedIndex si(std::move(inSeeds));
    std::vector<PacBio::Pancake::Int128t> results;
    int64_t n = si.GetSeeds(inKey, results);

    EXPECT_EQ(expected, results);
    EXPECT_EQ(static_cast<int64_t>(expected.size()), n);
}

TEST(SeedIndex, GetSeeds2)
{
    /*
     * Fetch the byte span of a single small block of 1 sequence.
     * It doesn't span more than 1 file or sequences.
    */

    const int32_t seqId = 123;
    const int32_t span = 28;
    std::vector<PacBio::Pancake::Int128t> inSeeds = {
        PacBio::Pancake::Seed::Encode(0, span, seqId, 0, false),
        PacBio::Pancake::Seed::Encode(123, span, seqId, 1, false),
        PacBio::Pancake::Seed::Encode(5, span, seqId, 2, false),
        PacBio::Pancake::Seed::Encode(7, span, seqId, 3, false),
        PacBio::Pancake::Seed::Encode(5, span, seqId, 4, false),
        PacBio::Pancake::Seed::Encode(0, span, seqId, 5, false),
        PacBio::Pancake::Seed::Encode(123, span, seqId, 6, false),
        PacBio::Pancake::Seed::Encode(5, span, seqId, 7, false),
        PacBio::Pancake::Seed::Encode(123, span, seqId, 8, false),
    };
    uint64_t inKey = 5;

    // Expected results.
    std::vector<PacBio::Pancake::Int128t> expected = {
        PacBio::Pancake::Seed::Encode(5, span, seqId, 2, false),
        PacBio::Pancake::Seed::Encode(5, span, seqId, 4, false),
        PacBio::Pancake::Seed::Encode(5, span, seqId, 7, false),
    };

    // Run the unit under test.
    PacBio::Pancake::SeedIndex si(std::move(inSeeds));
    std::vector<PacBio::Pancake::Int128t> results;
    int64_t n = si.GetSeeds(inKey, results);

    EXPECT_EQ(expected, results);
    EXPECT_EQ(static_cast<int64_t>(expected.size()), n);
}

TEST(SeedIndex, ComputeSeedSpanStats_EmptyInput)
{
    std::vector<PacBio::Pancake::Int128t> inSeeds = {};

    // Expected results.
    const int32_t expectedMinSeedSpan = 0;
    const int32_t expectedMaxSeedSpan = 0;
    const double expectedAvgSeedSpan = 0.0;

    // Run the unit under test.
    PacBio::Pancake::SeedIndex si(std::move(inSeeds));

    // Evaluate.
    EXPECT_EQ(expectedMinSeedSpan, si.GetMinSeedSpan());
    EXPECT_EQ(expectedMaxSeedSpan, si.GetMaxSeedSpan());
    EXPECT_EQ(expectedAvgSeedSpan, si.GetAvgSeedSpan());
}

TEST(SeedIndex, ComputeSeedSpanStats_SameSpan)
{
    const int32_t seqId = 123;
    const int32_t span = 15;
    std::vector<PacBio::Pancake::Int128t> inSeeds = {
        PacBio::Pancake::Seed::Encode(0, span, seqId, 0, false),
        PacBio::Pancake::Seed::Encode(0, span, seqId, 1, false),
        PacBio::Pancake::Seed::Encode(0, span, seqId, 2, false),
        PacBio::Pancake::Seed::Encode(0, span, seqId, 3, false),
        PacBio::Pancake::Seed::Encode(0, span, seqId, 4, false),
        PacBio::Pancake::Seed::Encode(0, span, seqId, 5, false),
    };

    // Expected results.
    const int32_t expectedMinSeedSpan = span;
    const int32_t expectedMaxSeedSpan = span;
    const double expectedAvgSeedSpan = span;

    // Run the unit under test.
    PacBio::Pancake::SeedIndex si(std::move(inSeeds));

    // Evaluate.
    EXPECT_EQ(expectedMinSeedSpan, si.GetMinSeedSpan());
    EXPECT_EQ(expectedMaxSeedSpan, si.GetMaxSeedSpan());
    EXPECT_EQ(expectedAvgSeedSpan, si.GetAvgSeedSpan());
}

TEST(SeedIndex, ComputeSeedSpanStats_MixOfDifferentSpans)
{
    const int32_t seqId = 123;
    const int32_t span1 = 15;
    const int32_t span2 = 19;
    const int32_t span3 = 28;
    std::vector<PacBio::Pancake::Int128t> inSeeds = {
        PacBio::Pancake::Seed::Encode(0, span1, seqId, 0, false),
        PacBio::Pancake::Seed::Encode(0, span1, seqId, 1, false),
        PacBio::Pancake::Seed::Encode(0, span3, seqId, 2, false),
        PacBio::Pancake::Seed::Encode(0, span2, seqId, 3, false),
        PacBio::Pancake::Seed::Encode(0, span3, seqId, 4, false),
        PacBio::Pancake::Seed::Encode(0, span1, seqId, 5, false),
    };

    // Expected results.
    const int32_t expectedMinSeedSpan = span1;
    const int32_t expectedMaxSeedSpan = span3;
    const double expectedAvgSeedSpan =
        static_cast<double>(span1 + span1 + span3 + span2 + span3 + span1) /
        static_cast<double>(inSeeds.size());

    // Run the unit under test.
    PacBio::Pancake::SeedIndex si(std::move(inSeeds));

    // Evaluate.
    EXPECT_EQ(expectedMinSeedSpan, si.GetMinSeedSpan());
    EXPECT_EQ(expectedMaxSeedSpan, si.GetMaxSeedSpan());
    EXPECT_EQ(expectedAvgSeedSpan, si.GetAvgSeedSpan());
}

TEST(SeedIndex, GetSeeds3NonexistentKey)
{
    /*
     * Fetch the byte span of a single small block of 1 sequence.
     * It doesn't span more than 1 file or sequences.
    */

    const int32_t seqId = 123;
    const int32_t span = 19;
    std::vector<PacBio::Pancake::Int128t> inSeeds = {
        PacBio::Pancake::Seed::Encode(0, span, seqId, 0, false),
        PacBio::Pancake::Seed::Encode(123, span, seqId, 1, false),
        PacBio::Pancake::Seed::Encode(5, span, seqId, 2, false),
        PacBio::Pancake::Seed::Encode(7, span, seqId, 3, false),
        PacBio::Pancake::Seed::Encode(5, span, seqId, 4, false),
        PacBio::Pancake::Seed::Encode(0, span, seqId, 5, false),
        PacBio::Pancake::Seed::Encode(123, span, seqId, 6, false),
        PacBio::Pancake::Seed::Encode(5, span, seqId, 7, false),
        PacBio::Pancake::Seed::Encode(123, span, seqId, 8, false),
    };
    uint64_t inKey = 1024;

    // Expected results.
    std::vector<PacBio::Pancake::Int128t> expected = {};

    // Run the unit under test.
    PacBio::Pancake::SeedIndex si(std::move(inSeeds));
    std::vector<PacBio::Pancake::Int128t> results;
    int64_t n = si.GetSeeds(inKey, results);

    EXPECT_EQ(expected, results);
    EXPECT_EQ(static_cast<int64_t>(expected.size()), n);
}

TEST(SeedIndex, CollectHitsEmptyQueryNonemptyTarget)
{
    /*
     * Tests collecting hits for a given vector of query keys and a maximum
     * allowed seed frequency.
     *
     * This particular test builds the SeedIndex with a set of seeds from
     * an imaginary target.
     * It uses the same set of seeds for the input to CollectHits.
     * For every input seed, we should get all other positions in the target.
     *
     * Here, the target seeds contains actual seeds.
     * Query seed set is empty.
     * There should be no hits.
    */
    const int32_t targetId = 0;
    const int32_t span = 17;
    std::vector<PacBio::Pancake::SeedRaw> targetSeeds = {
        PacBio::Pancake::Seed::Encode(0, span, targetId, 0, false),
        PacBio::Pancake::Seed::Encode(123, span, targetId, 1, false),
        PacBio::Pancake::Seed::Encode(5, span, targetId, 2, false),
        PacBio::Pancake::Seed::Encode(7, span, targetId, 3, false),
        PacBio::Pancake::Seed::Encode(5, span, targetId, 4, false),
        PacBio::Pancake::Seed::Encode(0, span, targetId, 5, false),
        PacBio::Pancake::Seed::Encode(123, span, targetId, 6, false),
        PacBio::Pancake::Seed::Encode(5, span, targetId, 7, false),
        PacBio::Pancake::Seed::Encode(123, span, targetId, 8, false),
    };
    const std::vector<PacBio::Pancake::SeedRaw> querySeeds = {};
    const int32_t queryLen = 38;

    // Expected results.
    const std::vector<PacBio::Pancake::SeedHit> expected = {};

    // Run the unit under test.
    PacBio::Pancake::SeedIndex si(std::move(targetSeeds));
    std::vector<PacBio::Pancake::SeedHit> results;
    si.CollectHits(querySeeds, queryLen, results, 100);

    // Evaluate.
    EXPECT_EQ(expected, results);
}

TEST(SeedIndex, CollectHitsNonemptyQueryEmptyTarget)
{
    /*
     * Tests collecting hits for a given vector of query keys and a maximum
     * allowed seed frequency.
     *
     * This particular test builds the SeedIndex with a set of seeds from
     * an imaginary target.
     * It uses the same set of seeds for the input to CollectHits.
     * For every input seed, we should get all other positions in the target.
     *
     * Target seed set is empty.
     * Query seed set contains seeds which should not match anything.
     * There should be no hits.
    */
    std::vector<PacBio::Pancake::SeedRaw> targetSeeds = {};
    const int32_t queryId = 0;
    const int32_t span = 15;
    const std::vector<PacBio::Pancake::SeedRaw> querySeeds = {
        PacBio::Pancake::Seed::Encode(0, span, queryId, 0, false),
        PacBio::Pancake::Seed::Encode(123, span, queryId, 1, false),
        PacBio::Pancake::Seed::Encode(5, span, queryId, 2, false),
        PacBio::Pancake::Seed::Encode(7, span, queryId, 3, false),
        PacBio::Pancake::Seed::Encode(5, span, queryId, 4, false),
        PacBio::Pancake::Seed::Encode(0, span, queryId, 5, false),
        PacBio::Pancake::Seed::Encode(123, span, queryId, 6, false),
        PacBio::Pancake::Seed::Encode(5, span, queryId, 7, false),
        PacBio::Pancake::Seed::Encode(123, span, queryId, 8, false),
    };
    const int32_t queryLen = 38;

    // Expected results.
    const std::vector<PacBio::Pancake::SeedHit> expected = {};

    // Run the unit under test.
    PacBio::Pancake::SeedIndex si(std::move(targetSeeds));
    std::vector<PacBio::Pancake::SeedHit> results;
    si.CollectHits(querySeeds, queryLen, results, 100);

    // Evaluate.
    EXPECT_EQ(expected, results);
}

TEST(SeedIndex, CollectHitsNonemptyQueryNonemptyTargetNoHits)
{
    /*
     * Tests collecting hits for a given vector of query keys and a maximum
     * allowed seed frequency.
     *
     * This particular test builds the SeedIndex with a set of seeds from
     * an imaginary target.
     * It uses the same set of seeds for the input to CollectHits.
     * For every input seed, we should get all other positions in the target.
     *
     * Target seed set not empty.
     * Query seed set contains seeds which should not match anything in the target set.
     * There should be no hits.
    */
    const int32_t span = 15;
    const int32_t targetId = 0;
    std::vector<PacBio::Pancake::SeedRaw> targetSeeds = {
        PacBio::Pancake::Seed::Encode(0, span, targetId, 0, false),
        PacBio::Pancake::Seed::Encode(7, span, targetId, 3, false),
        PacBio::Pancake::Seed::Encode(0, span, targetId, 5, false),
    };
    const int32_t queryId = 0;
    const std::vector<PacBio::Pancake::SeedRaw> querySeeds = {
        PacBio::Pancake::Seed::Encode(123, span, queryId, 1, false),
        PacBio::Pancake::Seed::Encode(5, span, queryId, 2, false),
        PacBio::Pancake::Seed::Encode(5, span, queryId, 4, false),
        PacBio::Pancake::Seed::Encode(5, span, queryId, 7, false),
    };
    const int32_t queryLen = 38;

    // Expected results.
    const std::vector<PacBio::Pancake::SeedHit> expected = {};

    // Run the unit under test.
    PacBio::Pancake::SeedIndex si(std::move(targetSeeds));
    std::vector<PacBio::Pancake::SeedHit> results;
    si.CollectHits(querySeeds, queryLen, results, 100);

    // Evaluate.
    EXPECT_EQ(expected, results);
}

TEST(SeedIndex, CollectHitsPerfectMatch)
{
    /*
     * Tests collecting hits for a given vector of query keys and a maximum
     * allowed seed frequency.
     *
     * This particular test builds the SeedIndex with a set of seeds from
     * an imaginary target.
     * It uses the same set of seeds for the input to CollectHits.
     * For every input seed, we should get all other positions in the target.
    */
    const int32_t span = 15;
    const int32_t targetId = 0;
    std::vector<PacBio::Pancake::SeedRaw> targetSeeds = {
        PacBio::Pancake::Seed::Encode(0, span, targetId, 0, false),
        PacBio::Pancake::Seed::Encode(123, span, targetId, 1, false),
        PacBio::Pancake::Seed::Encode(5, span, targetId, 2, false),
        PacBio::Pancake::Seed::Encode(7, span, targetId, 3, false),
        PacBio::Pancake::Seed::Encode(5, span, targetId, 4, false),
        PacBio::Pancake::Seed::Encode(0, span, targetId, 5, false),
        PacBio::Pancake::Seed::Encode(123, span, targetId, 6, false),
        PacBio::Pancake::Seed::Encode(5, span, targetId, 7, false),
        PacBio::Pancake::Seed::Encode(123, span, targetId, 8, false),
    };
    const std::vector<PacBio::Pancake::SeedRaw> querySeeds = targetSeeds;
    const int32_t queryLen = 38;

    // Expected results.
    const std::vector<PacBio::Pancake::SeedHit> expected = {
        // 0
        {targetId, false, 0, 0, span, span, 0},
        {targetId, false, 5, 0, span, span, 0},
        // 123
        {targetId, false, 1, 1, span, span, 0},
        {targetId, false, 6, 1, span, span, 0},
        {targetId, false, 8, 1, span, span, 0},
        /// 5
        {targetId, false, 2, 2, span, span, 0},
        {targetId, false, 4, 2, span, span, 0},
        {targetId, false, 7, 2, span, span, 0},
        /// 7
        {targetId, false, 3, 3, span, span, 0},
        // 5
        {targetId, false, 2, 4, span, span, 0},
        {targetId, false, 4, 4, span, span, 0},
        {targetId, false, 7, 4, span, span, 0},
        // 0
        {targetId, false, 0, 5, span, span, 0},
        {targetId, false, 5, 5, span, span, 0},
        // 123
        {targetId, false, 1, 6, span, span, 0},
        {targetId, false, 6, 6, span, span, 0},
        {targetId, false, 8, 6, span, span, 0},
        // 5
        {targetId, false, 2, 7, span, span, 0},
        {targetId, false, 4, 7, span, span, 0},
        {targetId, false, 7, 7, span, span, 0},
        // 123
        {targetId, false, 1, 8, span, span, 0},
        {targetId, false, 6, 8, span, span, 0},
        {targetId, false, 8, 8, span, span, 0},
    };

    // Run the unit under test.
    PacBio::Pancake::SeedIndex si(std::move(targetSeeds));
    std::vector<PacBio::Pancake::SeedHit> results;
    si.CollectHits(querySeeds, queryLen, results, 100);

    // Evaluate.
    EXPECT_EQ(expected, results);
}

TEST(SeedIndex, CollectHitsFrequencyThreshold)
{
    /*
     * Tests collecting hits for a given vector of query keys and a maximum
     * allowed seed frequency.
     *
     * Same as before, but test frequency cutoff.
     * Do not fetch hith with frequency > 2.
    */
    const int32_t span = 15;
    const int32_t targetId = 0;
    std::vector<PacBio::Pancake::SeedRaw> targetSeeds = {
        PacBio::Pancake::Seed::Encode(0, span, targetId, 0, false),
        PacBio::Pancake::Seed::Encode(123, span, targetId, 1, false),
        PacBio::Pancake::Seed::Encode(5, span, targetId, 2, false),
        PacBio::Pancake::Seed::Encode(7, span, targetId, 3, false),
        PacBio::Pancake::Seed::Encode(5, span, targetId, 4, false),
        PacBio::Pancake::Seed::Encode(0, span, targetId, 5, false),
        PacBio::Pancake::Seed::Encode(123, span, targetId, 6, false),
        PacBio::Pancake::Seed::Encode(5, span, targetId, 7, false),
        PacBio::Pancake::Seed::Encode(123, span, targetId, 8, false),
    };
    const std::vector<PacBio::Pancake::SeedRaw> querySeeds = targetSeeds;
    const int32_t queryLen = 38;
    int32_t freqCutoff = 2;

    // Expected results.
    const std::vector<PacBio::Pancake::SeedHit> expected = {
        // 0
        {targetId, false, 0, 0, span, span, 0},
        {targetId, false, 5, 0, span, span, 0},
        /// 7
        {targetId, false, 3, 3, span, span, 0},
        // 0
        {targetId, false, 0, 5, span, span, 0},
        {targetId, false, 5, 5, span, span, 0},
    };

    // Run the unit under test.
    PacBio::Pancake::SeedIndex si(std::move(targetSeeds));
    std::vector<PacBio::Pancake::SeedHit> results;
    si.CollectHits(querySeeds, queryLen, results, freqCutoff);

    // Evaluate.
    EXPECT_EQ(expected, results);
}

TEST(SeedIndex, CollectHitsReverseStrand)
{
    /*
     * Tests collecting hits for a given vector of query keys and a maximum
     * allowed seed frequency.
     *
     * Some of the hits are on the reverse strand of the query.
    */
    const int32_t k = 30;

    const int32_t targetId = 0;
    std::vector<PacBio::Pancake::SeedRaw> targetSeeds = {
        PacBio::Pancake::Seed::Encode(0, k, targetId, 0, false),
        PacBio::Pancake::Seed::Encode(123, k, targetId, 1, false),
        PacBio::Pancake::Seed::Encode(5, k, targetId, 2, false),
        PacBio::Pancake::Seed::Encode(7, k, targetId, 3, false),
        PacBio::Pancake::Seed::Encode(5, k, targetId, 4, true),
        PacBio::Pancake::Seed::Encode(0, k, targetId, 5, false),
        PacBio::Pancake::Seed::Encode(123, k, targetId, 6, false),
        PacBio::Pancake::Seed::Encode(5, k, targetId, 7, false),
        PacBio::Pancake::Seed::Encode(123, k, targetId, 8, false),
    };
    const int32_t queryId = 0;
    const std::vector<PacBio::Pancake::SeedRaw> querySeeds = {
        PacBio::Pancake::Seed::Encode(0, k, queryId, 0, false),
        PacBio::Pancake::Seed::Encode(123, k, queryId, 1, true),
        PacBio::Pancake::Seed::Encode(5, k, queryId, 2, true),
        PacBio::Pancake::Seed::Encode(7, k, queryId, 3, true),
        PacBio::Pancake::Seed::Encode(5, k, queryId, 4, false),
        PacBio::Pancake::Seed::Encode(0, k, queryId, 5, false),
        PacBio::Pancake::Seed::Encode(123, k, queryId, 6, false),
        PacBio::Pancake::Seed::Encode(5, k, queryId, 7, false),
        PacBio::Pancake::Seed::Encode(123, k, queryId, 8, false),
    };
    const int32_t queryLen = 38;
    int32_t freqCutoff = 0;

    // Expected results.
    const std::vector<PacBio::Pancake::SeedHit> expected = {

        // 0
        {targetId, false, 0, 0, 30, 30, 0},
        {targetId, false, 5, 0, 30, 30, 0},

        // 123
        {targetId, true, 1, 7, 30, 30, 0},
        {targetId, true, 6, 7, 30, 30, 0},
        {targetId, true, 8, 7, 30, 30, 0},

        /// 5
        {targetId, true, 2, 6, 30, 30, 0},
        {targetId, true, 7, 6, 30, 30, 0},
        {targetId, false, 4, 2, 30, 30, 0},

        /// 7
        {targetId, true, 3, 5, 30, 30, 0},

        // 5
        {targetId, false, 2, 4, 30, 30, 0},
        {targetId, false, 7, 4, 30, 30, 0},
        {targetId, true, 4, 4, 30, 30, 0},

        // 0
        {targetId, false, 0, 5, 30, 30, 0},
        {targetId, false, 5, 5, 30, 30, 0},

        // 123
        {targetId, false, 1, 6, 30, 30, 0},
        {targetId, false, 6, 6, 30, 30, 0},
        {targetId, false, 8, 6, 30, 30, 0},

        // 5
        {targetId, false, 2, 7, 30, 30, 0},
        {targetId, false, 7, 7, 30, 30, 0},
        {targetId, true, 4, 1, 30, 30, 0},

        // 123
        {targetId, false, 1, 8, 30, 30, 0},
        {targetId, false, 6, 8, 30, 30, 0},
        {targetId, false, 8, 8, 30, 30, 0},
    };

    // Run the unit under test.
    PacBio::Pancake::SeedIndex si(std::move(targetSeeds));
    std::vector<PacBio::Pancake::SeedHit> results;
    si.CollectHits(querySeeds, queryLen, results, freqCutoff);

    for (size_t i = 0; i < results.size(); ++i) {
        const auto& hit = results[i];
        // std::cerr << "[i = " << i << "] targetId = " << hit.targetId
        //           << ", targetRev = " << hit.targetRev << ", targetPos = " << hit.targetPos
        //           << ", flags = " << hit.flags << ", queryPos = " << hit.queryPos << "\n";
        std::cerr << "{targetId"
                  << ", " << (hit.targetRev ? "true" : "false") << ", " << hit.targetPos << ", "
                  << hit.queryPos << ", " << hit.targetSpan << ", " << hit.querySpan << ", "
                  << hit.flags << "},\n";
    }

    EXPECT_EQ(expected, results);
}

TEST(SeedIndex, CollectHitsLongSeedSpan)
{
    const int32_t targetId = 0;
    std::vector<PacBio::Pancake::SeedRaw> targetSeeds = {
        PacBio::Pancake::Seed::Encode(0, 164, targetId, 0, 0),
        PacBio::Pancake::Seed::Encode(1, 15, targetId, 1, 0),
        PacBio::Pancake::Seed::Encode(2, 128, targetId, 2, 0),
        PacBio::Pancake::Seed::Encode(3, 255, targetId, 3, 0),
    };
    const int32_t queryId = 0;
    const std::vector<PacBio::Pancake::SeedRaw> querySeeds = {
        // PacBio::Pancake::Seed::Encode(0, 164, queryId, 0, 0),
        // PacBio::Pancake::Seed::Encode(1, 15, queryId, 153, 0),
        // PacBio::Pancake::Seed::Encode(1, 15, queryId, 155, 0),
        // PacBio::Pancake::Seed::Encode(1, 15, queryId, 157, 0),

        PacBio::Pancake::Seed::Encode(0, 164, queryId, 0, 0),
        PacBio::Pancake::Seed::Encode(1, 15, queryId, 1, 0),
        PacBio::Pancake::Seed::Encode(2, 128, queryId, 2, 0),
        PacBio::Pancake::Seed::Encode(3, 255, queryId, 3, 0),
    };
    const int32_t queryLen = 173;

    // clang-format off
    // Expected results.
    const std::vector<PacBio::Pancake::SeedHit> expected = {
        {targetId, false, 0, 0, 164, 164, 0},
        {targetId, false, 1, 1, 15, 15, 0},
        {targetId, false, 2, 2, 128, 128, 0},
        {targetId, false, 3, 3, 255, 255, 0},
    };
    // clang-format on

    // Run the unit under test.
    PacBio::Pancake::SeedIndex si(std::move(targetSeeds));
    std::vector<PacBio::Pancake::SeedHit> results;
    si.CollectHits(querySeeds, queryLen, results, 100);

    for (size_t i = 0; i < results.size(); ++i) {
        const auto& hit = results[i];
        std::cerr << "{targetId"
                  << ", " << (hit.targetRev ? "true" : "false") << ", " << hit.targetPos << ", "
                  << hit.queryPos << ", " << hit.targetSpan << ", " << hit.querySpan << ", "
                  << hit.flags << "},\n";
    }

    // Evaluate.
    // This test is not very useful to check for positive/negative spans because
    // the SeedHit constructor will convert those positive values >=128 to a signed 8-bit form.
    EXPECT_EQ(expected, results);

    // This actually captures that something went wrong.
    for (const auto& result : results) {
        std::ostringstream oss;
        oss << result;
        SCOPED_TRACE("Testing hit: " + oss.str());
        EXPECT_TRUE(static_cast<int32_t>(result.querySpan) >= 0);
    }
}

TEST(SeedIndex, ComputeFrequencyStatsEmptyIndex)
{
    /*
     * Initialize a SeedIndex, and compute the seed statistics.
    */
    // Input values.
    std::vector<PacBio::Pancake::Int128t> inSeeds = {};
    double freqPercentile = 0.60;  // Large percentile, just because this is a small test.

    // Expected results.
    int64_t expectedFreqMax = 0;
    double expectedFreqAvg = 0;
    double expectedFreqMedian = 0;
    int64_t expectedFreqCutoff = 0;

    // Run the unit under test.
    PacBio::Pancake::SeedIndex si(std::move(inSeeds));
    int64_t resultFreqMax = 0;
    double resultFreqAvg = 0.0;
    double resultFreqMedian = 0.0;
    int64_t resultFreqCutoff = 0;
    si.ComputeFrequencyStats(freqPercentile, resultFreqMax, resultFreqAvg, resultFreqMedian,
                             resultFreqCutoff);

    // Evaluate.
    EXPECT_EQ(
        std::make_tuple(expectedFreqMax, expectedFreqAvg, expectedFreqMedian, expectedFreqCutoff),
        std::make_tuple(resultFreqMax, resultFreqAvg, resultFreqMedian, resultFreqCutoff));
}

TEST(SeedIndex, ComputeFrequencyStatsNormal)
{
    /*
     * Initialize a SeedIndex, and compute the seed statistics.
    */
    // Input values.
    const int32_t seqId = 123;
    const int32_t span = 15;
    std::vector<PacBio::Pancake::Int128t> inSeeds = {
        PacBio::Pancake::Seed::Encode(0, span, seqId, 0, false),
        PacBio::Pancake::Seed::Encode(123, span, seqId, 1, false),
        PacBio::Pancake::Seed::Encode(5, span, seqId, 2, false),
        PacBio::Pancake::Seed::Encode(7, span, seqId, 3, false),
        PacBio::Pancake::Seed::Encode(5, span, seqId, 4, false),
        PacBio::Pancake::Seed::Encode(0, span, seqId, 5, false),
        PacBio::Pancake::Seed::Encode(123, span, seqId, 6, false),
        PacBio::Pancake::Seed::Encode(5, span, seqId, 7, false),
        PacBio::Pancake::Seed::Encode(123, span, seqId, 8, false),
    };
    const double freqPercentile = 0.60;  // Large percentile, just because this is a small test.

    // Expected results.
    const int64_t expectedFreqMax = 3;
    const double expectedFreqAvg = 2.25;
    const double expectedFreqMedian = 2.5;
    const int64_t expectedFreqCutoff = 2;

    // Run the unit under test.
    const PacBio::Pancake::SeedIndex si(std::move(inSeeds));
    int64_t resultFreqMax = 0;
    double resultFreqAvg = 0.0;
    double resultFreqMedian = 0.0;
    int64_t resultFreqCutoff = 0;
    si.ComputeFrequencyStats(freqPercentile, resultFreqMax, resultFreqAvg, resultFreqMedian,
                             resultFreqCutoff);

    // Evaluate.
    EXPECT_EQ(
        std::make_tuple(expectedFreqMax, expectedFreqAvg, expectedFreqMedian, expectedFreqCutoff),
        std::make_tuple(resultFreqMax, resultFreqAvg, resultFreqMedian, resultFreqCutoff));
}

TEST(SeedIndex, ComputeFrequencyThresholdOutOfBounds)
{
    /*
     * This should throw, because freqPercentile should be in [0.0, 1.0].
    */
    // Input values.
    const int32_t seqId = 123;
    const int32_t span = 15;
    std::vector<PacBio::Pancake::Int128t> inSeeds = {
        PacBio::Pancake::Seed::Encode(0, span, seqId, 0, false),
        PacBio::Pancake::Seed::Encode(123, span, seqId, 1, false),
        PacBio::Pancake::Seed::Encode(5, span, seqId, 2, false),
        PacBio::Pancake::Seed::Encode(7, span, seqId, 3, false),
        PacBio::Pancake::Seed::Encode(5, span, seqId, 4, false),
        PacBio::Pancake::Seed::Encode(0, span, seqId, 5, false),
        PacBio::Pancake::Seed::Encode(123, span, seqId, 6, false),
        PacBio::Pancake::Seed::Encode(5, span, seqId, 7, false),
        PacBio::Pancake::Seed::Encode(123, span, seqId, 8, false),
    };

    // Run the unit under test.
    const PacBio::Pancake::SeedIndex si(std::move(inSeeds));
    int64_t resultFreqMax = 0;
    double resultFreqAvg = 0.0;
    double resultFreqMedian = 0.0;
    int64_t resultFreqCutoff = 0;

    // Evaluate.
    EXPECT_THROW(
        {
            si.ComputeFrequencyStats(1.5, resultFreqMax, resultFreqAvg, resultFreqMedian,
                                     resultFreqCutoff);
        },
        std::runtime_error);
    EXPECT_THROW(
        {
            si.ComputeFrequencyStats(-1.0, resultFreqMax, resultFreqAvg, resultFreqMedian,
                                     resultFreqCutoff);
        },
        std::runtime_error);
}

TEST(SeedIndex, ParsingSeedIndexCache_Stream_RoundTrip_Good)
{
    // Load the SeedDB cache.
    const std::string targetSeedDBString =
        R"(V	0.1.0
P	k=15,w=10,s=0,hpc=0,rc=1
F	0	reads.seeddb.0.seeds	2	992
S	0	read1-fwd	0	0	496	180	31
S	1	read4-rev	0	496	496	180	31
B	0	0	2	992
)";
    std::istringstream is(targetSeedDBString);
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> targetSeedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(is, "filename.seeddb");

    std::ostringstream results;
    results << *targetSeedDBCache;

    EXPECT_EQ(targetSeedDBString, results.str());
}

TEST(SeedIndex, ParsingSeedIndexCache_File_RoundTrip_Good)
{
    // Load the SeedDB cache.
    const std::string targetSeedDBString =
        R"(V	0.1.0
P	k=15,w=10,s=0,hpc=0,rc=1
F	0	reads.seeddb.0.seeds	2	992
S	0	read1-fwd	0	0	496	180	31
S	1	read4-rev	0	496	496	180	31
B	0	0	2	992
)";
    std::string seedDBFn = PacBio::PancakeTestsConfig::GeneratedData_Dir + "/test.seeddb";
    {
        std::ofstream ofs(seedDBFn);
        ofs << targetSeedDBString;
    }

    std::FILE* fpIn = std::fopen(seedDBFn.c_str(), "r");
    std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> targetSeedDBCache =
        PacBio::Pancake::LoadSeedDBIndexCache(fpIn, "test.seeddb");

    std::ostringstream results;
    results << *targetSeedDBCache;

    EXPECT_EQ(targetSeedDBString, results.str());
}

TEST(SeedIndex, ParsingSeedIndexCache_RoundTrip_Stream_ExtraWhitespaceAtTheEndShouldThrow)
{
    // Load the SeedDB cache.
    // NOTE: The whitespace on the last line is intentional!
    const std::string targetSeedDBString =
        R"(V	0.1.0
P	k=15,w=10,s=0,hpc=0,hpc_len=10,rc=1
F	0	reads.seeddb.0.seeds	2	992
S	0	read1-fwd	0	0	496	180	31
S	1	read4-rev	0	496	496	180	31
B	0	0	2	992
    )";
    std::istringstream is(targetSeedDBString);

    EXPECT_THROW(
        {
            std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> targetSeedDBCache =
                PacBio::Pancake::LoadSeedDBIndexCache(is, "filename.seeddb");
        },
        std::runtime_error);
}

TEST(SeedIndex, ParsingSeedIndexCache_File_ExtraWhitespaceAtTheEndShouldThrow)
{
    // Load the SeedDB cache.
    const std::string targetSeedDBString =
        R"(V	0.1.0
P	k=15,w=10,s=0,hpc=0,hpc_len=10,rc=1
F	0	reads.seeddb.0.seeds	2	992
S	0	read1-fwd	0	0	496	180	31
S	1	read4-rev	0	496	496	180	31
B	0	0	2	992
    )";
    std::string seedDBFn = PacBio::PancakeTestsConfig::GeneratedData_Dir + "/test.seeddb";
    {
        std::ofstream ofs(seedDBFn);
        ofs << targetSeedDBString;
    }

    std::FILE* fpIn = std::fopen(seedDBFn.c_str(), "r");
    EXPECT_THROW(
        {
            std::shared_ptr<PacBio::Pancake::SeedDBIndexCache> targetSeedDBCache =
                PacBio::Pancake::LoadSeedDBIndexCache(fpIn, "test.seeddb");
        },
        std::runtime_error);
}

TEST(SeedIndex, ComputeSeedHitHistogram_ArrayOfTests)
{
    // clang-format off
    struct TestData
    {
        const std::string testName;
        const std::vector<PacBio::Pancake::SeedRaw> targetSeeds;
        const std::vector<PacBio::Pancake::SeedRaw> querySeeds;
        const std::vector<std::pair<int64_t, int64_t>> expected;
    };

    const int32_t span = 15;
    const int32_t targetId = 0;
    const std::vector<TestData> testData = {
        {   // Test 1.
            "Empty input",
            {},
            {},
            {},
        },
        {   // Test 2.
            "Perfect match (query and target seeds are the same)",
            {   // Target seeds.
                PacBio::Pancake::Seed::Encode(0, span, targetId, 0, false),
                PacBio::Pancake::Seed::Encode(123, span, targetId, 1, false),
                PacBio::Pancake::Seed::Encode(5, span, targetId, 2, false),
                PacBio::Pancake::Seed::Encode(7, span, targetId, 3, false),
                PacBio::Pancake::Seed::Encode(5, span, targetId, 4, false),
                PacBio::Pancake::Seed::Encode(0, span, targetId, 5, false),
                PacBio::Pancake::Seed::Encode(123, span, targetId, 6, false),
                PacBio::Pancake::Seed::Encode(5, span, targetId, 7, false),
                PacBio::Pancake::Seed::Encode(123, span, targetId, 8, false),
            },
            {   // Query seeds.
                PacBio::Pancake::Seed::Encode(0, span, targetId, 0, false),
                PacBio::Pancake::Seed::Encode(123, span, targetId, 1, false),
                PacBio::Pancake::Seed::Encode(5, span, targetId, 2, false),
                PacBio::Pancake::Seed::Encode(7, span, targetId, 3, false),
                PacBio::Pancake::Seed::Encode(5, span, targetId, 4, false),
                PacBio::Pancake::Seed::Encode(0, span, targetId, 5, false),
                PacBio::Pancake::Seed::Encode(123, span, targetId, 6, false),
                PacBio::Pancake::Seed::Encode(5, span, targetId, 7, false),
                PacBio::Pancake::Seed::Encode(123, span, targetId, 8, false),
            },
            {   // Expected results.
                // Seed hit occurrence:
                //  "0":    2
                //  "123":  3
                //  "5":    3
                //  "7":    1
                //  "5":    3
                //  "0":    2
                //  "123":  3
                //  "5":    3
                //  "123":  3
                // Sorted order:
                {1, 1},
                {2, 2},
                {3, 6},
            },
        },
        {   // Test 3.
            "Query and target seed set are different.",
            {   // Target seeds.
                PacBio::Pancake::Seed::Encode(0, span, targetId, 0, false),
                PacBio::Pancake::Seed::Encode(123, span, targetId, 1, false),
                PacBio::Pancake::Seed::Encode(5, span, targetId, 2, false),
                PacBio::Pancake::Seed::Encode(7, span, targetId, 3, false),
                PacBio::Pancake::Seed::Encode(5, span, targetId, 4, false),
                PacBio::Pancake::Seed::Encode(0, span, targetId, 5, false),
                PacBio::Pancake::Seed::Encode(123, span, targetId, 6, false),
                PacBio::Pancake::Seed::Encode(5, span, targetId, 7, false),
                PacBio::Pancake::Seed::Encode(123, span, targetId, 8, false),
            },
            {   // Query seeds.
                PacBio::Pancake::Seed::Encode(0, span, targetId, 0, false),
                PacBio::Pancake::Seed::Encode(99999, span, targetId, 1, false),
                PacBio::Pancake::Seed::Encode(5, span, targetId, 2, false),
                PacBio::Pancake::Seed::Encode(7, span, targetId, 3, false),
                PacBio::Pancake::Seed::Encode(5, span, targetId, 4, false),
                PacBio::Pancake::Seed::Encode(0, span, targetId, 5, false),
                PacBio::Pancake::Seed::Encode(99999, span, targetId, 6, false),
                PacBio::Pancake::Seed::Encode(5, span, targetId, 7, false),
                PacBio::Pancake::Seed::Encode(123, span, targetId, 8, false),
            },
            {   // Expected results.
                //  "0":        2
                //  "99999":    0
                //  "5":        3
                //  "7":        1
                //  "5":        3
                //  "0":        2
                //  "99999":    0
                //  "5":        3
                //  "123":      3
                // Sorted order:
                {0, 2},
                {1, 1},
                {2, 2},
                {3, 4},
            },
        },
    };
    // clang-format on

    for (const auto& data : testData) {
        SCOPED_TRACE(data.testName);
        std::cerr << "testName = " << data.testName << "\n";

        // Construct the hash needed for this test.
        auto targetSeeds = data.targetSeeds;
        const PacBio::Pancake::SeedIndex si(std::move(targetSeeds));

        // Run unit under test.
        const std::vector<std::pair<int64_t, int64_t>> results =
            PacBio::Pancake::ComputeSeedHitHistogram(data.querySeeds, si.GetHash());

        EXPECT_EQ(data.expected, results);
    }
}

TEST(SeedIndex, CollectSeedHitsFromSequences)
{
    // clang-format off
    // Inputs.
    // Note: first 15 bases of the query are reverse complemented.
    const std::string target = "CATGGTGAGTCACCTCTGACTGAGAGTTTACTCACTTAGCCGCGTGTCCACTATTGCTGGGTAAGATCAGACCGTTATTCTCGACAGCGGAA";
    const std::string query =   "GAGGTGACTCACCAT"      "AGTTTACTCACTTAG"                     "AAGATCAGACCGTTA";
    const int32_t k = 15;
    const int32_t w = 1;
    const int32_t spacing = 0;
    const bool useRC = true;
    const bool useHomopolymerCompression = false;
    const double seedFreqPercentile = 0.0;
    const int64_t seedOccurrenceMin = 0;
    const int64_t seedOccurrenceMax = 0;
    const int64_t seedOccurrenceMaxMemory = 0;

    const FastaSequenceCachedStore queryStore({FastaSequenceCached("query", query.c_str(), query.size(), 0)});
    const FastaSequenceCachedStore targetStore({FastaSequenceCached("target", target.c_str(), target.size(), 0)});
    // clang-format on

    // Expected results.
    // Target coordinate is always FWD, query coordinate is strand specific.
    const std::vector<PacBio::Pancake::SeedHit> expected = {
        {0, true, 1, 30, 15, 15, 0},
        {0, false, 24, 15, 15, 15, 0},
        {0, false, 62, 30, 15, 15, 0},
    };

    // Unit under test.
    // The outter vector contains number of elements equal to the number of queries.
    const std::vector<std::vector<SeedHit>> hits = CollectSeedHitsFromSequences(
        queryStore, targetStore, k, w, spacing, useRC, useHomopolymerCompression,
        seedFreqPercentile, seedOccurrenceMin, seedOccurrenceMax, seedOccurrenceMaxMemory);

    // Debug print.
    for (size_t queryId = 0; queryId < hits.size(); ++queryId) {
        for (size_t j = 0; j < hits[queryId].size(); ++j) {
            std::cerr << "[queryId = " << queryId << "]\n";
            const auto& hit = hits[queryId][j];
            std::cerr << "{targetId"
                      << ", " << (hit.targetRev ? "true" : "false") << ", " << hit.targetPos << ", "
                      << hit.queryPos << ", " << static_cast<int32_t>(hit.targetSpan) << ", "
                      << static_cast<int32_t>(hit.querySpan) << ", " << hit.flags << "},\n";
            std::cerr << "-----\n";
        }
    }

    ASSERT_EQ(1, hits.size());
    EXPECT_EQ(expected, hits[0]);
}