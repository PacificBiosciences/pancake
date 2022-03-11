#include <gtest/gtest.h>

#include <istl/lcs.hpp>

#include <cstdint>
#include <cstring>
#include <functional>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

TEST(LCS, EmptyInput)
{
    std::vector<std::pair<int32_t, int32_t>> data = {};

    std::vector<std::pair<int32_t, int32_t>> result = istl::LCSk<int32_t>(data, 1);

    std::vector<std::pair<int32_t, int32_t>> expected = {};

    ASSERT_EQ(expected, result);
}

TEST(LCS, SinglePoint)
{
    std::vector<std::pair<int32_t, int32_t>> data = {
        {1, 1},
    };

    std::vector<std::pair<int32_t, int32_t>> result = istl::LCSk<int32_t>(data, 1);

    std::vector<std::pair<int32_t, int32_t>> expected = {
        {1, 1},
    };

    ASSERT_EQ(expected, result);
}

TEST(LCS, SimpleTest2D)
{
    std::vector<std::pair<int32_t, int32_t>> data = {
        {1, 1}, {2, 2}, {3, 1}, {4, 4}, {5, 5},
    };

    std::vector<std::pair<int32_t, int32_t>> result = istl::LCSk<int32_t>(data, 1);

    std::vector<std::pair<int32_t, int32_t>> expected = {
        {1, 1},
        {2, 2},
        {4, 4},
        {5, 5},
    };

    ASSERT_EQ(expected, result);
}

TEST(LCS, SimpleTest2D_2)
{
    /**
     * @brief This test won't work well with 1D LIS.
     */

    std::vector<std::pair<int32_t, int32_t>> data = {
        {1, 1}, {1, 4}, {3, 3}, {4, 4}, {5, 5}, {6, 6},
    };

    std::vector<std::pair<int32_t, int32_t>> result = istl::LCSk<int32_t>(data, 1);

    std::vector<std::pair<int32_t, int32_t>> expected = {
        {1, 1}, {3, 3}, {4, 4}, {5, 5}, {6, 6},
    };

    ASSERT_EQ(expected, result);
}

TEST(LCS, RealTest1)
{
    /**
     * @brief This test won't work well with 1D LIS. It would result in a shorter chain.
     */

    std::vector<std::pair<int32_t, int32_t>> data = {
        {4342, 24},   {4349, 31},   {4353, 38},   {4940, 1113}, {4975, 679},  {4983, 687},
        {5035, 947},  {5035, 1040}, {5035, 1101}, {5043, 1046}, {5043, 1108}, {5048, 1113},
        {5048, 1228}, {5065, 1299}, {5072, 846},  {5074, 848},  {5095, 947},  {5095, 1040},
        {5095, 1101}, {5095, 1331}, {5103, 1046}, {5103, 1108}, {5258, 980},  {5265, 988},
        {5299, 722},  {5304, 859},  {5339, 734},  {5494, 1116}, {5568, 915},  {5592, 1046},
        {5592, 1108}, {5627, 1299}, {5658, 1113}, {5658, 1228},
    };

    std::vector<std::pair<int32_t, int32_t>> result = istl::LCSk<int32_t>(data, 1);

    std::vector<std::pair<int32_t, int32_t>> expected = {
        {4342, 24},  {4349, 31},  {4353, 38},  {4975, 679}, {4983, 687},  {5072, 846},
        {5074, 848}, {5095, 947}, {5258, 980}, {5265, 988}, {5592, 1108}, {5627, 1299},
    };

    ASSERT_EQ(expected, result);
}

TEST(LCS, RealTest2_TestMultiplePointsWithSameXCoordinate)
{
    std::vector<std::pair<int32_t, int32_t>> data = {
        {6585, 38493}, {6585, 40327}, {29565, 11779}, {29577, 11791}, {29637, 11849},
    };

    std::vector<std::pair<int32_t, int32_t>> result = istl::LCSk<int32_t>(data, 1);

    std::vector<std::pair<int32_t, int32_t>> expected = {
        {29565, 11779},
        {29577, 11791},
        {29637, 11849},
    };

    ASSERT_EQ(expected, result);
}
