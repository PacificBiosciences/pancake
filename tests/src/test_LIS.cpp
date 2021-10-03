#include <gtest/gtest.h>

#include <lib/istl/lis.hpp>

#include <cstdint>
#include <cstring>
#include <functional>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

struct Point
{
    int32_t x = 0;
    int32_t y = 0;

    Point() : x(0), y(0) {}
    Point(int32_t _x, int32_t _y) : x(_x), y(_y) {}

    bool operator==(const Point& b) const { return x == b.x && y == b.y; }
};
inline std::ostream& operator<<(std::ostream& os, const Point& a)
{
    os << "x = " << a.x << ", y = " << a.y;
    return os;
}
std::function<bool(const Point& a, const Point& b)> ComparisonLIS = [](
    const Point& a, const Point& b) { return (a.x < b.x && a.y < b.y); };

TEST(LIS, EmptyInput)
{
    std::vector<Point> data = {};

    std::vector<Point> result = istl::LIS<Point>(data, 0, data.size(), ComparisonLIS);

    std::vector<Point> expected = {};

    ASSERT_EQ(expected, result);
}

TEST(LIS, SinglePoint)
{
    std::vector<Point> data = {
        {1, 1},
    };

    std::vector<Point> result = istl::LIS<Point>(data, 0, data.size(), ComparisonLIS);

    std::vector<Point> expected = {
        {1, 1},
    };

    ASSERT_EQ(expected, result);
}

TEST(LIS, SimpleTest1D)
{
    std::vector<int32_t> data = {10, 22, 9, 33, 21, 50, 41, 60};

    std::vector<int32_t> result = istl::LIS<int32_t>(data, 0, data.size());

    std::vector<int32_t> expected = {10, 22, 33, 41, 60};

    ASSERT_EQ(expected, result);
}

TEST(LIS, SimpleTest1D_2)
{
    std::vector<int32_t> data = {3, 7, 5};

    std::vector<int32_t> result = istl::LIS<int32_t>(data, 0, data.size());

    std::vector<int32_t> expected = {3, 5};

    ASSERT_EQ(expected, result);
}

TEST(LIS, SimpleTest2D)
{
    std::vector<Point> data = {
        {1, 1}, {2, 2}, {3, 1}, {4, 4}, {5, 5},
    };

    std::vector<Point> result = istl::LIS<Point>(data, 0, data.size(), ComparisonLIS);

    std::vector<Point> expected = {
        {1, 1}, {2, 2}, {4, 4}, {5, 5},
    };

    ASSERT_EQ(expected, result);
}

TEST(LIS, MoreComplexTest1D)
{
    std::vector<int32_t> data = {
        24,   31,   38,  1113, 679, 687,  947,  1040, 1101, 1046, 1108, 1113,
        1228, 1299, 846, 848,  947, 1040, 1101, 1331, 1046, 1108, 980,  988,
        722,  859,  734, 1116, 915, 1046, 1108, 1299, 1113, 1228,
    };

    std::vector<int32_t> result = istl::LIS<int32_t>(data, 0, data.size());

    std::vector<int32_t> expected = {24,  31,  38,  679,  687,  846,  848,
                                     947, 980, 988, 1046, 1108, 1113, 1228};

    ASSERT_EQ(expected, result);
}

TEST(LIS, RealTest1)
{
    std::vector<Point> data = {
        {4342, 24},   {4349, 31},   {4353, 38},   {4940, 1113}, {4975, 679},  {4983, 687},
        {5035, 947},  {5035, 1040}, {5035, 1101}, {5043, 1046}, {5043, 1108}, {5048, 1113},
        {5048, 1228}, {5065, 1299}, {5072, 846},  {5074, 848},  {5095, 947},  {5095, 1040},
        {5095, 1101}, {5095, 1331}, {5103, 1046}, {5103, 1108}, {5258, 980},  {5265, 988},
        {5299, 722},  {5304, 859},  {5339, 734},  {5494, 1116}, {5568, 915},  {5592, 1046},
        {5592, 1108}, {5627, 1299}, {5658, 1113}, {5658, 1228},
    };

    std::vector<Point> result = istl::LIS<Point>(data, 0, data.size(), ComparisonLIS);

    std::vector<Point> expected = {
        {4342, 24},  {4349, 31},  {4353, 38},  {4975, 679},  {4983, 687},  {5072, 846},
        {5074, 848}, {5304, 859}, {5568, 915}, {5592, 1108}, {5658, 1228},
    };

    for (size_t i = 0; i < result.size(); ++i) {
        std::cerr << "[i = " << i << "] " << result[i] << "\n";
    }

    ASSERT_EQ(expected, result);
}

TEST(LIS, RealTest2_TestMultiplePointsWithSameXCoordinate)
{
    std::vector<Point> data = {
        {6585, 38493}, {6585, 40327}, {29565, 11779}, {29577, 11791}, {29637, 11849},
    };

    std::vector<Point> result = istl::LIS<Point>(data, 0, data.size(), ComparisonLIS);

    std::vector<Point> expected = {
        {29565, 11779}, {29577, 11791}, {29637, 11849},
    };

    for (size_t i = 0; i < result.size(); ++i) {
        std::cerr << "[i = " << i << "] " << result[i] << "\n";
    }

    ASSERT_EQ(expected, result);
}
