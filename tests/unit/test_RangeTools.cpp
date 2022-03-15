#include <gtest/gtest.h>

#include <pancake/third-party/istl/range_tools.hpp>

#include <cstdint>
#include <cstring>
#include <string>
#include <tuple>
#include <vector>

// clang-format off
TEST(RangeTools, FindRange1) {
    const std::vector<int32_t> data = {1, 2, 3, 4, 5, 6, 7, 8, 9, 0};

    const auto result = istl::FindRanges<int32_t>(data);

    const std::vector<std::pair<size_t, size_t>> expected = {
                                            {0, 1},
                                            {1, 2},
                                            {2, 3},
                                            {3, 4},
                                            {4, 5},
                                            {5, 6},
                                            {6, 7},
                                            {7, 8},
                                            {8, 9},
                                            {9, 10},
    };

    ASSERT_EQ(result, expected);
}

TEST(RangeTools, FindRange2) {
    const std::vector<int32_t> data = {1, 2, 2, 3, 4, 5, 6, 7, 8, 9, 0};

    const auto result= istl::FindRanges<int32_t>(data);

    const  std::vector<std::pair<size_t, size_t>> expected = {
                                            {0, 1},
                                            {1, 3},
                                            {3, 4},
                                            {4, 5},
                                            {5, 6},
                                            {6, 7},
                                            {7, 8},
                                            {8, 9},
                                            {9, 10},
                                            {10, 11},
    };

    ASSERT_EQ(result, expected);
}

TEST(RangeTools, FindRange3) {
    const std::vector<int32_t> data = {1, 2, 2, 3, 3, 3, 3, 3, 4};

    const auto result= istl::FindRanges<int32_t>(data);

    const std::vector<std::pair<size_t, size_t>> expected = {
                                            {0, 1},
                                            {1, 3},
                                            {3, 8},
                                            {8, 9},
    };

    ASSERT_EQ(result, expected);
}

TEST(RangeTools, FindRange4) {
    const std::vector<int32_t> data = {1, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4};

    const auto result= istl::FindRanges<int32_t>(data);

    const std::vector<std::pair<size_t, size_t>> expected = {
                                            {0, 1},
                                            {1, 3},
                                            {3, 8},
                                            {8, 13},
    };

    ASSERT_EQ(result, expected);
}

TEST(RangeTools, FindRange5) {
    const std::vector<int32_t> data = {};

    const auto result= istl::FindRanges<int32_t>(data);

    const std::vector<std::pair<size_t, size_t>> expected = { };

    ASSERT_EQ(result, expected);
}

TEST(RangeTools, FindRange6) {
    const std::vector<int32_t> data = {1};

    const auto result= istl::FindRanges<int32_t>(data);

    const std::vector<std::pair<size_t, size_t>> expected = {
                                            {0, 1},
    };

    ASSERT_EQ(result, expected);
}

TEST(RangeTools, FindRangeInCustomClass1) {
    class TestDataClass {
        public:
        TestDataClass(int64_t _target_id, bool _target_rev)
                    : target_id(_target_id), target_rev(_target_rev) { }
        int64_t target_id;
        bool target_rev;
    };

    const std::vector<TestDataClass> data = {
                                    TestDataClass(0, false),
                                    TestDataClass(0, false),
                                    TestDataClass(0, false),
                                    TestDataClass(0, true),
                                    TestDataClass(0, true),
                                    TestDataClass(1, true),
                                    TestDataClass(1, true),
                                    TestDataClass(1, true),
                                    TestDataClass(1, true),
                                    TestDataClass(1, true),
                                    TestDataClass(17, false),
                                };

    const auto result = istl::FindRanges<TestDataClass>(data,
                                    [](const TestDataClass& a, const TestDataClass& b) {
                                        return (a.target_id == b.target_id && a.target_rev == b.target_rev);
                                    }
    );

    const std::vector<std::pair<size_t, size_t>> expected = {
                                            {0, 3},
                                            {3, 5},
                                            {5, 10},
                                            {10, 11},
    };

    ASSERT_EQ(result, expected);
}

TEST(RangeTools, FindRangeInCustomClass2) {
    const std::string data = "122233455555";

    const auto result = istl::FindRanges<char, std::string>(data, [](char a, char b) { return a == b; });

    const std::vector<std::pair<size_t, size_t>> expected = {
                                            {0, 1},
                                            {1, 4},
                                            {4, 6},
                                            {6, 7},
                                            {7, 12},
    };

    ASSERT_EQ(result, expected);
}

// clang-format on
