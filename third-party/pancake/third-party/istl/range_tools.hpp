/*
 *  range_tools.hpp
 *
 *  Created on: Sep 12, 2018
 *      Author: Ivan Sovic
 */

#include <functional>
#include <tuple>
#include <vector>

#ifndef SRC_ISTL_RANGE_TOOLS_HPP_
#define SRC_ISTL_RANGE_TOOLS_HPP_

namespace istl {

/* Returns a vector of pair numbers: [start, end> indexes within the data vector,
 * which demarcate the range satisfied by the comparison function.
 * The comparison function should return true if the two elements are supposed to be
 * groupped together, otherwise false.
 * The first template class Q is the atomic data type.
 * The second template class T is an iterable container of atomic types, std::vector by default.
*/
template <class Q, class T = std::vector<Q>>
std::vector<std::pair<size_t, size_t>> FindRanges(
    const T& data, std::function<bool(const Q& a, const Q& b)> comp_eq =
                       [](const Q& a, const Q& b) { return a == b; })
{

    std::vector<std::pair<size_t, size_t>> ret;
    if (data.empty()) {
        return ret;
    }
    size_t prevEnd = 0;
    for (size_t pos = 1; pos < data.size(); ++pos) {
        if (comp_eq(data[pos], data[pos - 1])) {
            continue;
        }
        ret.emplace_back(prevEnd, pos);
        prevEnd = pos;
    }
    if (prevEnd != data.size()) {
        ret.emplace_back(prevEnd, data.size());
    }
    return ret;
}

}  // namespace istl

#endif
