/**
 * Implementation of the Fenwick tree data structure,
 * supporting the following operations:
 *  a) get(x): gets the maximum element in the interval [0,x]
 *  b) update(x, v): sets the x-th element to v, if v
 *                   is bigger than the previously stored value
 *                   on position x, otherwise nothing happens
 *  Note: x >= 0
 * @author: Filip Pavetic (fpavetic@gmail.com)
 * https://github.com/fpavetic/lcskpp
 *
 *
    The MIT License (MIT)

    Copyright (c) 2014 Filip Pavetic

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
 */

#ifndef LCSKPP_FENWICK
#define LCSKPP_FENWICK

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <iostream>
#include <utility>
#include <vector>

namespace istl {

template <class T>
class FenwickMax
{
public:
    FenwickMax(size_t n) { elements_ = std::vector<T>(n + 1); }

    void update(size_t pos, const T& val)
    {
        ++pos;
        for (; pos < elements_.size(); pos += lobit(pos)) {
            elements_[pos] = std::max(elements_[pos], val);
        }
    }

    T get(size_t pos)
    {
        ++pos;
        T ret = T();
        for (; pos > 0; pos -= lobit(pos)) {
            ret = std::max(ret, elements_[pos]);
        }
        return ret;
    }

private:
    size_t lobit(const size_t& a) { return a & -a; }

private:
    std::vector<T> elements_;
};
}

#endif
