/*
 * TicToc.h
 *
 * This source file was obtained and adjusted from the Raptor
 * graph-based mapping tool codebase.
 *
 *  Created on: Nov 29, 2016
 *      Author: Ivan Sovic
 */

#ifndef PANCAKE_UTIL_TIC_TOC_HPP
#define PANCAKE_UTIL_TIC_TOC_HPP

#include <chrono>
#include <ctime>
#include <string>

class TicToc
{
public:
    TicToc();
    ~TicToc();

    void Start();
    void Stop();
    double GetSecs(bool current = false) const;
    double GetMillisecs(bool current = false) const;
    double GetMicrosecs(bool current = false) const;
    double GetNanosecs(bool current = false) const;

    std::string VerboseSecs(bool current = false) const;

private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start_;
    std::chrono::time_point<std::chrono::high_resolution_clock> end_;
};

#endif  // PANCAKE_UTIL_TIC_TOC_HPP
