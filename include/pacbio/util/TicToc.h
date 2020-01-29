/*
 * TicToc.h
 *
 * This source file was obtained and adjusted from the Raptor
 * graph-based mapping tool codebase.
 *
 *  Created on: Nov 29, 2016
 *      Author: Ivan Sovic
 */

#ifndef SRC_UTIL_TICTOC_H_
#define SRC_UTIL_TICTOC_H_

#include <chrono>
#include <ctime>

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

    double GetCpuSecs(bool current = false) const;
    double GetCpuMillisecs(bool current = false) const;
    double GetCpuMicrosecs(bool current = false) const;
    double GetCpuNanosecs(bool current = false) const;

private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start_;
    std::chrono::time_point<std::chrono::high_resolution_clock> end_;
    std::clock_t startCpu_;
    std::clock_t endCpu_;

    inline double GetCpuDuration_(bool current, double factor) const;
};

#endif /* SRC_UTIL_TICTOC_H_ */
