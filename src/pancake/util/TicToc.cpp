/*
 * TicToc.cpp
 *
 * This source file was obtained and adjusted from the Raptor
 * graph-based mapping tool codebase.
 *
 *  Created on: Nov 29, 2016
 *      Author: Ivan Sovic
 */

#include <pancake/util/TicToc.hpp>

#include <sstream>

TicToc::TicToc() : start_(), end_()
{
    Start();
    end_ = start_;
}

TicToc::~TicToc() = default;

void TicToc::Start() { start_ = std::chrono::high_resolution_clock::now(); }

void TicToc::Stop() { end_ = std::chrono::high_resolution_clock::now(); }

double TicToc::GetSecs(bool current) const
{
    auto end = (current) ? std::chrono::high_resolution_clock::now() : end_;
    double elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start_).count();
    return elapsed;
}

double TicToc::GetMillisecs(bool current) const
{
    auto end = (current) ? std::chrono::high_resolution_clock::now() : end_;
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start_).count();
    return elapsed;
}

double TicToc::GetMicrosecs(bool current) const
{
    auto end = (current) ? std::chrono::high_resolution_clock::now() : end_;
    double elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start_).count();
    return elapsed;
}

double TicToc::GetNanosecs(bool current) const
{
    auto end = (current) ? std::chrono::high_resolution_clock::now() : end_;
    double elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start_).count();
    return elapsed;
}

std::string TicToc::VerboseSecs(bool current) const
{
    std::ostringstream oss;
    oss << "Real: " << GetSecs(current) << " sec.";
    return oss.str();
}
