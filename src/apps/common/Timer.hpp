#pragma once

#include <chrono>
#include <string>

namespace fem
{
class Timer
{
public:
    using Duration = std::chrono::high_resolution_clock::duration;

public:
    Timer();

    void start();
    Duration lap();
    Duration getElapsedTime();

private:
    bool m_started;
    std::chrono::high_resolution_clock::time_point m_startTime;
    std::chrono::high_resolution_clock::time_point m_previousLapTime;
};

std::string minutesSecondsMilliseconds(Timer::Duration d);
} // namespace fem
