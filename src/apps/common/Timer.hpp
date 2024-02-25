#pragma once

#include <chrono>
#include <ostream>

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

struct MinutesSecondsMilliseconds
{
    Timer::Duration duration;
};

std::ostream& operator<<(std::ostream& out, MinutesSecondsMilliseconds msms);
} // namespace fem
