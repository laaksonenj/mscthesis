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
    Timer() = default;

    void start(const std::string& msg);
    void stop();

private:
    std::chrono::high_resolution_clock::time_point m_startTime;
};

std::string minutesSecondsMilliseconds(Timer::Duration d);
} // namespace fem
