#include "apps/common/Timer.hpp"

#include <cassert>
#include <iomanip>
#include <sstream>

namespace fem
{
Timer::Timer()
    : m_started(false)
{
}

void Timer::start()
{
    assert(!m_started);
    m_started = true;
    m_startTime = std::chrono::high_resolution_clock::now();
    m_previousLapTime = m_startTime;
}

Timer::Duration Timer::lap()
{
    assert(m_started);
    const auto currentTime = std::chrono::high_resolution_clock::now();
    const auto lapTime = currentTime - m_previousLapTime;
    m_previousLapTime = currentTime;
    return lapTime;
}

Timer::Duration Timer::getElapsedTime()
{
    assert(m_started);
    const auto currentTime = std::chrono::high_resolution_clock::now();
    return currentTime - m_startTime;
}

std::string minutesSecondsMilliseconds(Timer::Duration d)
{
    using namespace std::chrono;

    const minutes min = duration_cast<minutes>(d);
    d -= min;
    const seconds sec = duration_cast<seconds>(d);
    d -= sec;
    const milliseconds ms = duration_cast<milliseconds>(d);

    std::stringstream ss;
    ss << min.count() << "m";
    ss << sec.count() << "." << std::setfill('0') << std::setw(3) << ms.count() << "s";
    return ss.str();
}
} // namespace fem
