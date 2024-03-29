#include "apps/common/Timer.hpp"

#include <cassert>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace fem
{
void Timer::start(const std::string& msg)
{
    std::cout << msg << std::flush;
    m_startTime = std::chrono::high_resolution_clock::now();
}

void Timer::stop()
{
    const auto currentTime = std::chrono::high_resolution_clock::now();
    const auto elapsedTime = currentTime - m_startTime;
    std::cout << minutesSecondsMilliseconds(elapsedTime) << std::endl;
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
