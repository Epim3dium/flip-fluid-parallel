#ifndef EMP_TIME_HPP
#define EMP_TIME_HPP
#include <chrono>
#include <string>
#include <unordered_map>

class Stopwatch {
    double delta_time = -1.0;
    typedef  std::chrono::duration<double> duration_t;
    typedef std::chrono::time_point<std::chrono::high_resolution_clock> time_point_t;
    time_point_t m_start;
    time_point_t m_stop;
public:
    Stopwatch() {
        m_start = std::chrono::high_resolution_clock::now();
    }
    double getElapsedTime() {
        m_stop = std::chrono::high_resolution_clock::now();
        auto dur = m_stop - m_start;
        return static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(dur).count()) / 1e9;
    }
    double stop() {
        m_stop = std::chrono::high_resolution_clock::now();
        auto dur = m_stop - m_start;
        delta_time = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(dur).count()) / 1e9;
        return delta_time;
    }
    double restart() {
        m_stop = std::chrono::high_resolution_clock::now();
        auto dur = m_stop - m_start;
        double elapsed = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(dur).count()) / 1e9;
        m_start = std::chrono::high_resolution_clock::now();
        return elapsed;
    }
};

#endif// EMP_TIME_HPP
