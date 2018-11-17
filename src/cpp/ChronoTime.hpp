/**
    ChronoTime
    ChronoTime.hpp
    Purpose: Time functions

    @author J.M.G. Salmerón
    @version 1.0 11/10/2018
    
    Created on: 11/10/2018
*/

#ifndef CHRONOTIME_H_
#define CHRONOTIME_H_

#include <chrono>
#include <sstream>
#include <iomanip> // Format numbers

namespace chronotime {

    std::chrono::steady_clock::time_point temp_start_time;
    std::chrono::steady_clock::time_point temp_second_time;

    inline
    void start_time() {
        temp_start_time = std::chrono::steady_clock::now();
    }

    inline
    auto final_time() { // return auto -> c++14
        return std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - temp_start_time).count();
    }

    inline
    void start_second_time() {
        temp_second_time = std::chrono::steady_clock::now();
    }

    inline
    auto final_second_time() {
        return std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - temp_second_time).count();
    }

    inline
    std::string format_time(auto time) {
        int hour = 0;
        int min = 0;
        int sec = 0;

        hour = time / 3600;
        time = time % 3600;
        min = time / 60;
        time = time % 60;
        sec = time;

        std::stringstream ss;

        ss  << std::setw(2) << std::setfill('0') << hour
            << ":" << std::setw(2) << std::setfill('0') << min
            << ":" << std::setw(2) << std::setfill('0') << sec;

        return ss.str();
    }

    template<typename F, typename... Args>
    double funcTime(F func, Args&&... args){
        std::chrono::steady_clock::time_point t1;
        t1 = std::chrono::steady_clock::now();
        func(std::forward<Args>(args)...);
        return std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - temp_start_time).count();
    }

}

#endif /* CHRONOTIME_H_ */