#pragma once

#include "../include/utilities.hpp"

namespace ZrRender{

class clock{
public:
    clock() : start_time(), end_time(), is_running(false) {}

    void start() {
        start_time = std::chrono::high_resolution_clock::now();
        is_running = true;
    }

    void stop() {
        end_time = std::chrono::high_resolution_clock::now();
        is_running = false;
    }

    double elapsed() const {
        if (is_running) {
            auto current_time = std::chrono::high_resolution_clock::now();
            return std::chrono::duration<double>(current_time - start_time).count();
        } else {
            return std::chrono::duration<double>(end_time - start_time).count();
        }
    }
    std::tuple<int, int> elapsed_minutes_and_seconds() const {
        double total_seconds = elapsed(); // 获取总的秒数
        int minutes = static_cast<int>(total_seconds) / 60; // 计算分钟
        int seconds = static_cast<int>(total_seconds) % 60; // 计算秒数
        return std::make_tuple(minutes, seconds);
    }
private:
    std::chrono::high_resolution_clock::time_point start_time;
    std::chrono::high_resolution_clock::time_point end_time;
    bool is_running;
};    
}