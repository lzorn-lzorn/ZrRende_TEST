#pragma once
#include <cmath>
#include <limits>
#include <memory>
#include <vector>
#include <string>
#include <random>
#include <iostream>
#include <format>
#include <ostream>
#include <sstream>
#include <chrono>
#include <algorithm>
#include <functional>
#include <thread>
#include <mutex>
#include <array>
#include <map>
#include <variant>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../include/stb_image_write.h"

namespace ZrRender {
using std::jthread;
using std::make_shared;
using std::mutex;
using std::shared_ptr;
using std::unique_ptr;
static constexpr double epsilon  = std::numeric_limits<double>::epsilon();
static constexpr double infinity = std::numeric_limits<double>::infinity();
static constexpr double pi       = 3.14159265358979323846;

// 将给定的数字截取到[min,max]范围内
inline double clamp(double value, double min, double max) {
    return value < min ? min : (value > max ? max : value);
}

inline double degrees_to_radians(double degrees) {
    return degrees * pi / 180.0;
}

}  // namespace ZrRender
