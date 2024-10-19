#pragma once

#include "../include/vec3.hpp"

namespace ZrRender {

inline double random_double() {
    static std::uniform_real_distribution<double> distribution(0.0, 1.0);
    static std::mt19937                           generator;
    return distribution(generator);
}

inline double random_double(double min, double max) {
    return min + (max - min) * random_double();
}

inline vec3 random_vec() {
    return vec3(random_double(), random_double(), random_double());
}

inline vec3 random_unit_vec() {
    return normalize(random_vec());
}

inline vec3 random_vec(double min, double max) {
    return vec3(random_double(min, max), random_double(min, max), random_double(min, max));
}

inline int random_int(int min, int max) {
    return static_cast<int>(random_double(min, max + 1));
}

inline point3 random_in_unit_sphere() {
    point3 p{};
    do {
        p = random_vec(-1.0, 1.0);
    } while (p.length_squared() >= 1.0);
    return p;
}

inline vec3 random_unit_sphere() {
    return normalize(random_in_unit_sphere());
}

// 在半球内随机取点
inline vec3 random_in_hemisphere(const vec3 &normal) {
    vec3 in_unit_sphere = random_in_unit_sphere();
    // 判断该偏移量是否落入了下半球，如果落入下半球则偏移量应该取反
    if (dot(in_unit_sphere, normal) > 0.0)
        return in_unit_sphere;
    else
        return -in_unit_sphere;
}

// 生成单位圆盘内随机一点
inline vec3 random_in_unit_disk() {
    while (true) {
        auto p = vec3(random_double(-1, 1), random_double(-1, 1), 0);
        if (p.length_squared() >= 1)
            continue;
        return p;
    }
}

// 生成概率分布为cos(theta)/pi的随机方向
inline vec3 random_cosine_direction() {
    auto r1 = random_double();
    auto r2 = random_double();
    // 方向是单位向量，所以z坐标就是cos(theta)
    auto z = sqrt(1 - r2);

    auto phi = 2 * pi * r1;
    auto x   = cos(phi) * sqrt(r2);
    auto y   = sin(phi) * sqrt(r2);

    return vec3(x, y, z);
}

// 在球体外对球体随机采样
inline vec3 random_to_sphere(double radius, double distance_squared) {
    auto r1 = random_double();
    auto r2 = random_double();
    auto z  = 1 + r2 * (sqrt(1 - radius * radius / distance_squared) - 1);

    auto phi = 2 * pi * r1;
    auto x   = cos(phi) * sqrt(1 - z * z);
    auto y   = sin(phi) * sqrt(1 - z * z);

    return vec3(x, y, z);
}
};  // namespace ZrRender