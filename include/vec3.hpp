#pragma once

#include "../include/utilities.hpp"

namespace ZrRender {

class vec3 {
public:
    double x, y, z;

    vec3() : x(0.0), y(0.0), z(0.0) {
    }

    vec3(double x, double y, double z) : x(x), y(y), z(z) {
    }

    const double &operator[](int i) const {
        switch (i) {
            case 0:
                return x;
            case 1:
                return y;
            case 2:
                return z;
            default:
                throw std::runtime_error("from vec3->operator[]: index out of range");
        }
    }

    double &operator[](int i) {
        switch (i) {
            case 0:
                return x;
            case 1:
                return y;
            case 2:
                return z;
            default:
                throw std::runtime_error("from vec3->operator[]: index out of range");
        }
    }

    vec3 operator-() const {
        return vec3(-x, -y, -z);
    }

    vec3 operator=(const vec3 &v) {
        x = v.x;
        y = v.y;
        z = v.z;
        return *this;
    }

    vec3 operator=(vec3 &v) {
        x = v.x;
        y = v.y;
        z = v.z;
        return *this;
    }

    vec3 &operator+=(const vec3 &other) {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }

    vec3 &operator-=(const vec3 &other) {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        return *this;
    }

    vec3 &operator*=(double t) {
        x *= t;
        y *= t;
        z *= t;
        return *this;
    }

    vec3 &operator/=(double t) {
        if (std::fabs(t) < epsilon) {
            throw std::runtime_error("from vec3->operator/=: division by zero");
            return *this;
        }
        auto inv_t = 1.0 / t;

        x *= inv_t;
        y *= inv_t;
        z *= inv_t;
        return *this;
    }

    bool operator==(const vec3 &other) const {
        return x == other.x && y == other.y && z == other.z;
    }

    bool operator!=(const vec3 &other) const {
        return !(*this == other);
    }

    double length() const {
        return std::sqrt(x * x + y * y + z * z);
    }

    double length_squared() const {
        return x * x + y * y + z * z;
    }

    double length_squared() {
        return x * x + y * y + z * z;
    }

    vec3 normalize() const {
        double len = length();

        if (len < epsilon) {
            throw std::runtime_error("from vec3->normalize(): normalize failure.division by zero");
        }
        auto inv_len = 1.0 / len;
        return vec3(x * inv_len, y * inv_len, z * inv_len);
    }

    double dot(const vec3 &other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    vec3 cross(const vec3 &other) const {
        return vec3(y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x);
    }

public:
    bool near_zero() const {
        return (std::fabs(x) < epsilon) && (std::fabs(y) < epsilon) && (std::fabs(z) < epsilon);
    }
};

using point3 = vec3;
using color  = vec3;

inline std::ostream &operator<<(std::ostream &out, const vec3 &v) {
    // return out << "(" << v.x << ", " << v.y << ", " << v.z << ")";
    out << std::format("({}, {}, {})\n", v.x, v.y, v.z);
    return out;
}

inline double dot(const vec3 &a, const vec3 &b) {
    return a.dot(b);
}

inline vec3 cross(const vec3 &a, const vec3 &b) {
    return a.cross(b);
}

inline vec3 normalize(const vec3 &v) {
    return v.normalize();
}

inline double length(const vec3 &v) {
    return v.length();
}

inline vec3 operator+(const vec3 &a, const vec3 &b) {
    return vec3(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline vec3 operator+(const vec3 &a, double t) {
    return vec3(a.x + t, a.y + t, a.z + t);
}

inline vec3 operator+(double t, const vec3 &a) {
    return a + t;
}

inline vec3 operator-(const vec3 &a, const vec3 &b) {
    return vec3(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline vec3 operator-(const vec3 &a, double t) {
    return vec3(a.x - t, a.y - t, a.z - t);
}

inline vec3 operator*(const vec3 &a, const vec3 &b) {
    return vec3(a.x * b.x, a.y * b.y, a.z * b.z);
}

inline vec3 operator*(const vec3 &a, double t) {
    return vec3(a.x * t, a.y * t, a.z * t);
}

inline vec3 operator*(double t, const vec3 &a) {
    return a * t;
}

inline vec3 operator/(const vec3 &a, double t) {
    if (std::fabs(t) < epsilon) {
        throw std::runtime_error("from inline operator/: division by zero");
        // return vec3(0.0, 0.0, 0.0);
    }
    auto inv_t = 1.0 / t;
    return vec3(a.x * inv_t, a.y * inv_t, a.z * inv_t);
}
}  // namespace ZrRender

/*
if (abs(t))
*/