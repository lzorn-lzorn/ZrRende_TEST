#pragma once

#include "../include/vec3.hpp"

namespace ZrRender {
class ray {
public:
  ray() {}
  ray(const vec3 &origin, const vec3 &direction, double time = 0.0)
      : origin(origin), direction(direction), time(time) {}

  point3 at(double t) const {
    return origin + direction * t;
  }
  point3 origin;
  vec3 direction;
  double time;
};

inline vec3 reflect(const vec3 &v, const vec3 &n) {
  return v - 2 * dot(v, n) * n;
}

inline vec3 refract(const vec3& R, const vec3& n, double etai_over_etat){
  auto cos_theta = std::fmin(dot(-R, n), 1.0);
  vec3 r_out_perp = etai_over_etat * (R + cos_theta * n);
  vec3 r_out_parallel = -sqrt(std::fabs(1.0-r_out_perp.length_squared())) * n;
  return r_out_perp + r_out_parallel;
}
} // namespace ZrRender