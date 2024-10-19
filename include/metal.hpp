#pragma once
#include "../include/material.hpp"
#include "../include/object.hpp"
#include "../include/ray.hpp"

#include "../include/vec3.hpp"

namespace ZrRender {
// 金属材质类
class metal : public material {
public:
    metal(const color &a, const double f) : albedo(a), fuzz(f < 1 ? f : 1) {
    }

    virtual bool scatter(const ray &r_in, const hit_record &rec, color &attenuation, ray &scattered) const override {
        // 传入镜面反射函数都是单位向量
        vec3 reflected = reflect(normalize(r_in.direction), rec.normal);
        scattered      = ray(rec.p, reflected, r_in.time);
        attenuation    = albedo;
        return (dot(scattered.direction, rec.normal) > 0);
    }

    virtual bool scatter(ray &r_in, const hit_record &rec, scatter_record &srec) const override {
        vec3 reflected    = reflect(normalize(r_in.direction), rec.normal);
        srec.specular_ray = ray(rec.p, reflected + fuzz * random_in_unit_sphere());
        srec.attenuation  = albedo;
        srec.is_specular  = true;
        srec.pdf_ptr      = nullptr;
        return true;
    }

public:
    color  albedo;
    double fuzz;  // 反射率
};
}  // namespace ZrRender