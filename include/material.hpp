#pragma once

#include "../include/ray.hpp"
#include "../include/vec3.hpp"

namespace ZrRender {

struct hit_record;
struct pdf;

// 统一管理散射光线
struct scatter_record {
    ray             specular_ray;  // 散射光线
    bool            is_specular;   // 是否是镜面反射，金属或者电介质为true
    color           attenuation;   // 反射率
    shared_ptr<pdf> pdf_ptr;       // 散射光线pdf，如果是金属或者电介质就是空指针
};

class material {
public:
    virtual bool scatter(const ray &r_in, const hit_record &rec, color &attenuation, ray &scattered) const {
        return false;
    }

    virtual bool scatter(const ray &r_in, const hit_record &rec, color &attenuation, ray &scattered,
                         double &pdf) const {
        return false;
    }

    virtual bool scatter(ray &r_in, const hit_record &rec, scatter_record &srec) const {
        return false;
    }

    virtual double scattering_pdf(const ray &r_in, const hit_record &rec, const ray &scattered) const {
        return 0.0;
    }

    virtual color emitted(double u, double v, const point3 &p) const {
        return color(0, 0, 0);
    }

    virtual color emitted(const ray &r_in, const hit_record &rec, double u, double v, const point3 &p) const {
        return color(0, 0, 0);
    }
};

}  // namespace ZrRender
