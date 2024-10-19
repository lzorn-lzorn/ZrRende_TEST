#pragma once
#include "../include/material.hpp"
#include "../include/object.hpp"
#include "../include/random.hpp"
#include "../include/vec3.hpp"
#include "../include/ray.hpp"

namespace ZrRender {

class dielectric : public material {
public:
    dielectric(double index_of_refraction) : ir(index_of_refraction) {
    }

    virtual bool scatter(const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered) const override {
        // 电介质不吸收任何光，全部被折射或者反射，所以衰减系数恒为1
        attenuation = color(1.0, 1.0, 1.0);
        // 如果是正面，则是从空气进入介质，反之从介质折射出去
        auto   inv_ir           = 1.0 / ir;
        double refraction_ratio = rec.front_face ? (1.0 * inv_ir) : ir;
        // 入射光线记得单位化，折射函数传入的参数都是单位向量
        vec3 unit_direction = normalize(r_in.direction);
        // 计算是否发生折射
        double cos_theta = fmin(dot(-unit_direction, rec.normal), 1.0);
        double sin_theta = sqrt(1.0 - cos_theta * cos_theta);

        bool cannot_refract = refraction_ratio * sin_theta > 1.0;
        vec3 direction;

        if (cannot_refract || reflectance(cos_theta, refraction_ratio) > random_double())
            direction = reflect(unit_direction, rec.normal);
        else
            direction = refract(unit_direction, rec.normal, refraction_ratio);

        scattered = ray(rec.p, direction, r_in.time);
        return true;
    }

    virtual bool scatter(ray& r_in, const hit_record& rec, scatter_record& srec) const override {
        srec.is_specular = true;
        srec.pdf_ptr     = nullptr;
        // 电介质不吸收任何光，全部被折射或者反射，所以衰减系数恒为1
        srec.attenuation = color(1.0, 1.0, 1.0);
        // 如果是正面，则是从空气进入介质，反之从介质折射出去
        double refraction_ratio = rec.front_face ? (1.0 / ir) : ir;
        // 入射光线记得单位化，折射函数传入的参数都是单位向量
        vec3 unit_direction = normalize(r_in.direction);
        // 计算是否发生折射
        double cos_theta = fmin(dot(-unit_direction, rec.normal), 1.0);
        double sin_theta = sqrt(1.0 - cos_theta * cos_theta);

        bool cannot_refract = refraction_ratio * sin_theta > 1.0;
        vec3 direction;

        if (cannot_refract || reflectance(cos_theta, refraction_ratio) > random_double())
            direction = reflect(unit_direction, rec.normal);
        else
            direction = refract(unit_direction, rec.normal, refraction_ratio);

        srec.specular_ray = ray(rec.p, direction, r_in.time);
        return true;
    }

public:
    double ir;  // 介质的折射率

private:
    static double reflectance(double cosine, double ref_idx) {
        // 使用Schlick's近似计算菲涅尔项
        auto inv_ref_idx_temp = 1.0 / (1 + ref_idx);
        auto r0               = (1 - ref_idx) * inv_ref_idx_temp;
        r0 *= r0;
        return r0 + (1 - r0) * pow((1 - cosine), 5);
    }
};

}  // namespace ZrRender