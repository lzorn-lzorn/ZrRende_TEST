#pragma once

#include "../include/object.hpp"
#include "../include/random.hpp"
#include "../include/texture.hpp"
#include "../include/vec3.hpp"
#include "../include/onb.hpp"
#include "../include/pdf.hpp"

namespace ZrRender {
class lambertian : public material {
public:
    lambertian(const color& a) : albedo(make_shared<solid_color>(a)) {
    }

    lambertian(const shared_ptr<texture>& a) : albedo(a) {
    }

    virtual bool scatter(const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered) const override {
        // 这里省略了rec.p + rec.normal + random_unit_vector() - rec.p中的rec.p;
        auto scatter_direction = rec.normal + random_unit_sphere();
        // 如果散射方向为0，则取法线方向作为散射方向
        if (scatter_direction.near_zero()) {
            scatter_direction = rec.normal;
        }
        // 散射光线时刻和入射光线一样
        scattered   = ray(rec.p, scatter_direction, r_in.time);
        attenuation = albedo->value(rec.u, rec.v, rec.p);
        return true;
    }

    virtual bool scatter(const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered,
                         double& pdf) const override {
        // 构建法线空间的标准正交基
        onb uvw;
        uvw.build_from_w(rec.normal);
        // 得到法线空间下概率分布为cos(theta)/pi的随机方向
        auto direction = uvw.local(random_cosine_direction());
        scattered      = ray(rec.p, normalize(direction), r_in.time);
        attenuation    = albedo->value(rec.u, rec.v, rec.p);
        // 采样光线的概率密度，即cos(theta)/pi
        pdf = dot(uvw.w(), scattered.direction) / pi;
        return true;
    }

    // 采样散射光线
    virtual bool scatter(ray& r_in, const hit_record& rec, scatter_record& srec) const override {
        srec.is_specular = false;
        srec.attenuation = albedo->value(rec.u, rec.v, rec.p);
        // 采样光线的概率密度
        srec.pdf_ptr = make_shared<cosine_pdf>(rec.normal);
        return true;
    }

    double scattering_pdf(const ray& r_in, const hit_record& rec, const ray& scattered) const override {
        auto cosine = dot(rec.normal, normalize(scattered.direction));
        return cosine < 0 ? 0 : cosine / pi;
    }

public:
    shared_ptr<texture> albedo;
};

}  // namespace ZrRender