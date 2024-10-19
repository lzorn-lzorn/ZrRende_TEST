#pragma once

#include "material.hpp"
#include "texture.hpp"
#include "object.hpp"

namespace ZrRender {
// 各向同性材质
class isotropic : public material {
public:
    isotropic(color c) : albedo(make_shared<solid_color>(c)) {
    }

    isotropic(shared_ptr<texture> a) : albedo(a) {
    }

    virtual bool scatter(const ray &r_in, const hit_record &rec, color &attenuation, ray &scattered) const override {
        // 光线向各个方向等概率均匀散射
        scattered   = ray(rec.p, random_in_unit_sphere(), r_in.time);
        attenuation = albedo->value(rec.u, rec.v, rec.p);
        return true;
    }

public:
    shared_ptr<texture> albedo;
};
}  // namespace ZrRender