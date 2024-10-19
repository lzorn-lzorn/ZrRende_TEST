#pragma once

#include "../include/material.hpp"
#include "../include/texture.hpp"
#include "../include/object.hpp"

namespace ZrRender {

class diffuse_light : public material {

public:
    diffuse_light(shared_ptr<texture> a) : emit(a) {
    }

    diffuse_light(color c) : emit(make_shared<solid_color>(c)) {
    }

    virtual bool scatter(const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered) const override {
        return false;
    }

    virtual color emitted(double u, double v, const point3& p) const override {
        return emit->value(u, v, p);
    }

    // 只有正面发光
    virtual color emitted(const ray& r_in, const hit_record& rec, double u, double v, const point3& p) const override {
        if (rec.front_face)
            return emit->value(u, v, p);
        else
            return color(0, 0, 0);
    }

public:
    shared_ptr<texture> emit;
};
}  // namespace ZrRender