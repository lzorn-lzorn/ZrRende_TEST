#pragma once

#include "object.hpp"
#include "material.hpp"
#include "texture.hpp"
#include "isotropic.hpp"

namespace ZrRender {

class constant_medium : public object {
public:
    constant_medium(shared_ptr<object> b, double d, shared_ptr<texture> a)
        : boundary(b), neg_inv_density(-1 / d), phase_function(make_shared<isotropic>(a)) {
    }

    constant_medium(shared_ptr<object> b, double d, color c)
        : boundary(b), neg_inv_density(-1 / d), phase_function(make_shared<isotropic>(c)) {
    }

    virtual bool hit(const ray &r, double t_min, double t_max, hit_record &rec) const override {
        // 用于debug
        const bool enableDebug = false;
        const bool debugging   = enableDebug && random_double() < 0.00001;

        // 求光线和边界的两个交点
        hit_record rec1, rec2;

        if (!boundary->hit(r, -infinity, infinity, rec1))
            return false;

        if (!boundary->hit(r, rec1.t + 0.0001, infinity, rec2))
            return false;

        if (debugging)
            std::cerr << "\nt_min=" << rec1.t << ", t_max=" << rec2.t << '\n';

        if (rec1.t < t_min)
            rec1.t = t_min;
        if (rec2.t > t_max)
            rec2.t = t_max;

        if (rec1.t >= rec2.t)
            return false;

        if (rec1.t < 0)
            rec1.t = 0;

        // 光线在介质中的距离
        const auto ray_length               = r.direction.length();
        const auto distance_inside_boundary = (rec2.t - rec1.t) * ray_length;
        // 光线发生散射的距离，两个相乘的数都是小于1的负数，所以密度越大值越小
        const auto hit_distance = neg_inv_density * log(random_double());

        // 发生散射的距离大于光线在介质中的距离则没有发生散射，直接穿过介质
        if (hit_distance > distance_inside_boundary)
            return false;

        // 散射发生的位置
        rec.t = rec1.t + hit_distance / ray_length;
        rec.p = r.at(rec.t);

        if (debugging) {
            std::cerr << "hit_distance = " << hit_distance << '\n'
                      << "rec.t = " << rec.t << '\n'
                      << "rec.p = " << rec.p << '\n';
        }

        // 法线方向这些属性可以随便设置
        rec.normal     = vec3(1, 0, 0);
        rec.front_face = true;
        rec.mat_ptr    = phase_function;

        return true;
    }

    virtual bool bounding_box(double time0, double time1, aabb &output_box) const override {
        return boundary->bounding_box(time0, time1, output_box);
    }

public:
    shared_ptr<object>   boundary;         // 边界
    shared_ptr<material> phase_function;   // 各向同性材质，保证光线向各个方向等概率均匀散射
    double               neg_inv_density;  // 密度的负倒数
};

}  // namespace ZrRender