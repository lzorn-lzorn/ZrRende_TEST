#pragma once

// #include "../include/material.hpp"
#include "../include/object.hpp"
#include "../include/utilities.hpp"
#include "../include/aabb.hpp"
#include "../include/onb.hpp"
#include "../include/random.hpp"

namespace ZrRender {
class sphere : public object {
public:
    sphere() {
    }

    sphere(point3 cen0, point3 cen1, double _time0, double _time1, double r, std::shared_ptr<material> m)
        : center0(cen0), center1(cen1), time0(_time0), time1(_time1), radius(r), mat_ptr(m) {
    }

    // 重载虚函数
    virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override {
        vec3 oc     = r.origin - center(r.time);
        auto a      = r.direction.length_squared();
        auto half_b = dot(oc, r.direction);
        auto c      = oc.length_squared() - radius * radius;

        auto discriminant = half_b * half_b - a * c;
        if (discriminant < 0)
            return false;

        auto sqrtd = sqrt(discriminant);

        // 找到满足条件的最近的交点
        if (a < epsilon) {
            // std::cerr <<"ray direction is too small" << std::flush;
            throw std::runtime_error("from sphere->hit(): ray direction is too small");
        }

        auto inv_a = 1.0 / a;

        auto root = (-half_b - sqrtd) * inv_a;
        if (root < t_min || t_max < root) {
            root = (-half_b + sqrtd) * inv_a;
            if (root < t_min || t_max < root)
                return false;
        }
        // 记录该交点的相关信息
        rec.t = root;
        rec.p = r.at(rec.t);
        // 法线记得归一化
        auto inv_radius     = 1.0 / radius;
        vec3 outward_normal = (rec.p - center(r.time)) * inv_radius;
        // 判断交点在正面还是背面，并设置正确的法线方向
        rec.set_face_normal(r, outward_normal);
        // 计算纹理坐标
        get_sphere_uv(outward_normal, rec.u, rec.v);
        // 记录材质
        rec.mat_ptr = mat_ptr;
        return true;
    }

    virtual double pdf_value(const point3& o, const vec3& v) const override {
        hit_record rec;
        if (!this->hit(ray(o, v), 0.001, infinity, rec))
            return 0;

        auto cos_theta_max = sqrt(1 - radius * radius / (center0 - o).length_squared());
        auto solid_angle   = 2 * pi * (1 - cos_theta_max);

        return 1 / solid_angle;
    }

    virtual vec3 random(const point3& o) const override {
        vec3 direction        = center0 - o;
        auto distance_squared = direction.length_squared();
        onb  uvw;
        uvw.build_from_w(direction);
        return uvw.local(random_to_sphere(radius, distance_squared));
    }

    point3 center(double time) const {
        auto time_diff = time - time0;
        if (time_diff == 0)
            return center0;
        auto inv_time_diff = 1.0 / (time1 - time0);
        return center0 + (time_diff * inv_time_diff) * (center1 - center0);
    }

    virtual bool bounding_box(double time0, double time1, aabb& output_box) const override {
        aabb box0(center(time0) - vec3(radius, radius, radius), center(time0) + vec3(radius, radius, radius));
        aabb box1(center(time1) - vec3(radius, radius, radius), center(time1) + vec3(radius, radius, radius));

        output_box = surrounding_box(box0, box1);
        return true;
    }

public:
    point3                    center0, center1;
    double                    radius;
    double                    time0, time1;
    std::shared_ptr<material> mat_ptr;  // 材质

private:
    // 计算给定球面上的点p的纹理坐标，p是圆心在原点的单位球面上的坐标，一般用归一化的法线
    static void get_sphere_uv(const point3& p, double& u, double& v) {
        auto theta = std::acos(-p.y);
        auto phi   = std::atan2(-p.z, p.x) + pi;

        auto inv_pi = 1.0 / pi;
        u           = phi * inv_pi * 0.5;
        v           = theta * inv_pi;
    }
};

class xy_rect : public object {
public:
    xy_rect() {
    }

    xy_rect(double _x0, double _x1, double _y0, double _y1, double _k, std::shared_ptr<material> mat)
        : x0(_x0), x1(_x1), y0(_y0), y1(_y1), k(_k), mat_ptr(mat) {
    }

    virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override {
        auto t = (k - r.origin.z) / r.direction.z;
        if (t < t_min || t > t_max)
            return false;
        auto x = r.origin.x + t * r.direction.x;
        auto y = r.origin.y + t * r.direction.y;
        if (x < x0 || x > x1 || y < y0 || y > y1)
            return false;
        rec.u               = (x - x0) / (x1 - x0);
        rec.v               = (y - y0) / (y1 - y0);
        rec.t               = t;
        auto outward_normal = vec3(0, 0, 1);
        rec.set_face_normal(r, outward_normal);
        rec.mat_ptr = mat_ptr;
        rec.p       = r.at(t);
        return true;
    }

    virtual bool bounding_box(double time0, double time1, aabb& output_box) const override {
        // z 方向填充一个很小的长度，防止 BVH 划分出问题
        output_box = aabb(point3(x0, y0, k - 0.0001), point3(x1, y1, k + 0.0001));
        return true;
    }

public:
    double                    x0, x1, y0, y1, k;
    std::shared_ptr<material> mat_ptr;
};

class xz_rect : public object {
public:
    xz_rect() {
    }

    xz_rect(double _x0, double _x1, double _z0, double _z1, double _k, shared_ptr<material> mat)
        : x0(_x0), x1(_x1), z0(_z0), z1(_z1), k(_k), mat_ptr(mat) {};

    virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override {
        auto t = (k - r.origin.y) / r.direction.y;
        if (t < t_min || t > t_max)
            return false;
        auto x = r.origin.x + t * r.direction.x;
        auto z = r.origin.z + t * r.direction.z;
        if (x < x0 || x > x1 || z < z0 || z > z1)
            return false;
        auto inv_diff_x     = 1.0 / (x1 - x0);
        auto inv_diff_z     = 1.0 / (z1 - z0);
        rec.u               = (x - x0) * inv_diff_x;
        rec.v               = (z - z0) * inv_diff_z;
        rec.t               = t;
        auto outward_normal = vec3(0, 1, 0);
        rec.set_face_normal(r, outward_normal);
        rec.mat_ptr = mat_ptr;
        rec.p       = r.at(t);
        return true;
    }

    virtual bool bounding_box(double time0, double time1, aabb& output_box) const override {
        output_box = aabb(point3(x0, k - 0.0001, z0), point3(x1, k + 0.0001, z1));
        return true;
    }

    // 计算随机采样的pdf
    virtual double pdf_value(const point3& origin, const vec3& v) const override {
        hit_record rec;
        if (!this->hit(ray(origin, v), 0.001, infinity, rec))
            return 0;
        // 光源平面面积
        auto area = (x1 - x0) * (z1 - z0);
        // 光源采样点到着色点的距离平方
        auto distance_squared = rec.t * rec.t * v.length_squared();
        // 光线和光源平面法线cos
        auto cosine = fabs(dot(v, rec.normal) / v.length());
        // 概率密度
        return distance_squared / (cosine * area);
    }

    // 随机采样一点，作为随机采样的方向
    virtual vec3 random(const point3& origin) const override {
        auto random_point = point3(random_double(x0, x1), k, random_double(z0, z1));
        return random_point - origin;
    }

public:
    shared_ptr<material> mat_ptr;
    double               x0, x1, z0, z1, k;
};

class yz_rect : public object {
public:
    yz_rect() {
    }

    yz_rect(double _y0, double _y1, double _z0, double _z1, double _k, shared_ptr<material> mat)
        : y0(_y0), y1(_y1), z0(_z0), z1(_z1), k(_k), mat_ptr(mat) {};

    virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override {
        auto t = (k - r.origin.x) / r.direction.x;
        if (t < t_min || t > t_max)
            return false;
        auto y = r.origin.y + t * r.direction.y;
        auto z = r.origin.z + t * r.direction.z;
        if (y < y0 || y > y1 || z < z0 || z > z1)
            return false;
        auto inv_diff_y     = 1.0 / (y1 - y0);
        auto inv_diff_z     = 1.0 / (z1 - z0);
        rec.u               = (y - y0) * inv_diff_y;
        rec.v               = (z - z0) * inv_diff_z;
        rec.t               = t;
        auto outward_normal = vec3(1, 0, 0);
        rec.set_face_normal(r, outward_normal);
        rec.mat_ptr = mat_ptr;
        rec.p       = r.at(t);
        return true;
    }

    virtual bool bounding_box(double time0, double time1, aabb& output_box) const override {
        output_box = aabb(point3(k - 0.0001, y0, z0), point3(k + 0.0001, y1, z1));
        return true;
    }

    void print() const noexcept {
        std::cout << "yz_rect: y0=" << y0 << " y1=" << y1 << " z0=" << z0 << " z1=" << z1 << " k=" << k << std::endl;
    }

public:
    shared_ptr<material> mat_ptr;
    double               y0, y1, z0, z1, k;
};

class box : public object {
public:
    box() {
    }

    box(const point3& p0, const point3& p1, shared_ptr<material> ptr) {
        box_min = p0;
        box_max = p1;

        sides.add(make_shared<xy_rect>(p0.x, p1.x, p0.y, p1.y, p1.z, ptr));
        sides.add(make_shared<xy_rect>(p0.x, p1.x, p0.y, p1.y, p0.z, ptr));

        sides.add(make_shared<xz_rect>(p0.x, p1.x, p0.z, p1.z, p1.y, ptr));
        sides.add(make_shared<xz_rect>(p0.x, p1.x, p0.z, p1.z, p0.y, ptr));

        sides.add(make_shared<yz_rect>(p0.y, p1.y, p0.z, p1.z, p1.x, ptr));
        sides.add(make_shared<yz_rect>(p0.y, p1.y, p0.z, p1.z, p0.x, ptr));
    }

    virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override {
        return sides.hit(r, t_min, t_max, rec);
    }

    virtual bool bounding_box(double time0, double time1, aabb& output_box) const override {
        output_box = aabb(box_min, box_max);
        return true;
    }

public:
    point3       box_min;
    point3       box_max;
    object_group sides;
};

/*
这里将三种变换全部看作是对光线的反向变换，即物体的相对位置不变，光线的方向发生了变化。
*/
// 平移变换
class translate : public object {
public:
    translate(shared_ptr<object> p, const vec3& displacement) : ptr(p), offset(displacement) {
    }

    virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override {
        // 光线向反方向平移
        ray moved_r(r.origin - offset, r.direction, r.time);

        // 计算交点，这里计算出的交点是相对坐标，物体还在原本的地方
        if (!ptr->hit(moved_r, t_min, t_max, rec))
            return false;

        // 把物体和光线的交点加上偏移，得到平移后的物体和光线的交点在世界空间的绝对坐标
        // 这才相当于把物体移动了
        rec.p += offset;
        rec.set_face_normal(moved_r, rec.normal);

        return true;
    }

    virtual bool bounding_box(double time0, double time1, aabb& output_box) const override {
        if (!ptr->bounding_box(time0, time1, output_box))
            return false;

        output_box = aabb(output_box.min() + offset, output_box.max() + offset);

        return true;
    }

public:
    shared_ptr<object> ptr;
    vec3               offset;
};

class rotate_y : public object {
public:
    rotate_y(shared_ptr<object> p, double angle) : ptr(p) {
        auto radians = degrees_to_radians(angle);
        sin_theta    = sin(radians);
        cos_theta    = cos(radians);
        hasbox       = ptr->bounding_box(0, 1, bbox);

        point3 min(infinity, infinity, infinity);
        point3 max(-infinity, -infinity, -infinity);
        // 遍历bounding box的每个顶点，并进行变换
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                for (int k = 0; k < 2; k++) {
                    auto x = i * bbox.max().x + (1 - i) * bbox.min().x;
                    auto y = j * bbox.max().y + (1 - j) * bbox.min().y;
                    auto z = k * bbox.max().z + (1 - k) * bbox.min().z;

                    auto newx = cos_theta * x + sin_theta * z;
                    auto newz = -sin_theta * x + cos_theta * z;

                    vec3 tester(newx, y, newz);

                    for (int c = 0; c < 3; c++) {
                        min[c] = fmin(min[c], tester[c]);
                        max[c] = fmax(max[c], tester[c]);
                    }
                }
            }
        }

        bbox = aabb(min, max);
    }

    virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override {
        auto origin    = r.origin;
        auto direction = r.direction;

        // 光线向反方向旋转
        origin[0] = cos_theta * r.origin[0] - sin_theta * r.origin[2];
        origin[2] = sin_theta * r.origin[0] + cos_theta * r.origin[2];
        // 因为光线方向实际上是两个点的差，所以也可以直接应用变换矩阵
        direction[0] = cos_theta * r.direction[0] - sin_theta * r.direction[2];
        direction[2] = sin_theta * r.direction[0] + cos_theta * r.direction[2];

        ray rotated_r(origin, direction, r.time);

        // 得到的交点同样是相对的坐标
        if (!ptr->hit(rotated_r, t_min, t_max, rec))
            return false;

        auto& p      = rec.p;
        auto  normal = rec.normal;

        // 将交点进行旋转
        p[0] = cos_theta * rec.p[0] + sin_theta * rec.p[2];
        p[2] = -sin_theta * rec.p[0] + cos_theta * rec.p[2];
        // 法线也要旋转，法线变换应该用原变换矩阵的逆转置矩阵，旋转矩阵正交因此逆转置矩阵就是原矩阵
        normal[0] = cos_theta * rec.normal[0] + sin_theta * rec.normal[2];
        normal[2] = -sin_theta * rec.normal[0] + cos_theta * rec.normal[2];

        rec.p = p;
        rec.set_face_normal(rotated_r, normal);

        return true;
    }

    virtual bool bounding_box(double time0, double time1, aabb& output_box) const override {
        output_box = bbox;
        return hasbox;
    }

public:
    shared_ptr<object> ptr;
    double             sin_theta;
    double             cos_theta;
    bool               hasbox;
    aabb               bbox;
};

}  // namespace ZrRender