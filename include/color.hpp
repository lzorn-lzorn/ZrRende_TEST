#pragma once

#include "../include/ray.hpp"
#include "../include/vec3.hpp"
#include "../include/utilities.hpp"
#include "../include/object.hpp"
#include "../include/random.hpp"
#include "../include/BVH.hpp"
#include "../include/pdf.hpp"
#include "tools/log.h"

namespace ZrRender {
// 向数组中写入一个颜色，用于最后的图像输出，用到了指针的引用传递
// 注意输入的color是[0,1]范围的
inline void write_color(unsigned char *&p, color PixelColor, int SamplesPerPixel) {
    auto r = PixelColor.x;
    auto g = PixelColor.y;
    auto b = PixelColor.z;

    // 处理异常像素值
    if (r != r)
        r = 0.0;
    if (g != g)
        g = 0.0;
    if (b != b)
        b = 0.0;

    auto scale = 1.0 / SamplesPerPixel;

    // 线性空间到非线性空间的转换  gamma = 2.0 所以是平方根
    auto inv_gamma = 1.0 / 2.2;
    r              = pow(scale * r, 1.0 * inv_gamma);
    g              = pow(scale * g, 1.0 * inv_gamma);
    b              = pow(scale * b, 1.0 * inv_gamma);

    *p++ = (unsigned char)(256 * clamp(r, 0.0, 0.999));
    *p++ = (unsigned char)(256 * clamp(g, 0.0, 0.999));
    *p++ = (unsigned char)(256 * clamp(b, 0.0, 0.999));
}

// 计算光线的颜色，使用BVH加速
inline color ray_color_with_bvh(ray &r, const bvh_node &world, int depth, double RR) {
    hit_record rec{};

    if (world.hit(r, 0.001, infinity, rec)) {
        // 一定概率停止弹射
        if (random_double() >= RR)
            return color(0, 0, 0);
        // 根据物体材质得到光线传播方向和反射率
        ray   scattered;
        color attenuation;
        if (rec.mat_ptr->scatter(r, rec, attenuation, scattered))
            return attenuation * ray_color_with_bvh(scattered, world, depth - 1, RR) / RR;

        return color(0, 0, 0);
    }
    vec3 unit_direction = normalize(r.direction);
    auto t              = 0.5 * (unit_direction.y + 1.0);
    return (1.0 - t) * color(1.0, 1.0, 1.0) + t * color(0.5, 0.7, 1.0);
}

inline color ray_color_with_background(
    const ray &r, const color &background, const bvh_node &world, int depth, double RR) {
    hit_record rec{};

    if (depth <= 0 && random_double() >= RR)
        return color(0, 0, 0);

    if (!world.hit(r, 0.001, infinity, rec))
        return background;

    ray   scattered;
    color attenuation;
    color emitted = rec.mat_ptr->emitted(rec.u, rec.v, rec.p);

    if (!rec.mat_ptr->scatter(r, rec, attenuation, scattered))
        return emitted;

    return emitted + attenuation * ray_color_with_background(scattered, background, world, depth - 1, RR) / RR;
}

// 得到光线颜色
inline color ray_color_without_bvh(
    const ray &r, const color &background, const object_group &world, int depth, double RR) {
    hit_record rec;

    // 光线弹射指定次数后开始用RR算法终止递归
    if (depth < 0 && random_double() >= RR)
        return color(0, 0, 0);

    // 如果光线没有打到任何物体，返回背景颜色
    // 这里的t的下界设为0.001是为了防止一些光线弹射到物体上得到的t非常接近0，比如可能出现0.000001这样的值
    if (!world.hit(r, 0.001, infinity, rec))
        return background;

    // 根据物体材质得到光线传播方向、反射率及自发光颜色
    ray   scattered;
    color attenuation;
    color emitted = rec.mat_ptr->emitted(rec.u, rec.v, rec.p);
    // 对于光源，不会发生散射，返回光源颜色
    if (!rec.mat_ptr->scatter(r, rec, attenuation, scattered))
        return emitted;

    return emitted + attenuation * ray_color_without_bvh(scattered, background, world, depth - 1, RR) / RR;
}

// 得到光线颜色
inline color ray_color(
    ray &r, const color &background, const bvh_node &world, shared_ptr<object> lights, int depth, double RR) {
    hit_record rec;

    // 光线弹射指定次数后开始用RR算法终止递归
    if (depth < 0 && random_double() >= RR)
        return color(0, 0, 0);

    // 如果光线没有打到任何物体，返回背景颜色
    // 这里的t的下界设为0.001是为了防止一些光线弹射到物体上得到的t非常接近0，比如可能出现0.000001这样的值
    if (!world.hit(r, 0.001, infinity, rec)) {
        return background;
    }

    // 根据物体材质得到光线传播方向、反射率及自发光颜色
    scatter_record srec;
    color          emitted = rec.mat_ptr->emitted(r, rec, rec.u, rec.v, rec.p);

    // 对于光源，不会发生散射，返回光源颜色
    if (!rec.mat_ptr->scatter(r, rec, srec)) {
        return emitted;
    }

    // 如果是高光反射或折射，采用之前的渲染方程，隐式的使采样pdf和散射pdf保持一致
    if (srec.is_specular) {
        return srec.attenuation * ray_color(srec.specular_ray, background, world, lights, depth - 1, RR) / RR;
    }

    // 对光源采样的pdf
    auto light_ptr = make_shared<hittable_pdf>(lights, rec.p);

    // 混合pdf
    mixture_pdf p(light_ptr, srec.pdf_ptr);

    // 采样光线
    ray scattered = ray(rec.p, p.generate(), r.time);

    // 采样光线的pdf值
    auto pdf_val = p.value(scattered.direction);

    // 渲染方程
    return emitted
           + srec.attenuation * rec.mat_ptr->scattering_pdf(r, rec, scattered)
                 * ray_color(scattered, background, world, lights, depth - 1, RR) / pdf_val / RR;
}
}  // namespace ZrRender