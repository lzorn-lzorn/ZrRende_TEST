#pragma once

#include "../include/utilities.hpp"
#include "../include/vec3.hpp"
#include "../include/noise.hpp"
#define STB_IMAGE_IMPLEMENTATION
#include "../include/stb_image.h"

namespace ZrRender {

class texture {
public:
    virtual color value(double u, double v, const point3& p) const = 0;
};

// 纯色纹理
class solid_color : public texture {
public:
    solid_color() {
    }

    solid_color(const color& c) : color_value(c) {
    }

    solid_color(double red, double green, double blue) : solid_color(color(red, green, blue)) {
    }

    virtual color value(double u, double v, const point3& p) const override {
        return color_value;
    }

private:
    color color_value;
};

// 格子纹理
class checker_texture : public texture {
public:
    checker_texture() {
    }

    checker_texture(shared_ptr<texture> _even, shared_ptr<texture> _odd) : even(_even), odd(_odd) {
    }

    checker_texture(const color& c1, const color& c2)
        : even(make_shared<solid_color>(c1)), odd(make_shared<solid_color>(c2)) {
    }

    virtual color value(double u, double v, const point3& p) const override {
        auto sines = sin(10 * p.x) * sin(10 * p.y) * sin(10 * p.z);
        if (sines < 0)
            return odd->value(u, v, p);
        else
            return even->value(u, v, p);
    }

public:
    shared_ptr<texture> even;
    shared_ptr<texture> odd;
};

class perlin_noise_texture : public texture {
public:
    perlin_noise_texture() = default;

    perlin_noise_texture(double sc) : scale(sc) {
    }

    virtual color value(double u, double v, const point3& p) const override {
        return color(1, 1, 1) * 0.5 * (1 + sin(scale * p.z + 10 * noise.turb(p)));
    }

public:
    perlin noise;
    double scale;
};

class image_texture : public texture {
public:
    const static int bytes_per_pixel = 3;

    image_texture() : data(nullptr), width(0), height(0), bytes_per_scanline(0) {
    }

    image_texture(const char* filename) {
        auto components_per_pixel = bytes_per_pixel;

        data = stbi_load(filename, &width, &height, &components_per_pixel, components_per_pixel);

        if (!data) {
            std::cerr << "ERROR: Could not load texture image file '" << filename << "'.\n";
            width = height = 0;
        }

        bytes_per_scanline = bytes_per_pixel * width;
    }

    ~image_texture() {
        delete data;
    }

    virtual color value(double u, double v, const vec3& p) const override {
        // 没有图片的话返回一个固定的颜色，便于debug
        if (data == nullptr)
            return color(0, 1, 1);

        // 输入的纹理坐标截断在[0,1]区间
        u = clamp(u, 0.0, 1.0);
        // 由于图片是从上到下存储的，相当于原点在左上角，而纹理坐标原点在左下角，因此纵坐标要翻转一下
        v = 1.0 - clamp(v, 0.0, 1.0);
        // 纹理坐标映射到图片坐标
        auto i = static_cast<int>(u * width);
        auto j = static_cast<int>(v * height);
        if (i >= width)
            i = width - 1;
        if (j >= height)
            j = height - 1;
        // 我们返回的颜色都在[0,1]之间，因此要除以255
        const auto color_scale = 1.0 / 255.0;
        auto       pixel       = data + j * bytes_per_scanline + i * bytes_per_pixel;

        return color(color_scale * pixel[0], color_scale * pixel[1], color_scale * pixel[2]);
    }

private:
    unsigned char* data;
    int            width, height;
    int            bytes_per_scanline;
};
}  // namespace ZrRender