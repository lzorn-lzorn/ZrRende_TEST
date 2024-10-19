#pragma once

#include "BVH.hpp"
#include "camera.hpp"
#include "object.hpp"
#include "vec3.hpp"
#include "utilities.hpp"
#include <memory>
#include <string>
#include <variant>
#include "color.hpp"
#include "../include/tools/log.h"

namespace ZrRender {
class scene {
public:
    scene() {
    }

    scene &add(std::string _name) {
        name = _name;
        return *this;
    }

    scene &add(const char *_name) {
        this->name = std::string(_name);
        return *this;
    }

    scene &add(std::shared_ptr<object> object) {
        world.add(object);
        return *this;
    }

    scene &add(object_group objects) {
        for (auto obj : objects.group) {
            world.add(obj);
        }
        return *this;
    }

    scene &set_light(std::shared_ptr<object> light) {
        this->light = light;
        return *this;
    }

    bvh_node &build_bvh() {
        bvh_world = bvh_node(world, 0, 1);
        return bvh_world;
    }

public:
    bvh_node                bvh_world;
    std::shared_ptr<object> light;
    object_group            world;
    std::string             name;
};

struct camera_args {
    point3 lookfrom{};
    point3 lookat{};
    vec3   vup{};
    double vfov         = 40.0;
    double aspect_ratio = 16.0 / 9.0;
    double aperture     = 0.0;   // 光圈大小，光圈为0就是之前的针孔相机
    double focus_dist   = 10.0;  // 焦点距离，在焦点距离处的物体不会发生散焦模糊
    double _time0       = 0;     // 快门开启时间
    double _time1       = 0;     // 快门关闭时间

    camera build() const {
        return camera(lookfrom, lookat, vup, vfov, aspect_ratio, aperture, focus_dist, _time0, _time1);
    }
};
enum class tag {
    lookfrom,
    lookat,
    vup,
    vfov,
    aspect_ratio,
    aperture,
    focus_dist,
    time_range
};

class render {
public:
    render() = default;

    ~render() {
    }

    render &build(std::string &name) {
        detail_build(name);
        return *this;
    }

    render &build(const char *name) {
        detail_build(std::string(name));
        return *this;
    }

    render &rendering() {
        if (mutli_render_check()) {
            single_render();
            // mutli_render();
        } else {
            single_render();
        }
        return *this;
    }

    void make_pic() {
        std::string filepath = output_path + filename;
        stbi_write_png(filepath.c_str(), image_width, image_height, channel, picture_data.get(), 0);
    }

public:  // 设置相机参数
    render &set_camera(tag tag, double value) {
        detail_switch(tag, value);
        return *this;
    }

    render &set_camera(tag tag, vec3 &&value) {
        detail_switch(tag, value);
        return *this;
    }

    render &set_camera(tag tag, vec3 &value) {
        detail_switch(tag, value);
        return *this;
    }

private:
    template <typename _Ty>
    void detail_switch(tag tag, _Ty &&_value) {
        std::variant<vec3, double> value = _value;
        switch (tag) {
            case tag::lookfrom:
                args.lookfrom = std::get<vec3>(value);
                break;
            case tag::lookat:
                args.lookat = std::get<vec3>(value);
                break;
            case tag::vup:
                args.vup = std::get<vec3>(value);
                break;
            case tag::vfov:
                args.vfov = std::get<double>(value);
                break;
            case tag::aspect_ratio:
                args.aspect_ratio = std::get<double>(value);
                break;
            case tag::aperture:
                args.aperture = std::get<double>(value);
                break;
            case tag::focus_dist:
                args.focus_dist = std::get<double>(value);
                break;
            case tag::time_range:
                args._time0 = std::get<double>(value);
                args._time1 = std::get<double>(value);
                break;
            default:
                break;
        }
    }

    template <typename _String>
    void detail_build(_String &&name) {
        cam                  = args.build();
        cur_world            = worlds.at(name);
        bvh_world            = cur_world->build_bvh();
        image_height         = static_cast<int>(image_width / args.aspect_ratio);
        std::string filepath = output_path + filename;
        picture_data         = std::make_unique<unsigned char[]>(image_width * image_height * channel);
        //= (unsigned char *)malloc(image_width * image_height * channel);
    }

public:  // 设置场景参数
    render &add_scene(std::shared_ptr<scene> _scene) {
        worlds[_scene->name] = _scene;
        return *this;
    }

    render &set_background(color &color) {
        background = color;
        return *this;
    }

    render &set_background(color &&color) {
        background = color;
        return *this;
    }

public:  // 设置渲染参数
    render &set_bounce(int bounce) {
        min_bounce = bounce;
        return *this;
    }

    render &set_samples_per_pixel(int samples) {
        samples_per_pixel = samples;
        return *this;
    }

public:  // 设置输出图片参数
    render &set_image(int width) {
        image_width  = width;
        image_height = static_cast<int>(image_width / args.aspect_ratio);
        return *this;
    }

    render &set_channels(int channel) {
        channel = channel;
        return *this;
    }

    render &set_output_filename(std::string filename) {
        filename = filename;
        return *this;
    }

    render &set_output_filename(const char *filename) {
        this->filename = std::string(filename);
        return *this;
    }

    render &set_output_path(std::string path) {
        output_path = path;
        return *this;
    }

    render &set_output_path(const char *path) {
        output_path = std::string(path);
        return *this;
    }

private:
    void render_detail(std::shared_ptr<unsigned char[]> data, int image_height, int image_width) {
        unsigned char *data_ptr = data.get();
        for (int j = image_height - 1; j >= 0; --j) {
            // 标准错误流显示进度信息，单行刷新显示
#ifndef MUTLI_THREAD_RENDER
            std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
#endif
            for (int i = 0; i < image_width; ++i) {
                color pixel_color(0, 0, 0);
                for (int s = 0; s < samples_per_pixel; ++s) {
                    // 在像素内部随机采样
                    auto u = (i + random_double()) / (image_width - 1);
                    auto v = (j + random_double()) / (image_height - 1);
                    ray  r = cam.get_ray(u, v);
                    // pixel_color += ray_color_with_bvh(r, world_box,
                    // min_bounce, RR); pixel_color +=
                    // ray_color_with_background(r, background, world_box,
                    // min_bounce, RR); pixel_color += ray_color_without_bvh(r,
                    // background, world, min_bounce, RR);

                    pixel_color += ray_color(r, background, bvh_world, cur_world->light, min_bounce, RR);
                }
                write_color(data_ptr, pixel_color, samples_per_pixel);
            }
        }
        data_ptr = nullptr;
    }

    bool mutli_render_check() {
        int threads_num = std::thread::hardware_concurrency();
        if (image_width % threads_num == 0) {
            return true;
        } else {
            return false;
        }
    }

    bool mutli_render() {
        using Task = std::function<void(int, int)>;

        Task task = [&](int blocks, int number) {
            thread_local std::shared_ptr<unsigned char[]> local_data =
                std::make_unique<unsigned char[]>(image_height * blocks);  // 线程的本地存储

            render_detail(local_data, image_height, blocks);  // 将渲染任务的结果存入本地存储中

            int offset = number * image_height * blocks;
            std::copy(local_data.get(), local_data.get() + (image_height * blocks), picture_data.get() + offset);
            // 将线程的本地存取写入全局的数组

            local_data = nullptr;
        };

        int threads_num = std::thread::hardware_concurrency();
        int block_size  = image_width / threads_num;

        std::vector<jthread> threads(threads_num);

        for (int i = 0; i < 8; i++) {
            threads[i] = jthread(task, block_size, i);
        }
        return true;
    }

    void single_render() {
        render_detail(picture_data, image_height, image_width);
    }

private:
    camera                                        cam;
    camera_args                                   args;
    std::map<std::string, std::shared_ptr<scene>> worlds;
    std::shared_ptr<scene>                        cur_world;
    bvh_node                                      bvh_world;
    int                                           samples_per_pixel = 100;
    int                                           min_bounce        = 50;
    const double                                  RR                = 0.9;
    int                                           image_width       = 1080;
    int                                           image_height      = 1920;
    int                                           channel           = 3;
    std::string                                   filename          = "output.png";
    std::string                                   output_path       = "";
    std::shared_ptr<unsigned char[]>              picture_data      = nullptr;
    color                                         background        = color(0.0, 0.0, 0.0);
};
}  // namespace ZrRender