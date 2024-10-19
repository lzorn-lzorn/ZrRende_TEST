#include "../include/camera.hpp"
#include "../include/object.hpp"
#include "../include/ray.hpp"
#include "../include/shape.hpp"
#include "../include/utilities.hpp"
#include "../include/utilities.hpp"
#include "../include/vec3.hpp"
#include "../include/material.hpp"
#include "../include/lambertian.hpp"
#include "../include/metal.hpp"
#include "../include/color.hpp"
#include "../include/dielectric.hpp"
#include "../include/random.hpp"
#include "../include/tools/clock.hpp"
#include "../include/texture.hpp"
#include "../include/BVH.hpp"
#include "../include/light.hpp"
#include "../include/constant_medium.hpp"
#include "../include/texture.hpp"
#include "../include/render.hpp"
#include <memory>
#include "../include/tools/log.h"

using namespace ZrRender;

shared_ptr<scene> load_random_scene() {
    std::cout << "Loading random scene..." << std::endl;
    auto _scene          = std::make_shared<scene>();
    auto checker         = make_shared<checker_texture>(color(0.2, 0.3, 0.1), color(0.9, 0.9, 0.9));
    auto ground_material = make_shared<lambertian>(checker);
    _scene->add("random_scene")
        .add(make_shared<sphere>(point3(0, -1000, 0), point3(0, -1000, 0), 0.0, 1.0, 1000, ground_material));
    std::cout << "Step1 over..." << std::endl;
    for (int a = -11; a < 11; a++) {
        for (int b = -11; b < 11; b++) {
            auto   choose_mat = random_double();
            point3 center(a + 0.9 * random_double(), 0.2, b + 0.9 * random_double());

            if ((center - point3(4, 0.2, 0)).length() > 0.9) {
                std::shared_ptr<material> sphere_material;

                if (choose_mat < 0.8) {
                    // diffuse
                    auto albedo     = random_vec() * random_vec();
                    sphere_material = make_shared<lambertian>(albedo);
                    // 移动的小球
                    auto center2 = center + vec3(0, random_double(0, 0.5), 0);
                    _scene->add(make_shared<sphere>(center, center2, 0.0, 1.0, 0.2, sphere_material));
                } else if (choose_mat < 0.9) {
                    // metal
                    auto albedo     = random_vec(0.5, 1);
                    auto fuzz       = random_double(0, 0.5);
                    sphere_material = make_shared<metal>(albedo, fuzz);
                    _scene->add(make_shared<sphere>(center, center, 0.0, 1.0, 0.2, sphere_material));
                } else {
                    // glass
                    sphere_material = make_shared<dielectric>(1.5);
                    _scene->add(make_shared<sphere>(center, center, 0.0, 1.0, 0.2, sphere_material));
                }
            }
        }
    }
    _scene->build_bvh();
    std::cout << "Step2 over..." << std::endl;
    return _scene;
}

template <typename... Args>
std::string LOG(int level, int target, std::string format, Args &&...args) {
    std::string head_format = "[{},{},level:{}]:Message:[\n\t your_contexts \n]\n ";

    std::string tail_str = std::string(std::vformat(format, std::make_format_args(std::forward<Args>(args)...)));
    std::string str = std::format("[{},{},level:{}]:Message:[ ", __FILE_NAME__, __LINE__, level) + tail_str + "]\n";
    std::cout << "str:" << str << std::endl;
}

int main() {
    ZrRender::clock timer;
    timer.start();

    render render;
    render.add_scene(load_random_scene());
    std::cout << "Scene loaded." << std::endl;
    render.set_camera(ZrRender::tag::lookfrom, point3(13, 2, 3))
        .set_camera(ZrRender::tag::lookat, point3(0, 0, 0))
        .set_camera(ZrRender::tag::vfov, 20.0)
        .set_camera(ZrRender::tag::aperture, 0.1)
        .set_camera(ZrRender::tag::aspect_ratio, 16.0 / 9)
        .set_camera(ZrRender::tag::vup, vec3(0, 1, 0))
        .set_background(color(0.7, 0.8, 1.0))
        .set_output_filename("random_scene.png")
        .set_output_path("C:/Program User/code/ZrRender/output/")
        .set_bounce(100)
        .set_image(600);
    std::cout << "Render set done." << std::endl;
    render.build("random_scene").rendering().make_pic();
    std::cout << "Rendering done." << std::endl;
    timer.stop();
    auto [minutes, seconds] = timer.elapsed_minutes_and_seconds();
    std::cerr << "\nDone.Elapsed time: " << minutes << " minutes and " << seconds << " seconds." << std::endl;

    return 0;
}