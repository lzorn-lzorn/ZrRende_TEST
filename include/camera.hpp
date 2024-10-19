#pragma once

#include "../include/ray.hpp"
#include "../include/utilities.hpp"
#include "../include/random.hpp"
#include "vec3.hpp"

namespace ZrRender {
class camera;

class camera {
public:
    camera() {
    }

    camera(point3 lookfrom, point3 lookat, vec3 vup, double vfov, double aspect_ratio,
           double aperture,    // 光圈大小，光圈为0就是之前的针孔相机
           double focus_dist,  // 焦点距离，在焦点距离处的物体不会发生散焦模糊
           double _time0 = 0,  // 快门开启时间
           double _time1 = 0   // 快门关闭时间
    ) {
        // 通过视场角和宽高比计算视口大小，视场角是垂直方向的，单位是角度
        auto theta           = degrees_to_radians(vfov);
        auto h               = tan(theta / 2);
        auto viewport_height = 2.0 * h;
        auto viewport_width  = aspect_ratio * viewport_height;
        // 计算视口平面的三个坐标轴
        w = normalize(lookfrom - lookat);
        u = normalize(cross(vup, w));
        v = cross(w, u);
        // 相机位置
        origin = lookfrom;
        // 视口
        horizontal        = focus_dist * viewport_width * u;
        vertical          = focus_dist * viewport_height * v;
        lower_left_corner = origin - horizontal * 0.5 - vertical * 0.5 - focus_dist * w;
        // 镜头半径等于光圈大小的一半
        lens_radius = aperture * 0.5;
        // 快门时间
        time0 = _time0;
        time1 = _time1;
    }

    ray get_ray(double s, double t) const {
        // 圆盘内随机取一点作为偏移
        vec3 rd     = lens_radius * random_in_unit_disk();
        vec3 offset = u * rd.x + v * rd.y;
        // 从偏离远镜头的位置投射光线，模拟散焦
        return ray(origin + offset, lower_left_corner + s * horizontal + t * vertical - origin - offset,
                   random_double(time0, time1));
    }

private:
    point3 origin;
    point3 lower_left_corner;
    vec3   horizontal;
    vec3   vertical;
    vec3   u, v, w;
    double lens_radius;
    double time0, time1;
    //  camera_args args;
};

}  // namespace ZrRender