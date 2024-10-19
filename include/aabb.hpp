#pragma once

#include "../include/utilities.hpp"
#include "../include/ray.hpp"
#include "../include/vec3.hpp"

namespace ZrRender{

class aabb{
public:
    
    aabb(){}
    aabb(const point3& min, const point3& max):minslab(min), maxslab(max){};
    aabb(const aabb& other):minslab(other.minslab), maxslab(other.maxslab){};

    point3 min() const { return minslab; }
    point3 max() const { return maxslab; }

    point3 min()  { return minslab; }
    point3 max()  { return maxslab; }

    bool hit(const ray& r, double tmin, double tmax) const
    {
        for (int i = 0; i < 3; ++i)
        {
            auto inv_dir_i = 1.0 / r.direction[i];
            auto a = (min()[i] - r.origin[i]) * inv_dir_i;
            auto b = (max()[i] - r.origin[i]) * inv_dir_i;
            if(inv_dir_i < 0.0f){
                std::swap(a, b);
            }
            tmin = std::max(tmin, a);
            tmax = std::min(tmax, b);
            if (tmax <= tmin)
                return false;
        }
        return true;
    }
private:
    point3 minslab, maxslab;    
};    

    // 合并两个包围盒形成一个大的包围盒
aabb surrounding_box(aabb box0, aabb box1) {
    point3 small(fmin(box0.min().x, box1.min().x),
        fmin(box0.min().y, box1.min().y),
        fmin(box0.min().z, box1.min().z));

    point3 big(fmax(box0.max().x, box1.max().x),
        fmax(box0.max().y, box1.max().y),
        fmax(box0.max().z, box1.max().z));

    return aabb(small, big);
}
}