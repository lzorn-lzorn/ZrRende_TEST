#pragma once

#include "../include/utilities.hpp"
#include "../include/object.hpp"
#include "../include/random.hpp"

namespace ZrRender {
bool box_x_compare(const shared_ptr<object> a, const shared_ptr<object> b);
bool box_y_compare(const shared_ptr<object> a, const shared_ptr<object> b);
bool box_z_compare(const shared_ptr<object> a, const shared_ptr<object> b);

class bvh_node : public object {
public:
    bvh_node() {
    }

    bvh_node(const std::vector<std::shared_ptr<object>>& src_objects, size_t start, size_t end, double time0,
             double time1) {
        // 不要修改原来的物体列表
        auto objects = src_objects;
        // 随机选择一个划分方向，并传递相应的比较函数指针
        int  axis       = random_int(0, 2);
        auto comparator = (axis == 0) ? box_x_compare : (axis == 1) ? box_y_compare : box_z_compare;
        // 区间内的物体数量
        size_t object_span = end - start;
        // 只有一个物体，同时放到两个孩子节点
        if (object_span == 1) {
            left = right = objects[start];
        }
        // 有两个物体，左右孩子结点各一个
        else if (object_span == 2) {
            if (comparator(objects[start], objects[start + 1])) {
                left  = objects[start];
                right = objects[start + 1];
            } else {
                left  = objects[start + 1];
                right = objects[start];
            }
        }
        // 多个物体递归的进行划分
        else {
            // 先按照指定的维度进行排序
            std::sort(objects.begin() + start, objects.begin() + end, comparator);
            // 二分
            auto mid = start + object_span / 2;
            left     = make_shared<bvh_node>(objects, start, mid, time0, time1);
            right    = make_shared<bvh_node>(objects, mid, end, time0, time1);
        }

        // 当前节点的box由左右两个子树的box的合并而来
        aabb box_left, box_right;

        if (!left->bounding_box(time0, time1, box_left) || !right->bounding_box(time0, time1, box_right))
            std::cerr << "No bounding box in bvh_node constructor.\n";

        box = surrounding_box(box_left, box_right);
    }

    bvh_node(const object_group& list, double time0, double time1)
        : bvh_node(list.group, 0, list.group.size(), time0, time1) {
    }

    virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override {
        if (!box.hit(r, t_min, t_max))
            return false;

        bool hit_left = left->hit(r, t_min, t_max, rec);
        // 如果和左子树中有交点，要更新t_max，使得最终的交点是最近的交点
        bool hit_right = right->hit(r, t_min, hit_left ? rec.t : t_max, rec);

        return hit_left || hit_right;
    }

    virtual bool bounding_box(double time0, double time1, aabb& output_box) const override {
        output_box = box;
        return true;
    }

public:
    // 孩子节点可以是其他的BVH节点也可以是物体
    std::shared_ptr<object> left;
    std::shared_ptr<object> right;
    aabb                    box;
};

inline bool box_compare(const shared_ptr<object> a, const shared_ptr<object> b, int axis) {
    aabb box_a;
    aabb box_b;

    if (!a->bounding_box(0, 0, box_a) || !b->bounding_box(0, 0, box_b))
        std::cerr << "No bounding box in bvh_node constructor.\n";

    return box_a.min()[axis] < box_b.min()[axis];
}

bool box_x_compare(const shared_ptr<object> a, const shared_ptr<object> b) {
    return box_compare(a, b, 0);
}

bool box_y_compare(const shared_ptr<object> a, const shared_ptr<object> b) {
    return box_compare(a, b, 1);
}

bool box_z_compare(const shared_ptr<object> a, const shared_ptr<object> b) {
    return box_compare(a, b, 2);
}
}  // namespace ZrRender