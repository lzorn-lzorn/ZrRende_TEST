#pragma once

#include "../include/utilities.hpp"
#include "../include/random.hpp"

namespace ZrRender{

class perlin {
public:
    perlin() {
        // 生成随机梯度向量
        ranvec = new vec3[point_count];
        for (int i = 0; i < point_count; ++i) {
            ranvec[i] = normalize(random_vec(-1, 1));
        }
        // 生成三个维度的P[n]数组
        perm_x = perlin_generate_perm();
        perm_y = perlin_generate_perm();
        perm_z = perlin_generate_perm();
    }

    ~perlin() {
        delete[] ranvec;
        delete[] perm_x;
        delete[] perm_y;
        delete[] perm_z;
    }

    // Perlin 噪声函数，给定三维空间一点，返回噪声值
    double noise(const point3& p) const {
        // 三个维度坐标的小数部分相当于到一个单位立方体各个维度的距离
        auto u = p.x - floor(p.x);
        auto v = p.y - floor(p.y);
        auto w = p.z - floor(p.z);
        // 三个维度坐标的整数部分，用于到P[n]数组中找到查找表下标
        auto i = static_cast<int>(floor(p.x));
        auto j = static_cast<int>(floor(p.y));
        auto k = static_cast<int>(floor(p.z));
        // 存储离该点最近的八个晶格顶点的随机梯度向量
        vec3 c[2][2][2];
        // 到查找表中查找八个随机梯度向量
        for (int di = 0; di < 2; di++)
            for (int dj = 0; dj < 2; dj++)
                for (int dk = 0; dk < 2; dk++)
                    // 异或是不进位加法，相当于把三个P[n]数组中的下标加起来并限制在0~255之间，避免下标越界
                    c[di][dj][dk] = ranvec[
                        perm_x[(i + di) & 255] ^
                            perm_y[(j + dj) & 255] ^
                            perm_z[(k + dk) & 255]
                    ];
        // 用八个梯度向量和距离向量点乘并线性插值
        return perlin_interp(c, u, v, w);
    }
    // 多个noise组合形成turbulence
    double turb(const point3& p, int depth = 7) const {
        auto accum = 0.0;
        auto temp_p = p;
        auto weight = 1.0;

        for (int i = 0; i < depth; i++) {
            accum += weight * noise(temp_p);
            weight *= 0.5;
            temp_p *= 2;
        }

        return fabs(accum);
    }
private:
    static const int point_count = 256;
    vec3* ranvec;       // 存储随机梯度向量查找表
    int* perm_x;        // x维度的P[n]数组
    int* perm_y;        // y维度的P[n]数组
    int* perm_z;        // z维度的P[n]数组

    // 生成从 0 到 point_count - 1 的乱序数组
    // 思路是先生成顺序的数组，然后使用洗牌算法打乱
    static int* perlin_generate_perm() 
    {
        auto p = new int[point_count];

        for (int i = 0; i < perlin::point_count; i++)
            p[i] = i;

        permute(p, point_count);

        return p;
    }

    // 洗牌算法打乱数组，从后向前遍历数组，每次随机挑一个下标和当前下标交换
    static void permute(int* p, int n) 
    {
        for (int i = n - 1; i > 0; i--) 
        {
            int target = random_int(0, i);
            int tmp = p[i];
            p[i] = p[target];
            p[target] = tmp;
        }
    }

    // 将八个顶点的梯度向量和点p到八个顶点的距离向量点乘，并做三线性插值
    static double perlin_interp(vec3 c[2][2][2], double u, double v, double w) 
    {
        // 经过改进的线性插值系数
        auto uu = fade(u);
        auto vv = fade(v);
        auto ww = fade(w);

        auto accum = 0.0;

        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                for (int k = 0; k < 2; k++) 
                {
                    vec3 weight_v(u - i, v - j, w - k);
                    accum += (i * uu + (1 - i) * (1 - uu))
                        * (j * vv + (1 - j) * (1 - vv))
                        * (k * ww + (1 - k) * (1 - ww))
                        * dot(c[i][j][k], weight_v);
                }

        return accum;
    }

    // 缓和曲线，用于缓和线性插值系数
    static double fade(double x)
    {
        return x * x * x * (x * (x * 6 - 15) + 10);
    }
};

}