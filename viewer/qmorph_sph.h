#pragma once
#include <vector>
#include <map>
#include <glm/glm.hpp>
#include "cgal_mesh_generator_sph.h"

// 为了避免循环引用，只使用前向声明
class CGALMeshGenerator;

class Qmorph {
public:
    // 转换结果的结构体
    struct Result {
        std::vector<CGALMeshGenerator::Quad> quads;
        std::vector<CGALMeshGenerator::Triangle> remaining_triangles;
    };

    Qmorph() = default;

    // 核心函数：输入一个三角网格，返回一个四边形为主的网格
    Result run(const CGALMeshGenerator& delaunay_mesh);

private:
    // --- 内部数据结构 ---
    // 使用别名来简化代码
    using Tri_idx = size_t;
    using Vert_idx = unsigned int;
    using Edge = std::pair<Vert_idx, Vert_idx>;

    // 存储每个三角形的邻接信息
    struct Tri_adj {
        // 使用 SIZE_MAX 而不是 -1 来表示无效索引
        Tri_idx neighbors[3] = { SIZE_MAX, SIZE_MAX, SIZE_MAX };
        // 邻接边通过共享边连接
        Edge edges[3];
    };

    // --- 多阶段转换算法 ---
    void build_adjacency(const std::vector<CGALMeshGenerator::Triangle>& triangles);
    void initial_pair_merging(const std::vector<glm::vec2>& vertices, const std::vector<CGALMeshGenerator::Triangle>& triangles);
    void iterative_cleanup(const std::vector<glm::vec2>& vertices, std::vector<CGALMeshGenerator::Triangle>& triangles);
    // void final_smoothing(std::vector<glm::vec2>& vertices, const Boundary& boundary); // 注意：需要边界信息

    // --- 辅助函数和工具函数 ---
    float calculate_quad_quality(const glm::vec2& p0, const glm::vec2& p1, const glm::vec2& p2, const glm::vec2& p3);
    Edge make_sorted_edge(Vert_idx v1, Vert_idx v2);

    // --- 成员变量 ---
    std::vector<bool> merged_triangles_;
    std::vector<Tri_adj> adj_list_;
    std::map<Edge, std::vector<Tri_idx>> edge_to_tri_map_;
    Result result_;
};

