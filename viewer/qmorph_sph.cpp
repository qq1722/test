#include "qmorph_sph.h"
#include "cgal_mesh_generator_sph.h"
#include <glm/glm.hpp>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <set>
#include <cmath>

// 辅助函数：将弧度转换为角度
inline float radians_to_degrees(float radians) {
    return radians * 180.0f / 3.14159265358979323846f;
}

// --- 辅助函数：创建排序后的边，用于作为map的键 ---
Qmorph::Edge Qmorph::make_sorted_edge(Vert_idx v1, Vert_idx v2) {
    if (v1 > v2) std::swap(v1, v2);
    return { v1, v2 };
}

// --- 主执行函数 ---
Qmorph::Result Qmorph::run(const CGALMeshGenerator& delaunay_mesh) {
    const auto& initial_triangles = delaunay_mesh.get_triangles();
    // 创建一个可修改的顶点副本（用于平滑）
    auto vertices = delaunay_mesh.get_vertices();

    if (initial_triangles.empty()) {
        std::cout << "Qmorph: No triangles to convert." << std::endl;
        return {};
    }

    // --- Stage 0: 初始化 ---
    result_.quads.clear();
    result_.remaining_triangles.clear();
    merged_triangles_.assign(initial_triangles.size(), false);

    // --- Stage 1: 初始对合并 ---
    std::cout << "Qmorph Stage 1: Starting initial pair merging..." << std::endl;
    build_adjacency(initial_triangles);
    initial_pair_merging(vertices, initial_triangles);
    std::cout << "Qmorph Stage 1: Done. Quads created: " << result_.quads.size() << std::endl;

    // --- Stage 2: 迭代清理 ---
    // 这一步是可选的，取决于具体需求
    // iterative_cleanup(vertices, initial_triangles); 

    // --- Stage 3: 收集最终结果 ---
    for (size_t i = 0; i < initial_triangles.size(); ++i) {
        if (!merged_triangles_[i]) {
            result_.remaining_triangles.push_back(initial_triangles[i]);
        }
    }

    std::cout << "Qmorph Conversion complete: " << result_.quads.size() << " quads, "
        << result_.remaining_triangles.size() << " triangles remaining." << std::endl;

    return result_;
}


// --- 阶段 1: 构建邻接关系，建立三角形和共享边的映射 ---
void Qmorph::build_adjacency(const std::vector<CGALMeshGenerator::Triangle>& triangles) {
    edge_to_tri_map_.clear();
    adj_list_.assign(triangles.size(), Tri_adj{});

    for (Tri_idx i = 0; i < triangles.size(); ++i) {
        const auto& t = triangles[i];
        Vert_idx v[3] = { t.v0, t.v1, t.v2 };
        for (int j = 0; j < 3; ++j) {
            Edge edge = make_sorted_edge(v[j], v[(j + 1) % 3]);
            adj_list_[i].edges[j] = edge;
            edge_to_tri_map_[edge].push_back(i);
        }
    }

    for (Tri_idx i = 0; i < triangles.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
            const auto& edge = adj_list_[i].edges[j];
            const auto& tri_pair = edge_to_tri_map_[edge];
            if (tri_pair.size() == 2) {
                adj_list_[i].neighbors[j] = (tri_pair[0] == i) ? tri_pair[1] : tri_pair[0];
            }
            // else { neighbor remains SIZE_MAX, indicating a boundary edge }
        }
    }
}


// --- 阶段 2: 初始对合并，将相邻三角形对转换为四边形 ---
void Qmorph::initial_pair_merging(const std::vector<glm::vec2>& vertices, const std::vector<CGALMeshGenerator::Triangle>& triangles) {
    for (const auto& pair : edge_to_tri_map_) {
        const auto& edge = pair.first;
        const auto& tri_indices = pair.second;

        if (tri_indices.size() == 2) { // 内部边
            Tri_idx idx1 = tri_indices[0];
            Tri_idx idx2 = tri_indices[1];

            if (merged_triangles_[idx1] || merged_triangles_[idx2]) {
                continue; // 其中一个已经被合并
            }

            const auto& t1 = triangles[idx1];
            const auto& t2 = triangles[idx2];

            // 找到两个三角形的另外两个顶点
            Vert_idx other_v1 = (t1.v0 != edge.first && t1.v0 != edge.second) ? t1.v0 : ((t1.v1 != edge.first && t1.v1 != edge.second) ? t1.v1 : t1.v2);
            Vert_idx other_v2 = (t2.v0 != edge.first && t2.v0 != edge.second) ? t2.v0 : ((t2.v1 != edge.first && t2.v1 != edge.second) ? t2.v1 : t2.v2);

            // 计算质量
            float quality = calculate_quad_quality(vertices[other_v1], vertices[edge.first], vertices[other_v2], vertices[edge.second]);

            // 阈值可以根据需要调整，0.5是一个比较宽松的值
            if (quality > 0.5f) {
                result_.quads.push_back({ other_v1, edge.first, other_v2, edge.second });
                merged_triangles_[idx1] = true;
                merged_triangles_[idx2] = true;
            }
        }
    }
}


// --- 质量评估函数 (可以根据需要变得更复杂) ---
float Qmorph::calculate_quad_quality(const glm::vec2& p0, const glm::vec2& p1, const glm::vec2& p2, const glm::vec2& p3) {
    // 计算四个内角
    auto angle = [](const glm::vec2& prev, const glm::vec2& curr, const glm::vec2& next) {
        glm::vec2 v1 = glm::normalize(prev - curr);
        glm::vec2 v2 = glm::normalize(next - curr);
        return std::acos(glm::clamp(glm::dot(v1, v2), -1.0f, 1.0f));
        };

    float angles[4];
    angles[0] = angle(p3, p0, p1);
    angles[1] = angle(p0, p1, p2);
    angles[2] = angle(p1, p2, p3);
    angles[3] = angle(p2, p3, p0);

    float min_angle_deg = 180.0f, max_angle_deg = 0.0f;
    for (int i = 0; i < 4; ++i) {
        float deg = radians_to_degrees(angles[i]);
        if (std::isnan(deg)) return 0.0f; // 无效角度
        min_angle_deg = std::min(min_angle_deg, deg);
        max_angle_deg = std::max(max_angle_deg, deg);
    }

    // 一个简单的评估：角度越接近90度越好
    float deviation = 0.0f;
    for (int i = 0; i < 4; ++i) {
        deviation += std::abs(radians_to_degrees(angles[i]) - 90.0f);
    }

    // 返回一个0到1之间的值，1表示最好
    return 1.0f - (deviation / 360.0f);
}

