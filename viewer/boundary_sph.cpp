#include "boundary_sph.h"
#include <algorithm>
#include <cfloat>

// Helper functions for finding closest points
glm::vec2 closest_point_on_segment(const glm::vec2& p, const glm::vec2& a, const glm::vec2& b) {
    glm::vec2 ab = b - a;
    glm::vec2 ap = p - a;
    float proj = glm::dot(ap, ab);
    float ab_len_sq = glm::dot(ab, ab);
    if (ab_len_sq < 1e-9f) return a;
    float d = proj / ab_len_sq;
    if (d <= 0.0f) return a;
    if (d >= 1.0f) return b;
    return a + d * ab;
}

glm::vec2 closest_point_on_polygon(const glm::vec2& p, const std::vector<glm::vec2>& vertices) {
    if (vertices.empty()) return p;
    glm::vec2 closest_point = vertices[0];
    float min_dist_sq = FLT_MAX;
    for (size_t i = 0, j = vertices.size() - 1; i < vertices.size(); j = i++) {
        glm::vec2 closest_pt_on_edge = closest_point_on_segment(p, vertices[j], vertices[i]);
        float dist_sq = glm::dot(p - closest_pt_on_edge, p - closest_pt_on_edge);
        if (dist_sq < min_dist_sq) {
            min_dist_sq = dist_sq;
            closest_point = closest_pt_on_edge;
        }
    }
    return closest_point;
}

// 构造函数：支持外环和内环
// 如果只提供外环，holes使用默认值（空列表）
Boundary::Boundary(const std::vector<glm::vec2>& outer_loop, 
                   const std::vector<std::vector<glm::vec2>>& holes) 
    : vertices_(outer_loop), holes_(holes)
{
    calculate_aabb();
}

const std::vector<glm::vec2>& Boundary::get_vertices() const 
{
    return vertices_;
}

const glm::vec4& Boundary::get_aabb() const 
{
    return aabb_;
}

void Boundary::calculate_aabb() 
{
    if (vertices_.empty()) {
        aabb_ = glm::vec4(0.0f, 0.0f, 0.0f, 0.0f);
        return;
    }

    glm::vec2 min_coords = vertices_[0];
    glm::vec2 max_coords = vertices_[0];

    // 计算外环的AABB
    for (size_t i = 1; i < vertices_.size(); ++i) 
    {
        min_coords.x = std::min(min_coords.x, vertices_[i].x);
        min_coords.y = std::min(min_coords.y, vertices_[i].y);
        max_coords.x = std::max(max_coords.x, vertices_[i].x);
        max_coords.y = std::max(max_coords.y, vertices_[i].y);
    }
    
    // 内环不影响AABB（它们在外环内）
    aabb_ = glm::vec4(min_coords.x, min_coords.y, max_coords.x, max_coords.y);
}

// 辅助函数：判断点是否在多边形内（射线法）
bool Boundary::is_point_in_polygon(const glm::vec2& point, const std::vector<glm::vec2>& polygon) const
{
    if (polygon.empty()) {
        return false;
    }
    int num_vertices = polygon.size();
    bool inside = false;

    for (int i = 0, j = num_vertices - 1; i < num_vertices; j = i++) 
    {
        const auto& p1 = polygon[i];
        const auto& p2 = polygon[j];
        float dy = p2.y - p1.y;
        if (std::abs(dy) < 1e-9f) continue;

        if (((p1.y > point.y) != (p2.y > point.y)) &&
            (point.x < (p2.x - p1.x) * (point.y - p1.y) / dy + p1.x)) 
        {
            inside = !inside;
        }
    }
    return inside;
}

// 判断点是否在边界内：在外环内，但在所有内环外
bool Boundary::is_inside(const glm::vec2& point) const 
{
    // 首先检查是否在外环内
    if (!is_point_in_polygon(point, vertices_)) {
        return false;
    }
    
    // 然后检查是否在任何内环内（如果在内环内，则不在边界内）
    for (const auto& hole : holes_) {
        if (is_point_in_polygon(point, hole)) {
            return false;  // 点在洞内，不在边界内
        }
    }
    
    return true;  // 在外环内且不在任何洞内
}


// [新增] 实现计算最近点逻辑
glm::vec2 Boundary::get_closest_point(const glm::vec2& p) const {
    // 1. 先计算到外环的最近点
    glm::vec2 best_point = closest_point_on_polygon(p, vertices_);
    float min_dist_sq = glm::dot(p - best_point, p - best_point);

    // 2. 遍历所有内洞，看有没有更近的
    for (const auto& hole : holes_) {
        if (hole.empty()) continue;
        glm::vec2 hole_pt = closest_point_on_polygon(p, hole);
        float dist_sq = glm::dot(p - hole_pt, p - hole_pt);

        if (dist_sq < min_dist_sq) {
            min_dist_sq = dist_sq;
            best_point = hole_pt;
        }
    }

    return best_point;
}

// [新增] 核心算法：同时寻找最近点和切线
void Boundary::get_closest_point_and_tangent(const glm::vec2& p, glm::vec2& out_closest, glm::vec2& out_tangent) const
{
    float min_dist_sq = FLT_MAX;

    // 定义一个 Lambda 函数来处理单个环（外环或内洞）
    auto process_ring = [&](const std::vector<glm::vec2>& ring) {
        if (ring.empty()) return;
        int n = (int)ring.size();
        for (int i = 0; i < n; ++i) {
            glm::vec2 a = ring[i];
            glm::vec2 b = ring[(i + 1) % n]; // 循环连接到下一个点

            glm::vec2 ab = b - a;
            float len_sq = glm::dot(ab, ab);
            if (len_sq < 1e-9f) continue; // 忽略退化边

            // 计算点 p 在线段 ab 上的投影比例 t
            float t = glm::dot(p - a, ab) / len_sq;
            t = std::max(0.0f, std::min(1.0f, t)); // 限制在线段范围内

            glm::vec2 pt = a + t * ab;
            float dist_sq = glm::dot(p - pt, p - pt);

            if (dist_sq < min_dist_sq) {
                min_dist_sq = dist_sq;
                out_closest = pt;
                // 切线就是边的方向 (归一化)
                out_tangent = glm::normalize(ab);
            }
        }
        };

    // 1. 检查外环
    process_ring(vertices_);

    // 2. 检查所有内洞
    for (const auto& hole : holes_) {
        process_ring(hole);
    }
}

// get_outer_boundary() 的实现
const std::vector<glm::vec2>& Boundary::get_outer_boundary() const
{
    return vertices_;
}

