#include "viewer_sph.h"
#include "viewer.h"
#include <glm/glm.hpp>
#include <algorithm>
#include <cmath>
#include <unordered_set>
#include <map>
#include <set>
#include <cstdlib>
#include <cfloat>

#include "../xatlas/xatlas.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
// 暂时注释掉SPHMesh1.3-master的依赖，避免链接错误
// #include "../SPHMesh1.3-master/Boundary.h"
// #include "../SPHMesh1.3-master/Simulation2D.h"
// #include "../SPHMesh1.3-master/MeshGenerator2D.h"
// #include "../SPHMesh1.3-master/halfedge/HalfedgeBuilder.hpp"
// #include "../SPHMesh1.3-master/halfedge/HalfedgeDS.hpp"

namespace sph_integration {

SPHMeshIntegrator::SPHMeshIntegrator() 
    : simulation_ready_(false) {
}

SPHMeshIntegrator::~SPHMeshIntegrator() = default;

void SPHMeshIntegrator::createBoundaryFromChart(const ChartData& chart_data) {
    current_chart_ = chart_data;
    simulation_ready_ = true;
    particles_.clear();
    particle_positions_for_render_.clear();
    printf("Created boundary with %zu points\n", chart_data.boundary_points.size());
}

// 辅助函数：判断点是否在多边形内（射线法）
static bool is_point_in_polygon(const glm::vec2& point, const std::vector<glm::vec2>& polygon) {
    if (polygon.empty()) {
        return false;
    }
    int num_vertices = polygon.size();
    bool inside = false;

    for (int i = 0, j = num_vertices - 1; i < num_vertices; j = i++) {
        const auto& p1 = polygon[i];
        const auto& p2 = polygon[j];
        float dy = p2.y - p1.y;
        if (std::abs(dy) < 1e-9f) continue;

        if (((p1.y > point.y) != (p2.y > point.y)) &&
            (point.x < (p2.x - p1.x) * (point.y - p1.y) / dy + p1.x)) {
            inside = !inside;
        }
    }
    return inside;
}

// 存储边界环和洞的结构
struct BoundaryWithHoles {
    std::vector<glm::vec2> outer_loop;
    std::vector<std::vector<glm::vec2>> holes;
};

// 全局缓存：存储每个chart的内环（洞）列表
static std::vector<std::vector<std::vector<glm::vec2>>> g_chart_holes_cache;

// 从xatlas全局UV网格提取边界环的核心函数
// 按chart分组提取每个chart的边界
std::vector<std::vector<glm::vec2>> SPHMeshIntegrator::extractBoundaryLoopsFromAtlas() {
    std::vector<std::vector<glm::vec2>> all_boundary_loops;
    g_chart_holes_cache.clear();  // 清空缓存
    
    if (!atlasIsReady()) {
        return all_boundary_loops;
    }
    
    // 获取atlas数据 - 使用viewer.h中声明的函数
    const xatlas::Atlas* atlas_data = (const xatlas::Atlas*)atlasGetData();
    if (!atlas_data || atlas_data->meshCount == 0) {
        return all_boundary_loops;
    }
    
    const xatlas::Mesh& mesh = atlas_data->meshes[0]; // 完整的输出网格
    
    printf("Mesh has %u charts\n", mesh.chartCount);
    
    // 为每个chart分别提取边界
    for (uint32_t chartIdx = 0; chartIdx < mesh.chartCount; ++chartIdx) {
        const xatlas::Chart& chart = mesh.chartArray[chartIdx];
        
        // 使用viewer的边界检测算法：基于UV坐标的边检测
        // 1. 先将chart中所有边的正向存入map
        // 2. 然后检查每条边的反向是否在map中
        // 3. 如果反向不在map中，说明这是边界边
        
        // 使用UV坐标作为边的键（使用glm::vec2比较）
        struct UVEdgeKey {
            glm::vec2 p0, p1;
            bool operator<(const UVEdgeKey& other) const {
                if (p0.x != other.p0.x) return p0.x < other.p0.x;
                if (p0.y != other.p0.y) return p0.y < other.p0.y;
                if (p1.x != other.p1.x) return p1.x < other.p1.x;
                return p1.y < other.p1.y;
            }
        };
        
        std::set<UVEdgeKey> forward_edges;  // 存储所有正向边
        
        // 第一步：将所有正向边存入set
        for (uint32_t f = 0; f < chart.faceCount; ++f) {
            uint32_t faceIdx = chart.faceArray[f];
            int v0 = mesh.indexArray[faceIdx * 3 + 0];
            int v1 = mesh.indexArray[faceIdx * 3 + 1];
            int v2 = mesh.indexArray[faceIdx * 3 + 2];
            
            const xatlas::Vertex& v0_uv = mesh.vertexArray[v0];
            const xatlas::Vertex& v1_uv = mesh.vertexArray[v1];
            const xatlas::Vertex& v2_uv = mesh.vertexArray[v2];
            
            UVEdgeKey e0, e1, e2;
            e0.p0 = glm::vec2(v0_uv.uv[0], v0_uv.uv[1]);
            e0.p1 = glm::vec2(v1_uv.uv[0], v1_uv.uv[1]);
            e1.p0 = glm::vec2(v1_uv.uv[0], v1_uv.uv[1]);
            e1.p1 = glm::vec2(v2_uv.uv[0], v2_uv.uv[1]);
            e2.p0 = glm::vec2(v2_uv.uv[0], v2_uv.uv[1]);
            e2.p1 = glm::vec2(v0_uv.uv[0], v0_uv.uv[1]);
            
            forward_edges.insert(e0);
            forward_edges.insert(e1);
            forward_edges.insert(e2);
        }
        
        // 第二步：找到边界边（反向边不在forward_edges中的边）
        std::vector<std::pair<int, int>> chart_boundary_edges;
        std::vector<std::pair<glm::vec2, glm::vec2>> chart_boundary_uv_edges;  // 存储UV坐标
        
        for (uint32_t f = 0; f < chart.faceCount; ++f) {
            uint32_t faceIdx = chart.faceArray[f];
            int v0 = mesh.indexArray[faceIdx * 3 + 0];
            int v1 = mesh.indexArray[faceIdx * 3 + 1];
            int v2 = mesh.indexArray[faceIdx * 3 + 2];
            
            const xatlas::Vertex& v0_uv = mesh.vertexArray[v0];
            const xatlas::Vertex& v1_uv = mesh.vertexArray[v1];
            const xatlas::Vertex& v2_uv = mesh.vertexArray[v2];
            
            // 检查三条边的反向
            UVEdgeKey rev_e0, rev_e1, rev_e2;
            rev_e0.p0 = glm::vec2(v1_uv.uv[0], v1_uv.uv[1]);
            rev_e0.p1 = glm::vec2(v0_uv.uv[0], v0_uv.uv[1]);
            rev_e1.p0 = glm::vec2(v2_uv.uv[0], v2_uv.uv[1]);
            rev_e1.p1 = glm::vec2(v1_uv.uv[0], v1_uv.uv[1]);
            rev_e2.p0 = glm::vec2(v0_uv.uv[0], v0_uv.uv[1]);
            rev_e2.p1 = glm::vec2(v2_uv.uv[0], v2_uv.uv[1]);
            
            if (forward_edges.find(rev_e0) == forward_edges.end()) {
                chart_boundary_edges.push_back({std::min(v0, v1), std::max(v0, v1)});
                chart_boundary_uv_edges.push_back({glm::vec2(v0_uv.uv[0], v0_uv.uv[1]), 
                                                   glm::vec2(v1_uv.uv[0], v1_uv.uv[1])});
            }
            if (forward_edges.find(rev_e1) == forward_edges.end()) {
                chart_boundary_edges.push_back({std::min(v1, v2), std::max(v1, v2)});
                chart_boundary_uv_edges.push_back({glm::vec2(v1_uv.uv[0], v1_uv.uv[1]), 
                                                   glm::vec2(v2_uv.uv[0], v2_uv.uv[1])});
            }
            if (forward_edges.find(rev_e2) == forward_edges.end()) {
                chart_boundary_edges.push_back({std::min(v2, v0), std::max(v2, v0)});
                chart_boundary_uv_edges.push_back({glm::vec2(v2_uv.uv[0], v2_uv.uv[1]), 
                                                   glm::vec2(v0_uv.uv[0], v0_uv.uv[1])});
            }
        }
        
        printf("Chart %u: Found %zu boundary edges\n", chartIdx, chart_boundary_edges.size());
        
        // 构建这个chart的所有边界环（可能包括外环和内环/洞）
        // 使用UV坐标直接构建循环，更精确
        std::vector<std::vector<glm::vec2>> chart_loops;
        std::set<size_t> used_edge_indices;  // 使用索引而不是边对
        
        // 辅助函数：计算多边形有向面积（鞋带公式）
        auto calculate_signed_area = [](const std::vector<glm::vec2>& loop) -> float {
            if (loop.size() < 3) return 0.0f;
            float area = 0.0f;
            for (size_t i = 0; i < loop.size(); ++i) {
                size_t j = (i + 1) % loop.size();
                area += loop[i].x * loop[j].y;
                area -= loop[j].x * loop[i].y;
            }
            return area * 0.5f;
        };
        
        // 辅助函数：检查两个UV点是否接近（用于连接循环）
        auto uv_points_close = [](const glm::vec2& p1, const glm::vec2& p2, float threshold = 1e-5f) -> bool {
            return glm::length(p1 - p2) < threshold;
        };
        
        // 使用UV坐标构建循环
        for (size_t start_idx = 0; start_idx < chart_boundary_uv_edges.size(); ++start_idx) {
            if (used_edge_indices.find(start_idx) != used_edge_indices.end()) continue;
            
            std::vector<glm::vec2> current_loop;
            std::set<size_t> current_used_indices;
            
            // 从当前边界边开始构建环
            const auto& start_edge = chart_boundary_uv_edges[start_idx];
            current_loop.push_back(start_edge.first);
            current_loop.push_back(start_edge.second);
            current_used_indices.insert(start_idx);
            used_edge_indices.insert(start_idx);
            
            // 尝试继续构建环
            glm::vec2 current_point = start_edge.second;
            glm::vec2 start_point = start_edge.first;
            bool found_next = true;
            int max_attempts = 1000; // 防止无限循环
            int attempts = 0;
            
            while (found_next && current_loop.size() < 1000 && attempts < max_attempts) {
                found_next = false;
                attempts++;
                
                // 检查循环是否已闭合（当前点接近起始点）
                if (uv_points_close(current_point, start_point) && current_loop.size() >= 2) {
                    // 循环已闭合，确保最后一个点等于第一个点
                    if (!uv_points_close(current_loop.back(), current_loop.front())) {
                        current_loop.push_back(current_loop.front());
                    }
                    break;  // 循环已闭合，退出
                }
                
                // 查找连接到当前点的边
                for (size_t i = 0; i < chart_boundary_uv_edges.size(); ++i) {
                    if (current_used_indices.find(i) != current_used_indices.end()) continue;
                    
                    const auto& edge = chart_boundary_uv_edges[i];
                    // 检查边的起点是否与当前点接近
                    if (uv_points_close(edge.first, current_point)) {
                        current_loop.push_back(edge.second);
                        current_used_indices.insert(i);
                        used_edge_indices.insert(i);
                        current_point = edge.second;
                        found_next = true;
                        break;
                    } else if (uv_points_close(edge.second, current_point)) {
                        // 边的方向相反
                        current_loop.push_back(edge.first);
                        current_used_indices.insert(i);
                        used_edge_indices.insert(i);
                        current_point = edge.first;
                        found_next = true;
                        break;
                    }
                }
            }
            
            // 只添加有效的闭合环（至少4个点，且是闭合的）
            // 注意：至少4个点是因为需要有3个不同的顶点（起点和终点相同）
            if (current_loop.size() >= 4) {
                // 检查循环是否闭合（第一个点和最后一个点应该相同或非常接近）
                const glm::vec2& first = current_loop.front();
                const glm::vec2& last = current_loop.back();
                float dist = glm::length(first - last);
                
                // 如果循环未闭合，尝试闭合它（如果当前点接近起始点）
                if (dist > 1e-5f && uv_points_close(current_point, start_point)) {
                    // 循环在点层面闭合了，但UV坐标可能略有不同，添加起点来显式闭合
                    current_loop.push_back(first);
                }
                
                // 再次检查是否闭合
                const glm::vec2& final_first = current_loop.front();
                const glm::vec2& final_last = current_loop.back();
                float final_dist = glm::length(final_first - final_last);
                
                if (final_dist < 1e-5f || current_loop.size() >= 4) {
                    float area = std::abs(calculate_signed_area(current_loop));
                    
                    // 过滤策略：只保留较大的循环（外环）
                    // 对于泰迪这样的复杂模型，真正的chart边界应该是最大的循环
                    // 我们暂时保留所有循环，然后在后面选择最大的作为外环
                    chart_loops.push_back(current_loop);
                    printf("  Found loop for chart %u with %zu points, area: %.2f\n", 
                           chartIdx, current_loop.size(), area);
                }
            }
        }
        
        // 对于每个chart，找到外环（面积最大的）和所有内环（洞）
        if (!chart_loops.empty()) {
            // 按面积排序，找到外环和所有内环
            std::vector<std::pair<float, std::vector<glm::vec2>*>> loops_with_area;
            for (auto& loop : chart_loops) {
                float area = std::abs(calculate_signed_area(loop));
                loops_with_area.push_back({area, &loop});
            }
            
            // 按面积从大到小排序
            std::sort(loops_with_area.begin(), loops_with_area.end(), 
                     [](const auto& a, const auto& b) { return a.first > b.first; });
            
            // 外环是面积最大的
            std::vector<glm::vec2> outer_loop = *loops_with_area[0].second;
            std::vector<std::vector<glm::vec2>> holes;
            
            // 所有其他环都是内环（洞）
            printf("  Checking %zu loops for holes (outer loop area: %.2f)\n", 
                   loops_with_area.size() - 1, loops_with_area[0].first);
            
            for (size_t i = 1; i < loops_with_area.size(); ++i) {
                // 检查内环是否在外环内（简单检查：内环的中心点是否在外环内）
                const auto& inner_loop = *loops_with_area[i].second;
                if (inner_loop.empty()) {
                    printf("  Warning: Loop %zu for chart %u is empty, skipping\n", i, chartIdx);
                    continue;
                }
                
                // 计算内环的中心点
                glm::vec2 center(0.0f);
                for (const auto& p : inner_loop) {
                    center += p;
                }
                center /= float(inner_loop.size());
                
                // 检查中心点是否在外环内
                bool center_in_outer = is_point_in_polygon(center, outer_loop);
                printf("  Loop %zu for chart %u: center (%.3f, %.3f), area: %.2f, center_in_outer: %d\n",
                       i, chartIdx, center.x, center.y, loops_with_area[i].first, center_in_outer ? 1 : 0);
                
                if (center_in_outer) {
                    holes.push_back(inner_loop);
                    printf("  -> Added as hole #%zu for chart %u with %zu points\n", 
                           holes.size() - 1, chartIdx, inner_loop.size());
                } else {
                    printf("  -> Skipping loop %zu (center outside outer loop)\n", i);
                }
            }
            
            // 存储外环（兼容现有接口）
            all_boundary_loops.push_back(outer_loop);
            
            // 存储内环信息到全局缓存（索引对应chart索引）
            g_chart_holes_cache.push_back(holes);
            
            printf("  Chart %u: outer loop with %zu points, %zu holes\n", 
                   chartIdx, outer_loop.size(), holes.size());
            
            // 调试：打印边界框
            float minX = outer_loop.front().x, maxX = outer_loop.front().x;
            float minY = outer_loop.front().y, maxY = outer_loop.front().y;
            for (const auto& p : outer_loop) {
                minX = std::min(minX, p.x); maxX = std::max(maxX, p.x);
                minY = std::min(minY, p.y); maxY = std::max(maxY, p.y);
            }
            printf("    Bounding box: (%.1f,%.1f) to (%.1f,%.1f), size: %.1f x %.1f\n",
                   minX, minY, maxX, maxY, maxX - minX, maxY - minY);
        } else {
            printf("  Warning: Chart %u has no valid boundary loops!\n", chartIdx);
        }
    }
    
    printf("Total boundary loops found: %zu (should be %u)\n", all_boundary_loops.size(), mesh.chartCount);
    return all_boundary_loops;
}

// 获取指定chart的内环（洞）列表
std::vector<std::vector<glm::vec2>> SPHMeshIntegrator::getChartHoles(uint32_t chartIndex) {
    if (chartIndex < g_chart_holes_cache.size()) {
        return g_chart_holes_cache[chartIndex];
    }
    return {};
}

// 更新指定chart的内环（洞）列表（用于归一化后更新）
void SPHMeshIntegrator::updateChartHoles(uint32_t chartIndex, const std::vector<std::vector<glm::vec2>>& holes) {
    if (chartIndex < g_chart_holes_cache.size()) {
        g_chart_holes_cache[chartIndex] = holes;
    }
}

void SPHMeshIntegrator::initializeParticles() {
    particles_.clear();
    
    printf("initializeParticles: boundary_points.size() = %zu\n", current_chart_.boundary_points.size());
    
    if (current_chart_.boundary_points.empty()) {
        printf("initializeParticles: boundary_points is empty, returning\n");
        return;
    }
    
    // 计算边界框
    float minX = current_chart_.boundary_points[0].x;
    float maxX = current_chart_.boundary_points[0].x;
    float minY = current_chart_.boundary_points[0].y;
    float maxY = current_chart_.boundary_points[0].y;
    
    for (const auto& p : current_chart_.boundary_points) {
        minX = std::min(minX, p.x);
        maxX = std::max(maxX, p.x);
        minY = std::min(minY, p.y);
        maxY = std::max(maxY, p.y);
    }
    
    printf("initializeParticles: bounding box (%.1f, %.1f) to (%.1f, %.1f), size %.1f x %.1f\n",
           minX, minY, maxX, maxY, maxX - minX, maxY - minY);
    
    // 检测坐标是否归一化（在[0,1]范围内）
    bool isNormalized = (maxX <= 1.1f && maxY <= 1.1f && minX >= -0.1f && minY >= -0.1f);
    
    // 根据坐标系统调整平滑半径
    float effective_radius = smoothing_radius_;
    if (isNormalized) {
        // 归一化坐标，使用小的平滑半径（相对于1.0）
        effective_radius = 0.015f; // 大约1.5%的边界框大小
    } else {
        // 像素坐标，使用像素单位
        effective_radius = smoothing_radius_;
    }
    
    // 在边界内生成均匀分布的粒子
    const float spacing = effective_radius * 0.5f;
    const float jitter = spacing * 0.1f;
    
    printf("initializeParticles: isNormalized=%d, effective_radius=%.4f, spacing=%.4f, jitter=%.4f\n",
           isNormalized ? 1 : 0, effective_radius, spacing, jitter);
    
    int total_attempts = 0;
    int inside_count = 0;
    
    // 确保循环能执行：如果spacing太大，至少尝试一些点
    float startY = minY + spacing * 0.5f;
    float endY = maxY - spacing * 0.5f;
    float startX = minX + spacing * 0.5f;
    float endX = maxX - spacing * 0.5f;
    
    // 如果spacing太大导致没有点，使用更小的步长
    float actual_spacing = spacing;
    float actual_jitter = jitter;
    if (startY > endY || startX > endX) {
        // 使用边界框大小的1/20作为间距
        actual_spacing = std::min(maxX - minX, maxY - minY) / 20.0f;
        actual_jitter = actual_spacing * 0.1f;
        startY = minY + actual_spacing * 0.5f;
        endY = maxY - actual_spacing * 0.5f;
        startX = minX + actual_spacing * 0.5f;
        endX = maxX - actual_spacing * 0.5f;
        printf("initializeParticles: spacing too large, using adaptive spacing=%.4f\n", actual_spacing);
    }
    
    for (float y = startY; y <= endY; y += actual_spacing) {
        for (float x = startX; x <= endX; x += actual_spacing) {
            total_attempts++;
            glm::vec2 pos(
                x + ((float)rand() / RAND_MAX - 0.5f) * actual_jitter,
                y + ((float)rand() / RAND_MAX - 0.5f) * actual_jitter
            );
            
            if (isPointInsideBoundary(pos)) {
                inside_count++;
                Particle p;
                p.position = pos;
                p.velocity = glm::vec2(0.0f);
                p.force = glm::vec2(0.0f);
                p.density = rest_density_;
                p.pressure = 0.0f;
                particles_.push_back(p);
            }
        }
    }
    
    printf("initializeParticles: total_attempts=%d, inside_count=%d, final_particles=%zu\n",
           total_attempts, inside_count, particles_.size());
    
    particle_positions_for_render_.resize(particles_.size());
    // 初始化渲染位置
    for (size_t i = 0; i < particles_.size(); i++) {
        particle_positions_for_render_[i] = particles_[i].position;
    }
    if (!particles_.empty()) {
        printf("Initialized %zu particles, first at (%.1f, %.1f), last at (%.1f, %.1f)\n", 
               particles_.size(),
               particles_[0].position.x, particles_[0].position.y,
               particles_.back().position.x, particles_.back().position.y);
    } else {
        printf("Warning: initializeParticles generated 0 particles!\n");
    }
}

bool SPHMeshIntegrator::isPointInsideBoundary(const glm::vec2& point) const {
    if (current_chart_.boundary_points.empty()) {
        return false;
    }
    
    // 使用射线法判断点是否在多边形内
    bool inside = false;
    size_t n = current_chart_.boundary_points.size();
    
    // 处理边界情况：如果边界点太少，直接返回false
    if (n < 3) {
        return false;
    }
    
    for (size_t i = 0, j = n - 1; i < n; j = i++) {
        const glm::vec2& pi = current_chart_.boundary_points[i];
        const glm::vec2& pj = current_chart_.boundary_points[j];
        
        // 检查除零
        if (fabs(pj.y - pi.y) < 1e-6f) {
            continue;
        }
        
        if (((pi.y > point.y) != (pj.y > point.y)) &&
            (point.x < (pj.x - pi.x) * (point.y - pi.y) / (pj.y - pi.y) + pi.x)) {
            inside = !inside;
        }
    }
    
    return inside;
}

float SPHMeshIntegrator::distanceToBoundary(const glm::vec2& point) const {
    if (current_chart_.boundary_points.empty()) {
        return 1000.0f;
    }
    
    float minDist = 1000.0f;
    size_t n = current_chart_.boundary_points.size();
    
    for (size_t i = 0; i < n; i++) {
        const glm::vec2& p1 = current_chart_.boundary_points[i];
        const glm::vec2& p2 = current_chart_.boundary_points[(i + 1) % n];
        
        glm::vec2 edge = p2 - p1;
        float edgeLen = glm::length(edge);
        if (edgeLen < 1e-6f) continue;
        
        glm::vec2 toPoint = point - p1;
        float t = glm::clamp(glm::dot(toPoint, edge) / (edgeLen * edgeLen), 0.0f, 1.0f);
        glm::vec2 closest = p1 + edge * t;
        
        float dist = glm::length(point - closest);
        minDist = std::min(minDist, dist);
    }
    
    return minDist;
}

float SPHMeshIntegrator::sphKernel(float distance, float smoothingRadius) const {
    float q = distance / smoothingRadius;
    if (q >= 1.0f) return 0.0f;
    
    float factor = 1.0f - q;
    return factor * factor * (2.0f * M_PI * smoothingRadius * smoothingRadius);
}

glm::vec2 SPHMeshIntegrator::sphKernelGradient(const glm::vec2& r, float smoothingRadius) const {
    float dist = glm::length(r);
    if (dist < 1e-6f) return glm::vec2(0.0f);
    
    float q = dist / smoothingRadius;
    if (q >= 1.0f) return glm::vec2(0.0f);
    
    float factor = 1.0f - q;
    float gradientMag = -2.0f * factor / smoothingRadius;
    float coeff = gradientMag * (2.0f * M_PI * smoothingRadius * smoothingRadius);
    return (r / dist) * coeff;
}

void SPHMeshIntegrator::computeDensity() {
    const float mass = 1.0f;
    
    for (size_t i = 0; i < particles_.size(); i++) {
        float density = 0.0f;
        
        for (size_t j = 0; j < particles_.size(); j++) {
            glm::vec2 r = particles_[i].position - particles_[j].position;
            float dist = glm::length(r);
            density += mass * sphKernel(dist, smoothing_radius_);
        }
        
        particles_[i].density = density;
    }
}

void SPHMeshIntegrator::computePressure() {
    for (auto& p : particles_) {
        p.pressure = gas_constant_ * (p.density - rest_density_);
    }
}

void SPHMeshIntegrator::computeForces() {
    const float mass = 1.0f;
    
    for (size_t i = 0; i < particles_.size(); i++) {
        glm::vec2 pressureForce(0.0f);
        glm::vec2 viscosityForce(0.0f);
        
        for (size_t j = 0; j < particles_.size(); j++) {
            if (i == j) continue;
            
            glm::vec2 r = particles_[i].position - particles_[j].position;
            float dist = glm::length(r);
            
            if (dist < smoothing_radius_ && dist > 1e-6f) {
                // 压力力
                float avgPressure = (particles_[i].pressure + particles_[j].pressure) * 0.5f;
                glm::vec2 gradW = sphKernelGradient(r, smoothing_radius_);
                pressureForce -= mass * avgPressure / particles_[j].density * gradW;
                
                // 粘性力
                glm::vec2 velDiff = particles_[j].velocity - particles_[i].velocity;
                viscosityForce += mass * viscosity_ * velDiff / particles_[j].density * sphKernel(dist, smoothing_radius_);
            }
        }
        
        particles_[i].force = pressureForce + viscosityForce;
    }
}

void SPHMeshIntegrator::updateParticles() {
    for (auto& p : particles_) {
        p.velocity += p.force * time_step_;
        p.velocity *= damping_;
        p.position += p.velocity * time_step_;
    }
}

void SPHMeshIntegrator::handleBoundaries() {
    for (auto& p : particles_) {
        if (!isPointInsideBoundary(p.position)) {
            float dist = distanceToBoundary(p.position);
            glm::vec2 normal(0.0f);
            
            // 找到最近的边界点
            float minDist = 1000.0f;
            size_t n = current_chart_.boundary_points.size();
            for (size_t i = 0; i < n; i++) {
                const glm::vec2& p1 = current_chart_.boundary_points[i];
                const glm::vec2& p2 = current_chart_.boundary_points[(i + 1) % n];
                
                glm::vec2 edge = p2 - p1;
                float edgeLen = glm::length(edge);
                if (edgeLen < 1e-6f) continue;
                
                glm::vec2 toPoint = p.position - p1;
                float t = glm::clamp(glm::dot(toPoint, edge) / (edgeLen * edgeLen), 0.0f, 1.0f);
                glm::vec2 closest = p1 + edge * t;
                
                float d = glm::length(p.position - closest);
                if (d < minDist) {
                    minDist = d;
                    normal = glm::normalize(p.position - closest);
                }
            }
            
            // 将粒子推回边界内
            p.position += normal * (minDist + smoothing_radius_ * 0.1f);
            p.velocity -= 2.0f * glm::dot(p.velocity, normal) * normal;
            p.velocity *= 0.5f; // 阻尼
        }
    }
}

void SPHMeshIntegrator::stepSimulation() {
    if (!simulation_ready_ || particles_.empty()) {
        printf("stepSimulation: not ready or no particles (ready=%d, particles=%zu)\n", 
               simulation_ready_, particles_.size());
        return;
    }
    
    computeDensity();
    computePressure();
    computeForces();
    updateParticles();
    handleBoundaries();
    
    // 更新渲染位置
    for (size_t i = 0; i < particles_.size(); i++) {
        particle_positions_for_render_[i] = particles_[i].position;
    }
}

void SPHMeshIntegrator::runSimulation(int steps) {
    if (!simulation_ready_ || particles_.empty()) {
        if (particles_.empty() && simulation_ready_) {
            initializeParticles();
        }
        if (particles_.empty()) {
            return;
        }
    }
    
    for (int i = 0; i < steps; i++) {
        stepSimulation();
    }
    
    printf("Ran simulation for %d steps, %zu particles\n", steps, particles_.size());
}

const std::vector<glm::vec2>& SPHMeshIntegrator::getParticlePositions() const {
    return particle_positions_for_render_;
}

void SPHMeshIntegrator::generateQuadMesh() {
    if (!simulation_ready_) {
        return;
    }
    
    // 暂时简化实现，避免依赖MeshGenerator2D
    printf("Quad mesh generation not implemented yet\n");
}

const std::vector<glm::vec2>& SPHMeshIntegrator::getMeshVertices() const {
    static std::vector<glm::vec2> empty;
    // 暂时返回空结果，避免依赖MeshGenerator2D
    return empty;
}

const std::vector<sph_integration::Quad>& SPHMeshIntegrator::getMeshQuads() const {
    static std::vector<sph_integration::Quad> empty;
    // 暂时返回空结果，避免依赖MeshGenerator2D
    return empty;
}

void SPHMeshIntegrator::mapTo3DModel(const std::vector<glm::vec3>& original_vertices,
                                   const std::vector<glm::vec2>& original_uvs,
                                   std::vector<glm::vec3>& output_vertices,
                                   std::vector<glm::vec2>& output_uvs) {
    // 暂时简化实现，避免依赖MeshGenerator2D
    printf("3D mapping not implemented yet\n");
    output_vertices.clear();
    output_uvs.clear();
}

} // namespace sph_integration
