#include "cgal_mesh_generator_sph.h"
#include <iostream>
#include <map>
#include <vector>

#ifdef HAVE_CGAL
void CGALMeshGenerator::generate_mesh(const std::vector<Simulation2D::Particle>& particles, const Boundary& boundary) {
    vertices_.clear();
    triangles_.clear();
    quads_.clear();

    const auto& boundary_vertices = boundary.get_vertices();
    if (particles.empty() || boundary_vertices.empty()) {
        std::cerr << "CGALMeshGenerator: Empty particles or boundary." << std::endl;
        return;
    }

    CDT cdt;

    // --- 步骤 1: 使用所有粒子位置作为点集（包含边界粒子），并做轻度去重 ---
    struct QuantHash {
        size_t operator()(const std::pair<int, int>& v) const noexcept {
            return (size_t(uint32_t(v.first)) << 32) ^ size_t(uint32_t(v.second));
        }
    };
    auto quantize = [](const glm::vec2& p) -> std::pair<int,int> {
        // 量化到1e-6 精度，避免重复点
        const double scale = 1e6;
        int qx = (int)llround(p.x * scale);
        int qy = (int)llround(p.y * scale);
        return { qx, qy };
    };

    std::unordered_set<std::pair<int,int>, QuantHash> seen;
    std::vector<Point> all_points;
    all_points.reserve(particles.size());
    for (const auto& p : particles) {
        auto q = quantize(p.position);
        if (seen.insert(q).second) {
            all_points.emplace_back(p.position.x, p.position.y);
        }
    }

    if (all_points.empty()) {
        std::cerr << "Warning: No valid particle positions to generate mesh from." << std::endl;
        return;
    }

    cdt.insert(all_points.begin(), all_points.end());

    // --- 步骤 2 & 3: 筛选并提取内部三角形 (使用边界判断法) ---
    // 创建map将CGAL的内部顶点映射到我们自己的顶点索引
    std::map<CDT::Vertex_handle, unsigned int> vertex_map;

    for (auto face_it = cdt.finite_faces_begin(); face_it != cdt.finite_faces_end(); ++face_it) {
        // a) 计算三角形的重心
        Point p0 = face_it->vertex(0)->point();
        Point p1 = face_it->vertex(1)->point();
        Point p2 = face_it->vertex(2)->point();
        glm::vec2 centroid((p0.x() + p1.x() + p2.x()) / 3.0f, (p0.y() + p1.y() + p2.y()) / 3.0f);

        // b) 使用边界判断三角形是否在有效区域内
        if (boundary.is_inside(centroid)) {
            unsigned int v_indices[3];
            for (int i = 0; i < 3; ++i) {
                CDT::Vertex_handle vh = face_it->vertex(i);

                // c) 如果是新顶点，添加到我们的顶点列表并记录映射
                if (vertex_map.find(vh) == vertex_map.end()) {
                    vertex_map[vh] = (unsigned int)vertices_.size();
                    vertices_.emplace_back(vh->point().x(), vh->point().y());
                }
                v_indices[i] = vertex_map[vh];
            }
            // d) 添加筛选后的三角形
            triangles_.push_back({ v_indices[0], v_indices[1], v_indices[2] });
        }
    }

    std::cout << "CGAL generated (Clip Method): " << vertices_.size() << " vertices, " << triangles_.size() << " triangles." << std::endl;
}
#else
void CGALMeshGenerator::generate_mesh(const std::vector<Simulation2D::Particle>& particles, const Boundary& boundary) {
    vertices_.clear();
    triangles_.clear();
    quads_.clear();
    std::cerr << "CGALMeshGenerator: CGAL library not available. Please install CGAL or define HAVE_CGAL." << std::endl;
}
#endif

