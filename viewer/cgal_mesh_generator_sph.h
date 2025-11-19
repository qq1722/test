#pragma once
#include <vector>
#include <glm/glm.hpp>
#include "simulation2d_sph.h"
#include "boundary_sph.h"

// 自动检测是否可用 CGAL（若未显式定义 HAVE_CGAL）
#if !defined(HAVE_CGAL)
#  if defined(__has_include)
#    if __has_include(<CGAL/Exact_predicates_inexact_constructions_kernel.h>)
#      define HAVE_CGAL 1
#    endif
#  endif
#endif

// 检查是否有CGAL库
#ifdef HAVE_CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>

// 面信息结构
struct FaceInfo2 {
    bool in_domain() const { return nesting_level % 2 == 1; }
    int nesting_level = 0;
};

// CGAL 内核
using K = CGAL::Exact_predicates_inexact_constructions_kernel;

// 创建一个带约束的面基类：可以存储约束边的信息，也可以存储自定义面信息
using Cfb = CGAL::Constrained_triangulation_face_base_2<K>;
using Fbb = CGAL::Triangulation_face_base_with_info_2<FaceInfo2, K, Cfb>;

// 顶点数据结构
using Vb = CGAL::Triangulation_vertex_base_2<K>;
using Tds = CGAL::Triangulation_data_structure_2<Vb, Fbb>;

// 最终的约束Delaunay三角剖分数据结构
using CDT = CGAL::Constrained_Delaunay_triangulation_2<K, Tds>;
using Point = CDT::Point;
#endif

class CGALMeshGenerator {
public:
    struct Triangle {
        unsigned int v0, v1, v2;
    };
    struct Quad {
        unsigned int v0, v1, v2, v3;
    };

    CGALMeshGenerator() = default;
    ~CGALMeshGenerator() = default;

    // 从粒子和边界生成Delaunay三角网格
    void generate_mesh(const std::vector<Simulation2D::Particle>& particles, const Boundary& boundary);

    const std::vector<glm::vec2>& get_vertices() const { return vertices_; }
    const std::vector<Triangle>& get_triangles() const { return triangles_; }
    const std::vector<Quad>& get_quads() const { return quads_; }

private:
    std::vector<glm::vec2> vertices_;
    std::vector<Triangle> triangles_;
    std::vector<Quad> quads_;
};

