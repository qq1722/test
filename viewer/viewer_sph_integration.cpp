#include "viewer_sph.h"
#include "boundary_sph.h"
#include "simulation2d_sph.h"
#include "cgal_mesh_generator_sph.h"
#include "qmorph_sph.h"
#include "viewer.h"
#include "../thirdparty/imgui/imgui.h"
#include <vector>
#include <memory>
#include <limits>
#include <cmath>
#include <algorithm>
#include <glm/glm.hpp>
#include "../xatlas/xatlas.h"

#if BX_PLATFORM_WINDOWS
#include <direct.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#endif
void sphRunSimulationAll(int steps);
bool sphHasAnySimulation();
void sphSetShowAllTriOverlay(bool enabled);
void sphSetShowParticlesIn3D(bool enabled);
static void ensure_all_containers_sized();
static void ensureExportDir();
void updateBoundaryLoopsCache();
namespace {
    std::unique_ptr<sph_integration::SPHMeshIntegrator> g_sph_integrator;
    std::unique_ptr<Boundary> g_current_boundary;
    std::unique_ptr<Simulation2D> g_sph_simulation;

    // 多图表并行支持
    std::vector<std::unique_ptr<Boundary>> g_boundaries_all;                 // 与图表索引对齐
    std::vector<std::unique_ptr<Simulation2D>> g_simulations_all;            // 与图表索引对齐
    std::vector<glm::vec2> g_all_particles_cache;                            // 汇总粒子位置缓存（用于3D渲染）
    bool g_all_particles_cache_dirty = true;

    bool g_sph_enabled = false;
    uint32_t g_current_chart_index = 0;
    int g_simulation_steps = 100;
    bool g_auto_simulation = false;
    
    // 缓存边界环数据，避免重复计算
    std::vector<std::vector<glm::vec2>> g_cached_boundary_loops;
    bool g_boundary_loops_cached = false;

    // 坐标选项
    bool g_normalize_uv_01 = true;       // 归一化到[0,1]（按整张atlas尺寸）
    bool g_chart_local_01 = false;       // 归一化到图表局部[0,1]
    
    // 网格生成相关
    std::unique_ptr<CGALMeshGenerator> g_delaunay_generator;
    std::unique_ptr<Qmorph> g_qmorph;
    Qmorph::Result g_quad_mesh_result;
    bool g_tri_mesh_generated = false;
    bool g_quad_mesh_generated = false;

    // 视图切换：true 显示网格，false 显示粒子
    bool g_view_mesh = false;
    // 三维粒子渲染开关
    bool g_show_particles_in_3d = true;

    // 用户意图：若用户点击过“Generate Quad Mesh”，则切换图表时也自动尝试生成四边形网格
    bool g_request_quad_on_switch = false;

    // 3D 全局网格叠加控制
    bool g_show_all_tri_overlay = false;
    bool g_show_all_quad_overlay = false;

    // 所有图表的网格缓存（UV空间）
    std::vector<std::vector<glm::vec2>> g_all_tri_vertices;
    std::vector<std::vector<glm::uvec3>> g_all_tri_indices;
    std::vector<std::vector<glm::vec2>> g_all_quad_vertices;
    std::vector<std::vector<glm::uvec4>> g_all_quad_indices;
    std::vector<std::vector<glm::uvec3>> g_all_quad_extra_tris;

    struct CombinedTriMesh3D {
        std::vector<glm::vec3> positions;
        std::vector<glm::uvec3> triangles;
    };

    struct CombinedQuadMesh3D {
        std::vector<glm::vec3> positions;
        std::vector<glm::uvec4> quads;
        std::vector<glm::uvec3> triangles;
    };

    struct MeshQualitySummary {
        size_t faceCount = 0;
        size_t validFaces = 0;
        size_t degenerateFaces = 0;
        float minArea = std::numeric_limits<float>::max();
        float maxArea = 0.0f;
        float minEdge = std::numeric_limits<float>::max();
        float maxEdge = 0.0f;
        float minAngle = 180.0f;
        float maxAngle = 0.0f;
        float minAspect = std::numeric_limits<float>::max();
        float maxAspect = 0.0f;
    };

    CombinedTriMesh3D g_combined_tri_mesh3d;
    CombinedQuadMesh3D g_combined_quad_mesh3d;
    MeshQualitySummary g_tri_quality_summary;
    MeshQualitySummary g_quad_quality_summary;
    bool g_tri_mesh3d_valid = false;
    bool g_quad_mesh3d_valid = false;

    constexpr float kQualityEpsilon = 1e-6f;

    static void invalidateCombinedTriMesh()
    {
        g_combined_tri_mesh3d.positions.clear();
        g_combined_tri_mesh3d.triangles.clear();
        g_tri_mesh3d_valid = false;
        g_tri_quality_summary = {};
    }

    static void invalidateCombinedQuadMesh()
    {
        g_combined_quad_mesh3d.positions.clear();
        g_combined_quad_mesh3d.quads.clear();
        g_quad_mesh3d_valid = false;
        g_quad_quality_summary = {};
    }

    static void ensureBoundaryAndSimulation(size_t chartIdx)
    {
        if (chartIdx >= g_cached_boundary_loops.size())
            return;
        if (!g_boundaries_all[chartIdx]) {
            auto holes = sph_integration::SPHMeshIntegrator::getChartHoles((uint32_t)chartIdx);
            g_boundaries_all[chartIdx] = std::make_unique<Boundary>(g_cached_boundary_loops[chartIdx], holes);
        }
        if (!g_simulations_all[chartIdx] && g_boundaries_all[chartIdx]) {
            g_simulations_all[chartIdx] = std::make_unique<Simulation2D>(*g_boundaries_all[chartIdx]);
            g_all_particles_cache_dirty = true;
        }
    }

    static void ensureTriMeshCache()
    {
        if (!atlasIsReady())
            return;
        if (!g_boundary_loops_cached)
            updateBoundaryLoopsCache();
            ensure_all_containers_sized();
        bool generatedAny = false;
        for (size_t i = 0; i < g_cached_boundary_loops.size(); ++i) {
            ensureBoundaryAndSimulation(i);
            if (i >= g_simulations_all.size())
                continue;
            if (!g_simulations_all[i] || !g_boundaries_all[i])
                continue;
            if (!g_all_tri_indices[i].empty())
                continue;
            CGALMeshGenerator gen;
            gen.generate_mesh(g_simulations_all[i]->get_particles(), *g_boundaries_all[i]);
            g_all_tri_vertices[i] = gen.get_vertices();
            std::vector<glm::uvec3> triIdx;
            const auto& tris = gen.get_triangles();
            triIdx.reserve(tris.size());
            for (const auto& t : tris)
                triIdx.emplace_back(t.v0, t.v1, t.v2);
            g_all_tri_indices[i] = std::move(triIdx);
            generatedAny = true;
        }
        if (generatedAny)
            invalidateCombinedTriMesh();
    }

    static void ensureQuadMeshCache()
    {
        ensureTriMeshCache();
        ensure_all_containers_sized();
        bool generatedAny = false;
        for (size_t i = 0; i < g_cached_boundary_loops.size(); ++i) {
            ensureBoundaryAndSimulation(i);
            if (i >= g_simulations_all.size())
                continue;
            if (!g_simulations_all[i] || !g_boundaries_all[i])
                continue;
            if (!g_all_quad_indices[i].empty())
                continue;

            CGALMeshGenerator gen;
            gen.generate_mesh(g_simulations_all[i]->get_particles(), *g_boundaries_all[i]);
            Qmorph qm;
            auto res = qm.run(gen);
            g_all_quad_vertices[i] = gen.get_vertices();
            std::vector<glm::uvec4> quadIdx;
            quadIdx.reserve(res.quads.size());
            for (const auto& q : res.quads)
                quadIdx.emplace_back(q.v0, q.v1, q.v2, q.v3);
            g_all_quad_indices[i] = std::move(quadIdx);
            std::vector<glm::uvec3> remainIdx;
            remainIdx.reserve(res.remaining_triangles.size());
            for (const auto& t : res.remaining_triangles)
                remainIdx.emplace_back(t.v0, t.v1, t.v2);
            g_all_quad_extra_tris[i] = std::move(remainIdx);
            generatedAny = true;
        }
        if (generatedAny)
            invalidateCombinedQuadMesh();
    }

    static bool mapUvTo3DPosition(const glm::vec2& uv, glm::vec3& outPos)
    {
        float uvCoords[2] = { uv.x, uv.y };
        bx::Vec3 pos;
        if (!modelMapUvTo3D(uvCoords, &pos))
            return false;
        outPos = glm::vec3(pos.x, pos.y, pos.z);
        return true;
    }

    static MeshQualitySummary analyzeTriMeshQuality(const CombinedTriMesh3D& mesh);
    static MeshQualitySummary analyzeQuadMeshQuality(const CombinedQuadMesh3D& mesh);

    static bool rebuildCombinedTriMesh3D()
    {
        g_combined_tri_mesh3d.positions.clear();
        g_combined_tri_mesh3d.triangles.clear();
        g_tri_mesh3d_valid = false;
        g_tri_quality_summary = {};

        if (!modelIsLoaded()) {
            printf("Cannot merge triangle mesh: model not loaded.\n");
            return false;
        }
        if (!atlasIsReady()) {
            printf("Cannot merge triangle mesh: atlas not ready.\n");
            return false;
        }
        ensureTriMeshCache();
        size_t vertexFailures = 0;
        size_t skippedFaces = 0;
        size_t chartsUsed = 0;

        for (size_t chart = 0; chart < g_all_tri_vertices.size(); ++chart) {
            if (chart >= g_all_tri_indices.size())
                break;
            const auto& uvVerts = g_all_tri_vertices[chart];
            const auto& tris = g_all_tri_indices[chart];
            if (uvVerts.empty() || tris.empty())
                continue;
            chartsUsed++;

            std::vector<int> localToGlobal(uvVerts.size(), -1);
            for (size_t v = 0; v < uvVerts.size(); ++v) {
                glm::vec3 pos3d;
                if (!mapUvTo3DPosition(uvVerts[v], pos3d)) {
                    vertexFailures++;
                    continue;
                }
                localToGlobal[v] = (int)g_combined_tri_mesh3d.positions.size();
                g_combined_tri_mesh3d.positions.push_back(pos3d);
            }

            for (const auto& tri : tris) {
                if (tri.x >= localToGlobal.size() || tri.y >= localToGlobal.size() || tri.z >= localToGlobal.size()) {
                    skippedFaces++;
                    continue;
                }
                int ia = localToGlobal[tri.x];
                int ib = localToGlobal[tri.y];
                int ic = localToGlobal[tri.z];
                if (ia < 0 || ib < 0 || ic < 0) {
                    skippedFaces++;
                    continue;
                }
                g_combined_tri_mesh3d.triangles.emplace_back((uint32_t)ia, (uint32_t)ib, (uint32_t)ic);
            }
        }

        if (!g_combined_tri_mesh3d.triangles.empty()) {
            g_tri_quality_summary = analyzeTriMeshQuality(g_combined_tri_mesh3d);
            g_tri_mesh3d_valid = true;
            printf("Combined triangle mesh ready: %zu charts, %zu vertices, %zu faces (skipped faces: %zu, invalid vertices: %zu)\n",
                   chartsUsed,
                   g_combined_tri_mesh3d.positions.size(),
                   g_combined_tri_mesh3d.triangles.size(),
                   skippedFaces,
                   vertexFailures);
        } else {
            printf("Combined triangle mesh is empty (skipped faces: %zu, invalid vertices: %zu)\n", skippedFaces, vertexFailures);
        }
        return g_tri_mesh3d_valid;
    }

    static bool rebuildCombinedQuadMesh3D()
    {
        g_combined_quad_mesh3d.positions.clear();
        g_combined_quad_mesh3d.quads.clear();
        g_quad_mesh3d_valid = false;
        g_quad_quality_summary = {};

        if (!modelIsLoaded()) {
            printf("Cannot merge quad mesh: model not loaded.\n");
            return false;
        }
        if (!atlasIsReady()) {
            printf("Cannot merge quad mesh: atlas not ready.\n");
            return false;
        }
        ensureQuadMeshCache();
        size_t vertexFailures = 0;
        size_t skippedFaces = 0;
        size_t chartsUsed = 0;

        for (size_t chart = 0; chart < g_all_quad_vertices.size(); ++chart) {
            if (chart >= g_all_quad_indices.size())
                break;
            const auto& uvVerts = g_all_quad_vertices[chart];
            const auto& quads = g_all_quad_indices[chart];
            const std::vector<glm::uvec3>* extra = nullptr;
            if (chart < g_all_quad_extra_tris.size())
                extra = &g_all_quad_extra_tris[chart];
            if (uvVerts.empty() || (quads.empty() && (!extra || extra->empty())))
                continue;
            chartsUsed++;

            std::vector<int> localToGlobal(uvVerts.size(), -1);
            for (size_t v = 0; v < uvVerts.size(); ++v) {
                glm::vec3 pos3d;
                if (!mapUvTo3DPosition(uvVerts[v], pos3d)) {
                    vertexFailures++;
                    continue;
                }
                localToGlobal[v] = (int)g_combined_quad_mesh3d.positions.size();
                g_combined_quad_mesh3d.positions.push_back(pos3d);
            }

            for (const auto& quad : quads) {
                const uint32_t idx[4] = { quad.x, quad.y, quad.z, quad.w };
                bool valid = true;
                uint32_t mapped[4];
                for (int k = 0; k < 4; ++k) {
                    if (idx[k] >= localToGlobal.size()) {
                        valid = false;
                        break;
                    }
                    int mappedIdx = localToGlobal[idx[k]];
                    if (mappedIdx < 0) {
                        valid = false;
                        break;
                    }
                    mapped[k] = (uint32_t)mappedIdx;
                }
                if (!valid) {
                    skippedFaces++;
                    continue;
                }
                g_combined_quad_mesh3d.quads.emplace_back(mapped[0], mapped[1], mapped[2], mapped[3]);
            }

            if (extra) {
                for (const auto& tri : *extra) {
                    const uint32_t idx[3] = { tri.x, tri.y, tri.z };
                    bool valid = true;
                    uint32_t mapped[3];
                    for (int k = 0; k < 3; ++k) {
                        if (idx[k] >= localToGlobal.size()) {
                            valid = false;
                            break;
                        }
                        int mappedIdx = localToGlobal[idx[k]];
                        if (mappedIdx < 0) {
                            valid = false;
                            break;
                        }
                        mapped[k] = (uint32_t)mappedIdx;
                    }
                    if (!valid) {
                        skippedFaces++;
                        continue;
                    }
                    g_combined_quad_mesh3d.triangles.emplace_back(mapped[0], mapped[1], mapped[2]);
                }
            }
        }

        size_t totalFaces = g_combined_quad_mesh3d.quads.size() + g_combined_quad_mesh3d.triangles.size();
        if (totalFaces > 0) {
            g_quad_quality_summary = analyzeQuadMeshQuality(g_combined_quad_mesh3d);
            g_quad_mesh3d_valid = true;
            printf("Combined quad mesh ready: %zu charts, %zu vertices, %zu quads + %zu tris (skipped faces: %zu, invalid vertices: %zu)\n",
                   chartsUsed,
                   g_combined_quad_mesh3d.positions.size(),
                   g_combined_quad_mesh3d.quads.size(),
                   g_combined_quad_mesh3d.triangles.size(),
                   skippedFaces,
                   vertexFailures);
        } else {
            printf("Combined quad mesh is empty (skipped faces: %zu, invalid vertices: %zu)\n", skippedFaces, vertexFailures);
        }

        return g_quad_mesh3d_valid;
    }

    static float safeAngleDeg(const glm::vec3& a, const glm::vec3& b)
    {
        float lenA = glm::length(a);
        float lenB = glm::length(b);
        if (lenA <= kQualityEpsilon || lenB <= kQualityEpsilon)
            return 0.0f;
        float cosValue = glm::dot(a, b) / (lenA * lenB);
        cosValue = glm::clamp(cosValue, -1.0f, 1.0f);
        return std::acos(cosValue) * 57.2957795f;
    }

    static MeshQualitySummary analyzeTriMeshQuality(const CombinedTriMesh3D& mesh)
    {
        MeshQualitySummary summary;
        summary.faceCount = mesh.triangles.size();
        for (const auto& tri : mesh.triangles) {
            if (tri.x >= mesh.positions.size() ||
                tri.y >= mesh.positions.size() ||
                tri.z >= mesh.positions.size()) {
                summary.degenerateFaces++;
                continue;
            }
            const glm::vec3& a = mesh.positions[tri.x];
            const glm::vec3& b = mesh.positions[tri.y];
            const glm::vec3& c = mesh.positions[tri.z];
            glm::vec3 ab = b - a;
            glm::vec3 ac = c - a;
            glm::vec3 bc = c - b;

            float lenAB = glm::length(ab);
            float lenAC = glm::length(ac);
            float lenBC = glm::length(bc);
            summary.minEdge = std::min(summary.minEdge, std::min(lenAB, std::min(lenAC, lenBC)));
            summary.maxEdge = std::max(summary.maxEdge, std::max(lenAB, std::max(lenAC, lenBC)));

            float area = 0.5f * glm::length(glm::cross(ab, ac));
            if (area <= kQualityEpsilon) {
                summary.degenerateFaces++;
                continue;
            }

            summary.validFaces++;
            summary.minArea = std::min(summary.minArea, area);
            summary.maxArea = std::max(summary.maxArea, area);

            float aspect = (lenAB * lenAB + lenAC * lenAC + lenBC * lenBC) / (4.0f * std::sqrt(3.0f) * area);
            summary.minAspect = std::min(summary.minAspect, aspect);
            summary.maxAspect = std::max(summary.maxAspect, aspect);

            float angleA = safeAngleDeg(ab, ac);
            float angleB = safeAngleDeg(-ab, bc);
            float angleC = safeAngleDeg(-ac, -bc);
            summary.minAngle = std::min(summary.minAngle, std::min(angleA, std::min(angleB, angleC)));
            summary.maxAngle = std::max(summary.maxAngle, std::max(angleA, std::max(angleB, angleC)));
        }

        if (summary.validFaces == 0) {
            summary.minArea = 0.0f;
            summary.minEdge = 0.0f;
            summary.minAngle = 0.0f;
            summary.minAspect = 0.0f;
        }
        return summary;
    }

    static MeshQualitySummary analyzeQuadMeshQuality(const CombinedQuadMesh3D& mesh)
    {
        MeshQualitySummary summary;
        summary.faceCount = mesh.quads.size() + mesh.triangles.size();
        for (const auto& quad : mesh.quads) {
            if (quad.x >= mesh.positions.size() ||
                quad.y >= mesh.positions.size() ||
                quad.z >= mesh.positions.size() ||
                quad.w >= mesh.positions.size()) {
                summary.degenerateFaces++;
                continue;
            }
            const glm::vec3& p0 = mesh.positions[quad.x];
            const glm::vec3& p1 = mesh.positions[quad.y];
            const glm::vec3& p2 = mesh.positions[quad.z];
            const glm::vec3& p3 = mesh.positions[quad.w];

            glm::vec3 edges[4] = { p1 - p0, p2 - p1, p3 - p2, p0 - p3 };
            float edgeLens[4];
            for (int i = 0; i < 4; ++i) {
                edgeLens[i] = glm::length(edges[i]);
                summary.minEdge = std::min(summary.minEdge, edgeLens[i]);
                summary.maxEdge = std::max(summary.maxEdge, edgeLens[i]);
            }

            float area = 0.5f * glm::length(glm::cross(p1 - p0, p2 - p0));
            area += 0.5f * glm::length(glm::cross(p3 - p0, p2 - p0));
            if (area <= kQualityEpsilon) {
                summary.degenerateFaces++;
                continue;
            }

            summary.validFaces++;
            summary.minArea = std::min(summary.minArea, area);
            summary.maxArea = std::max(summary.maxArea, area);

            float minEdgeQuad = std::min(std::min(edgeLens[0], edgeLens[1]), std::min(edgeLens[2], edgeLens[3]));
            float maxEdgeQuad = std::max(std::max(edgeLens[0], edgeLens[1]), std::max(edgeLens[2], edgeLens[3]));
            float aspect = (minEdgeQuad <= kQualityEpsilon) ? std::numeric_limits<float>::infinity() : (maxEdgeQuad / minEdgeQuad);
            summary.minAspect = std::min(summary.minAspect, aspect);
            summary.maxAspect = std::max(summary.maxAspect, aspect);

            glm::vec3 verts[4] = { p0, p1, p2, p3 };
            for (int i = 0; i < 4; ++i) {
                const glm::vec3& curr = verts[i];
                const glm::vec3& prev = verts[(i + 3) % 4];
                const glm::vec3& next = verts[(i + 1) % 4];
                float angle = safeAngleDeg(prev - curr, next - curr);
                summary.minAngle = std::min(summary.minAngle, angle);
                summary.maxAngle = std::max(summary.maxAngle, angle);
            }
        }

        for (const auto& tri : mesh.triangles) {
            if (tri.x >= mesh.positions.size() ||
                tri.y >= mesh.positions.size() ||
                tri.z >= mesh.positions.size()) {
                summary.degenerateFaces++;
                continue;
            }
            const glm::vec3& a = mesh.positions[tri.x];
            const glm::vec3& b = mesh.positions[tri.y];
            const glm::vec3& c = mesh.positions[tri.z];
            glm::vec3 ab = b - a;
            glm::vec3 ac = c - a;
            glm::vec3 bc = c - b;

            float lenAB = glm::length(ab);
            float lenAC = glm::length(ac);
            float lenBC = glm::length(bc);
            summary.minEdge = std::min(summary.minEdge, std::min(lenAB, std::min(lenAC, lenBC)));
            summary.maxEdge = std::max(summary.maxEdge, std::max(lenAB, std::max(lenAC, lenBC)));

            float area = 0.5f * glm::length(glm::cross(ab, ac));
            if (area <= kQualityEpsilon) {
                summary.degenerateFaces++;
                continue;
            }

            summary.validFaces++;
            summary.minArea = std::min(summary.minArea, area);
            summary.maxArea = std::max(summary.maxArea, area);

            float aspect = (lenAB * lenAB + lenAC * lenAC + lenBC * lenBC) / (4.0f * std::sqrt(3.0f) * area);
            summary.minAspect = std::min(summary.minAspect, aspect);
            summary.maxAspect = std::max(summary.maxAspect, aspect);

            float angleA = safeAngleDeg(ab, ac);
            float angleB = safeAngleDeg(-ab, bc);
            float angleC = safeAngleDeg(-ac, -bc);
            summary.minAngle = std::min(summary.minAngle, std::min(angleA, std::min(angleB, angleC)));
            summary.maxAngle = std::max(summary.maxAngle, std::max(angleA, std::max(angleB, angleC)));
        }

        if (summary.validFaces == 0) {
            summary.minArea = 0.0f;
            summary.minEdge = 0.0f;
            summary.minAngle = 0.0f;
            summary.minAspect = 0.0f;
        }
        return summary;
    }

    static bool exportTriMeshObj(const CombinedTriMesh3D& mesh, const char* filepath)
    {
        if (mesh.positions.empty() || mesh.triangles.empty())
            return false;
        FILE* f = nullptr;
#if _MSC_VER
        fopen_s(&f, filepath, "w");
#else
        f = fopen(filepath, "w");
#endif
        if (!f) {
            printf("Failed to open '%s' for writing\n", filepath);
            return false;
        }
        for (const auto& v : mesh.positions)
            fprintf(f, "v %.9f %.9f %.9f\n", v.x, v.y, v.z);
        for (const auto& tri : mesh.triangles)
            fprintf(f, "f %u %u %u\n", tri.x + 1, tri.y + 1, tri.z + 1);
        fclose(f);
        printf("Exported combined triangle mesh to '%s'\n", filepath);
        return true;
    }

    static bool exportQuadMeshObj(const CombinedQuadMesh3D& mesh, const char* filepath)
    {
        if (mesh.positions.empty() || (mesh.quads.empty() && mesh.triangles.empty()))
            return false;
        FILE* f = nullptr;
#if _MSC_VER
        fopen_s(&f, filepath, "w");
#else
        f = fopen(filepath, "w");
#endif
        if (!f) {
            printf("Failed to open '%s' for writing\n", filepath);
            return false;
        }
        for (const auto& v : mesh.positions)
            fprintf(f, "v %.9f %.9f %.9f\n", v.x, v.y, v.z);
        for (const auto& quad : mesh.quads)
            fprintf(f, "f %u %u %u %u\n", quad.x + 1, quad.y + 1, quad.z + 1, quad.w + 1);
        for (const auto& tri : mesh.triangles)
            fprintf(f, "f %u %u %u\n", tri.x + 1, tri.y + 1, tri.z + 1);
        fclose(f);
        printf("Exported combined quad mesh to '%s'\n", filepath);
        return true;
    }
}

static const Simulation2D* get_current_simulation()
{
    if (g_boundary_loops_cached && g_current_chart_index < g_simulations_all.size()) {
        if (g_simulations_all[g_current_chart_index]) return g_simulations_all[g_current_chart_index].get();
    }
    return g_sph_simulation.get();
}

static const Boundary* get_current_boundary()
{
    if (g_boundary_loops_cached && g_current_chart_index < g_boundaries_all.size()) {
        if (g_boundaries_all[g_current_chart_index]) return g_boundaries_all[g_current_chart_index].get();
    }
    return g_current_boundary.get();
}

void sphInit() {
    g_sph_integrator = std::make_unique<sph_integration::SPHMeshIntegrator>();
}

void sphShutdown() {
    g_sph_integrator.reset();
    g_current_boundary.reset();
    g_sph_simulation.reset();
    g_boundaries_all.clear();
    g_simulations_all.clear();
    g_all_particles_cache.clear();
    g_all_particles_cache_dirty = true;
    g_delaunay_generator.reset();
    g_qmorph.reset();
    g_cached_boundary_loops.clear();
    g_boundary_loops_cached = false;
    g_tri_mesh_generated = false;
    g_quad_mesh_generated = false;
    g_view_mesh = false;
}


void updateBoundaryLoopsCache() {
    if (!g_sph_integrator) return;
    
    g_cached_boundary_loops = g_sph_integrator->extractBoundaryLoopsFromAtlas();

    // 可选：UV归一化（同时归一化外环和内环）
    const xatlas::Atlas* atlas = (const xatlas::Atlas*)atlasGetData();
    if (atlas && g_normalize_uv_01) {
        const float invW = atlas->width > 0 ? 1.0f / float(atlas->width) : 1.0f;
        const float invH = atlas->height > 0 ? 1.0f / float(atlas->height) : 1.0f;
        
        // 归一化外环
        for (auto& loop : g_cached_boundary_loops) {
            for (auto& p : loop) {
                p.x *= invW;
                p.y *= invH;
            }
        }
        
        // 归一化内环（洞）
        for (uint32_t chartIdx = 0; chartIdx < g_cached_boundary_loops.size(); ++chartIdx) {
            std::vector<std::vector<glm::vec2>> holes = 
                sph_integration::SPHMeshIntegrator::getChartHoles(chartIdx);
            for (auto& hole : holes) {
                for (auto& p : hole) {
                    p.x *= invW;
                    p.y *= invH;
                }
            }
            // 更新归一化后的内环到缓存
            if (!holes.empty()) {
                sph_integration::SPHMeshIntegrator::updateChartHoles(chartIdx, holes);
            }
        }
    }

    // 可选：图表局部归一化到[0,1]
    if (g_chart_local_01) {
        for (size_t chartIdx = 0; chartIdx < g_cached_boundary_loops.size(); ++chartIdx) {
            auto& loop = g_cached_boundary_loops[chartIdx];
            if (loop.empty()) continue;
            
            glm::vec2 mn = loop[0], mx = loop[0];
            for (const auto& p : loop) {
                mn.x = std::min(mn.x, p.x); mn.y = std::min(mn.y, p.y);
                mx.x = std::max(mx.x, p.x); mx.y = std::max(mx.y, p.y);
            }
            glm::vec2 size = mx - mn;
            if (size.x == 0.0f) size.x = 1.0f;
            if (size.y == 0.0f) size.y = 1.0f;
            
            // 归一化外环
            for (auto& p : loop) {
                p = (p - mn) / size;
            }
            
            // 归一化内环（使用相同的变换）
            std::vector<std::vector<glm::vec2>> holes = 
                sph_integration::SPHMeshIntegrator::getChartHoles((uint32_t)chartIdx);
            for (auto& hole : holes) {
                for (auto& p : hole) {
                    p = (p - mn) / size;
                }
            }
            // 更新归一化后的内环到缓存
            if (!holes.empty()) {
                sph_integration::SPHMeshIntegrator::updateChartHoles((uint32_t)chartIdx, holes);
            }
        }
    }
    g_boundary_loops_cached = true;
    printf("Updated boundary loops cache: %zu loops found\n", g_cached_boundary_loops.size());
}

static void ensure_all_containers_sized()
{
    const size_t n = g_cached_boundary_loops.size();
    if (g_boundaries_all.size() < n) g_boundaries_all.resize(n);
    if (g_simulations_all.size() < n) g_simulations_all.resize(n);
    if (g_all_tri_vertices.size() < n) g_all_tri_vertices.resize(n);
    if (g_all_tri_indices.size() < n) g_all_tri_indices.resize(n);
    if (g_all_quad_vertices.size() < n) g_all_quad_vertices.resize(n);
    if (g_all_quad_indices.size() < n) g_all_quad_indices.resize(n);
    if (g_all_quad_extra_tris.size() < n) g_all_quad_extra_tris.resize(n);
}

void sphProcessChart(uint32_t chartIndex) {
    if (!g_sph_integrator || !atlasIsReady()) {
        return;
    }
    
    g_current_chart_index = chartIndex;
    
    // 使用缓存的边界环，避免重复计算
    if (!g_boundary_loops_cached) {
        updateBoundaryLoopsCache();
    }

    if (chartIndex < g_cached_boundary_loops.size()) {
        // 获取内环（洞）信息
        std::vector<std::vector<glm::vec2>> holes = 
            sph_integration::SPHMeshIntegrator::getChartHoles(chartIndex);
        
        // 创建Boundary对象（包含外环和内环）
        // 第二个参数有默认值，可以直接传入
        g_current_boundary = std::make_unique<Boundary>(g_cached_boundary_loops[chartIndex], holes);
        
        // 创建Simulation2D（会自动初始化粒子）
        // 这会根据grid_cell_size自动调整粒子数量
        g_sph_simulation = std::make_unique<Simulation2D>(*g_current_boundary);
        
        // 不再使用旧的integrator，只使用新的Simulation2D
        // (旧代码保留但不调用，以防需要)
        // sph_integration::SPHMeshIntegrator::ChartData chart_data;
        // chart_data.chart_index = chartIndex;
        // chart_data.boundary_points = g_cached_boundary_loops[chartIndex];
        // g_sph_integrator->createBoundaryFromChart(chart_data);
        // g_sph_integrator->initializeParticles();
        
        printf("Processed chart %d with %zu boundary points, %zu particles\n", 
               chartIndex, g_cached_boundary_loops[chartIndex].size(),
               g_sph_simulation ? g_sph_simulation->get_particles().size() : 0);

        // 同步到批量容器
        ensure_all_containers_sized();
        g_boundaries_all[chartIndex] = std::make_unique<Boundary>(*g_current_boundary);
        g_simulations_all[chartIndex] = std::make_unique<Simulation2D>(*g_boundaries_all[chartIndex]);
        g_all_particles_cache_dirty = true;
    } else {
        printf("Chart index %d out of range (found %zu boundary loops)\n", chartIndex, g_cached_boundary_loops.size());
    }
}

void sphRunSimulation() {
    if (!g_sph_integrator) {
        return;
    }
    
    g_sph_integrator->runSimulation(g_simulation_steps);
}

// 生成Delaunay三角网格
void sphGenerateTriMesh() {
    printf("=== Button: Generate Tri Mesh clicked ===\n");
    // 明确用户只请求三角网格
    g_request_quad_on_switch = false;
    const Simulation2D* sim = get_current_simulation();
    const Boundary* boundary = get_current_boundary();
    if (!sim || !boundary) {
        printf("ERROR: Simulation or boundary not ready (current chart index %u)!\n", g_current_chart_index);
        return;
    }
    
    // 创建Delaunay网格生成器
    g_delaunay_generator = std::make_unique<CGALMeshGenerator>();
    
    // 从粒子生成三角网格
    const auto& particles = sim->get_particles();
    g_delaunay_generator->generate_mesh(particles, *boundary);
    
    g_tri_mesh_generated = true;
    g_quad_mesh_generated = false;  // 重置四边形网格状态
    
    const auto& vertices = g_delaunay_generator->get_vertices();
    const auto& triangles = g_delaunay_generator->get_triangles();
    printf("Tri mesh generated: %zu vertices, %zu triangles\n", vertices.size(), triangles.size());

    // 缓存到全局（当前图表）
    ensure_all_containers_sized();
    if (g_current_chart_index < g_all_tri_vertices.size()) {
        g_all_tri_vertices[g_current_chart_index] = vertices;
        std::vector<glm::uvec3> triIdx;
        triIdx.reserve(triangles.size());
        for (const auto &t : triangles) triIdx.emplace_back(t.v0, t.v1, t.v2);
        g_all_tri_indices[g_current_chart_index]  = std::move(triIdx);
    }
    invalidateCombinedTriMesh();
    invalidateCombinedQuadMesh();
}

// 合并三角网格为四边形网格
void sphGenerateQuadMesh() {
    printf("=== Button: Generate Quad Mesh clicked ===\n");
    // 用户明确请求四边形网格，后续切换图表时也保留这一意图
    g_request_quad_on_switch = true;
    if (!g_tri_mesh_generated || !g_delaunay_generator) {
        printf("ERROR: Tri mesh not generated yet! Please generate tri mesh first.\n");
        return;
    }
    
    // 创建Qmorph转换器
    g_qmorph = std::make_unique<Qmorph>();
    
    // 将三角网格转换为四边形网格
    g_quad_mesh_result = g_qmorph->run(*g_delaunay_generator);
    
    g_quad_mesh_generated = true;
    
    printf("Quad mesh generated: %zu quads, %zu remaining triangles\n", 
           g_quad_mesh_result.quads.size(), 
           g_quad_mesh_result.remaining_triangles.size());

    // 缓存到全局（当前图表）
    ensure_all_containers_sized();
    if (g_current_chart_index < g_all_quad_vertices.size()) {
        g_all_quad_vertices[g_current_chart_index] = g_delaunay_generator->get_vertices();
        std::vector<glm::uvec4> quadIdx;
        quadIdx.reserve(g_quad_mesh_result.quads.size());
        for (const auto &q : g_quad_mesh_result.quads) quadIdx.emplace_back(q.v0, q.v1, q.v2, q.v3);
        g_all_quad_indices[g_current_chart_index]  = std::move(quadIdx);
        std::vector<glm::uvec3> remainIdx;
        remainIdx.reserve(g_quad_mesh_result.remaining_triangles.size());
        for (const auto& t : g_quad_mesh_result.remaining_triangles)
            remainIdx.emplace_back(t.v0, t.v1, t.v2);
        g_all_quad_extra_tris[g_current_chart_index] = std::move(remainIdx);
    }
    invalidateCombinedQuadMesh();
}

// 保持向后兼容
void sphGenerateMesh() {
    sphGenerateQuadMesh();
}

void sphShowGuiOptions() {
    if (!atlasIsReady()) {
        return;
    }
    
    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Spacing();
    ImGui::Text("SPH Mesh Generator");
    ImGui::Spacing();
    
    
    ImGui::Checkbox("Enable SPH Mesh Generation", &g_sph_enabled);
    
    if (g_sph_enabled) {
        ImGui::Spacing();
        
       
        ImGui::Text("Chart Selection:");
        
        
        if (!g_boundary_loops_cached) {
            updateBoundaryLoopsCache();
        }
        
        int max_charts = (int)g_cached_boundary_loops.size() - 1;
        if (max_charts < 0) max_charts = 0;
        
        ImGui::Text("Found %zu boundary loops", g_cached_boundary_loops.size());
        uint32_t prev_index = g_current_chart_index;
        ImGui::SliderInt("Boundary Loop Index", (int*)&g_current_chart_index, 0, max_charts);
        if (g_boundary_loops_cached && g_current_chart_index != prev_index) {
            ensure_all_containers_sized();
            
            if (g_current_chart_index < g_boundaries_all.size() && g_boundaries_all[g_current_chart_index]) {
                g_current_boundary = std::make_unique<Boundary>(*g_boundaries_all[g_current_chart_index]);
            }
           
            if (g_view_mesh) {
                sphGenerateTriMesh();
                if (g_tri_mesh_generated && g_request_quad_on_switch) {
                    sphGenerateQuadMesh();
                }
            }
        }
        ImGui::Checkbox("Normalize UV to [0,1] (atlas)", &g_normalize_uv_01);
        ImGui::Checkbox("Chart-local [0,1]", &g_chart_local_01);
        
        if (ImGui::Button("Refresh Boundary Loops")) {
            printf("=== Button: Refresh Boundary Loops clicked ===\n");
            updateBoundaryLoopsCache();
        }
        
        if (ImGui::Button("Process Boundary Loop")) {
            printf("=== Button: Process Boundary Loop clicked (chart %d) ===\n", g_current_chart_index);
            sphProcessChart(g_current_chart_index);
        }
        
        if (ImGui::Button("Initialize Particles")) {
            printf("=== Button: Initialize Particles clicked ===\n");
            if (g_current_boundary) {
                printf("Creating new Simulation2D...\n");
                g_sph_simulation = std::make_unique<Simulation2D>(*g_current_boundary);
                printf("Simulation2D created with %zu particles\n", g_sph_simulation->get_particles().size());
            } else {
                printf("ERROR: g_current_boundary is null! Call Process Boundary Loop first.\n");
            }
        }
        
        if (ImGui::Button("Process All Boundary Loops")) {
            printf("=== Button: Process All Boundary Loops clicked ===\n");
            sphProcessAllBoundaryLoops();
        }
        
        ImGui::Spacing();
        
       
        ImGui::Text("Simulation Control:");
        ImGui::SliderInt("Simulation Steps", &g_simulation_steps, 10, 1000);
        ImGui::Checkbox("Auto Simulation", &g_auto_simulation);
        
        if (g_auto_simulation && g_sph_simulation) {
            
            sphStepSimulationAll();
            if (g_sph_simulation) g_sph_simulation->step();
        }
        
        if (ImGui::Button("Run Simulation")) {
            printf("=== Button: Run Simulation clicked ===\n");
            if (sphHasAnySimulation()) {
                sphRunSimulationAll(g_simulation_steps);
                printf("Run all simulations completed\n");
            }
            if (g_sph_simulation) {
                for (int i = 0; i < g_simulation_steps; i++) g_sph_simulation->step();
                printf("Run current simulation completed\n");
            }
        }
        
        if (ImGui::Button("Step Simulation")) {
            printf("=== Button: Step Simulation clicked ===\n");
            sphStepSimulationAll();
            if (g_sph_simulation) g_sph_simulation->step();
        }
        
        ImGui::Spacing();
        
        
        ImGui::Text("Mesh Generation:");
        if (ImGui::Button("Generate Tri Mesh")) {
            sphGenerateTriMesh();
        }
        if (g_tri_mesh_generated) {
            ImGui::SameLine();
            ImGui::TextColored(ImVec4(0.0f, 1.0f, 0.0f, 1.0f), " [Done]");
            if (g_delaunay_generator) {
                const auto& vertices = g_delaunay_generator->get_vertices();
                const auto& triangles = g_delaunay_generator->get_triangles();
                ImGui::Text("  Tri mesh: %zu vertices, %zu triangles", vertices.size(), triangles.size());
            }
        }
        
        if (ImGui::Button("Generate Quad Mesh")) {
            sphGenerateQuadMesh();
        }
        if (g_quad_mesh_generated) {
            ImGui::SameLine();
            ImGui::TextColored(ImVec4(0.0f, 1.0f, 0.0f, 1.0f), " [Done]");
            ImGui::Text("  Quad mesh: %zu quads, %zu remaining triangles", 
                       g_quad_mesh_result.quads.size(), 
                       g_quad_mesh_result.remaining_triangles.size());
        }

       
        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Spacing();
        {
            bool triOverlay = g_show_all_tri_overlay;
            bool quadOverlay = g_show_all_quad_overlay;
            if (ImGui::Checkbox("Show All Tri Mesh (3D)", &triOverlay)) {
                sphSetShowAllTriOverlay(triOverlay);
            }
            if (ImGui::Checkbox("Show All Quad Mesh (3D)", &quadOverlay)) {
                sphSetShowAllQuadOverlay(quadOverlay);
            }
            bool showParticles3D = g_show_particles_in_3d;
            if (ImGui::Checkbox("Show Particles in 3D", &showParticles3D)) {
                sphSetShowParticlesIn3D(showParticles3D);
            }
        }

        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Spacing();
        ImGui::Text("Surface Mesh Merge & Export");
        if (ImGui::Button("Merge Tri Mesh (3D)")) {
            rebuildCombinedTriMesh3D();
        }
        if (g_tri_mesh3d_valid) {
            ImGui::Text("Tri Mesh: %zu faces / %zu vertices (degenerate %zu)", g_tri_quality_summary.faceCount, g_combined_tri_mesh3d.positions.size(), g_tri_quality_summary.degenerateFaces);
            if (g_tri_quality_summary.validFaces > 0) {
                ImGui::Text("  Area range: %.6f - %.6f", g_tri_quality_summary.minArea, g_tri_quality_summary.maxArea);
                ImGui::Text("  Edge length range: %.6f - %.6f", g_tri_quality_summary.minEdge, g_tri_quality_summary.maxEdge);
                ImGui::Text("  Angle range: %.2f° - %.2f°", g_tri_quality_summary.minAngle, g_tri_quality_summary.maxAngle);
                ImGui::Text("  Aspect ratio range: %.3f - %.3f", g_tri_quality_summary.minAspect, g_tri_quality_summary.maxAspect);
            } else {
                ImGui::Text("  Quality: no valid faces");
            }
            if (ImGui::Button("Export Tri Mesh (OBJ)")) {
                ensureExportDir();
                char path[512];
                const char* base = modelGetFileBaseName();
                snprintf(path, sizeof(path), "exportdata/%s_combined_tri.obj", base && base[0] ? base : "model");
                exportTriMeshObj(g_combined_tri_mesh3d, path);
            }
        } else {
            ImGui::Text("Tri Mesh: not merged yet");
        }
        ImGui::Spacing();
        if (ImGui::Button("Merge Quad Mesh (3D)")) {
            rebuildCombinedQuadMesh3D();
        }
        if (g_quad_mesh3d_valid) {
            ImGui::Text("Quad Mesh: %zu quads + %zu tris / %zu vertices (degenerate %zu)",
                g_combined_quad_mesh3d.quads.size(),
                g_combined_quad_mesh3d.triangles.size(),
                g_combined_quad_mesh3d.positions.size(),
                g_quad_quality_summary.degenerateFaces);
            if (g_quad_quality_summary.validFaces > 0) {
                ImGui::Text("  Area range: %.6f - %.6f", g_quad_quality_summary.minArea, g_quad_quality_summary.maxArea);
                ImGui::Text("  Edge length range: %.6f - %.6f", g_quad_quality_summary.minEdge, g_quad_quality_summary.maxEdge);
                ImGui::Text("  Angle range: %.2f° - %.2f°", g_quad_quality_summary.minAngle, g_quad_quality_summary.maxAngle);
                ImGui::Text("  Aspect ratio range: %.3f - %.3f", g_quad_quality_summary.minAspect, g_quad_quality_summary.maxAspect);
            } else {
                ImGui::Text("  Quality: no valid faces");
            }
            if (ImGui::Button("Export Quad Mesh (OBJ)")) {
                ensureExportDir();
                char path[512];
                const char* base = modelGetFileBaseName();
                snprintf(path, sizeof(path), "exportdata/%s_combined_quad.obj", base && base[0] ? base : "model");
                exportQuadMeshObj(g_combined_quad_mesh3d, path);
            }
        } else {
            ImGui::Text("Quad Mesh: not merged yet");
        }

        ImGui::Checkbox("View Mesh (V)", &g_view_mesh);
        
        ImGui::Spacing();
        
       
        if (!g_cached_boundary_loops.empty() && g_current_chart_index < g_cached_boundary_loops.size()) {
            const auto& current_loop = g_cached_boundary_loops[g_current_chart_index];
            ImGui::Text("Current Boundary Loop:");
            ImGui::Text("Points: %zu", current_loop.size());
            
            if (ImGui::CollapsingHeader("Boundary Points")) {
                for (size_t i = 0; i < current_loop.size() && i < 10; i++) {
                    ImGui::Text("  [%zu]: (%.3f, %.3f)", i, current_loop[i].x, current_loop[i].y);
                }
                if (current_loop.size() > 10) {
                    ImGui::Text("  ... and %zu more points", current_loop.size() - 10);
                }
            }
        }
        
       
        if (g_sph_simulation) {
            const auto& particles = g_sph_simulation->get_particles();
            ImGui::Text("Particles: %zu", particles.size());
            if (!particles.empty()) {
                ImGui::Text("Kinetic Energy: %.6f", g_sph_simulation->get_kinetic_energy());
            }
        }
    }
}


const std::vector<glm::vec2>* sphGetCurrentBoundaryLoop() {
    if (!g_boundary_loops_cached || g_cached_boundary_loops.empty()) {
        return nullptr;
    }
    if (g_current_chart_index >= g_cached_boundary_loops.size()) {
        return nullptr;
    }
    return &g_cached_boundary_loops[g_current_chart_index];
}

uint32_t sphGetCurrentBoundaryLoopIndex() {
    return g_current_chart_index;
}

// 静态变量：存储当前chart的内环列表（用于显示）
static std::vector<std::vector<glm::vec2>> g_current_boundary_holes_cache;

// 获取当前选择的边界环的内环（洞）列表（用于在Atlas窗口中显示）
const std::vector<std::vector<glm::vec2>>* sphGetCurrentBoundaryHoles() {
    if (!g_boundary_loops_cached || g_cached_boundary_loops.empty()) {
        return nullptr;
    }
    if (g_current_chart_index >= g_cached_boundary_loops.size()) {
        return nullptr;
    }
    // 从缓存获取内环
    g_current_boundary_holes_cache = sph_integration::SPHMeshIntegrator::getChartHoles(g_current_chart_index);
    return &g_current_boundary_holes_cache;
}

// 获取粒子位置（用于在Atlas窗口中显示）
const std::vector<glm::vec2>* sphGetParticlePositions() {
    // 网格视图打开时不再返回粒子，避免重复渲染
    if (g_view_mesh) return nullptr;
    // 优先返回“当前选择图表”的批量模拟粒子
    if (g_boundary_loops_cached && g_current_chart_index < g_simulations_all.size()) {
        auto &sim = g_simulations_all[g_current_chart_index];
        if (sim) {
            const auto& positions = sim->get_particle_positions();
            if (!positions.empty()) return &positions;
        }
    }
    // 其次返回当前单独的模拟
    if (g_sph_simulation) {
        const auto& positions = g_sph_simulation->get_particle_positions();
        if (!positions.empty()) {
            return &positions;
        }
    }
    return nullptr;
}

void sphProcessAllBoundaryLoops()
{
    if (!atlasIsReady()) return;
    if (!g_boundary_loops_cached) updateBoundaryLoopsCache();
    ensure_all_containers_sized();
    for (size_t i = 0; i < g_cached_boundary_loops.size(); i++) {
        // 外环与内环
        std::vector<std::vector<glm::vec2>> holes =
            sph_integration::SPHMeshIntegrator::getChartHoles((uint32_t)i);
        g_boundaries_all[i] = std::make_unique<Boundary>(g_cached_boundary_loops[i], holes);
        g_simulations_all[i] = std::make_unique<Simulation2D>(*g_boundaries_all[i]);
    }
    g_all_particles_cache_dirty = true;
    printf("Processed ALL boundary loops: %zu simulations ready\n", g_simulations_all.size());
}

void sphRunSimulationAll(int steps)
{
    if (g_simulations_all.empty()) return;
    for (auto &sim : g_simulations_all) {
        if (!sim) continue;
        for (int i = 0; i < steps; ++i) sim->step();
    }
    g_all_particles_cache_dirty = true;
}

void sphStepSimulationAll()
{
    if (g_simulations_all.empty()) return;
    for (auto &sim : g_simulations_all) {
        if (!sim) continue;
        sim->step();
    }
    g_all_particles_cache_dirty = true;
}

bool sphHasAnySimulation()
{
    for (auto &sim : g_simulations_all) if (sim) return true;
    return false;
}

static void rebuild_all_particles_cache_if_needed()
{
    if (!g_all_particles_cache_dirty) return;
    g_all_particles_cache.clear();
    for (auto &sim : g_simulations_all) {
        if (!sim) continue;
        const auto &pp = sim->get_particle_positions();
        g_all_particles_cache.insert(g_all_particles_cache.end(), pp.begin(), pp.end());
    }
    g_all_particles_cache_dirty = false;
}

const std::vector<glm::vec2>* sphGetAllParticlePositions()
{
    if (g_simulations_all.empty()) return nullptr;
    rebuild_all_particles_cache_if_needed();
    return &g_all_particles_cache;
}

void sphSetShowParticlesIn3D(bool enabled)
{
    g_show_particles_in_3d = enabled;
}

bool sphIsShowParticlesIn3D()
{
    return g_show_particles_in_3d;
}

static void ensureExportDir()
{
#if BX_PLATFORM_WINDOWS
	_mkdir("exportdata");
#else
	mkdir("exportdata", 0755);
#endif
}

void sphExportAllBoundaryLoopsNormalized()
{
    if (!atlasIsReady()) return;
    if (!g_boundary_loops_cached) updateBoundaryLoopsCache();
    ensureExportDir();
    const char *modelBase = modelGetFileBaseName();
    // 获取 atlas 尺寸以便在未归一化时进行归一化
    const xatlas::Atlas* atlas = (const xatlas::Atlas*)atlasGetData();
    const float invW = atlas && atlas->width  > 0 ? 1.0f / float(atlas->width)  : 1.0f;
    const float invH = atlas && atlas->height > 0 ? 1.0f / float(atlas->height) : 1.0f;
    for (size_t i = 0; i < g_cached_boundary_loops.size(); ++i) {
        // 外环
        char filename[512];
        snprintf(filename, sizeof(filename), "exportdata/%s_chart_%zu_boundary.txt", modelBase, i);
        FILE *f = nullptr;
#if _MSC_VER
        fopen_s(&f, filename, "w");
#else
        f = fopen(filename, "w");
#endif
        if (!f) continue;
        const auto &loop = g_cached_boundary_loops[i];
        for (const auto &p : loop) {
            float x = p.x;
            float y = p.y;
            if (!g_normalize_uv_01) {
                x = p.x * invW;
                y = p.y * invH;
            }
            fprintf(f, "%.6f %.6f\n", x, y);
        }
        fclose(f);
        // 内环（洞），若有则分别输出为独立文件
        std::vector<std::vector<glm::vec2>> holes = sph_integration::SPHMeshIntegrator::getChartHoles((uint32_t)i);
        for (size_t h = 0; h < holes.size(); ++h) {
            snprintf(filename, sizeof(filename), "exportdata/%s_chart_%zu_hole_%zu.txt", modelBase, i, h);
#if _MSC_VER
            fopen_s(&f, filename, "w");
#else
            f = fopen(filename, "w");
#endif
            if (!f) continue;
            for (const auto &hp : holes[h]) {
                float x = hp.x;
                float y = hp.y;
                if (!g_normalize_uv_01) {
                    x = hp.x * invW;
                    y = hp.y * invH;
                }
                fprintf(f, "%.6f %.6f\n", x, y);
            }
            fclose(f);
        }
    }
}


void sphSetShowAllTriOverlay(bool enabled)
{
    g_show_all_tri_overlay = enabled;
    if (enabled)
        ensureTriMeshCache();
}
void sphSetShowAllQuadOverlay(bool enabled)
{
    g_show_all_quad_overlay = enabled;
    if (enabled) {
        g_show_all_tri_overlay = true;
        ensureQuadMeshCache();
    }
}
bool sphIsShowAllTriOverlay() { return g_show_all_tri_overlay; }
bool sphIsShowAllQuadOverlay() { return g_show_all_quad_overlay; }

const std::vector<std::vector<glm::vec2>>& sphGetAllTriMeshVertices() { return g_all_tri_vertices; }
const std::vector<std::vector<glm::uvec3>>& sphGetAllTriMeshIndices() { return g_all_tri_indices; }
const std::vector<std::vector<glm::vec2>>& sphGetAllQuadMeshVertices() { return g_all_quad_vertices; }
const std::vector<std::vector<glm::uvec4>>& sphGetAllQuadMeshIndices() { return g_all_quad_indices; }
const std::vector<std::vector<glm::uvec3>>& sphGetAllQuadExtraTriangles() { return g_all_quad_extra_tris; }
// 视图切换
void sphToggleMeshView() {
    g_view_mesh = !g_view_mesh;
}

bool sphIsMeshView() {
    return g_view_mesh;
}

// 输出绘制用三角网格
void sphGetTriMeshForDraw(std::vector<glm::vec2>& outVertices,
                          std::vector<glm::uvec3>& outTriangles) {
    outVertices.clear();
    outTriangles.clear();
    if (!g_delaunay_generator || !g_tri_mesh_generated) return;
    const auto& verts = g_delaunay_generator->get_vertices();
    const auto& tris  = g_delaunay_generator->get_triangles();
    outVertices = verts;
    outTriangles.reserve(tris.size());
    for (const auto& t : tris) {
        outTriangles.emplace_back(t.v0, t.v1, t.v2);
    }
}

// 输出绘制用四边形网格
void sphGetQuadMeshForDraw(std::vector<glm::vec2>& outVertices,
                           std::vector<glm::uvec4>& outQuads) {
    outVertices.clear();
    outQuads.clear();
    if (!g_quad_mesh_generated || !g_delaunay_generator) return;
    const auto& verts = g_delaunay_generator->get_vertices();
    outVertices = verts;
    outQuads.reserve(g_quad_mesh_result.quads.size());
    for (const auto& q : g_quad_mesh_result.quads) {
        outQuads.emplace_back(q.v0, q.v1, q.v2, q.v3);
    }
}
