#pragma once
#include <vector>
#include <glm/glm.hpp>
#include <memory>

// 暂时注释掉SPHMesh1.3-master的依赖，避免链接错误
// class Boundary;
// class Simulation2D;

// 前向声明MeshGenerator2D结构
namespace sph_integration {
    struct Quad {
        unsigned int v0, v1, v2, v3;
    };
    struct Triangle {
        unsigned int v0, v1, v2;
    };
}

namespace sph_integration {

// SPH粒子数据结构
struct Particle {
    glm::vec2 position;
    glm::vec2 velocity;
    glm::vec2 force;
    float density;
    float pressure;
};

// SPH网格生成器集成类
class SPHMeshIntegrator {
public:
    struct ChartData {
        std::vector<glm::vec2> boundary_points;
        std::vector<glm::vec2> uv_coords;
        std::vector<uint32_t> face_indices;
        uint32_t chart_index;
    };

    SPHMeshIntegrator();
    ~SPHMeshIntegrator();

    // 从xatlas图表数据创建SPH边界
    void createBoundaryFromChart(const ChartData& chart_data);
    
    // 从xatlas全局UV网格提取所有边界环（核心功能）
    std::vector<std::vector<glm::vec2>> extractBoundaryLoopsFromAtlas();
    
    // 获取指定chart的内环（洞）列表
    static std::vector<std::vector<glm::vec2>> getChartHoles(uint32_t chartIndex);
    
    // 更新指定chart的内环（洞）列表（用于归一化后更新）
    static void updateChartHoles(uint32_t chartIndex, const std::vector<std::vector<glm::vec2>>& holes);
    
    // 初始化粒子（在边界内生成）
    void initializeParticles();
    
    // 运行SPH粒子模拟
    void runSimulation(int steps = 100);
    
    // 单步模拟（用于实时显示）
    void stepSimulation();
    
    // 生成四边形网格
    void generateQuadMesh();
    
    // 获取粒子位置（用于显示）
    const std::vector<Particle>& getParticles() const { return particles_; }
    const std::vector<glm::vec2>& getParticlePositions() const;
   
    const std::vector<glm::vec2>& getMeshVertices() const;
    const std::vector<Quad>& getMeshQuads() const;
    
    // 将2D网格映射回3D模型
    void mapTo3DModel(const std::vector<glm::vec3>& original_vertices,
                     const std::vector<glm::vec2>& original_uvs,
                     std::vector<glm::vec3>& output_vertices,
                     std::vector<glm::vec2>& output_uvs);

private:
    // 检查点是否在边界内
    bool isPointInsideBoundary(const glm::vec2& point) const;
    
    // 计算点到边界的最短距离
    float distanceToBoundary(const glm::vec2& point) const;
    
    // SPH核函数
    float sphKernel(float distance, float smoothingRadius) const;
    glm::vec2 sphKernelGradient(const glm::vec2& r, float smoothingRadius) const;
    
    // 计算粒子密度
    void computeDensity();
    
    // 计算粒子压力
    void computePressure();
    
    // 计算粒子受力
    void computeForces();
    
    // 更新粒子位置
    void updateParticles();
    
    // 处理边界碰撞
    void handleBoundaries();
    
    ChartData current_chart_;
    bool simulation_ready_;
    
    // SPH粒子
    std::vector<Particle> particles_;
    std::vector<glm::vec2> particle_positions_for_render_;
    
    // SPH参数
    float smoothing_radius_ = 3.0f;  // 调整为像素单位（原来的0.1对像素坐标太小）
    float rest_density_ = 0.0f;
    float gas_constant_ = 0.0f;
    float viscosity_ = 0.0f;
    float time_step_ = 0.001f;
    float damping_ = 0.99f;
};

} // namespace sph_integration
