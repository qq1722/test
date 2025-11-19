#pragma once
#include <vector>
#include <glm/glm.hpp>

class Boundary 
{
public:
    // 构造函数：支持外环和多个内环（洞）
    // 如果只提供外环，holes使用默认值（空列表）
    Boundary(const std::vector<glm::vec2>& outer_loop, 
             const std::vector<std::vector<glm::vec2>>& holes = {});
    
    // 判断一个点是否在边界内（考虑洞）
    // 点必须在外环内，但在所有内环外
    bool is_inside(const glm::vec2& point) const;
    
    // 获取外环顶点
    const std::vector<glm::vec2>& get_vertices() const;
    const std::vector<glm::vec2>& get_outer_boundary() const; // <-- [修改]
    
    // 获取内环（洞）列表
    const std::vector<std::vector<glm::vec2>>& get_holes() const { return holes_; }
    
    // 获取边界框（AABB）
    const glm::vec4& get_aabb() const;
    // --- [新增] 核心修复功能 ---
    // 计算点 p 到整个边界系统（外环或任意内洞）的最近点
    glm::vec2 get_closest_point(const glm::vec2& p) const;
    // [新增] 获取最近点 *以及* 该处的切线单位向量
    void get_closest_point_and_tangent(const glm::vec2& p, glm::vec2& out_closest, glm::vec2& out_tangent) const;

private:
    std::vector<glm::vec2> vertices_;  // 外环
    std::vector<std::vector<glm::vec2>> holes_;  // 内环（洞）列表
    glm::vec4 aabb_;
    
    // 辅助函数：判断点是否在多边形内（射线法）
    bool is_point_in_polygon(const glm::vec2& point, const std::vector<glm::vec2>& polygon) const;
    
    void calculate_aabb();
};

