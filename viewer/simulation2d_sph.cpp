#include "simulation2d_sph.h"
#include "boundary_sph.h"
#include "background_grid_sph.h"
#include <random>
#include <algorithm>
#include <iostream>

constexpr float PI = 3.1415926535f;

glm::vec2 closest_point_on_polygon(const glm::vec2& p, const std::vector<glm::vec2>& vertices);

glm::vec2 Simulation2D::transform_to_local(const glm::vec2& vec, const glm::mat2& rot_matrix) const {
    return glm::transpose(rot_matrix) * vec;
}
//Simulation2D::Simulation2D(const Boundary& boundary) : boundary_(boundary) {
//    const auto& aabb = boundary.get_aabb();
//    float domain_width = aabb.z - aabb.x;
//    // 背景网格的分辨率决定了 h_t 场的精度，但不直接决定粒子数量
//    float grid_cell_size = domain_width / 20.0f;
//    grid_ = std::make_unique<BackgroundGrid>(boundary, grid_cell_size);
//
//    // 初始化粒子 (现在调用的是新的四叉树版本)
//    initialize_particles(boundary);
//}
Simulation2D::Simulation2D(const Boundary& boundary) : boundary_(boundary) {
    const auto& aabb = boundary.get_aabb();
    float domain_width = aabb.z - aabb.x;
    // 背景网格的分辨率决定了 h_t 场的精度，但不直接决定粒子数量
    float grid_cell_size = domain_width / 20.0f;
    grid_ = std::make_unique<BackgroundGrid>(boundary, grid_cell_size);
    // 初始化粒子 (现在调用的是新的四叉树版本)
    initialize_particles(boundary);
}
//void Simulation2D::compute_forces() {
//    for (auto& p : particles_) { p.force = glm::vec2(0.0f); }
//
//    for (int i = 0; i < num_particles_; ++i) {
//        for (int j = i + 1; j < num_particles_; ++j) {
//            glm::vec2 diff_global = particles_[i].position - particles_[j].position;
//            float h_avg = (particles_[i].smoothing_h + particles_[j].smoothing_h) * 0.5f;
//
//            // 转换到粒子 i 的局部坐标系
//            glm::vec2 diff_local_i = transform_to_local(diff_global, particles_[i].rotation);
//            float r_inf = l_inf_norm(diff_local_i);
//
//            if (r_inf < 2.0f * h_avg) { // 2.0 是紧支集半径
//                float q = r_inf / h_avg;
//                if (q > 1e-6) {
//                    float rho_t_i = particles_[i].target_density;
//                    float rho_t_j = particles_[j].target_density;
//                    float P_term = (stiffness_ / (rho_t_i * rho_t_i)) + (stiffness_ / (rho_t_j * rho_t_j));
//                    float W_grad_mag = wendland_c6_kernel_derivative(q, h_avg);
//
//                    // ============================================================
//                    // [核心修复] 使用 L-infinity 的真实梯度方向
//                    // ============================================================
//                    glm::vec2 grad_direction(0.0f);
//                    float abs_x = std::abs(diff_local_i.x);
//                    float abs_y = std::abs(diff_local_i.y);
//
//                    // 容差，处理对角线情况 (此时 x 和 y 差不多大)
//                    // 如果在该容差内，可以分配给两个方向，或者平滑过渡
//                    float epsilon = 1e-5f;
//
//                    if (abs_x > abs_y + epsilon) {
//                        // X 轴主导：力只在 X 方向
//                        grad_direction = glm::vec2((diff_local_i.x > 0) ? 1.0f : -1.0f, 0.0f);
//                    }
//                    else if (abs_y > abs_x + epsilon) {
//                        // Y 轴主导：力只在 Y 方向
//                        grad_direction = glm::vec2(0.0f, (diff_local_i.y > 0) ? 1.0f : -1.0f);
//                    }
//                    else {
//                        // 对角线情况：两个方向都有 (避免除以零或突变)
//                        // 这种情况下 L_inf 的梯度是多值的，我们取归一化向量作为近似
//                        grad_direction = glm::normalize(diff_local_i);
//                    }
//
//                    // [原来的错误代码] 导致了不稳定性
//                    // glm::vec2 normalized_diff_local = diff_local_i / r_inf; 
//
//                    // 计算局部力 (注意：力沿着 grad_direction)
//                    glm::vec2 force_local = -mass_ * mass_ * P_term * W_grad_mag * grad_direction;
//
//                    // 转换回全局坐标系
//                    glm::vec2 force_global = particles_[i].rotation * force_local;
//
//                    particles_[i].force += force_global;
//                    particles_[j].force -= force_global;
//                }
//            }
//        }
//    }
//}
void Simulation2D::compute_forces() {
    for (auto& p : particles_) { p.force = glm::vec2(0.0f); }

    for (int i = 0; i < num_particles_; ++i) {
        for (int j = i + 1; j < num_particles_; ++j) {
            glm::vec2 diff_global = particles_[i].position - particles_[j].position;
            float h_avg = (particles_[i].smoothing_h + particles_[j].smoothing_h) * 0.5f;

            // 转换到粒子 i 的局部坐标系
            glm::vec2 diff_local_i = transform_to_local(diff_global, particles_[i].rotation);
            float r_inf = l_inf_norm(diff_local_i);

            if (r_inf < 2.0f * h_avg) { // 2.0 是紧支集半径
                float q = r_inf / h_avg;
                if (q > 1e-6) {
                    float rho_t_i = particles_[i].target_density;
                    float rho_t_j = particles_[j].target_density;
                    float P_term = (stiffness_ / (rho_t_i * rho_t_i)) + (stiffness_ / (rho_t_j * rho_t_j));
                    float W_grad_mag = wendland_c6_kernel_derivative(q, h_avg);

                    // ============================================================
                    // [核心修复] 使用 L-infinity 的真实梯度方向
                    // ============================================================
                    glm::vec2 grad_direction(0.0f);
                    float abs_x = std::abs(diff_local_i.x);
                    float abs_y = std::abs(diff_local_i.y);

                    // 容差，处理对角线情况 (此时 x 和 y 差不多大)
                    // 如果在该容差内，可以分配给两个方向，或者平滑过渡
                    float epsilon = 1e-5f;

                    if (abs_x > abs_y + epsilon) {
                        // X 轴主导：力只在 X 方向
                        grad_direction = glm::vec2((diff_local_i.x > 0) ? 1.0f : -1.0f, 0.0f);
                    }
                    else if (abs_y > abs_x + epsilon) {
                        // Y 轴主导：力只在 Y 方向
                        grad_direction = glm::vec2(0.0f, (diff_local_i.y > 0) ? 1.0f : -1.0f);
                    }
                    else {
                        // 对角线情况：两个方向都有 (避免除以零或突变)
                        // 这种情况下 L_inf 的梯度是多值的，我们取归一化向量作为近似
                        grad_direction = glm::normalize(diff_local_i);
                    }

                    // [原来的错误代码] 导致了不稳定性
                    // glm::vec2 normalized_diff_local = diff_local_i / r_inf; 

                    // 计算局部力 (注意：力沿着 grad_direction)
                    glm::vec2 force_local = -mass_ * mass_ * P_term * W_grad_mag * grad_direction;

                    // 转换回全局坐标系
                    glm::vec2 force_global = particles_[i].rotation * force_local;

                    particles_[i].force += force_global;
                    particles_[j].force -= force_global;
                }
            }
        }
    }
}

// --- [修改] 位置更新后，更新粒子的方向 ---
void Simulation2D::update_positions() {
    for (auto& p : particles_) {
        p.velocity += (p.force / mass_) * time_step_;
        p.velocity *= damping_;
        p.position += p.velocity * time_step_;

        // 从背景网格获取每个粒子的目标尺寸
        p.smoothing_h = grid_->get_target_size(p.position);
        p.target_density = 1.0f / (p.smoothing_h * p.smoothing_h);

        // 关键：更新粒子的旋转矩阵以对齐方向
        glm::vec2 target_dir = grid_->get_target_direction(p.position);
        glm::vec2 current_dir = p.rotation[0]; // 局部X轴

        // 使用平滑插值旋转到目标方向，防止突变
        glm::vec2 new_dir = glm::normalize(current_dir + (target_dir - current_dir) * 0.1f);

        p.rotation[0] = new_dir;
        p.rotation[1] = glm::vec2(-new_dir.y, new_dir.x); // 正交基
    }
}
//// [修改] 初始化入口
//void Simulation2D::initialize_particles(const Boundary& boundary) {
//    particles_.clear();
//    std::cout << "Initializing particles: Hybrid Method (Paper Boundary + Cartesian Interior)..." << std::endl;
//
//    // --- A. 生成边界粒子 (论文算法) ---
//    initialize_boundary_particles(boundary.get_outer_boundary());
//    for (const auto& hole : boundary.get_holes()) {
//        initialize_boundary_particles(hole);
//    }
//    std::cout << "  Boundary particles generated." << std::endl;
//
//    // --- B. 生成域内粒子 (你的笛卡尔算法) ---
//    // 构建覆盖全域的根正方形
//    const glm::vec4& aabb = boundary.get_aabb();
//    float width = aabb.z - aabb.x;
//    float height = aabb.w - aabb.y;
//    float max_dim = std::max(width, height);
//    glm::vec2 center_aabb = { (aabb.x + aabb.z) * 0.5f, (aabb.y + aabb.w) * 0.5f };
//
//    float root_size = max_dim * 1.2f; // 稍微扩大以覆盖边界
//    glm::vec2 root_min = center_aabb - glm::vec2(root_size * 0.5f);
//    glm::vec2 root_max = center_aabb + glm::vec2(root_size * 0.5f);
//
//    recursive_spawn_particles(root_min, root_max, boundary);
//    std::cout << "  Interior Cartesian particles generated." << std::endl;
//
//    num_particles_ = particles_.size();
//    positions_for_render_.resize(num_particles_);
//    std::cout << "Total particles: " << num_particles_ << std::endl;
//}
// [修改] 初始化入口
void Simulation2D::initialize_particles(const Boundary& boundary) {
    particles_.clear();
    std::cout << "Initializing particles: Hybrid Method (Paper Boundary + Cartesian Interior)..." << std::endl;

    // --- A. 生成边界粒子 (论文算法) ---
    initialize_boundary_particles(boundary.get_outer_boundary());
    for (const auto& hole : boundary.get_holes()) {
        initialize_boundary_particles(hole);
    }
    std::cout << "  Boundary particles generated." << std::endl;

    // --- B. 生成域内粒子 (四叉树递归算法) ---
    // 构建覆盖全域的根正方形
    const glm::vec4& aabb = boundary.get_aabb();
    float width = aabb.z - aabb.x;
    float height = aabb.w - aabb.y;
    float max_dim = std::max(width, height);
    glm::vec2 center_aabb = { (aabb.x + aabb.z) * 0.5f, (aabb.y + aabb.w) * 0.5f };

    float root_size = max_dim * 1.2f; // 稍微扩大以覆盖边界
    glm::vec2 root_min = center_aabb - glm::vec2(root_size * 0.5f);
    glm::vec2 root_max = center_aabb + glm::vec2(root_size * 0.5f);

    recursive_spawn_particles(root_min, root_max, boundary);
    std::cout << "  Interior Cartesian particles generated." << std::endl;

    num_particles_ = particles_.size();
    positions_for_render_.resize(num_particles_);
    std::cout << "Total particles: " << num_particles_ << std::endl;
}

float Simulation2D::wendland_c6_kernel(float q, float h) {
    if (q >= 0.0f && q < 2.0f) {
        float term = 1.0f - q / 2.0f;
        float term_sq = term * term;
        float term4 = term_sq * term_sq;
        float alpha_d = (78.0f / (28.0f * PI * h * h));
        return alpha_d * term4 * term4 * (4.0f * q * q * q + 6.25f * q * q + 4.0f * q + 1.0f);
    }
    return 0.0f;
}

float Simulation2D::wendland_c6_kernel_derivative(float q, float h) {
    if (q > 1e-6f && q < 2.0f) {
        float term = 1.0f - q / 2.0f;
        float term_sq = term * term;
        float term3 = term_sq * term;
        float term4 = term_sq * term_sq;
        float term7 = term3 * term4;
        float alpha_d = (78.0f / (28.0f * PI * h * h));
        return alpha_d * term7 * (-10.0f * q * q * q - 10.25f * q * q - 2.0f * q) / h;
    }
    return 0.0f;
}

float Simulation2D::l_inf_norm(const glm::vec2& v) const {
    return std::max(std::abs(v.x), std::abs(v.y));
}
void Simulation2D::handle_boundaries(const Boundary& boundary) {
    for (auto& p : particles_) {
        // 如果粒子出界（无论是在最外层外面，还是在内洞里面）
        if (!boundary.is_inside(p.position)) {

            glm::vec2 closest_pt;
            glm::vec2 tangent;

            // 1. 获取最近的边界点和对应的切线方向
            boundary.get_closest_point_and_tangent(p.position, closest_pt, tangent);

            // 2. 【位置修正】：强行将粒子"吸附"到边界线上
            p.position = closest_pt;

            // 3. 【速度修正】：实现滑动 (Sliding)
            // 将速度投影到切线方向。
            // 数学原理：v_new = (v_old · tangent) * tangent
            // 这样就去掉了垂直于边界的分量，只保留沿边界跑的分量。
            float v_dot_t = glm::dot(p.velocity, tangent);
            p.velocity = v_dot_t * tangent;

            // (可选) 如果你希望粒子在边界上移动时有一些"摩擦力"而慢慢停下
            // 可以乘一个阻尼系数，比如 0.9f。如果不乘，就是光滑滑动。
            // p.velocity *= 0.95f; 
        }
    }
}

void Simulation2D::step() {
    if (num_particles_ == 0) return;
    compute_forces();
    update_positions();
    handle_boundaries(boundary_);
    for (int i = 0; i < num_particles_; ++i) {
        positions_for_render_[i] = particles_[i].position;
    }
}

const std::vector<glm::vec2>& Simulation2D::get_particle_positions() const {
    return positions_for_render_;
}

float Simulation2D::get_kinetic_energy() const {
    float total_energy = 0.0f;
    for (const auto& p : particles_) {
        total_energy += 0.5f * mass_ * glm::dot(p.velocity, p.velocity);
    }
    return total_energy;
}



// ==========================================
// 1. [论文算法] 边界粒子生成 (Algorithm 1)
// ==========================================
void Simulation2D::initialize_boundary_particles(const std::vector<glm::vec2>& loop) {
    if (loop.size() < 2) return;

    float Q = 0.0f; // 累加器

    for (size_t i = 0; i < loop.size(); ++i) {
        glm::vec2 p0 = loop[i];
        glm::vec2 p1 = loop[(i + 1) % loop.size()];

        // 获取端点目标尺寸
        float h0 = grid_->get_target_size(p0);
        float h1 = grid_->get_target_size(p1);
        float h_bar = 0.5f * (h0 + h1);

        // 计算边长和质量度量
        float L = glm::distance(p0, p1);
        float m = L / h_bar;
        Q += m;

        int n = static_cast<int>(std::floor(Q));

        if (n > 0) {
            glm::vec2 dir = p1 - p0;
            for (int k = 0; k < n; ++k) {
                // 均匀插值
                float t = (float)k / (float)n;
                glm::vec2 pos = p0 + t * dir;

                float h_t = grid_->get_target_size(pos);
                Particle p;
                p.position = pos;
                p.smoothing_h = h_t;
                p.target_density = 1.0f / (h_t * h_t);
                p.is_boundary = true; // 【关键】标记为边界粒子
                p.velocity = glm::vec2(0.0f); // 边界粒子不动

                particles_.push_back(p);
            }
            Q -= n;
        }
    }
}

// ==========================================
// 2. [算法] 域内笛卡尔粒子生成 (四叉树递归)
// ==========================================
void Simulation2D::recursive_spawn_particles(glm::vec2 min_pt, glm::vec2 max_pt, const Boundary& boundary) {
    glm::vec2 center = (min_pt + max_pt) * 0.5f;
    float current_cell_size = max_pt.x - min_pt.x;
    float h_target = grid_->get_target_size(center);

    // 防止无限递归的最小尺寸
    float min_allowed_h = grid_->get_cell_size() * 0.2f;

    // 如果当前格子比目标尺寸大，继续分裂（笛卡尔加密）
    if (current_cell_size > std::max(h_target, min_allowed_h)) {
        recursive_spawn_particles(min_pt, center, boundary); // 左下
        recursive_spawn_particles({ center.x, min_pt.y }, { max_pt.x, center.y }, boundary); // 右下
        recursive_spawn_particles({ min_pt.x, center.y }, { center.x, max_pt.y }, boundary); // 左上
        recursive_spawn_particles(center, max_pt, boundary); // 右上
    }
    else {
        // 叶子节点：生成粒子
        // 【关键检查】只有在边界内部才生成流体粒子
        // 并且最好离边界有一点点距离，防止和边界粒子重叠太厉害
        if (boundary.is_inside(center)) {
            Particle p;
            p.position = center; // 完美的笛卡尔中心点
            p.velocity = glm::vec2(0.0f);
            p.force = glm::vec2(0.0f);
            p.smoothing_h = h_target;
            p.target_density = 1.0f / (h_target * h_target);
            p.is_boundary = false; // 【关键】标记为域内流体粒子

            // 初始旋转矩阵对齐坐标轴
            p.rotation = glm::mat2(1.0f);

            particles_.push_back(p);
        }
    }
}
