#include "background_grid_sph.h"
#include "boundary_sph.h"
#include <algorithm>
#include <vector>
#include <cfloat>

// Forward declaration
glm::vec2 closest_point_on_polygon(const glm::vec2& p, const std::vector<glm::vec2>& vertices);

BackgroundGrid::BackgroundGrid(const Boundary& boundary, float grid_cell_size) {
    cell_size_ = grid_cell_size;
    const auto& aabb = boundary.get_aabb();
    min_coords_ = { aabb.x, aabb.y };
    width_ = static_cast<int>((aabb.z - aabb.x) / cell_size_) + 3;
    height_ = static_cast<int>((aabb.w - aabb.y) / cell_size_) + 3;
    target_size_field_.resize(width_ * height_);
    target_direction_field_.resize(width_ * height_, { 1.0f, 0.0f });

    compute_fields(boundary);
}

void BackgroundGrid::compute_fields(const Boundary& boundary) {
    const auto& boundary_vertices = boundary.get_vertices();
    if (boundary_vertices.size() < 2) return;

    std::vector<float> boundary_target_sizes(boundary_vertices.size());
    float h_min = cell_size_ * 0.5f;
    float h_max = cell_size_ * 2.0f;
    for (size_t i = 0; i < boundary_vertices.size(); ++i) {
        const auto& p_prev = boundary_vertices[(i + boundary_vertices.size() - 1) % boundary_vertices.size()];
        const auto& p_curr = boundary_vertices[i];
        const auto& p_next = boundary_vertices[(i + 1) % boundary_vertices.size()];
        glm::vec2 v1 = glm::normalize(p_curr - p_prev);
        glm::vec2 v2 = glm::normalize(p_next - p_curr);
        float dot_product = glm::dot(v1, v2);
        dot_product = std::max(-1.0f, std::min(1.0f, dot_product));
        float angle = std::acos(dot_product);
        float t = angle / 3.1415926535f;
        boundary_target_sizes[i] = glm::mix(h_min, h_max, t * t);
    }

    std::vector<float> sdf(width_ * height_, FLT_MAX);
    for (int y = 0; y < height_; ++y) {
        for (int x = 0; x < width_; ++x) {
            glm::vec2 grid_pos = min_coords_ + glm::vec2(x * cell_size_, y * cell_size_);
            float dist_to_boundary = glm::distance(grid_pos, closest_point_on_polygon(grid_pos, boundary_vertices));
            sdf[y * width_ + x] = boundary.is_inside(grid_pos) ? dist_to_boundary : -dist_to_boundary;

            float influence_radius = h_max * 5.0f;
            float t = std::min(dist_to_boundary / influence_radius, 1.0f);
            target_size_field_[y * width_ + x] = glm::mix(h_min, h_max, t * t);
        }
    }

    for (int y = 1; y < height_ - 1; ++y) {
        for (int x = 1; x < width_ - 1; ++x) {
            float grad_x = (sdf[y * width_ + (x + 1)] - sdf[y * width_ + (x - 1)]) / (2.0f * cell_size_);
            float grad_y = (sdf[(y + 1) * width_ + x] - sdf[(y - 1) * width_ + x]) / (2.0f * cell_size_);
            glm::vec2 grad = { grad_x, grad_y };
            if (glm::length(grad) > 1e-6f) {
                glm::vec2 tangent = { -grad.y, grad.x };
                target_direction_field_[y * width_ + x] = glm::normalize(tangent);
            }
        }
    }
}

float BackgroundGrid::get_target_size(const glm::vec2& pos) const {
    glm::vec2 local_pos = (pos - min_coords_) / cell_size_;
    int x0 = static_cast<int>(local_pos.x);
    int y0 = static_cast<int>(local_pos.y);
    x0 = std::max(0, std::min(x0, width_ - 2));
    y0 = std::max(0, std::min(y0, height_ - 2));
    int x1 = x0 + 1;
    int y1 = y0 + 1;
    float tx = local_pos.x - x0;
    float ty = local_pos.y - y0;
    float s00 = target_size_field_[y0 * width_ + x0];
    float s10 = target_size_field_[y0 * width_ + x1];
    float s01 = target_size_field_[y1 * width_ + x0];
    float s11 = target_size_field_[y1 * width_ + x1];
    float s_y0 = glm::mix(s00, s10, tx);
    float s_y1 = glm::mix(s01, s11, tx);
    return glm::mix(s_y0, s_y1, ty);
}

glm::vec2 BackgroundGrid::get_target_direction(const glm::vec2& pos) const {
    glm::vec2 local_pos = (pos - min_coords_) / cell_size_;
    int x0 = static_cast<int>(local_pos.x);
    int y0 = static_cast<int>(local_pos.y);
    x0 = std::max(0, std::min(x0, width_ - 2));
    y0 = std::max(0, std::min(y0, height_ - 2));
    int x1 = x0 + 1;
    int y1 = y0 + 1;
    float tx = local_pos.x - x0;
    float ty = local_pos.y - y0;
    glm::vec2 d00 = target_direction_field_[y0 * width_ + x0];
    glm::vec2 d10 = target_direction_field_[y0 * width_ + x1];
    glm::vec2 d01 = target_direction_field_[y1 * width_ + x0];
    glm::vec2 d11 = target_direction_field_[y1 * width_ + x1];
    glm::vec2 d_y0 = glm::mix(d00, d10, tx);
    glm::vec2 d_y1 = glm::mix(d01, d11, tx);
    return glm::normalize(glm::mix(d_y0, d_y1, ty));
}

