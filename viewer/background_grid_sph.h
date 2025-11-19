#pragma once
#include <vector>
#include <glm/glm.hpp>

class Boundary;

class BackgroundGrid {
public:
    BackgroundGrid(const Boundary& boundary, float grid_cell_size);

    float get_target_size(const glm::vec2& pos) const;
    glm::vec2 get_target_direction(const glm::vec2& pos) const;

    int get_width() const { return width_; }
    int get_height() const { return height_; }
    float get_cell_size() const { return cell_size_; }
    glm::vec2 get_min_coords() const { return min_coords_; }
    const std::vector<float>& get_target_size_field() const { return target_size_field_; }

private:
    void compute_fields(const Boundary& boundary);

    glm::vec2 min_coords_;
    float cell_size_;
    int width_, height_;
    std::vector<float> target_size_field_;
    std::vector<glm::vec2> target_direction_field_;
};

