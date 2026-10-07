#pragma once

#include "MarchingCubes.hpp"

#include <algorithm>
#include <cstddef>
#include <string>
#include <vector>

#include <pmp/visualization/gl.h>
#include <pmp/visualization/mesh_viewer.h>
#include <pmp/visualization/shader.h>

#include <imgui.h>

namespace iso
{

class MyViewer : public pmp::MeshViewer
{
public:
    MyViewer(const char* title, int w, int h)
        : pmp::MeshViewer(title, w, h)
    {
    }

    void load_mesh(const char* path) override
    {
        pmp::MeshViewer::load_mesh(path);
    }

    int run()
    {
        return pmp::MeshViewer::run();
    }

    // Sets the grid displayed in the viewer.
    // This never changes or rebuilds the loaded surface mesh.
    void set_voxel_grid(const pmp::Point& min,
                        const pmp::Point& max,
                        int nx,
                        int ny,
                        int nz)
    {
        grid_min_ = min;
        grid_max_ = max;

        nx_ = clamp_resolution(nx);
        ny_ = clamp_resolution(ny);
        nz_ = clamp_resolution(nz);

        pending_nx_ = static_cast<int>(nx_);
        pending_ny_ = static_cast<int>(ny_);
        pending_nz_ = static_cast<int>(nz_);

        build_voxel_grid();
    }

protected:
    void process_imgui() override
    {
        pmp::MeshViewer::process_imgui();

        ImGui::Separator();
        ImGui::TextUnformatted("Marching Cubes");

        ImGui::Checkbox("Show Voxel Grid", &show_voxel_grid_);

        ImGui::Spacing();
        ImGui::TextUnformatted("Displayed Sampling Grid");

        ImGui::PushItemWidth(100.0f);

        ImGui::InputInt("Nx", &pending_nx_);
        ImGui::InputInt("Ny", &pending_ny_);
        ImGui::InputInt("Nz", &pending_nz_);

        ImGui::PopItemWidth();

        pending_nx_ = std::clamp(
            pending_nx_, 2, kMaxResolution);
        pending_ny_ = std::clamp(
            pending_ny_, 2, kMaxResolution);
        pending_nz_ = std::clamp(
            pending_nz_, 2, kMaxResolution);

        if (ImGui::Button("Apply Grid"))
        {
            apply_displayed_grid();
        }

        ImGui::SameLine();

        if (ImGui::Button("Use Same N"))
        {
            pending_ny_ = pending_nx_;
            pending_nz_ = pending_nx_;
        }

        ImGui::Text(
            "Current grid: %zu x %zu x %zu",
            nx_,
            ny_,
            nz_);

        ImGui::Text(
            "Voxel cells: %zu x %zu x %zu",
            nx_ - 1,
            ny_ - 1,
            nz_ - 1);

        const double dx =
            (nx_ > 1)
                ? (grid_max_[0] - grid_min_[0]) /
                      static_cast<double>(nx_ - 1)
                : 0.0;

        const double dy =
            (ny_ > 1)
                ? (grid_max_[1] - grid_min_[1]) /
                      static_cast<double>(ny_ - 1)
                : 0.0;

        const double dz =
            (nz_ > 1)
                ? (grid_max_[2] - grid_min_[2]) /
                      static_cast<double>(nz_ - 1)
                : 0.0;

        ImGui::Text(
            "Voxel spacing: %.6g, %.6g, %.6g",
            dx,
            dy,
            dz);
    }

    void draw(const std::string& draw_mode) override
    {
        pmp::MeshViewer::draw(draw_mode);

        if (!show_voxel_grid_ || grid_indices_.empty())
            return;

        if (!voxel_grid_initialized_)
            initialize_voxel_grid_rendering();

        voxel_shader_.use();

        voxel_shader_.set_uniform(
            "projection_matrix",
            projection_matrix_);

        voxel_shader_.set_uniform(
            "modelview_matrix",
            modelview_matrix_);

        glEnable(GL_BLEND);
        glBlendFunc(
            GL_SRC_ALPHA,
            GL_ONE_MINUS_SRC_ALPHA);

        glDepthMask(GL_FALSE);

        glBindVertexArray(grid_vao_);

        glDrawElements(
            GL_LINES,
            static_cast<GLsizei>(grid_indices_.size()),
            GL_UNSIGNED_INT,
            nullptr);

        glBindVertexArray(0);

        glDepthMask(GL_TRUE);
        glDisable(GL_BLEND);

        voxel_shader_.disable();
    }

private:
    static constexpr int kMaxResolution = 513;

    static std::size_t clamp_resolution(int value)
    {
        return std::clamp(
            static_cast<std::size_t>(std::max(2, value)),
            std::size_t(2),
            static_cast<std::size_t>(kMaxResolution));
    }

    void apply_displayed_grid()
    {
        nx_ = clamp_resolution(pending_nx_);
        ny_ = clamp_resolution(pending_ny_);
        nz_ = clamp_resolution(pending_nz_);

        pending_nx_ = static_cast<int>(nx_);
        pending_ny_ = static_cast<int>(ny_);
        pending_nz_ = static_cast<int>(nz_);

        build_voxel_grid();
    }

    struct LineVertex
    {
        float x;
        float y;
        float z;
    };

    void build_voxel_grid()
    {
        grid_vertices_.clear();
        grid_indices_.clear();

        voxel_grid_initialized_ = false;

        if (nx_ < 2 || ny_ < 2 || nz_ < 2)
            return;

        const double dx =
            (grid_max_[0] - grid_min_[0]) /
            static_cast<double>(nx_ - 1);

        const double dy =
            (grid_max_[1] - grid_min_[1]) /
            static_cast<double>(ny_ - 1);

        const double dz =
            (grid_max_[2] - grid_min_[2]) /
            static_cast<double>(nz_ - 1);

        auto add_vertex =
            [&](double x, double y, double z) -> unsigned int
        {
            grid_vertices_.push_back({
                static_cast<float>(x),
                static_cast<float>(y),
                static_cast<float>(z)
            });

            return static_cast<unsigned int>(
                grid_vertices_.size() - 1);
        };

        auto add_line =
            [&](double ax, double ay, double az,
                double bx, double by, double bz)
        {
            const unsigned int a =
                add_vertex(ax, ay, az);

            const unsigned int b =
                add_vertex(bx, by, bz);

            grid_indices_.push_back(a);
            grid_indices_.push_back(b);
        };

        const double xmin = grid_min_[0];
        const double xmax = grid_max_[0];

        const double ymin = grid_min_[1];
        const double ymax = grid_max_[1];

        const double zmin = grid_min_[2];
        const double zmax = grid_max_[2];

        // X-parallel edges.
        for (std::size_t j = 0; j < ny_; ++j)
        {
            const double y =
                ymin + static_cast<double>(j) * dy;

            for (std::size_t k = 0; k < nz_; ++k)
            {
                const double z =
                    zmin + static_cast<double>(k) * dz;

                add_line(
                    xmin, y, z,
                    xmax, y, z);
            }
        }

        // Y-parallel edges.
        for (std::size_t i = 0; i < nx_; ++i)
        {
            const double x =
                xmin + static_cast<double>(i) * dx;

            for (std::size_t k = 0; k < nz_; ++k)
            {
                const double z =
                    zmin + static_cast<double>(k) * dz;

                add_line(
                    x, ymin, z,
                    x, ymax, z);
            }
        }

        // Z-parallel edges.
        for (std::size_t i = 0; i < nx_; ++i)
        {
            const double x =
                xmin + static_cast<double>(i) * dx;

            for (std::size_t j = 0; j < ny_; ++j)
            {
                const double y =
                    ymin + static_cast<double>(j) * dy;

                add_line(
                    x, y, zmin,
                    x, y, zmax);
            }
        }

        if (grid_vao_ == 0)
        {
            initialize_voxel_grid_rendering();
        }
        else
        {
            update_voxel_grid_buffers();
        }
    }

    void initialize_voxel_grid_rendering()
    {
        const char* vertex_shader_source = R"GLSL(
            #version 330 core

            layout(location = 0) in vec3 vertex;

            uniform mat4 projection_matrix;
            uniform mat4 modelview_matrix;

            void main()
            {
                gl_Position =
                    projection_matrix *
                    modelview_matrix *
                    vec4(vertex, 1.0);
            }
        )GLSL";

        const char* fragment_shader_source = R"GLSL(
            #version 330 core

            out vec4 frag_color;

            void main()
            {
                frag_color =
                    vec4(0.2, 0.2, 0.2, 0.35);
            }
        )GLSL";

        voxel_shader_.source(
            vertex_shader_source,
            fragment_shader_source);

        glGenVertexArrays(1, &grid_vao_);
        glGenBuffers(1, &grid_vbo_);
        glGenBuffers(1, &grid_ebo_);

        glBindVertexArray(grid_vao_);

        glBindBuffer(
            GL_ARRAY_BUFFER,
            grid_vbo_);

        glBufferData(
            GL_ARRAY_BUFFER,
            static_cast<GLsizeiptr>(
                grid_vertices_.size() *
                sizeof(LineVertex)),
            grid_vertices_.data(),
            GL_DYNAMIC_DRAW);

        glVertexAttribPointer(
            0,
            3,
            GL_FLOAT,
            GL_FALSE,
            sizeof(LineVertex),
            reinterpret_cast<void*>(0));

        glEnableVertexAttribArray(0);

        glBindBuffer(
            GL_ELEMENT_ARRAY_BUFFER,
            grid_ebo_);

        glBufferData(
            GL_ELEMENT_ARRAY_BUFFER,
            static_cast<GLsizeiptr>(
                grid_indices_.size() *
                sizeof(unsigned int)),
            grid_indices_.data(),
            GL_DYNAMIC_DRAW);

        glBindVertexArray(0);

        voxel_grid_initialized_ = true;
    }

    void update_voxel_grid_buffers()
    {
        if (grid_vao_ == 0)
            return;

        glBindVertexArray(grid_vao_);

        glBindBuffer(
            GL_ARRAY_BUFFER,
            grid_vbo_);

        glBufferData(
            GL_ARRAY_BUFFER,
            static_cast<GLsizeiptr>(
                grid_vertices_.size() *
                sizeof(LineVertex)),
            grid_vertices_.data(),
            GL_DYNAMIC_DRAW);

        glBindBuffer(
            GL_ELEMENT_ARRAY_BUFFER,
            grid_ebo_);

        glBufferData(
            GL_ELEMENT_ARRAY_BUFFER,
            static_cast<GLsizeiptr>(
                grid_indices_.size() *
                sizeof(unsigned int)),
            grid_indices_.data(),
            GL_DYNAMIC_DRAW);

        glBindVertexArray(0);

        voxel_grid_initialized_ = true;
    }

private:
    // Bounds of the original MC sampling domain.
    // These stay fixed while the displayed sampling resolution changes.
    pmp::Point grid_min_{0.0f, 0.0f, 0.0f};
    pmp::Point grid_max_{1.0f, 1.0f, 1.0f};

    // Currently displayed sampling-grid resolution.
    std::size_t nx_{2};
    std::size_t ny_{2};
    std::size_t nz_{2};

    // Values currently entered in the ImGui controls.
    int pending_nx_{2};
    int pending_ny_{2};
    int pending_nz_{2};

    bool show_voxel_grid_{false};
    bool voxel_grid_initialized_{false};

    std::vector<LineVertex> grid_vertices_;
    std::vector<unsigned int> grid_indices_;

    GLuint grid_vao_{0};
    GLuint grid_vbo_{0};
    GLuint grid_ebo_{0};

    pmp::Shader voxel_shader_;
};

} // namespace iso