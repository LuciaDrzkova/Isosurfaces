#pragma once

#include "MarchingCubes.hpp"

#include <GL/glew.h>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <chrono>
#include <utility>
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
    using ScalarField = std::function<double(double, double, double)>;

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

    // Configure the exact sampling grid used by the initial MC extraction.
    // This does not rebuild the current mesh; the first mesh was already built
    // by main.cpp. Pressing "Apply Sampling Grid" in the viewer rebuilds it.
    void set_marching_cubes(
        ScalarField scalar_field,
        const iso::mc::Bounds& bounds,
        std::size_t nx,
        std::size_t ny,
        std::size_t nz,
        double isovalue = 0.0)
    {
        scalar_field_ = std::move(scalar_field);
        mc_bounds_ = bounds;
        mc_isovalue_ = isovalue;

        nx_ = clamp_resolution(nx);
        ny_ = clamp_resolution(ny);
        nz_ = clamp_resolution(nz);

        pending_nx_ = static_cast<int>(nx_);
        pending_ny_ = static_cast<int>(ny_);
        pending_nz_ = static_cast<int>(nz_);

        build_voxel_grid();
    }

    // Configure only the grid visualization. This is used by the local
    // unfolding viewer, where the surface itself is not regenerated here.
    void set_voxel_grid(const pmp::Point& min,
                        const pmp::Point& max,
                        int nx,
                        int ny,
                        int nz)
    {
        mc_bounds_.min = {
            static_cast<double>(min[0]),
            static_cast<double>(min[1]),
            static_cast<double>(min[2])};
        mc_bounds_.max = {
            static_cast<double>(max[0]),
            static_cast<double>(max[1]),
            static_cast<double>(max[2])};

        scalar_field_ = {};

        nx_ = clamp_resolution(static_cast<std::size_t>(std::max(2, nx)));
        ny_ = clamp_resolution(static_cast<std::size_t>(std::max(2, ny)));
        nz_ = clamp_resolution(static_cast<std::size_t>(std::max(2, nz)));

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

        ImGui::Text("Sampling grid: %zu x %zu x %zu",
                    nx_, ny_, nz_);

        ImGui::Text("Voxel cells: %zu x %zu x %zu",
                    nx_ - 1,
                    ny_ - 1,
                    nz_ - 1);

        const double dx =
            (nx_ > 1)
                ? (mc_bounds_.max[0] - mc_bounds_.min[0]) /
                      static_cast<double>(nx_ - 1)
                : 0.0;
        const double dy =
            (ny_ > 1)
                ? (mc_bounds_.max[1] - mc_bounds_.min[1]) /
                      static_cast<double>(ny_ - 1)
                : 0.0;
        const double dz =
            (nz_ > 1)
                ? (mc_bounds_.max[2] - mc_bounds_.min[2]) /
                      static_cast<double>(nz_ - 1)
                : 0.0;

        ImGui::Text("Voxel spacing: %.6g, %.6g, %.6g", dx, dy, dz);

        if (scalar_field_)
        {
            ImGui::Spacing();
            ImGui::TextUnformatted("Change MC sampling resolution");

            ImGui::PushItemWidth(110.0f);
            ImGui::InputInt("Nx", &pending_nx_);
            ImGui::InputInt("Ny", &pending_ny_);
            ImGui::InputInt("Nz", &pending_nz_);
            ImGui::PopItemWidth();

            pending_nx_ = std::clamp(pending_nx_, 2, kMaxResolution);
            pending_ny_ = std::clamp(pending_ny_, 2, kMaxResolution);
            pending_nz_ = std::clamp(pending_nz_, 2, kMaxResolution);

            if (ImGui::Button("Apply Sampling Grid"))
                rebuild_marching_cubes();

            ImGui::SameLine();
            if (ImGui::Button("Use Same N"))
            {
                pending_ny_ = pending_nx_;
                pending_nz_ = pending_nx_;
            }

            if (last_rebuild_ms_ >= 0.0)
            {
                ImGui::Text("Last MC rebuild: %.3f ms", last_rebuild_ms_);
                ImGui::Text("Surface: %zu vertices, %zu triangles",
                            last_vertex_count_, last_triangle_count_);
            }

            if (!error_message_.empty())
                ImGui::TextWrapped("Error: %s", error_message_.c_str());
        }
        else
        {
            ImGui::TextUnformatted("Surface regeneration is unavailable for this viewer.");
        }
    }

    void draw(const std::string& draw_mode) override
    {
        pmp::MeshViewer::draw(draw_mode);

        if (!show_voxel_grid_ || grid_indices_.empty())
            return;

        if (!voxel_grid_initialized_)
            initialize_voxel_grid_rendering();

        voxel_shader_.use();
        voxel_shader_.set_uniform("projection_matrix", projection_matrix_);
        voxel_shader_.set_uniform("modelview_matrix", modelview_matrix_);

        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
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

    static std::size_t clamp_resolution(std::size_t value)
    {
        return std::clamp<std::size_t>(value, 2, kMaxResolution);
    }

    void rebuild_marching_cubes()
    {
        if (!scalar_field_)
            return;

        error_message_.clear();

        const std::size_t new_nx = clamp_resolution(
            static_cast<std::size_t>(pending_nx_));
        const std::size_t new_ny = clamp_resolution(
            static_cast<std::size_t>(pending_ny_));
        const std::size_t new_nz = clamp_resolution(
            static_cast<std::size_t>(pending_nz_));

        const auto start = std::chrono::steady_clock::now();

        try
        {
            iso::mc::Options options;
            options.nx = new_nx;
            options.ny = new_ny;
            options.nz = new_nz;
            options.isovalue = mc_isovalue_;

            const iso::mc::Mesh extracted = iso::mc::extract(
                scalar_field_, mc_bounds_, options);

            mesh_.clear();

            std::vector<pmp::Vertex> vertices;
            vertices.reserve(extracted.vertices.size());

            for (const auto& p : extracted.vertices)
            {
                vertices.push_back(mesh_.add_vertex(
                    pmp::Point(
                        static_cast<float>(p[0]),
                        static_cast<float>(p[1]),
                        static_cast<float>(p[2]))));
            }

            for (const auto& triangle : extracted.triangles)
            {
                mesh_.add_triangle(
                    vertices[triangle[0]],
                    vertices[triangle[1]],
                    vertices[triangle[2]]);
            }

            nx_ = new_nx;
            ny_ = new_ny;
            nz_ = new_nz;

            pending_nx_ = static_cast<int>(nx_);
            pending_ny_ = static_cast<int>(ny_);
            pending_nz_ = static_cast<int>(nz_);

            build_voxel_grid();
            update_mesh();

            last_vertex_count_ = extracted.vertices.size();
            last_triangle_count_ = extracted.triangles.size();
        }
        catch (const std::exception& e)
        {
            error_message_ = e.what();
        }

        const auto end = std::chrono::steady_clock::now();
        last_rebuild_ms_ =
            std::chrono::duration<double, std::milli>(end - start).count();
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

        if (nx_ < 2 || ny_ < 2 || nz_ < 2)
        {
            voxel_grid_initialized_ = false;
            return;
        }

        const double dx =
            (mc_bounds_.max[0] - mc_bounds_.min[0]) /
            static_cast<double>(nx_ - 1);
        const double dy =
            (mc_bounds_.max[1] - mc_bounds_.min[1]) /
            static_cast<double>(ny_ - 1);
        const double dz =
            (mc_bounds_.max[2] - mc_bounds_.min[2]) /
            static_cast<double>(nz_ - 1);

        auto add_vertex = [&](double x, double y, double z) -> unsigned int
        {
            grid_vertices_.push_back({
                static_cast<float>(x),
                static_cast<float>(y),
                static_cast<float>(z)});
            return static_cast<unsigned int>(grid_vertices_.size() - 1);
        };

        auto add_line = [&](double ax, double ay, double az,
                            double bx, double by, double bz)
        {
            const unsigned int ia = add_vertex(ax, ay, az);
            const unsigned int ib = add_vertex(bx, by, bz);
            grid_indices_.push_back(ia);
            grid_indices_.push_back(ib);
        };

        const double xmin = mc_bounds_.min[0];
        const double xmax = mc_bounds_.max[0];
        const double ymin = mc_bounds_.min[1];
        const double ymax = mc_bounds_.max[1];
        const double zmin = mc_bounds_.min[2];
        const double zmax = mc_bounds_.max[2];

        // X-parallel grid edges.
        for (std::size_t j = 0; j < ny_; ++j)
        {
            const double y = ymin + static_cast<double>(j) * dy;
            for (std::size_t k = 0; k < nz_; ++k)
            {
                const double z = zmin + static_cast<double>(k) * dz;
                add_line(xmin, y, z, xmax, y, z);
            }
        }

        // Y-parallel grid edges.
        for (std::size_t i = 0; i < nx_; ++i)
        {
            const double x = xmin + static_cast<double>(i) * dx;
            for (std::size_t k = 0; k < nz_; ++k)
            {
                const double z = zmin + static_cast<double>(k) * dz;
                add_line(x, ymin, z, x, ymax, z);
            }
        }

        // Z-parallel grid edges.
        for (std::size_t i = 0; i < nx_; ++i)
        {
            const double x = xmin + static_cast<double>(i) * dx;
            for (std::size_t j = 0; j < ny_; ++j)
            {
                const double y = ymin + static_cast<double>(j) * dy;
                add_line(x, y, zmin, x, y, zmax);
            }
        }

        voxel_grid_initialized_ = false;
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
                frag_color = vec4(0.2, 0.2, 0.2, 0.35);
            }
        )GLSL";

        voxel_shader_.source(vertex_shader_source, fragment_shader_source);

        glGenVertexArrays(1, &grid_vao_);
        glGenBuffers(1, &grid_vbo_);
        glGenBuffers(1, &grid_ebo_);

        glBindVertexArray(grid_vao_);

        glBindBuffer(GL_ARRAY_BUFFER, grid_vbo_);
        glBufferData(
            GL_ARRAY_BUFFER,
            static_cast<GLsizeiptr>(
                grid_vertices_.size() * sizeof(LineVertex)),
            grid_vertices_.data(),
            GL_STATIC_DRAW);

        glVertexAttribPointer(
            0,
            3,
            GL_FLOAT,
            GL_FALSE,
            sizeof(LineVertex),
            reinterpret_cast<void*>(offsetof(LineVertex, x)));

        glEnableVertexAttribArray(0);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, grid_ebo_);
        glBufferData(
            GL_ELEMENT_ARRAY_BUFFER,
            static_cast<GLsizeiptr>(
                grid_indices_.size() * sizeof(unsigned int)),
            grid_indices_.data(),
            GL_STATIC_DRAW);

        glBindVertexArray(0);

        voxel_grid_initialized_ = true;
    }

private:
    ScalarField scalar_field_;
    iso::mc::Bounds mc_bounds_{{-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0}};
    double mc_isovalue_{0.0};

    std::size_t nx_{2};
    std::size_t ny_{2};
    std::size_t nz_{2};

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

    double last_rebuild_ms_{-1.0};
    std::size_t last_vertex_count_{0};
    std::size_t last_triangle_count_{0};
    std::string error_message_;
};

} // namespace iso
