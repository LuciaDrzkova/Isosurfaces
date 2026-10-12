#include "ImplicitSurface.hpp"
#include "LocalMarchingCubes.hpp"
#include "MarchingCubes.hpp"
#include "SingularityDetector.hpp"
#include "app/MyViewer.h"

#include <GL/glew.h>
#include <pmp/io/io.h>
#include <pmp/surface_mesh.h>
#include <imgui.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{
struct CornerSample
{
    std::size_t i = 0, j = 0, k = 0;
    double x = 0.0, y = 0.0, z = 0.0, value = 0.0;
};

std::unique_ptr<ImplicitSurface> make_surface(const std::string& name, double c)
{
    if (name == "cone") return std::make_unique<ParameterizedConeQuadric>(c);
    if (name == "elliptic-cone") return std::make_unique<EllipticConeSurface>(c);
    if (name == "rotated-cone") return std::make_unique<RotatedEllipticConeSurface>(c);
    if (name == "x-cone") return std::make_unique<XAxisEllipticConeSurface>(c);
    if (name == "y-cone") return std::make_unique<YAxisEllipticConeSurface>(c);
    if (name == "tilted-cone") return std::make_unique<TiltedEllipticConeSurface>(c);
    if (name == "crossing-planes") return std::make_unique<IntersectingPlanesSurface>(c);
    if (name == "quartic") return std::make_unique<QuarticSaddleSurface>(c, 0.05);
    if (name == "sphere") return std::make_unique<SphereImplicitSurface>(1.0);
    throw std::invalid_argument("Unknown surface. Use cone, elliptic-cone, rotated-cone, x-cone, y-cone, tilted-cone, crossing-planes, quartic or sphere.");
}

std::vector<CornerSample> sample_grid(const ImplicitSurface& surface,
                                      const iso::mc::Bounds& bounds,
                                      std::size_t nx, std::size_t ny, std::size_t nz)
{
    if (nx < 2 || ny < 2 || nz < 2)
        throw std::invalid_argument("A 3D grid must have at least two sample points on each axis.");
    std::vector<CornerSample> samples;
    samples.reserve(nx * ny * nz);
    const double dx = (bounds.max[0] - bounds.min[0]) / static_cast<double>(nx - 1);
    const double dy = (bounds.max[1] - bounds.min[1]) / static_cast<double>(ny - 1);
    const double dz = (bounds.max[2] - bounds.min[2]) / static_cast<double>(nz - 1);
    for (std::size_t k = 0; k < nz; ++k)
        for (std::size_t j = 0; j < ny; ++j)
            for (std::size_t i = 0; i < nx; ++i)
            {
                CornerSample s;
                s.i = i; s.j = j; s.k = k;
                s.x = bounds.min[0] + static_cast<double>(i) * dx;
                s.y = bounds.min[1] + static_cast<double>(j) * dy;
                s.z = bounds.min[2] + static_cast<double>(k) * dz;
                s.value = surface.eval(s.x, s.y, s.z);
                samples.push_back(s);
            }
    return samples;
}

void write_samples_csv(const std::filesystem::path& path,
                       const std::vector<CornerSample>& samples)
{
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot write scalar sample CSV: " + path.string());
    out << "i,j,k,x,y,z,scalar_value,isovalue,inside_or_outside\n" << std::setprecision(17);
    for (const auto& s : samples)
        out << s.i << ',' << s.j << ',' << s.k << ',' << s.x << ',' << s.y << ',' << s.z << ','
            << s.value << ",0," << (s.value < 0.0 ? "negative" : (s.value > 0.0 ? "positive" : "on_isosurface")) << '\n';
}

// Split face fans that only touch at a vertex. This preserves the intended
// topological separation around singular points when converting to PMP.
pmp::SurfaceMesh to_pmp_mesh(const iso::mc::Mesh& source)
{
    using Index = iso::mc::Index;
    using Edge = std::pair<Index, Index>;
    pmp::SurfaceMesh mesh;
    const std::size_t face_count = source.triangles.size();
    if (!face_count) return mesh;
    std::vector<std::size_t> parent(face_count * 3), rank(face_count * 3, 0);
    for (std::size_t i = 0; i < parent.size(); ++i) parent[i] = i;
    auto root = [&](std::size_t x) {
        while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
        return x;
    };
    auto unite = [&](std::size_t a, std::size_t b) {
        a = root(a); b = root(b); if (a == b) return;
        if (rank[a] < rank[b]) std::swap(a, b);
        parent[b] = a; if (rank[a] == rank[b]) ++rank[a];
    };
    auto edge_key = [](Index a, Index b) -> Edge { if (a > b) std::swap(a, b); return {a,b}; };
    std::map<Edge, std::vector<std::size_t>> edge_faces;
    for (std::size_t f = 0; f < face_count; ++f) {
        const auto& t = source.triangles[f];
        edge_faces[edge_key(t[0],t[1])].push_back(f);
        edge_faces[edge_key(t[1],t[2])].push_back(f);
        edge_faces[edge_key(t[2],t[0])].push_back(f);
    }
    auto corner = [&](std::size_t face, Index v) {
        const auto& t = source.triangles[face];
        for (std::size_t c = 0; c < 3; ++c) if (t[c] == v) return 3 * face + c;
        throw std::runtime_error("Invalid triangle corner while converting patch mesh");
    };
    for (const auto& item : edge_faces) {
        const auto& e = item.first; const auto& fs = item.second;
        if (fs.size() == 2) {
            unite(corner(fs[0], e.first), corner(fs[1], e.first));
            unite(corner(fs[0], e.second), corner(fs[1], e.second));
        }
    }
    std::map<std::size_t, pmp::Vertex> vertices;
    for (std::size_t f = 0; f < face_count; ++f) {
        std::array<pmp::Vertex,3> tri;
        for (std::size_t c = 0; c < 3; ++c) {
            const std::size_t node = 3*f+c, r = root(node);
            auto it = vertices.find(r);
            if (it == vertices.end()) {
                const auto& p = source.vertices[source.triangles[f][c]];
                auto v = mesh.add_vertex(pmp::Point(float(p[0]), float(p[1]), float(p[2])));
                it = vertices.emplace(r, v).first;
            }
            tri[c] = it->second;
        }
        mesh.add_triangle(tri[0], tri[1], tri[2]);
    }
    return mesh;
}

bool flag(int argc, char** argv, const std::string& wanted)
{ for (int i=1;i<argc;++i) if (argv[i] == wanted) return true; return false; }
std::string value(int argc, char** argv, const std::string& key, const std::string& fallback)
{ for (int i=1;i+1<argc;++i) if (argv[i] == key) return argv[i+1]; return fallback; }

class PatchViewer final : public iso::MyViewer
{
public:
    PatchViewer(const char* title, int w, int h,
                const ImplicitSurface& surface,
                iso::mc::Mesh original_patch_mesh,
                std::vector<CornerSample> initial_samples,
                iso::mc::Bounds bounds,
                std::size_t nx, std::size_t ny, std::size_t nz,
                iso::mc::Point singularity_center,
                std::filesystem::path output_dir,
                std::string surface_name)
        : iso::MyViewer(title, w, h), surface_(&surface),
          original_patch_mesh_(std::move(original_patch_mesh)),
          active_mesh_(original_patch_mesh_),
          samples_(initial_samples), original_samples_(std::move(initial_samples)),
          bounds_(bounds), singularity_center_(singularity_center),
          original_nx_(nx), original_ny_(ny), original_nz_(nz),
          nx_(nx), ny_(ny), nz_(nz), output_dir_(std::move(output_dir)),
          surface_name_(std::move(surface_name))
    {
        if (nx_ < 2 || ny_ < 2 || nz_ < 2 || samples_.empty())
            throw std::invalid_argument("The local patch grid must have at least 2 samples per axis.");
        const auto max_cells = std::max({nx_ - 1, ny_ - 1, nz_ - 1});
        pending_grid_cells_ = static_cast<int>(std::clamp<std::size_t>(max_cells, 2, kMaxGridCells));
        selected_sample_index_ = nearest_sample_to_singularity(samples_);
        build_voxel_grid();
    }

    ~PatchViewer() override
    {
        // This derived destructor runs before the viewer's OpenGL context is torn down.
        if (grid_ebo_) glDeleteBuffers(1, &grid_ebo_);
        if (grid_vbo_) glDeleteBuffers(1, &grid_vbo_);
        if (grid_vao_) glDeleteVertexArrays(1, &grid_vao_);
        if (marker_vbo_) glDeleteBuffers(1, &marker_vbo_);
        if (marker_vao_) glDeleteVertexArrays(1, &marker_vao_);
    }

protected:
    void process_imgui() override
    {
        // Call PMP directly so the generic whole-grid MC controls from MyViewer
        // do not imply that they rebuild the proposed local voxel patch.
        pmp::MeshViewer::process_imgui();
        ImGui::Separator();
        ImGui::TextUnformatted("Local patch visualization");
        ImGui::Text("Surface: %s", surface_name_.c_str());
        ImGui::Text("Displayed mode: %s", mode_ == DisplayMode::LocalVoxelRefinement
                    ? "LocalMarchingCubes research patch"
                    : "Classic Marching Cubes preview (comparison only)");
        ImGui::Text("Fixed bounds: [%.5g, %.5g, %.5g] - [%.5g, %.5g, %.5g]",
                    bounds_.min[0], bounds_.min[1], bounds_.min[2],
                    bounds_.max[0], bounds_.max[1], bounds_.max[2]);
        ImGui::Text("Sample vertices: %zu x %zu x %zu (%zu values)", nx_, ny_, nz_, samples_.size());
        ImGui::Checkbox("Show local voxel grid", &show_voxel_grid_);
        ImGui::Checkbox("Highlight selected corner", &show_selected_corner_);

        ImGui::Spacing();
        ImGui::TextUnformatted("Change preview grid (uniform on X/Y/Z)");
        ImGui::TextWrapped("N means cells per axis; N+1 scalar samples are evaluated on each axis. The patch bounds stay fixed.");
        ImGui::SetNextItemWidth(120.0f);
        ImGui::InputInt("N cells per axis", &pending_grid_cells_);
        pending_grid_cells_ = std::clamp(pending_grid_cells_, 2, kMaxGridCells);
        if (ImGui::Button("Apply classic MC preview")) apply_classic_preview();
        ImGui::SameLine();
        if (ImGui::Button("Restore LocalMarchingCubes patch")) restore_original_patch();
        ImGui::TextWrapped("The classic MC preview is a separate comparison mode. It does not replace or claim to rerun LocalMarchingCubes.");
        if (last_action_ms_ >= 0.0)
            ImGui::Text("Last action: %.2f ms", last_action_ms_);
        if (!action_message_.empty()) ImGui::TextWrapped("%s", action_message_.c_str());
        if (!error_message_.empty()) ImGui::TextColored(ImVec4(1.0f, 0.35f, 0.3f, 1.0f), "%s", error_message_.c_str());

        ImGui::Spacing();
        ImGui::Separator();
        ImGui::TextUnformatted("Corner scalar values");
        if (samples_.size() <= kMaxAllLabelSamples)
            ImGui::Checkbox("Show labels on all grid corners", &show_all_labels_);
        else {
            show_all_labels_ = false;
            ImGui::TextWrapped("For readability, all-corner labels are available only up to %zu samples. Select a table row to highlight one corner.", kMaxAllLabelSamples);
        }

        if (selected_sample_index_ < samples_.size()) {
            const auto& s = samples_[selected_sample_index_];
            ImGui::Text("Selected corner (%zu, %zu, %zu): f = %.12g", s.i, s.j, s.k, s.value);
            ImGui::Text("Position: (%.9g, %.9g, %.9g)", s.x, s.y, s.z);
        }

        if (ImGui::BeginChild("ScalarSamples", ImVec2(0, 265), true)) {
            if (ImGui::BeginTable("ScalarSampleTable", 5,
                                  ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg |
                                  ImGuiTableFlags_ScrollY | ImGuiTableFlags_Resizable |
                                  ImGuiTableFlags_SizingFixedFit)) {
                ImGui::TableSetupScrollFreeze(0, 1);
                ImGui::TableSetupColumn("Grid corner (i,j,k)");
                ImGui::TableSetupColumn("x");
                ImGui::TableSetupColumn("y");
                ImGui::TableSetupColumn("z");
                ImGui::TableSetupColumn("f(x,y,z)");
                ImGui::TableHeadersRow();
                ImGuiListClipper clipper;
                clipper.Begin(static_cast<int>(samples_.size()));
                while (clipper.Step()) {
                    for (int row = clipper.DisplayStart; row < clipper.DisplayEnd; ++row) {
                        const auto index = static_cast<std::size_t>(row);
                        const auto& s = samples_[index];
                        ImGui::TableNextRow();
                        ImGui::TableSetColumnIndex(0);
                        char label[96];
                        std::snprintf(label, sizeof(label), "(%zu,%zu,%zu)##corner_%zu", s.i, s.j, s.k, index);
                        const bool selected = index == selected_sample_index_;
                        if (ImGui::Selectable(label, selected, ImGuiSelectableFlags_SpanAllColumns))
                            selected_sample_index_ = index;
                        ImGui::TableSetColumnIndex(1); ImGui::Text("%.7g", s.x);
                        ImGui::TableSetColumnIndex(2); ImGui::Text("%.7g", s.y);
                        ImGui::TableSetColumnIndex(3); ImGui::Text("%.7g", s.z);
                        ImGui::TableSetColumnIndex(4);
                        if (std::abs(s.value) <= 1e-12)
                            ImGui::TextColored(ImVec4(0.2f, 0.9f, 0.3f, 1.0f), "%.7g", s.value);
                        else ImGui::Text("%.7g", s.value);
                    }
                }
                ImGui::EndTable();
            }
        }
        ImGui::EndChild();

        // Draw annotations after the table UI is built. They are anchored to the
        // projected 3D sample positions, not to table rows or arbitrary screen spots.
        draw_corner_annotations();
    }

    void draw(const std::string& draw_mode) override
    {
        pmp::MeshViewer::draw(draw_mode);
        if (show_voxel_grid_ && !grid_indices_.empty()) draw_voxel_grid();
        if (show_selected_corner_ && selected_sample_index_ < samples_.size())
            draw_selected_corner_marker();
    }

private:
    enum class DisplayMode { LocalVoxelRefinement, ClassicMCPreview };
    struct LineVertex { float x, y, z; };
    static constexpr int kMaxGridCells = 64;
    static constexpr std::size_t kMaxAllLabelSamples = 729; // up to 9^3 corners

    std::size_t nearest_sample_to_singularity(const std::vector<CornerSample>& samples) const
    {
        std::size_t best_index = 0;
        double best_distance = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < samples.size(); ++i) {
            const double dx = samples[i].x - singularity_center_[0];
            const double dy = samples[i].y - singularity_center_[1];
            const double dz = samples[i].z - singularity_center_[2];
            const double d2 = dx*dx + dy*dy + dz*dz;
            if (d2 < best_distance) { best_distance = d2; best_index = i; }
        }
        return best_index;
    }

    void install_raw_mesh(const iso::mc::Mesh& raw_mesh)
    {
        pmp::SurfaceMesh converted = to_pmp_mesh(raw_mesh);
        mesh_ = converted;
        update_mesh();
    }

    void apply_classic_preview()
    {
        const auto start = std::chrono::steady_clock::now();
        error_message_.clear();
        action_message_.clear();
        const std::size_t cells = static_cast<std::size_t>(std::clamp(pending_grid_cells_, 2, kMaxGridCells));
        const std::size_t samples_per_axis = cells + 1;
        std::cout << "[PatchViewer] Re-extracting classic MC preview: " << cells << " cells/axis, "
                  << samples_per_axis << " samples/axis, fixed bounds.\n" << std::flush;
        try {
            iso::mc::Options options;
            options.nx = samples_per_axis;
            options.ny = samples_per_axis;
            options.nz = samples_per_axis;
            options.isovalue = 0.0;
            const auto extracted = iso::mc::extract(
                [&](double x, double y, double z) { return surface_->eval(x, y, z); },
                bounds_, options);
            if (extracted.triangles.empty())
                throw std::runtime_error("Classic MC preview has no triangles at this resolution. Try a larger N; the original LocalMarchingCubes patch is preserved.");
            auto new_samples = sample_grid(*surface_, bounds_, samples_per_axis, samples_per_axis, samples_per_axis);
            pmp::SurfaceMesh converted = to_pmp_mesh(extracted);

            const auto mesh_path = output_dir_ / ("local_patch_preview_classic_mc_" + surface_name_ + "_N" + std::to_string(cells) + ".obj");
            const auto csv_path = output_dir_ / ("local_patch_preview_samples_" + surface_name_ + "_N" + std::to_string(cells) + ".csv");
            pmp::write(converted, mesh_path.string());
            write_samples_csv(csv_path, new_samples);

            active_mesh_ = extracted;
            samples_ = std::move(new_samples);
            nx_ = ny_ = nz_ = samples_per_axis;
            mesh_ = converted;
            update_mesh();
            mode_ = DisplayMode::ClassicMCPreview;
            selected_sample_index_ = nearest_sample_to_singularity(samples_);
            build_voxel_grid();
            action_message_ = "Classic MC preview applied. The original LocalMarchingCubes mesh remains available with Restore. Saved: " + mesh_path.filename().string() + " and " + csv_path.filename().string();
            std::cout << "[PatchViewer] Preview ready: " << extracted.vertices.size() << " vertices, "
                      << extracted.triangles.size() << " triangles, " << samples_.size()
                      << " scalar samples.\n[PatchViewer] OBJ: " << mesh_path
                      << "\n[PatchViewer] CSV: " << csv_path << "\n" << std::flush;
        } catch (const std::exception& e) {
            error_message_ = e.what();
            std::cerr << "[PatchViewer] Preview failed: " << error_message_ << '\n' << std::flush;
        }
        const auto end = std::chrono::steady_clock::now();
        last_action_ms_ = std::chrono::duration<double, std::milli>(end - start).count();
    }

    void restore_original_patch()
    {
        const auto start = std::chrono::steady_clock::now();
        error_message_.clear();
        action_message_.clear();
        active_mesh_ = original_patch_mesh_;
        samples_ = original_samples_;
        nx_ = original_nx_; ny_ = original_ny_; nz_ = original_nz_;
        install_raw_mesh(original_patch_mesh_);
        mode_ = DisplayMode::LocalVoxelRefinement;
        selected_sample_index_ = nearest_sample_to_singularity(samples_);
        build_voxel_grid();
        action_message_ = "Restored the original LocalMarchingCubes research patch. Its mesh and grid resolution are unchanged.";
        std::cout << "[PatchViewer] Restored original LocalMarchingCubes patch ("
                  << active_mesh_.vertices.size() << " vertices, " << active_mesh_.triangles.size()
                  << " triangles).\n" << std::flush;
        const auto end = std::chrono::steady_clock::now();
        last_action_ms_ = std::chrono::duration<double, std::milli>(end - start).count();
    }

    void build_voxel_grid()
    {
        grid_vertices_.clear();
        grid_indices_.clear();
        if (nx_ < 2 || ny_ < 2 || nz_ < 2) { grid_gpu_data_dirty_ = true; return; }
        const double dx = (bounds_.max[0] - bounds_.min[0]) / static_cast<double>(nx_ - 1);
        const double dy = (bounds_.max[1] - bounds_.min[1]) / static_cast<double>(ny_ - 1);
        const double dz = (bounds_.max[2] - bounds_.min[2]) / static_cast<double>(nz_ - 1);
        auto add_vertex = [&](double x, double y, double z) -> unsigned int {
            grid_vertices_.push_back({static_cast<float>(x), static_cast<float>(y), static_cast<float>(z)});
            return static_cast<unsigned int>(grid_vertices_.size() - 1);
        };
        auto add_line = [&](double ax, double ay, double az, double bx, double by, double bz) {
            const unsigned int ia = add_vertex(ax, ay, az), ib = add_vertex(bx, by, bz);
            grid_indices_.push_back(ia); grid_indices_.push_back(ib);
        };
        const double xmin = bounds_.min[0], xmax = bounds_.max[0];
        const double ymin = bounds_.min[1], ymax = bounds_.max[1];
        const double zmin = bounds_.min[2], zmax = bounds_.max[2];
        for (std::size_t j = 0; j < ny_; ++j) {
            const double y = ymin + static_cast<double>(j) * dy;
            for (std::size_t k = 0; k < nz_; ++k) {
                const double z = zmin + static_cast<double>(k) * dz;
                add_line(xmin, y, z, xmax, y, z);
            }
        }
        for (std::size_t i = 0; i < nx_; ++i) {
            const double x = xmin + static_cast<double>(i) * dx;
            for (std::size_t k = 0; k < nz_; ++k) {
                const double z = zmin + static_cast<double>(k) * dz;
                add_line(x, ymin, z, x, ymax, z);
            }
        }
        for (std::size_t i = 0; i < nx_; ++i) {
            const double x = xmin + static_cast<double>(i) * dx;
            for (std::size_t j = 0; j < ny_; ++j) {
                const double y = ymin + static_cast<double>(j) * dy;
                add_line(x, y, zmin, x, y, zmax);
            }
        }
        grid_gpu_data_dirty_ = true;
    }

    void initialize_voxel_grid_rendering()
    {
        if (!grid_shader_initialized_) {
            const char* vertex_shader_source = R"GLSL(
                #version 330 core
                layout(location = 0) in vec3 vertex;
                uniform mat4 projection_matrix;
                uniform mat4 modelview_matrix;
                void main() { gl_Position = projection_matrix * modelview_matrix * vec4(vertex, 1.0); }
            )GLSL";
            const char* fragment_shader_source = R"GLSL(
                #version 330 core
                out vec4 frag_color;
                void main() { frag_color = vec4(0.2, 0.2, 0.2, 0.42); }
            )GLSL";
            voxel_shader_.source(vertex_shader_source, fragment_shader_source);
            grid_shader_initialized_ = true;
        }
        if (!grid_vao_) glGenVertexArrays(1, &grid_vao_);
        if (!grid_vbo_) glGenBuffers(1, &grid_vbo_);
        if (!grid_ebo_) glGenBuffers(1, &grid_ebo_);
        glBindVertexArray(grid_vao_);
        glBindBuffer(GL_ARRAY_BUFFER, grid_vbo_);
        glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(grid_vertices_.size() * sizeof(LineVertex)),
                     grid_vertices_.data(), GL_DYNAMIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(LineVertex), nullptr);
        glEnableVertexAttribArray(0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, grid_ebo_);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, static_cast<GLsizeiptr>(grid_indices_.size() * sizeof(unsigned int)),
                     grid_indices_.data(), GL_DYNAMIC_DRAW);
        glBindVertexArray(0);
        grid_gpu_data_dirty_ = false;
    }

    void draw_voxel_grid()
    {
        if (!grid_shader_initialized_ || !grid_vao_ || grid_gpu_data_dirty_)
            initialize_voxel_grid_rendering();
        if (grid_indices_.empty()) return;
        voxel_shader_.use();
        voxel_shader_.set_uniform("projection_matrix", projection_matrix_);
        voxel_shader_.set_uniform("modelview_matrix", modelview_matrix_);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glDepthMask(GL_FALSE);
        glBindVertexArray(grid_vao_);
        glDrawElements(GL_LINES, static_cast<GLsizei>(grid_indices_.size()), GL_UNSIGNED_INT, nullptr);
        glBindVertexArray(0);
        glDepthMask(GL_TRUE);
        glDisable(GL_BLEND);
        voxel_shader_.disable();
    }

    void initialize_marker_rendering()
    {
        if (marker_shader_initialized_) return;
        const char* vertex_shader_source = R"GLSL(
            #version 330 core
            layout(location = 0) in vec3 vertex;
            uniform mat4 projection_matrix;
            uniform mat4 modelview_matrix;
            uniform float point_size;
            void main() {
                gl_Position = projection_matrix * modelview_matrix * vec4(vertex, 1.0);
                gl_PointSize = point_size;
            }
        )GLSL";
        const char* fragment_shader_source = R"GLSL(
            #version 330 core
            out vec4 frag_color;
            void main() {
                vec2 d = gl_PointCoord - vec2(0.5);
                if (dot(d, d) > 0.25) discard;
                frag_color = vec4(1.0, 0.12, 0.02, 1.0);
            }
        )GLSL";
        marker_shader_.source(vertex_shader_source, fragment_shader_source);
        glGenVertexArrays(1, &marker_vao_);
        glGenBuffers(1, &marker_vbo_);
        glBindVertexArray(marker_vao_);
        glBindBuffer(GL_ARRAY_BUFFER, marker_vbo_);
        glBufferData(GL_ARRAY_BUFFER, sizeof(LineVertex), nullptr, GL_DYNAMIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(LineVertex), nullptr);
        glEnableVertexAttribArray(0);
        glBindVertexArray(0);
        marker_shader_initialized_ = true;
    }

    void draw_selected_corner_marker()
    {
        if (selected_sample_index_ >= samples_.size()) return;
        initialize_marker_rendering();
        const auto& s = samples_[selected_sample_index_];
        const LineVertex p{static_cast<float>(s.x), static_cast<float>(s.y), static_cast<float>(s.z)};
        glBindBuffer(GL_ARRAY_BUFFER, marker_vbo_);
        glBufferData(GL_ARRAY_BUFFER, sizeof(LineVertex), &p, GL_DYNAMIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);

        const GLboolean depth_was_enabled = glIsEnabled(GL_DEPTH_TEST);
        const GLboolean point_size_was_enabled = glIsEnabled(GL_PROGRAM_POINT_SIZE);
        glDisable(GL_DEPTH_TEST); // ensure a selected corner remains identifiable through the mesh
        glEnable(GL_PROGRAM_POINT_SIZE);
        marker_shader_.use();
        marker_shader_.set_uniform("projection_matrix", projection_matrix_);
        marker_shader_.set_uniform("modelview_matrix", modelview_matrix_);
        marker_shader_.set_uniform("point_size", 16.0f);
        glBindVertexArray(marker_vao_);
        glDrawArrays(GL_POINTS, 0, 1);
        glBindVertexArray(0);
        marker_shader_.disable();
        if (depth_was_enabled) glEnable(GL_DEPTH_TEST); else glDisable(GL_DEPTH_TEST);
        if (!point_size_was_enabled) glDisable(GL_PROGRAM_POINT_SIZE);
    }

    bool project_to_screen(const CornerSample& s, ImVec2& screen) const
    {
        const pmp::vec4 world(static_cast<float>(s.x), static_cast<float>(s.y), static_cast<float>(s.z), 1.0f);
        const pmp::vec4 clip = projection_matrix_ * modelview_matrix_ * world;
        if (clip[3] <= 1e-6f) return false;
        const float nx = clip[0] / clip[3];
        const float ny = clip[1] / clip[3];
        const float nz = clip[2] / clip[3];
        if (nx < -1.05f || nx > 1.05f || ny < -1.05f || ny > 1.05f || nz < -1.0f || nz > 1.0f)
            return false;
        const ImVec2 display = ImGui::GetIO().DisplaySize;
        screen = ImVec2((nx * 0.5f + 0.5f) * display.x,
                        (1.0f - (ny * 0.5f + 0.5f)) * display.y);
        return true;
    }

    static ImU32 scalar_color(double value)
    {
        if (std::abs(value) <= 1e-12) return IM_COL32(80, 255, 120, 255);
        return value < 0.0 ? IM_COL32(100, 180, 255, 255) : IM_COL32(255, 190, 70, 255);
    }

    void draw_corner_annotations()
    {
        ImDrawList* draw_list = ImGui::GetForegroundDrawList();
        if (show_all_labels_ && samples_.size() <= kMaxAllLabelSamples) {
            for (std::size_t i = 0; i < samples_.size(); ++i) {
                ImVec2 pos;
                if (!project_to_screen(samples_[i], pos)) continue;
                const bool selected = i == selected_sample_index_;
                draw_list->AddCircleFilled(pos, selected ? 4.5f : 2.5f,
                                           selected ? IM_COL32(255, 70, 40, 255) : scalar_color(samples_[i].value));
                char text[48];
                std::snprintf(text, sizeof(text), "%.3g", samples_[i].value);
                draw_list->AddText(ImVec2(pos.x + 4.0f, pos.y - 7.0f), scalar_color(samples_[i].value), text);
            }
        }
        if (show_selected_corner_ && selected_sample_index_ < samples_.size()) {
            const auto& s = samples_[selected_sample_index_];
            ImVec2 pos;
            if (project_to_screen(s, pos)) {
                const ImVec2 label_pos(pos.x + 14.0f, pos.y - 18.0f);
                char text[160];
                std::snprintf(text, sizeof(text), "(%zu,%zu,%zu)  f=%.7g", s.i, s.j, s.k, s.value);
                const ImVec2 text_size = ImGui::CalcTextSize(text);
                draw_list->AddLine(pos, ImVec2(label_pos.x, label_pos.y + text_size.y * 0.5f), IM_COL32(255, 75, 35, 255), 2.0f);
                draw_list->AddCircle(pos, 8.0f, IM_COL32(255, 75, 35, 255), 0, 2.0f);
                draw_list->AddRectFilled(ImVec2(label_pos.x - 3.0f, label_pos.y - 2.0f),
                                         ImVec2(label_pos.x + text_size.x + 5.0f, label_pos.y + text_size.y + 2.0f),
                                         IM_COL32(20, 20, 20, 225), 3.0f);
                draw_list->AddText(label_pos, IM_COL32(255, 255, 255, 255), text);
            }
        }
    }

    const ImplicitSurface* surface_;
    iso::mc::Mesh original_patch_mesh_;
    iso::mc::Mesh active_mesh_;
    std::vector<CornerSample> samples_;
    std::vector<CornerSample> original_samples_;
    iso::mc::Bounds bounds_;
    iso::mc::Point singularity_center_;
    std::size_t original_nx_, original_ny_, original_nz_;
    std::size_t nx_, ny_, nz_;
    std::filesystem::path output_dir_;
    std::string surface_name_;
    DisplayMode mode_{DisplayMode::LocalVoxelRefinement};
    std::size_t selected_sample_index_{0};
    int pending_grid_cells_{8};
    bool show_voxel_grid_{true};
    bool show_selected_corner_{true};
    bool show_all_labels_{false};

    std::vector<LineVertex> grid_vertices_;
    std::vector<unsigned int> grid_indices_;
    GLuint grid_vao_{0}, grid_vbo_{0}, grid_ebo_{0};
    GLuint marker_vao_{0}, marker_vbo_{0};
    pmp::Shader voxel_shader_;
    pmp::Shader marker_shader_;
    bool grid_shader_initialized_{false};
    bool grid_gpu_data_dirty_{true};
    bool marker_shader_initialized_{false};
    double last_action_ms_{-1.0};
    std::string action_message_;
    std::string error_message_;
};
} // namespace

int main(int argc, char** argv)
{
    try {
        if (flag(argc, argv, "--help")) {
            std::cout << "Local patch-only interactive viewer\n"
                      << "Usage: IsosurfacesPatchViewer [c] [resolution] [extent] [radius] [level] [--surface NAME] [--output DIR] [--no-gui]\n"
                      << "In the GUI, select a table row to highlight that exact 3D grid corner and its scalar value.\n"
                      << "Apply classic MC preview with N cells per axis; patch bounds stay fixed. Restore returns to the original LocalMarchingCubes patch.\n";
            return 0;
        }
        const double c = argc > 1 ? std::stod(argv[1]) : 0.0;
        const std::size_t resolution = argc > 2 ? std::stoul(argv[2]) : 65;
        const double extent = argc > 3 ? std::stod(argv[3]) : 2.0;
        const double radius = argc > 4 ? std::stod(argv[4]) : 0.30;
        const std::size_t level = argc > 5 ? std::stoul(argv[5]) : 2;
        const std::string surface_name = value(argc, argv, "--surface", "cone");
        const std::filesystem::path output_dir = value(argc, argv, "--output", "outputs/local_patch");
        if (resolution < 2 || !(extent > 0.0) || !(radius > 0.0) || level == 0)
            throw std::invalid_argument("resolution >= 2, extent > 0, radius > 0, and level >= 1 are required");
        std::filesystem::create_directories(output_dir);
        std::cout << "[1/5] Creating surface: " << surface_name << "\n" << std::flush;
        const auto surface = make_surface(surface_name, c);
        const iso::mc::Bounds bounds{{-extent,-extent,-extent},{extent,extent,extent}};
        iso::mc::Options global_options; global_options.nx=resolution; global_options.ny=resolution; global_options.nz=resolution;
        std::cout << "[2/5] Extracting coarse mesh to locate singularity...\n" << std::flush;
        const auto global = iso::mc::extract([&](double x,double y,double z){return surface->eval(x,y,z);}, bounds, global_options);
        const SingularityDetector detector(*surface, 1e-2, 1e-3);
        iso::mc::Point center{0,0,0}; bool found=false; double best=std::numeric_limits<double>::infinity();
        for (const auto& p : global.vertices) {
            const pmp::Point q{float(p[0]), float(p[1]), float(p[2])};
            const auto r = detector.classifyPoint(q);
            if (r.type == SingularityType::Regular || std::abs(r.function_value) > 1e-5) continue;
            if (!found || r.gradient_norm < best) { found=true; best=r.gradient_norm; center=p; }
        }
        if (!found) throw std::runtime_error("No on-isosurface singularity found. Try c=0 for the singular examples.");
        std::cout << "    singularity at (" << center[0] << ',' << center[1] << ',' << center[2] << ")\n" << std::flush;
        LocalMarchingCubes extractor; LocalMarchingCubes::Options opts;
        opts.radius=radius; opts.level=level; opts.padding_cells=1; opts.sampling_bounds=bounds; opts.base_resolution=resolution; opts.isovalue=0.0;
        std::cout << "[3/5] Extracting local voxel patch...\n" << std::flush;
        const auto patch = extractor.extract(*surface, center, opts);
        if (patch.mesh.triangles.empty()) throw std::runtime_error("Local patch has no triangles; increase radius or resolution.");
        std::cout << "    patch grid=" << patch.resolution_x << 'x' << patch.resolution_y << 'x' << patch.resolution_z
                  << ", vertices=" << patch.mesh.vertices.size() << ", triangles=" << patch.mesh.triangles.size() << "\n" << std::flush;
        std::cout << "[4/5] Sampling original local grid and writing CSV/mesh...\n" << std::flush;
        const auto samples = sample_grid(*surface, patch.local_bounds, patch.resolution_x, patch.resolution_y, patch.resolution_z);
        const auto csv_path = output_dir / ("local_patch_samples_" + surface_name + ".csv");
        write_samples_csv(csv_path, samples);
        const auto mesh_path = output_dir / ("local_patch_" + surface_name + ".obj");
        auto pmp_mesh = to_pmp_mesh(patch.mesh);
        pmp::write(pmp_mesh, mesh_path.string());
        std::cout << "    samples=" << samples.size() << "\n    CSV=" << csv_path << "\n    OBJ=" << mesh_path << "\n" << std::flush;
        if (flag(argc, argv, "--no-gui")) return 0;
        std::cout << "[5/5] Opening local patch viewer. Select a row to highlight its grid corner.\n" << std::flush;
        PatchViewer viewer("Local Patch + Scalar Values at Grid Corners", 1440, 960,
                           *surface, patch.mesh, samples, patch.local_bounds,
                           patch.resolution_x, patch.resolution_y, patch.resolution_z,
                           center, output_dir, surface_name);
        viewer.load_mesh(mesh_path.string().c_str());
        return viewer.run();
    } catch (const std::exception& e) {
        std::cerr << "Local patch viewer error: " << e.what() << '\n';
        return 1;
    }
}
