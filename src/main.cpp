#include "ImplicitSurface.hpp"
#include "SingularityDetector.hpp"
#include "MarchingCubes.hpp"
#include "LocalMarchingCubes.hpp"
#include "LocalUnfolder.hpp"
#include "app/MyViewer.h"

#include <pmp/io/io.h>
#include <pmp/surface_mesh.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace
{

const char* singularity_type_name(SingularityType type)
{
    switch (type)
    {
    case SingularityType::Regular:
        return "Regular";
    case SingularityType::NonDegenerateSingular:
        return "NonDegenerateSingular";
    case SingularityType::DegenerateSingular:
        return "DegenerateSingular";
    }
    return "Unknown";
}

void write_colored_ply(
    const pmp::SurfaceMesh& mesh,
    const std::vector<pmp::Color>& vertex_colors,
    const std::string& filename)
{
    if (vertex_colors.size() != mesh.n_vertices())
        throw std::runtime_error("Vertex/color count mismatch while writing PLY");

    const auto points = mesh.get_vertex_property<pmp::Point>("v:point");
    std::ofstream out(filename);
    if (!out)
        throw std::runtime_error("Could not open PLY output file: " + filename);

    out << "ply\n"
        << "format ascii 1.0\n"
        << "element vertex " << mesh.n_vertices() << "\n"
        << "property float x\n"
        << "property float y\n"
        << "property float z\n"
        << "property uchar red\n"
        << "property uchar green\n"
        << "property uchar blue\n"
        << "element face " << mesh.n_faces() << "\n"
        << "property list uchar int vertex_indices\n"
        << "end_header\n";

    for (auto v : mesh.vertices())
    {
        const pmp::Point p = points[v];
        const pmp::Color c = vertex_colors[v.idx()];
        const int r = static_cast<int>(std::clamp(c[0], 0.0f, 1.0f) * 255.0f);
        const int g = static_cast<int>(std::clamp(c[1], 0.0f, 1.0f) * 255.0f);
        const int b = static_cast<int>(std::clamp(c[2], 0.0f, 1.0f) * 255.0f);
        out << p[0] << ' ' << p[1] << ' ' << p[2] << ' '
            << r << ' ' << g << ' ' << b << '\n';
    }

    for (auto f : mesh.faces())
    {
        std::vector<int> indices;
        for (auto v : mesh.vertices(f))
            indices.push_back(v.idx());

        out << indices.size();
        for (int index : indices)
            out << ' ' << index;
        out << '\n';
    }
}

pmp::SurfaceMesh to_pmp_mesh(const iso::mc::Mesh& source)
{
    pmp::SurfaceMesh mesh;
    std::vector<pmp::Vertex> vertices;
    vertices.reserve(source.vertices.size());

    for (const auto& p : source.vertices)
    {
        vertices.push_back(mesh.add_vertex(
            pmp::Point(
                static_cast<float>(p[0]),
                static_cast<float>(p[1]),
                static_cast<float>(p[2]))));
    }

    for (const auto& triangle : source.triangles)
    {
        mesh.add_triangle(
            vertices[triangle[0]],
            vertices[triangle[1]],
            vertices[triangle[2]]);
    }

    return mesh;
}

using EdgeKey = std::pair<iso::mc::Index, iso::mc::Index>;

EdgeKey make_edge(iso::mc::Index a, iso::mc::Index b)
{
    if (a > b)
        std::swap(a, b);
    return {a, b};
}

struct EdgeStatistics
{
    std::size_t boundary_edges = 0;
    std::size_t boundary_vertices = 0;
    std::size_t nonmanifold_edges = 0;
};

EdgeStatistics compute_edge_statistics(const iso::mc::Mesh& mesh)
{
    std::map<EdgeKey, std::size_t> uses;
    for (const auto& t : mesh.triangles)
    {
        ++uses[make_edge(t[0], t[1])];
        ++uses[make_edge(t[1], t[2])];
        ++uses[make_edge(t[2], t[0])];
    }

    EdgeStatistics result;
    std::set<iso::mc::Index> boundary_vertices;

    for (const auto& [edge, count] : uses)
    {
        if (count == 1)
        {
            ++result.boundary_edges;
            boundary_vertices.insert(edge.first);
            boundary_vertices.insert(edge.second);
        }
        else if (count > 2)
        {
            ++result.nonmanifold_edges;
        }
    }

    result.boundary_vertices = boundary_vertices.size();
    return result;
}

std::vector<int> compute_connected_components(
    const pmp::SurfaceMesh& mesh,
    std::vector<std::size_t>& component_sizes)
{
    std::vector<int> component(mesh.n_vertices(), -1);
    int component_count = 0;

    for (auto start : mesh.vertices())
    {
        if (component[start.idx()] != -1)
            continue;

        const int id = component_count++;
        std::queue<pmp::Vertex> queue;
        queue.push(start);
        component[start.idx()] = id;
        std::size_t size = 0;

        while (!queue.empty())
        {
            const pmp::Vertex v = queue.front();
            queue.pop();
            ++size;

            for (auto vv : mesh.vertices(v))
            {
                if (component[vv.idx()] == -1)
                {
                    component[vv.idx()] = id;
                    queue.push(vv);
                }
            }
        }

        component_sizes.push_back(size);
    }

    return component;
}

struct CoincidentStatistics
{
    std::size_t groups = 0;
    std::size_t extra_vertices = 0;
};

CoincidentStatistics compute_coincident_statistics(const pmp::SurfaceMesh& mesh)
{
    using PositionKey = std::tuple<float, float, float>;
    std::map<PositionKey, std::vector<int>> groups;

    for (auto v : mesh.vertices())
    {
        const pmp::Point p = mesh.position(v);
        groups[{p[0], p[1], p[2]}].push_back(v.idx());
    }

    CoincidentStatistics result;
    for (const auto& [key, ids] : groups)
    {
        (void)key;
        if (ids.size() > 1)
        {
            ++result.groups;
            result.extra_vertices += ids.size() - 1;
        }
    }
    return result;
}

struct ClassificationCounts
{
    std::size_t regular = 0;
    std::size_t nondegenerate = 0;
    std::size_t degenerate = 0;
};

struct GeometricSingularity
{
    pmp::Point center{0.0f, 0.0f, 0.0f};
    SingularityType type = SingularityType::NonDegenerateSingular;
    std::vector<int> representative_vertices;
    double hessian_determinant = 0.0;
    Eigen::Vector3d eigenvalues = Eigen::Vector3d::Zero();
};

ClassificationCounts classify_mesh_counts(
    const pmp::SurfaceMesh& mesh,
    const SingularityDetector& detector)
{
    ClassificationCounts result;
    for (auto v : mesh.vertices())
    {
        const auto classification = detector.classifyPoint(mesh.position(v));
        switch (classification.type)
        {
        case SingularityType::Regular:
            ++result.regular;
            break;
        case SingularityType::NonDegenerateSingular:
            ++result.nondegenerate;
            break;
        case SingularityType::DegenerateSingular:
            ++result.degenerate;
            break;
        }
    }
    return result;
}

std::vector<GeometricSingularity> group_geometric_singularities(
    const pmp::SurfaceMesh& mesh,
    const std::vector<SingularityClassificationResult>& results)
{
    constexpr double tolerance = 1e-8;
    std::vector<int> singular_vertices;

    for (auto v : mesh.vertices())
    {
        if (results[v.idx()].type != SingularityType::Regular)
            singular_vertices.push_back(v.idx());
    }

    std::vector<bool> assigned(singular_vertices.size(), false);
    std::vector<GeometricSingularity> groups;

    for (std::size_t i = 0; i < singular_vertices.size(); ++i)
    {
        if (assigned[i])
            continue;

        const int seed_id = singular_vertices[i];
        const pmp::Point seed = mesh.position(pmp::Vertex(seed_id));

        GeometricSingularity group;
        group.type = results[seed_id].type;
        group.hessian_determinant = results[seed_id].hessian_determinant;
        group.eigenvalues = results[seed_id].eigenvalues;
        group.representative_vertices.push_back(seed_id);
        assigned[i] = true;

        for (std::size_t j = i + 1; j < singular_vertices.size(); ++j)
        {
            if (assigned[j])
                continue;

            const int candidate_id = singular_vertices[j];
            const pmp::Point candidate =
                mesh.position(pmp::Vertex(candidate_id));

            const double dx =
                static_cast<double>(candidate[0] - seed[0]);
            const double dy =
                static_cast<double>(candidate[1] - seed[1]);
            const double dz =
                static_cast<double>(candidate[2] - seed[2]);

            if (dx * dx + dy * dy + dz * dz <= tolerance * tolerance)
            {
                group.representative_vertices.push_back(candidate_id);
                assigned[j] = true;
            }
        }

        pmp::Point center(0.0f, 0.0f, 0.0f);
        for (int id : group.representative_vertices)
            center += mesh.position(pmp::Vertex(id));

        center /= static_cast<float>(
            group.representative_vertices.size());

        group.center = center;
        groups.push_back(std::move(group));
    }

    return groups;
}

std::vector<pmp::Color> make_local_colors(
    const iso::mc::Mesh& mesh,
    const iso::mc::Point& center,
    double radius)
{
    std::vector<pmp::Color> colors(
        mesh.vertices.size(),
        pmp::Color(0.6f, 0.6f, 0.6f));

    const double radius2 = radius * radius;

    for (std::size_t i = 0; i < mesh.vertices.size(); ++i)
    {
        const auto& p = mesh.vertices[i];

        const double dx = p[0] - center[0];
        const double dy = p[1] - center[1];
        const double dz = p[2] - center[2];
        const double d2 = dx * dx + dy * dy + dz * dz;

        if (d2 <= radius2)
            colors[i] = pmp::Color(0.2f, 0.6f, 1.0f);

        if (d2 <= 1e-16)
            colors[i] = pmp::Color(1.0f, 0.85f, 0.0f);
    }

    return colors;
}

std::vector<pmp::Color> make_uniform_colors(
    std::size_t count,
    const pmp::Color& color)
{
    return std::vector<pmp::Color>(count, color);
}

void print_scalar_statistics(
    const std::string& label,
    const LocalUnfolder::ScalarStatistics& scalar)
{
    std::cout << '\n' << label << ":\n"
              << "  Refinement factor : "
              << scalar.refinement_factor << "x\n"
              << "  Voxel spacing     : "
              << scalar.spacing << '\n'
              << "  Selected voxels   : "
              << scalar.selected_voxels << '\n'
              << "  Corner values     : "
              << scalar.corner_value_samples << '\n'
              << "  Scalar range      : ["
              << scalar.min_value << ", "
              << scalar.max_value << "]\n"
              << "  Scalar mean       : "
              << scalar.mean_value << '\n';
}

std::string format_parameter(double value)
{
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(6) << value;
    return stream.str();
}

bool has_no_gui_flag(int argc, char** argv)
{
    for (int i = 1; i < argc; ++i)
    {
        if (std::string(argv[i]) == "--no-gui")
            return true;
    }
    return false;
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        const double c_val =
            argc > 1 ? std::stod(argv[1]) : 0.0;

        const std::size_t resolution =
            argc > 2
                ? static_cast<std::size_t>(std::stoul(argv[2]))
                : 129;

        const double extent =
            argc > 3 ? std::stod(argv[3]) : 2.0;

        const double local_radius =
            argc > 4 ? std::stod(argv[4]) : 0.3;

        const std::size_t local_levels =
            argc > 5
                ? static_cast<std::size_t>(std::stoul(argv[5]))
                : 2;

        const bool no_gui = has_no_gui_flag(argc, argv);

        if (resolution < 2)
            throw std::invalid_argument(
                "resolution must be at least 2");

        if (!(extent > 0.0))
            throw std::invalid_argument(
                "extent must be positive");

        if (!(local_radius > 0.0))
            throw std::invalid_argument(
                "local radius must be positive");

        ParameterizedConeQuadric quadric(c_val);
        SingularityDetector detector(
            quadric,
            1e-2,
            1e-3);

        iso::mc::Bounds bounds{
            {-extent, -extent, -extent},
            {extent, extent, extent}};

        iso::mc::Options mc_options;
        mc_options.nx = resolution;
        mc_options.ny = resolution;
        mc_options.nz = resolution;
        mc_options.isovalue = 0.0;

        std::cout
            << "========================================================\n"
            << "  MARCHING CUBES ISOSURFACE GENERATION\n"
            << "========================================================\n\n"
            << "Parameter c      : " << c_val << '\n'
            << "Grid resolution  : " << resolution << "^3\n"
            << "Sampling bounds  : [-"
            << extent << ", " << extent << "]^3\n"
            << "Local radius     : " << local_radius << '\n'
            << "Local levels     : " << local_levels << '\n'
            << "GUI              : "
            << (no_gui ? "disabled" : "enabled")
            << "\n\n";

        const auto mc_start =
            std::chrono::steady_clock::now();

        const iso::mc::Mesh extracted =
            iso::mc::extract(
                [&](double x, double y, double z)
                {
                    return quadric.eval(x, y, z);
                },
                bounds,
                mc_options);

        const auto mc_end =
            std::chrono::steady_clock::now();

        const double mc_ms =
            std::chrono::duration<double, std::milli>(
                mc_end - mc_start).count();

        std::cout
            << "Marching Cubes produced "
            << extracted.vertices.size()
            << " vertices and "
            << extracted.triangles.size()
            << " triangles.\n"
            << "Marching Cubes time : "
            << mc_ms << " ms\n";

        pmp::SurfaceMesh mesh =
            to_pmp_mesh(extracted);

        std::vector<std::size_t> component_sizes;

        const std::vector<int> component =
            compute_connected_components(
                mesh,
                component_sizes);

        std::vector<SingularityClassificationResult>
            classifications(mesh.n_vertices());

        std::vector<pmp::Color> colors(
            mesh.n_vertices(),
            pmp::Color(0.6f, 0.6f, 0.6f));

        ClassificationCounts counts;

        std::cout
            << "\n========================================================\n"
            << "  SINGULARITY DIAGNOSTICS\n"
            << "========================================================\n";

        for (auto v : mesh.vertices())
        {
            const auto result =
                detector.classifyPoint(mesh.position(v));

            classifications[v.idx()] = result;

            switch (result.type)
            {
            case SingularityType::Regular:
                ++counts.regular;
                break;

            case SingularityType::NonDegenerateSingular:
                ++counts.nondegenerate;
                colors[v.idx()] =
                    pmp::Color(1.0f, 0.85f, 0.0f);
                break;

            case SingularityType::DegenerateSingular:
                ++counts.degenerate;
                colors[v.idx()] =
                    pmp::Color(1.0f, 0.0f, 0.0f);
                break;
            }

            if (result.type ==
                SingularityType::Regular)
            {
                continue;
            }

            std::cout
                << std::setprecision(12)
                << "\nCandidate vertex "
                << v.idx()
                << " position=("
                << result.position[0] << ", "
                << result.position[1] << ", "
                << result.position[2] << ")\n"
                << "  f(p)       = "
                << result.function_value << '\n'
                << "  |grad f|   = "
                << result.gradient_norm << '\n'
                << "  det(H)     = "
                << result.hessian_determinant << '\n'
                << "  eigenvalues= ("
                << result.eigenvalues[0] << ", "
                << result.eigenvalues[1] << ", "
                << result.eigenvalues[2] << ")\n"
                << "  class      = "
                << singularity_type_name(result.type)
                << '\n';
        }

        const auto geometric_singularities =
            group_geometric_singularities(
                mesh,
                classifications);

        const EdgeStatistics base_edges =
            compute_edge_statistics(extracted);

        const CoincidentStatistics coincident =
            compute_coincident_statistics(mesh);

        std::size_t pmp_boundary_edges = 0;
        std::size_t nonmanifold_vertices = 0;

        for (auto e : mesh.edges())
        {
            if (mesh.is_boundary(e))
                ++pmp_boundary_edges;
        }

        for (auto v : mesh.vertices())
        {
            if (!mesh.is_manifold(v))
                ++nonmanifold_vertices;
        }

        std::cout
            << "\n========================================================\n"
            << "  MESH TOPOLOGY + SINGULARITY REPORT\n"
            << "========================================================\n"
            << "Vertices            : "
            << mesh.n_vertices() << '\n'
            << "Edges               : "
            << mesh.n_edges() << '\n'
            << "Triangles           : "
            << mesh.n_faces() << '\n'
            << "Global surface boundary edges   : "
            << base_edges.boundary_edges << '\n'
            << "Global surface boundary vertices: "
            << base_edges.boundary_vertices << '\n'
            << "PMP boundary edges  : "
            << pmp_boundary_edges << '\n'
            << "Non-manifold edges  : "
            << base_edges.nonmanifold_edges << '\n'
            << "Non-manifold vertices: "
            << nonmanifold_vertices << '\n'
            << "Connected components : "
            << component_sizes.size() << '\n';

        for (std::size_t i = 0;
             i < component_sizes.size();
             ++i)
        {
            std::cout
                << "  Component " << i
                << " vertices: "
                << component_sizes[i] << '\n';
        }

        std::cout
            << "Coincident vertex groups: "
            << coincident.groups << '\n'
            << "Coincident extra vertices: "
            << coincident.extra_vertices << '\n'
            << "Regular vertices     : "
            << counts.regular << '\n'
            << "Nondegenerate singular vertices: "
            << counts.nondegenerate << '\n'
            << "Degenerate singular vertices    : "
            << counts.degenerate << '\n'
            << "Geometric singularity groups    : "
            << geometric_singularities.size()
            << '\n';

        for (std::size_t i = 0;
             i < geometric_singularities.size();
             ++i)
        {
            const auto& s =
                geometric_singularities[i];

            std::cout
                << "  Singularity " << i
                << " center=("
                << s.center[0] << ", "
                << s.center[1] << ", "
                << s.center[2] << ")"
                << " type="
                << singularity_type_name(s.type)
                << " mesh_representatives="
                << s.representative_vertices.size()
                << " det(H)="
                << s.hessian_determinant
                << " eigenvalues=("
                << s.eigenvalues[0] << ", "
                << s.eigenvalues[1] << ", "
                << s.eigenvalues[2] << ")\n";
        }

        std::cout
            << "========================================================\n";

        if (!geometric_singularities.empty())
        {
            const auto& singularity =
                geometric_singularities.front();

            const iso::mc::Point center{
                static_cast<double>(
                    singularity.center[0]),
                static_cast<double>(
                    singularity.center[1]),
                static_cast<double>(
                    singularity.center[2])};

            std::cout
                << "\nSingular mesh representatives:\n";

            for (int id :
                 singularity.representative_vertices)
            {
                const pmp::Vertex v(id);

                std::cout
                    << "  Vertex " << id
                    << " -> component "
                    << component[id]
                    << ", valence "
                    << mesh.valence(v)
                    << '\n';
            }

            LocalUnfolder::Options local_options;
            local_options.radius = local_radius;
            local_options.levels = local_levels;
            local_options.region_mode =
                LocalUnfolder::RegionMode::TopologicalBfs;
            local_options.sampling_bounds = bounds;
            local_options.base_resolution = resolution;
            local_options.isovalue =
                mc_options.isovalue;

            LocalUnfolder local_unfolder;

            const auto local_start =
                std::chrono::steady_clock::now();

            const auto local_result =
                local_unfolder.refine(
                    extracted,
                    quadric,
                    center,
                    local_options);

            const auto local_end =
                std::chrono::steady_clock::now();

            const double local_ms =
                std::chrono::duration<double, std::milli>(
                    local_end - local_start).count();

            std::cout
                << "\n========================================================\n"
                << "  LOCAL CONFORMING REFINEMENT BASELINE\n"
                << "========================================================\n"
                << "Input vertices      : "
                << local_result.input_vertices << '\n'
                << "Input triangles     : "
                << local_result.input_triangles << '\n'
                << "Sphere faces        : "
                << local_result.sphere_faces_in_region << '\n'
                << "Topological-BFS faces: "
                << local_result.topological_faces_in_region
                << '\n'
                << "Selected faces       : "
                << local_result.faces_in_region << '\n'
                << "Region vertices      : "
                << local_result.vertices_in_region << '\n'
                << "Requested levels    : "
                << local_levels << '\n'
                << "Executed levels     : "
                << local_result.levels.size() << '\n';

            print_scalar_statistics(
                "Original selected voxel scalar values",
                local_result.initial_scalar);

            for (const auto& level :
                 local_result.levels)
            {
                std::cout
                    << "\nLevel " << level.level
                    << ":\n"
                    << "  Selected faces     : "
                    << level.selected_faces << '\n'
                    << "  Sphere faces       : "
                    << level.selected_faces_sphere << '\n'
                    << "  Topological-BFS    : "
                    << level.selected_faces_topological
                    << '\n'
                    << "  Selection differs  : "
                    << (level.selection_changed
                            ? "YES"
                            : "NO")
                    << '\n'
                    << "  New vertices       : "
                    << level.new_vertices << '\n'
                    << "  New triangles      : "
                    << level.new_triangles << '\n'
                    << "  Output vertices    : "
                    << level.output_vertices << '\n'
                    << "  Output triangles   : "
                    << level.output_triangles << '\n'
                    << "  Projection failures: "
                    << level.projection_failures << '\n'
                    << "  Refinement interface edges    : "
                    << level.interface_edges << '\n'
                    << "  Refinement interface vertices : "
                    << level.interface_vertices << '\n'
                    << "  Global surface boundary edges : "
                    << level.global_boundary_edges << '\n'
                    << "  Global surface boundary verts : "
                    << level.global_boundary_vertices << '\n';

                print_scalar_statistics(
                    "  Local scalar values",
                    level.scalar);
            }

            std::cout
                << "\nLocal conforming refinement time: "
                << local_ms << " ms\n"
                << "Final local vertices: "
                << local_result.output_vertices << '\n'
                << "Final local triangles: "
                << local_result.output_triangles << '\n'
                << "Final refinement interface edges   : "
                << local_result.interface_edges << '\n'
                << "Final refinement interface vertices: "
                << local_result.interface_vertices << '\n'
                << "Final global surface boundary edges   : "
                << local_result.global_boundary_edges
                << '\n'
                << "Final global surface boundary vertices: "
                << local_result.global_boundary_vertices
                << '\n';

            std::cout
                << "\n========================================================\n"
                << "  LOCAL VOXEL MARCHING CUBES\n"
                << "========================================================\n";

            LocalMarchingCubes local_mc;

            std::filesystem::create_directories(
                "outputs");

            const std::string c_string =
                format_parameter(c_val);

            for (std::size_t level = 1;
                 level <= local_levels;
                 ++level)
            {
                LocalMarchingCubes::Options patch_options;
                patch_options.radius = local_radius;
                patch_options.level = level;
                patch_options.padding_cells = 1;
                patch_options.sampling_bounds = bounds;
                patch_options.base_resolution = resolution;
                patch_options.isovalue =
                    mc_options.isovalue;

                const auto patch =
                    local_mc.extract(
                        quadric,
                        center,
                        patch_options);

                std::cout
                    << "\nLevel " << level
                    << ":\n"
                    << "  Local MC resolution : "
                    << patch.resolution_x << 'x'
                    << patch.resolution_y << 'x'
                    << patch.resolution_z << '\n'
                    << "  Refinement factor   : "
                    << patch.refinement_factor
                    << "x\n"
                    << "  Base spacing        : "
                    << patch.base_spacing << '\n'
                    << "  Refined spacing     : "
                    << patch.refined_spacing << '\n'
                    << "  Cells near region   : "
                    << patch.selected_cells << '\n'
                    << "  Local MC vertices   : "
                    << patch.mesh.vertices.size()
                    << '\n'
                    << "  Local MC triangles   : "
                    << patch.mesh.triangles.size()
                    << '\n'
                    << "  Patch boundary edges: "
                    << patch.boundary_edges << '\n'
                    << "  Patch boundary verts: "
                    << patch.boundary_vertices << '\n'
                    << "  Non-manifold edges  : "
                    << patch.nonmanifold_edges << '\n'
                    << "  Local MC time       : "
                    << patch.extraction_time_ms
                    << " ms\n";

                const std::filesystem::path patch_ply =
                    std::filesystem::path("outputs") /
                    ("local_mc_patch_c_" +
                     c_string +
                     "_l" +
                     std::to_string(level) +
                     ".ply");

                const std::filesystem::path patch_obj =
                    std::filesystem::path("outputs") /
                    ("local_mc_patch_c_" +
                     c_string +
                     "_l" +
                     std::to_string(level) +
                     ".obj");

                pmp::SurfaceMesh patch_mesh =
                    to_pmp_mesh(patch.mesh);

                const auto patch_colors =
                    make_uniform_colors(
                        patch.mesh.vertices.size(),
                        pmp::Color(
                            0.3f,
                            0.8f,
                            1.0f));

                write_colored_ply(
                    patch_mesh,
                    patch_colors,
                    patch_ply.string());

                pmp::write(
                    patch_mesh,
                    patch_obj.string());

                std::cout
                    << "  Saved local MC PLY : "
                    << patch_ply.string()
                    << '\n'
                    << "  Saved local MC OBJ : "
                    << patch_obj.string()
                    << '\n';
            }

            pmp::SurfaceMesh local_mesh =
                to_pmp_mesh(local_result.mesh);

            const ClassificationCounts local_counts =
                classify_mesh_counts(
                    local_mesh,
                    detector);

            std::vector<std::size_t> local_components;

            compute_connected_components(
                local_mesh,
                local_components);

            const std::filesystem::path local_ply =
                std::filesystem::path("outputs") /
                ("local_unfolded_c_" +
                 c_string +
                 "_l" +
                 std::to_string(local_levels) +
                 ".ply");

            const std::filesystem::path local_obj =
                std::filesystem::path("outputs") /
                ("local_unfolded_c_" +
                 c_string +
                 "_l" +
                 std::to_string(local_levels) +
                 ".obj");

            const auto local_colors =
                make_local_colors(
                    local_result.mesh,
                    center,
                    local_radius);

            write_colored_ply(
                local_mesh,
                local_colors,
                local_ply.string());

            pmp::write(
                local_mesh,
                local_obj.string());

            std::cout
                << "\nLocal refined topology:\n"
                << "  Vertices            : "
                << local_mesh.n_vertices()
                << '\n'
                << "  Edges               : "
                << local_mesh.n_edges()
                << '\n'
                << "  Triangles           : "
                << local_mesh.n_faces()
                << '\n'
                << "  Connected components : "
                << local_components.size()
                << '\n'
                << "  Regular vertices     : "
                << local_counts.regular
                << '\n'
                << "  Nondegenerate singular vertices: "
                << local_counts.nondegenerate
                << '\n'
                << "  Degenerate singular vertices    : "
                << local_counts.degenerate
                << '\n'
                << "\nSaved local PLY: "
                << local_ply.string()
                << '\n'
                << "Saved local OBJ: "
                << local_obj.string()
                << '\n';

            if (no_gui)
                return 0;

            iso::MyViewer window(
                "Local Conforming Refinement",
                1024,
                768);

            window.load_mesh(
                local_obj.string().c_str());

            window.set_voxel_grid(
                pmp::Point(
                    static_cast<float>(bounds.min[0]),
                    static_cast<float>(bounds.min[1]),
                    static_cast<float>(bounds.min[2])),
                pmp::Point(
                    static_cast<float>(bounds.max[0]),
                    static_cast<float>(bounds.max[1]),
                    static_cast<float>(bounds.max[2])),
                resolution,
                resolution,
                resolution);

            return window.run();
        }

        const std::string c_string =
            format_parameter(c_val);

        const std::filesystem::path output_dir =
            "outputs";

        std::filesystem::create_directories(
            output_dir);

        const std::filesystem::path out_ply =
            output_dir /
            ("colored_result_c_" +
             c_string +
             ".ply");

        const std::filesystem::path out_obj =
            output_dir /
            ("viewer_mesh_c_" +
             c_string +
             ".obj");

        write_colored_ply(
            mesh,
            colors,
            out_ply.string());

        pmp::write(
            mesh,
            out_obj.string());

        std::cout
            << "\nSaved colored PLY: "
            << out_ply.string()
            << '\n'
            << "Saved viewer OBJ: "
            << out_obj.string()
            << '\n';

        if (no_gui)
            return 0;

        iso::MyViewer window(
            "Marching Cubes / Singularity Viewer",
            1024,
            768);

        window.load_mesh(
            out_obj.string().c_str());

        window.set_voxel_grid(
            pmp::Point(
                static_cast<float>(bounds.min[0]),
                static_cast<float>(bounds.min[1]),
                static_cast<float>(bounds.min[2])),
            pmp::Point(
                static_cast<float>(bounds.max[0]),
                static_cast<float>(bounds.max[1]),
                static_cast<float>(bounds.max[2])),
            resolution,
            resolution,
            resolution);

        return window.run();
    }
    catch (const std::exception& e)
    {
        std::cerr
            << "\nERROR: "
            << e.what()
            << '\n';

        return 1;
    }
}