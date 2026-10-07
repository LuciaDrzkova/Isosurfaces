#include "ImplicitSurface.hpp"
#include "LocalUnfolder.hpp"
#include "MarchingCubes.hpp"
#include "SingularityDetector.hpp"
#include "app/MyViewer.h"

#include <pmp/io/io.h>
#include <pmp/surface_mesh.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <utility>
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

const char* region_mode_name(LocalUnfolder::RegionMode mode)
{
    return mode == LocalUnfolder::RegionMode::Sphere
        ? "sphere"
        : "topological-BFS";
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
        indices.reserve(3);
        for (auto v : mesh.vertices(f))
            indices.push_back(v.idx());

        out << indices.size();
        for (int index : indices)
            out << ' ' << index;
        out << '\n';
    }
}

std::string format_parameter(double value)
{
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(6) << value;
    return stream.str();
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

    std::set<iso::mc::Index> boundary_vertices;
    EdgeStatistics result;
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

struct GeometricSingularity
{
    pmp::Point center{0.0f, 0.0f, 0.0f};
    SingularityType type = SingularityType::NonDegenerateSingular;
    std::vector<int> representative_vertices;
    double hessian_determinant = 0.0;
    Eigen::Vector3d eigenvalues = Eigen::Vector3d::Zero();
};

std::vector<GeometricSingularity> group_geometric_singularities(
    const pmp::SurfaceMesh& mesh,
    const std::vector<SingularityClassificationResult>& results,
    double position_tolerance = 1e-8)
{
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
            const pmp::Point candidate = mesh.position(pmp::Vertex(candidate_id));
            const double dx = static_cast<double>(candidate[0] - seed[0]);
            const double dy = static_cast<double>(candidate[1] - seed[1]);
            const double dz = static_cast<double>(candidate[2] - seed[2]);
            const double d = std::sqrt(dx * dx + dy * dy + dz * dz);

            if (d <= position_tolerance)
            {
                group.representative_vertices.push_back(candidate_id);
                assigned[j] = true;
            }
        }

        pmp::Point center(0.0f, 0.0f, 0.0f);
        for (int id : group.representative_vertices)
            center += mesh.position(pmp::Vertex(id));
        center /= static_cast<float>(group.representative_vertices.size());
        group.center = center;
        groups.push_back(std::move(group));
    }

    return groups;
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

void print_topology_report(
    const iso::mc::Mesh& extracted,
    const pmp::SurfaceMesh& mesh,
    const std::vector<std::size_t>& component_sizes,
    const std::vector<GeometricSingularity>& singularities,
    std::size_t regular_count,
    std::size_t nondegenerate_count,
    std::size_t degenerate_count)
{
    const EdgeStatistics edges = compute_edge_statistics(extracted);

    std::size_t pmp_boundary_edges = 0;
    for (auto e : mesh.edges())
    {
        if (mesh.is_boundary(e))
            ++pmp_boundary_edges;
    }

    std::size_t nonmanifold_vertices = 0;
    for (auto v : mesh.vertices())
    {
        if (!mesh.is_manifold(v))
            ++nonmanifold_vertices;
    }

    const CoincidentStatistics coincident = compute_coincident_statistics(mesh);

    std::cout
        << "\n========================================================\n"
        << "  MESH TOPOLOGY + SINGULARITY REPORT\n"
        << "========================================================\n"
        << "Vertices            : " << mesh.n_vertices() << '\n'
        << "Edges               : " << mesh.n_edges() << '\n'
        << "Triangles           : " << mesh.n_faces() << '\n'
        << "Global surface boundary edges   : " << edges.boundary_edges << '\n'
        << "Global surface boundary vertices: " << edges.boundary_vertices << '\n'
        << "PMP boundary edges  : " << pmp_boundary_edges << '\n'
        << "Non-manifold edges  : " << edges.nonmanifold_edges << '\n'
        << "Non-manifold vertices: " << nonmanifold_vertices << '\n'
        << "Connected components : " << component_sizes.size() << '\n';

    for (std::size_t i = 0; i < component_sizes.size(); ++i)
        std::cout << "  Component " << i << " vertices: " << component_sizes[i] << '\n';

    std::cout
        << "Coincident vertex groups: " << coincident.groups << '\n'
        << "Coincident extra vertices: " << coincident.extra_vertices << '\n'
        << "Regular vertices     : " << regular_count << '\n'
        << "Nondegenerate singular vertices: " << nondegenerate_count << '\n'
        << "Degenerate singular vertices    : " << degenerate_count << '\n'
        << "Geometric singularity groups    : " << singularities.size() << '\n';

    for (std::size_t i = 0; i < singularities.size(); ++i)
    {
        const auto& s = singularities[i];
        std::cout
            << "  Singularity " << i
            << " center=(" << s.center[0] << ", " << s.center[1] << ", " << s.center[2] << ")"
            << " type=" << singularity_type_name(s.type)
            << " mesh_representatives=" << s.representative_vertices.size()
            << " det(H)=" << s.hessian_determinant
            << " eigenvalues=(" << s.eigenvalues[0] << ", "
            << s.eigenvalues[1] << ", " << s.eigenvalues[2] << ")\n";
    }

    std::cout << "========================================================\n";
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
        const double c_val = argc > 1 ? std::stod(argv[1]) : 0.0;
        const std::size_t resolution = argc > 2
            ? static_cast<std::size_t>(std::stoul(argv[2]))
            : 129;
        const double extent = argc > 3 ? std::stod(argv[3]) : 2.0;
        const double local_radius = argc > 4 ? std::stod(argv[4]) : 0.3;
        const std::size_t local_levels = argc > 5
            ? static_cast<std::size_t>(std::stoul(argv[5]))
            : 2;
        const bool no_gui = has_no_gui_flag(argc, argv);

        if (resolution < 2)
            throw std::invalid_argument("resolution must be at least 2");
        if (!(extent > 0.0))
            throw std::invalid_argument("extent must be positive");
        if (!(local_radius > 0.0))
            throw std::invalid_argument("local radius must be positive");

        ParameterizedConeQuadric quadric(c_val);
        SingularityDetector detector(quadric, 1e-2, 1e-3);

        iso::mc::Bounds bounds{
            {-extent, -extent, -extent},
            { extent,  extent, extent}
        };

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
            << "Sampling bounds  : [-" << extent << ", " << extent << "]^3\n"
            << "Local radius     : " << local_radius << '\n'
            << "Local levels     : " << local_levels << '\n'
            << "GUI               : " << (no_gui ? "disabled" : "enabled") << "\n\n";

        const auto mc_start = std::chrono::steady_clock::now();
        const iso::mc::Mesh extracted = iso::mc::extract(
            [&](double x, double y, double z) { return quadric.eval(x, y, z); },
            bounds,
            mc_options);
        const auto mc_end = std::chrono::steady_clock::now();
        const double mc_ms = std::chrono::duration<double, std::milli>(
            mc_end - mc_start).count();

        std::cout
            << "Marching Cubes produced " << extracted.vertices.size()
            << " vertices and " << extracted.triangles.size() << " triangles.\n"
            << "Marching Cubes time : " << mc_ms << " ms\n";

        pmp::SurfaceMesh mesh = to_pmp_mesh(extracted);
        std::vector<std::size_t> component_sizes;
        const std::vector<int> component =
            compute_connected_components(mesh, component_sizes);

        std::vector<pmp::Color> colors(
            mesh.n_vertices(),
            pmp::Color(0.6f, 0.6f, 0.6f));
        std::vector<SingularityClassificationResult> results(mesh.n_vertices());

        std::size_t regular_count = 0;
        std::size_t nondegenerate_count = 0;
        std::size_t degenerate_count = 0;
        std::vector<pmp::Vertex> singular_vertices;

        std::cout
            << "\n========================================================\n"
            << "  SINGULARITY DIAGNOSTICS\n"
            << "========================================================\n";

        for (auto v : mesh.vertices())
        {
            const pmp::Point p = mesh.position(v);
            const auto result = detector.classifyPoint(p);
            results[v.idx()] = result;

            if (result.type == SingularityType::Regular)
            {
                ++regular_count;
                continue;
            }

            singular_vertices.push_back(v);
            std::cout << std::setprecision(12)
                      << "\nCandidate vertex " << v.idx()
                      << " position=(" << result.position[0] << ", "
                      << result.position[1] << ", " << result.position[2] << ")\n"
                      << "  f(p)       = " << result.function_value << '\n'
                      << "  |grad f|   = " << result.gradient_norm << '\n'
                      << "  det(H)     = " << result.hessian_determinant << '\n'
                      << "  eigenvalues= (" << result.eigenvalues[0] << ", "
                      << result.eigenvalues[1] << ", " << result.eigenvalues[2] << ")\n"
                      << "  class      = " << singularity_type_name(result.type) << '\n';

            if (result.type == SingularityType::NonDegenerateSingular)
            {
                ++nondegenerate_count;
                colors[v.idx()] = pmp::Color(1.0f, 0.85f, 0.0f);
            }
            else
            {
                ++degenerate_count;
                colors[v.idx()] = pmp::Color(1.0f, 0.0f, 0.0f);
            }
        }

        const auto geometric_singularities = group_geometric_singularities(
            mesh, results);

        print_topology_report(
            extracted,
            mesh,
            component_sizes,
            geometric_singularities,
            regular_count,
            nondegenerate_count,
            degenerate_count);

        // Preserve the useful c=0 apex diagnostic while explicitly keeping the
        // coincident apex copies as separate topology vertices/components.
        if (!singular_vertices.empty())
        {
            std::cout << "\nSingular mesh representatives:\n";
            for (auto v : singular_vertices)
            {
                std::cout << "  Vertex " << v.idx()
                          << " -> component " << component[v.idx()]
                          << ", valence " << mesh.valence(v) << '\n';
            }
        }

        if (!geometric_singularities.empty())
        {
            const auto& singularity = geometric_singularities.front();
            iso::mc::Point center{
                static_cast<double>(singularity.center[0]),
                static_cast<double>(singularity.center[1]),
                static_cast<double>(singularity.center[2])};

            LocalUnfolder::Options local_options;
            local_options.radius = local_radius;
            local_options.levels = local_levels;
            local_options.region_mode = LocalUnfolder::RegionMode::TopologicalBfs;

            LocalUnfolder local_unfolder;
            const auto local_start = std::chrono::steady_clock::now();
            const auto local_result = local_unfolder.refine(
                extracted,
                quadric,
                center,
                local_options);
            const auto local_end = std::chrono::steady_clock::now();
            const double local_ms = std::chrono::duration<double, std::milli>(
                local_end - local_start).count();

            std::cout
                << "\n========================================================\n"
                << "  LOCAL UNFOLDING / REFINEMENT\n"
                << "========================================================\n"
                << "Selection mode       : " << region_mode_name(local_options.region_mode) << '\n'
                << "Input vertices       : " << local_result.input_vertices << '\n'
                << "Input triangles      : " << local_result.input_triangles << '\n'
                << "Sphere faces         : " << local_result.sphere_faces_in_region << '\n'
                << "Topological-BFS faces: " << local_result.topological_faces_in_region << '\n'
                << "Selected faces       : " << local_result.faces_in_region << '\n'
                << "Region vertices      : " << local_result.vertices_in_region << '\n'
                << "Refinement interface edges      : " << local_result.interface_edges << '\n'
                << "Refinement interface vertices   : " << local_result.interface_vertices << '\n'
                << "\n";

            for (const auto& level : local_result.levels)
            {
                std::cout
                    << "Level " << level.level << ":\n"
                    << "  Selected faces      : " << level.selected_faces << '\n'
                    << "  Sphere faces        : " << level.selected_faces_sphere << '\n'
                    << "  Topological-BFS     : " << level.selected_faces_topological << '\n'
                    << "  Selection differs   : " << (level.selection_changed ? "YES" : "NO") << '\n'
                    << "  Split edges         : " << level.split_edges << '\n'
                    << "  Interface edges     : " << level.interface_edges << '\n'
                    << "  Interface vertices  : " << level.interface_vertices << '\n'
                    << "  Global surface boundary edges    : " << level.global_boundary_edges << '\n'
                    << "  Global surface boundary vertices : " << level.global_boundary_vertices << '\n'
                    << "  New vertices        : " << level.new_vertices << '\n'
                    << "  New triangles       : " << level.new_triangles << '\n'
                    << "  Projection failures : " << level.projection_failures << '\n';
            }

            std::cout
                << "\nLocal refinement time: " << local_ms << " ms\n"
                << "Final local vertices : " << local_result.output_vertices << '\n'
                << "Final local triangles: " << local_result.output_triangles << '\n';

            pmp::SurfaceMesh local_mesh = to_pmp_mesh(local_result.mesh);
            std::vector<std::size_t> local_component_sizes;
            compute_connected_components(local_mesh, local_component_sizes);
            const EdgeStatistics local_edges = compute_edge_statistics(local_result.mesh);
            std::cout
                << "Final global surface boundary edges    : " << local_edges.boundary_edges << '\n'
                << "Final global surface boundary vertices : " << local_edges.boundary_vertices << '\n'
                << "Final refinement interface edges       : " << local_result.interface_edges << '\n'
                << "Final refinement interface vertices    : " << local_result.interface_vertices << '\n';

            std::filesystem::create_directories("outputs");
            const std::string c_string = format_parameter(c_val);
            const std::filesystem::path local_ply =
                std::filesystem::path("outputs") /
                ("local_unfolded_c_" + c_string + "_l" +
                 std::to_string(local_levels) + ".ply");
            const std::filesystem::path local_obj =
                std::filesystem::path("outputs") /
                ("local_unfolded_c_" + c_string + "_l" +
                 std::to_string(local_levels) + ".obj");

            const auto local_colors = make_local_colors(
                local_result.mesh, center, local_radius);
            write_colored_ply(local_mesh, local_colors, local_ply.string());
            pmp::write(local_mesh, local_obj.string());

            std::cout
                << "Saved local PLY: " << local_ply.string() << '\n'
                << "Saved local OBJ: " << local_obj.string() << '\n';

            if (no_gui)
                return 0;

            iso::MyViewer window("Local Unfolding / Refinement", 1024, 768);
            window.load_mesh(local_obj.string().c_str());
            return window.run();
        }

        std::filesystem::create_directories("outputs");
        const std::string c_string = format_parameter(c_val);
        const std::filesystem::path out_ply =
            std::filesystem::path("outputs") /
            ("colored_result_c_" + c_string + ".ply");
        const std::filesystem::path out_obj =
            std::filesystem::path("outputs") /
            ("viewer_mesh_c_" + c_string + ".obj");

        write_colored_ply(mesh, colors, out_ply.string());
        pmp::write(mesh, out_obj.string());
        std::cout
            << "\nSaved colored PLY: " << out_ply.string() << '\n'
            << "Saved viewer OBJ: " << out_obj.string() << '\n';

        if (no_gui)
            return 0;

        iso::MyViewer window("Marching Cubes / Singularity Viewer", 1024, 768);
        window.load_mesh(out_obj.string().c_str());
        return window.run();
    }
    catch (const std::exception& e)
    {
        std::cerr << "\nERROR: " << e.what() << '\n';
        return 1;
    }
}
