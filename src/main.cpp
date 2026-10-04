#include "ImplicitSurface.hpp"
#include "SingularityDetector.hpp"
#include "MarchingCubes.hpp"
#include "LocalUnfolder.hpp"
#include "app/MyViewer.h"

#include <pmp/io/io.h>
#include <pmp/surface_mesh.h>

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

const char* singularity_type_name(
    SingularityType type)
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
    auto points =
        mesh.get_vertex_property<pmp::Point>(
            "v:point");

    std::ofstream out(filename);

    if (!out.is_open())
    {
        throw std::runtime_error(
            "Could not open PLY output file: " +
            filename);
    }

    out << "ply\n";
    out << "format ascii 1.0\n";

    out << "element vertex "
        << mesh.n_vertices()
        << "\n";

    out << "property float x\n";
    out << "property float y\n";
    out << "property float z\n";

    out << "property uchar red\n";
    out << "property uchar green\n";
    out << "property uchar blue\n";

    out << "element face "
        << mesh.n_faces()
        << "\n";

    out << "property list uchar int vertex_indices\n";
    out << "end_header\n";

    for (auto v :
         mesh.vertices())
    {
        const pmp::Point p =
            points[v];

        const pmp::Color c =
            vertex_colors[v.idx()];

        const int r =
            static_cast<int>(
                std::clamp(
                    c[0],
                    0.0f,
                    1.0f) *
                255.0f);

        const int g =
            static_cast<int>(
                std::clamp(
                    c[1],
                    0.0f,
                    1.0f) *
                255.0f);

        const int b =
            static_cast<int>(
                std::clamp(
                    c[2],
                    0.0f,
                    1.0f) *
                255.0f);

        out
            << p[0] << ' '
            << p[1] << ' '
            << p[2] << ' '
            << r << ' '
            << g << ' '
            << b
            << '\n';
    }

    for (auto f :
         mesh.faces())
    {
        std::vector<int> indices;
        indices.reserve(3);

        for (auto v :
             mesh.vertices(f))
        {
            indices.push_back(
                v.idx());
        }

        out << indices.size();

        for (int index :
             indices)
        {
            out
                << ' '
                << index;
        }

        out << '\n';
    }
}

std::string format_parameter(
    double value)
{
    std::ostringstream stream;

    stream
        << std::fixed
        << std::setprecision(6)
        << value;

    return stream.str();
}

using EdgeKey =
    std::pair<
        iso::mc::Index,
        iso::mc::Index>;

EdgeKey make_edge(
    iso::mc::Index a,
    iso::mc::Index b)
{
    if (a > b)
        std::swap(a, b);

    return {a, b};
}

struct EdgeStatistics
{
    std::size_t boundary_edges = 0;
    std::size_t nonmanifold_edges = 0;
};

EdgeStatistics compute_edge_statistics(
    const iso::mc::Mesh& extracted)
{
    std::map<
        EdgeKey,
        std::size_t> edge_use;

    for (const auto& triangle :
         extracted.triangles)
    {
        ++edge_use[
            make_edge(
                triangle[0],
                triangle[1])];

        ++edge_use[
            make_edge(
                triangle[1],
                triangle[2])];

        ++edge_use[
            make_edge(
                triangle[2],
                triangle[0])];
    }

    EdgeStatistics statistics;

    for (const auto& entry :
         edge_use)
    {
        const std::size_t uses =
            entry.second;

        if (uses == 1)
            ++statistics.boundary_edges;
        else if (uses > 2)
            ++statistics.nonmanifold_edges;
    }

    return statistics;
}

std::vector<int> compute_connected_components(
    const pmp::SurfaceMesh& mesh,
    std::vector<std::size_t>& component_sizes)
{
    std::vector<int> component(
        mesh.n_vertices(),
        -1);

    int component_count = 0;

    for (auto start :
         mesh.vertices())
    {
        const int start_id =
            start.idx();

        if (component[start_id] != -1)
            continue;

        const int current_component =
            component_count++;

        std::queue<pmp::Vertex> queue;

        queue.push(start);

        component[start_id] =
            current_component;

        std::size_t component_size = 0;

        while (!queue.empty())
        {
            const pmp::Vertex v =
                queue.front();

            queue.pop();

            ++component_size;

            for (auto vv :
                 mesh.vertices(v))
            {
                const int vv_id =
                    vv.idx();

                if (component[vv_id] == -1)
                {
                    component[vv_id] =
                        current_component;

                    queue.push(vv);
                }
            }
        }

        component_sizes.push_back(
            component_size);
    }

    return component;
}

struct CoincidentStatistics
{
    std::size_t groups = 0;
    std::size_t extra_vertices = 0;
};

CoincidentStatistics compute_coincident_statistics(
    const pmp::SurfaceMesh& mesh)
{
    using PositionKey =
        std::tuple<
            float,
            float,
            float>;

    std::map<
        PositionKey,
        std::vector<int>> groups;

    for (auto v :
         mesh.vertices())
    {
        const pmp::Point p =
            mesh.position(v);

        groups[
            PositionKey(
                p[0],
                p[1],
                p[2])]
            .push_back(
                v.idx());
    }

    CoincidentStatistics statistics;

    for (const auto& entry :
         groups)
    {
        const auto& ids =
            entry.second;

        if (ids.size() > 1)
        {
            ++statistics.groups;

            statistics.extra_vertices +=
                ids.size() - 1;
        }
    }

    return statistics;
}

struct GeometricSingularity
{
    pmp::Point center{
        0.0f,
        0.0f,
        0.0f
    };

    SingularityType type =
        SingularityType::NonDegenerateSingular;

    std::vector<int> representative_vertices;

    double hessian_determinant = 0.0;

    Eigen::Vector3d eigenvalues =
        Eigen::Vector3d::Zero();
};

std::vector<GeometricSingularity>
group_geometric_singularities(
    const pmp::SurfaceMesh& mesh,
    const std::vector<
        SingularityClassificationResult>& results)
{
    constexpr double position_tolerance =
        1e-8;

    std::vector<int> singular_vertices;

    for (auto v :
         mesh.vertices())
    {
        if (results[v.idx()].type !=
            SingularityType::Regular)
        {
            singular_vertices.push_back(
                v.idx());
        }
    }

    std::vector<bool> assigned(
        singular_vertices.size(),
        false);

    std::vector<
        GeometricSingularity> groups;

    for (std::size_t i = 0;
         i < singular_vertices.size();
         ++i)
    {
        if (assigned[i])
            continue;

        const int seed_id =
            singular_vertices[i];

        const pmp::Point seed =
            mesh.position(
                pmp::Vertex(seed_id));

        GeometricSingularity group;

        group.type =
            results[seed_id].type;

        group.hessian_determinant =
            results[seed_id]
                .hessian_determinant;

        group.eigenvalues =
            results[seed_id]
                .eigenvalues;

        group.representative_vertices
            .push_back(
                seed_id);

        assigned[i] = true;

        for (std::size_t j = i + 1;
             j < singular_vertices.size();
             ++j)
        {
            if (assigned[j])
                continue;

            const int candidate_id =
                singular_vertices[j];

            const pmp::Point candidate =
                mesh.position(
                    pmp::Vertex(candidate_id));

            const double dx =
                static_cast<double>(
                    candidate[0] -
                    seed[0]);

            const double dy =
                static_cast<double>(
                    candidate[1] -
                    seed[1]);

            const double dz =
                static_cast<double>(
                    candidate[2] -
                    seed[2]);

            const double distance =
                std::sqrt(
                    dx * dx +
                    dy * dy +
                    dz * dz);

            if (distance <=
                position_tolerance)
            {
                group.representative_vertices
                    .push_back(
                        candidate_id);

                assigned[j] = true;
            }
        }

        pmp::Point center(
            0.0f,
            0.0f,
            0.0f);

        for (int id :
             group.representative_vertices)
        {
            center +=
                mesh.position(
                    pmp::Vertex(id));
        }

        center /=
            static_cast<float>(
                group.representative_vertices
                    .size());

        group.center =
            center;

        groups.push_back(
            std::move(group));
    }

    return groups;
}

void print_topology_report(
    const iso::mc::Mesh& extracted,
    const pmp::SurfaceMesh& mesh,
    const std::vector<
        std::size_t>& component_sizes)
{
    std::cout
        << "\n========================================================\n"
        << "  MESH TOPOLOGY DIAGNOSTICS\n"
        << "========================================================\n";

    std::cout
        << "Vertices           : "
        << mesh.n_vertices()
        << '\n';

    std::cout
        << "Edges              : "
        << mesh.n_edges()
        << '\n';

    std::cout
        << "Triangles          : "
        << mesh.n_faces()
        << '\n';

    const EdgeStatistics
        edge_statistics =
            compute_edge_statistics(
                extracted);

    std::size_t pmp_boundary_edges = 0;

    for (auto e :
         mesh.edges())
    {
        if (mesh.is_boundary(e))
            ++pmp_boundary_edges;
    }

    std::size_t nonmanifold_vertices = 0;

    for (auto v :
         mesh.vertices())
    {
        if (!mesh.is_manifold(v))
            ++nonmanifold_vertices;
    }

    const CoincidentStatistics
        coincident =
            compute_coincident_statistics(
                mesh);

    std::cout
        << "\nBoundary edges      : "
        << edge_statistics.boundary_edges
        << '\n';

    std::cout
        << "PMP boundary edges  : "
        << pmp_boundary_edges
        << '\n';

    std::cout
        << "Non-manifold edges  : "
        << edge_statistics.nonmanifold_edges
        << '\n';

    std::cout
        << "Non-manifold vertices: "
        << nonmanifold_vertices
        << '\n';

    std::cout
        << "Connected components : "
        << component_sizes.size()
        << '\n';

    for (std::size_t i = 0;
         i < component_sizes.size();
         ++i)
    {
        std::cout
            << "  Component "
            << i
            << " vertices: "
            << component_sizes[i]
            << '\n';
    }

    std::cout
        << "Coincident vertex groups: "
        << coincident.groups
        << '\n';

    std::cout
        << "Coincident extra vertices: "
        << coincident.extra_vertices
        << '\n';

    std::cout
        << "========================================================\n";
}

pmp::SurfaceMesh to_pmp_mesh(
    const iso::mc::Mesh& source)
{
    pmp::SurfaceMesh mesh;

    std::vector<pmp::Vertex> vertices;

    vertices.reserve(
        source.vertices.size());

    for (const auto& p :
         source.vertices)
    {
        vertices.push_back(
            mesh.add_vertex(
                pmp::Point(
                    static_cast<float>(
                        p[0]),
                    static_cast<float>(
                        p[1]),
                    static_cast<float>(
                        p[2]))));
    }

    for (const auto& triangle :
         source.triangles)
    {
        mesh.add_triangle(
            vertices[triangle[0]],
            vertices[triangle[1]],
            vertices[triangle[2]]);
    }

    return mesh;
}

std::vector<pmp::Color>
make_local_colors(
    const iso::mc::Mesh& mesh,
    const iso::mc::Point& center,
    double radius)
{
    std::vector<pmp::Color> colors(
        mesh.vertices.size(),
        pmp::Color(
            0.6f,
            0.6f,
            0.6f));

    const double radius_squared =
        radius * radius;

    for (std::size_t i = 0;
         i < mesh.vertices.size();
         ++i)
    {
        const auto& p =
            mesh.vertices[i];

        const double dx =
            p[0] - center[0];

        const double dy =
            p[1] - center[1];

        const double dz =
            p[2] - center[2];

        const double d2 =
            dx * dx +
            dy * dy +
            dz * dz;

        if (d2 <= radius_squared)
        {
            colors[i] =
                pmp::Color(
                    0.2f,
                    0.6f,
                    1.0f);
        }

        if (std::sqrt(d2) <= 1e-8)
        {
            colors[i] =
                pmp::Color(
                    1.0f,
                    0.85f,
                    0.0f);
        }
    }

    return colors;
}

int main(
    int argc,
    char** argv)
{
    try
    {
        const double c_val =
            (argc > 1)
                ? std::stod(argv[1])
                : 0.0;

        const std::size_t resolution =
            (argc > 2)
                ? static_cast<std::size_t>(
                      std::stoul(argv[2]))
                : 129;

        const double extent =
            (argc > 3)
                ? std::stod(argv[3])
                : 2.0;

        const double local_radius =
            (argc > 4)
                ? std::stod(argv[4])
                : 0.3;

        const std::size_t local_levels =
            (argc > 5)
                ? static_cast<std::size_t>(
                      std::stoul(argv[5]))
                : 2;

        if (resolution < 2)
        {
            std::cerr
                << "ERROR: resolution must be at least 2.\n";

            return 1;
        }

        if (extent <= 0.0)
        {
            std::cerr
                << "ERROR: extent must be positive.\n";

            return 1;
        }

        if (local_radius <= 0.0)
        {
            std::cerr
                << "ERROR: local radius must be positive.\n";

            return 1;
        }

        ParameterizedConeQuadric quadric(
            c_val);

        SingularityDetector detector(
            quadric,
            1e-2,
            1e-3);

        iso::mc::Bounds bounds{
            {-extent, -extent, -extent},
            { extent,  extent, extent}
        };

        iso::mc::Options mc_options;

        mc_options.nx =
            resolution;

        mc_options.ny =
            resolution;

        mc_options.nz =
            resolution;

        mc_options.isovalue =
            0.0;

        std::cout
            << "========================================================\n"
            << "  MARCHING CUBES ISOSURFACE GENERATION\n"
            << "========================================================\n\n";

        std::cout
            << "Parameter c      : "
            << c_val
            << '\n';

        std::cout
            << "Grid resolution  : "
            << resolution
            << "^3\n";

        std::cout
            << "Sampling bounds  : [-"
            << extent
            << ", "
            << extent
            << "]^3\n";

        std::cout
            << "Local radius     : "
            << local_radius
            << '\n';

        std::cout
            << "Local levels     : "
            << local_levels
            << "\n\n";

        std::cout
            << "Sampling implicit surface on "
            << resolution
            << "^3 grid...\n";

        const auto mc_start =
            std::chrono::steady_clock::now();

        const iso::mc::Mesh extracted =
            iso::mc::extract(
                [&](double x,
                    double y,
                    double z)
                {
                    return quadric.eval(
                        x,
                        y,
                        z);
                },
                bounds,
                mc_options);

        const auto mc_end =
            std::chrono::steady_clock::now();

        const double mc_ms =
            std::chrono::duration<double,
                                  std::milli>(
                mc_end - mc_start)
                .count();

        std::cout
            << "Marching Cubes produced "
            << extracted.vertices.size()
            << " vertices and "
            << extracted.triangles.size()
            << " triangles.\n";

        std::cout
            << "Marching Cubes time : "
            << mc_ms
            << " ms\n";

        pmp::SurfaceMesh mesh =
            to_pmp_mesh(extracted);

        std::vector<std::size_t>
            component_sizes;

        const std::vector<int>
            component =
                compute_connected_components(
                    mesh,
                    component_sizes);

        print_topology_report(
            extracted,
            mesh,
            component_sizes);

        std::vector<pmp::Vertex>
            apex_vertices;

        for (auto v :
             mesh.vertices())
        {
            const pmp::Point p =
                mesh.position(v);

            if (p[0] == 0.0f &&
                p[1] == 0.0f &&
                p[2] == 0.0f)
            {
                apex_vertices.push_back(v);
            }
        }

        std::cout
            << "\nApex vertices at (0,0,0): "
            << apex_vertices.size()
            << '\n';

        for (auto v :
             apex_vertices)
        {
            std::cout
                << "  Vertex "
                << v.idx()
                << " -> component "
                << component[v.idx()]
                << ", valence "
                << mesh.valence(v)
                << '\n';

            std::cout
                << "    neighbors:";

            for (auto vv :
                 mesh.vertices(v))
            {
                const pmp::Point q =
                    mesh.position(vv);

                std::cout
                    << " "
                    << vv.idx()
                    << "("
                    << q[0]
                    << ","
                    << q[1]
                    << ","
                    << q[2]
                    << ")";
            }

            std::cout << '\n';
        }

        if (apex_vertices.size() > 1)
        {
            bool separate_components = false;

            for (std::size_t i = 0;
                 i < apex_vertices.size();
                 ++i)
            {
                for (std::size_t j = i + 1;
                     j < apex_vertices.size();
                     ++j)
                {
                    if (component[
                            apex_vertices[i].idx()]
                        !=
                        component[
                            apex_vertices[j].idx()])
                    {
                        separate_components =
                            true;
                    }
                }
            }

            std::cout
                << "\nApex vertices belong to different "
                   "components: "
                << (separate_components
                        ? "YES"
                        : "NO")
                << '\n';
        }

        std::vector<pmp::Color> colors(
            mesh.n_vertices());

        std::vector<
            SingularityClassificationResult>
            results(
                mesh.n_vertices());

        int regular_count = 0;
        int nondegenerate_count = 0;
        int degenerate_count = 0;

        std::cout
            << "\n========================================================\n"
            << "  SINGULARITY DIAGNOSTICS\n"
            << "========================================================\n";

        for (auto v :
             mesh.vertices())
        {
            const pmp::Point p =
                mesh.position(v);

            const SingularityClassificationResult
                result =
                    detector.classifyPoint(p);

            results[v.idx()] =
                result;

            if (result.type ==
                SingularityType::Regular)
            {
                colors[v.idx()] =
                    pmp::Color(
                        0.6f,
                        0.6f,
                        0.6f);

                ++regular_count;

                continue;
            }

            std::cout
                << "\n--- Singular candidate ---\n"
                << std::setprecision(12);

            std::cout
                << "Vertex index: "
                << v.idx()
                << '\n';

            std::cout
                << "Position: ("
                << result.position[0]
                << ", "
                << result.position[1]
                << ", "
                << result.position[2]
                << ")\n";

            std::cout
                << "f(p): "
                << result.function_value
                << '\n';

            std::cout
                << "Gradient: ("
                << result.gradient[0]
                << ", "
                << result.gradient[1]
                << ", "
                << result.gradient[2]
                << ")\n";

            std::cout
                << "|Gradient|: "
                << result.gradient_norm
                << '\n';

            std::cout
                << "Hessian determinant: "
                << result.hessian_determinant
                << '\n';

            std::cout
                << "Hessian eigenvalues: ("
                << result.eigenvalues[0]
                << ", "
                << result.eigenvalues[1]
                << ", "
                << result.eigenvalues[2]
                << ")\n";

            std::cout
                << "Classification: "
                << singularity_type_name(
                       result.type)
                << '\n';

            if (result.type ==
                SingularityType::NonDegenerateSingular)
            {
                colors[v.idx()] =
                    pmp::Color(
                        1.0f,
                        0.85f,
                        0.0f);

                ++nondegenerate_count;
            }
            else
            {
                colors[v.idx()] =
                    pmp::Color(
                        1.0f,
                        0.0f,
                        0.0f);

                ++degenerate_count;
            }
        }

        const std::vector<
            GeometricSingularity>
            geometric_singularities =
                group_geometric_singularities(
                    mesh,
                    results);

        std::cout
            << "\n--- Geometric Singularity Groups ---\n";

        std::cout
            << "Geometric singularities: "
            << geometric_singularities.size()
            << '\n';

        for (std::size_t i = 0;
             i < geometric_singularities.size();
             ++i)
        {
            const auto& singularity =
                geometric_singularities[i];

            std::cout
                << "\nSingularity #"
                << i
                << '\n';

            std::cout
                << "  Center: ("
                << singularity.center[0]
                << ", "
                << singularity.center[1]
                << ", "
                << singularity.center[2]
                << ")\n";

            std::cout
                << "  Type: "
                << singularity_type_name(
                       singularity.type)
                << '\n';

            std::cout
                << "  Mesh representatives: {";

            for (std::size_t j = 0;
                 j <
                 singularity
                     .representative_vertices
                     .size();
                 ++j)
            {
                if (j > 0)
                    std::cout << ", ";

                std::cout
                    << singularity
                           .representative_vertices[j];
            }

            std::cout
                << "}\n";

            std::cout
                << "  Hessian determinant: "
                << singularity.hessian_determinant
                << '\n';

            std::cout
                << "  Hessian eigenvalues: ("
                << singularity.eigenvalues[0]
                << ", "
                << singularity.eigenvalues[1]
                << ", "
                << singularity.eigenvalues[2]
                << ")\n";

            std::cout
                << "  Coincident mesh representatives: "
                << (singularity
                            .representative_vertices
                            .size() >
                        1
                        ? "YES"
                        : "NO")
                << '\n';
        }

        if (!geometric_singularities.empty())
        {
            const auto& singularity =
                geometric_singularities.front();

            iso::mc::Point center{
                static_cast<double>(
                    singularity.center[0]),
                static_cast<double>(
                    singularity.center[1]),
                static_cast<double>(
                    singularity.center[2])
            };

            int vertices_in_radius = 0;

            for (const auto& p :
                 extracted.vertices)
            {
                const double dx =
                    p[0] - center[0];

                const double dy =
                    p[1] - center[1];

                const double dz =
                    p[2] - center[2];

                const double d2 =
                    dx * dx +
                    dy * dy +
                    dz * dz;

                if (d2 <=
                    local_radius *
                    local_radius)
                {
                    ++vertices_in_radius;
                }
            }

            std::cout
                << "\n--- Localization Results ---\n";

            std::cout
                << "Singularity Center Estimated at: ("
                << center[0]
                << ", "
                << center[1]
                << ", "
                << center[2]
                << ")\n";

            std::cout
                << "Local Unfolding Region (Radius "
                << local_radius
                << ") contains "
                << vertices_in_radius
                << " vertices.\n";

            // ---------------------------------------------------------------
            // Local refinement.
            // ---------------------------------------------------------------
            LocalUnfolder local_unfolder;

            LocalUnfolder::Options
                local_options;

            local_options.radius =
                local_radius;

            local_options.levels =
                local_levels;

            const auto local_start =
                std::chrono::steady_clock::now();

            const LocalUnfolder::Result
                local_result =
                    local_unfolder.refine(
                        extracted,
                        quadric,
                        center,
                        local_options);

            const auto local_end =
                std::chrono::steady_clock::now();

            const double local_ms =
                std::chrono::duration<
                    double,
                    std::milli>(
                    local_end -
                    local_start)
                    .count();

            std::cout
                << "\n========================================================\n"
                << "  LOCAL UNFOLDING / REFINEMENT\n"
                << "========================================================\n";

            std::cout
                << "Input vertices      : "
                << local_result.input_vertices
                << '\n';

            std::cout
                << "Input triangles     : "
                << local_result.input_triangles
                << '\n';

            std::cout
                << "Region vertices     : "
                << local_result.vertices_in_region
                << '\n';

            std::cout
                << "Region faces        : "
                << local_result.faces_in_region
                << '\n';

            std::cout
                << "Refinement levels   : "
                << local_result.levels.size()
                << '\n';

            for (const auto& level :
                 local_result.levels)
            {
                std::cout
                    << "\nLevel "
                    << level.level
                    << ":\n";

                std::cout
                    << "  Selected faces    : "
                    << level.selected_faces
                    << '\n';

                std::cout
                    << "  New vertices      : "
                    << level.new_vertices
                    << '\n';

                std::cout
                    << "  New triangles     : "
                    << level.new_triangles
                    << '\n';

                std::cout
                    << "  Output vertices   : "
                    << level.output_vertices
                    << '\n';

                std::cout
                    << "  Output triangles  : "
                    << level.output_triangles
                    << '\n';

                std::cout
                    << "  Projection failures: "
                    << level.projection_failures
                    << '\n';
            }

            std::cout
                << "\nLocal refinement time: "
                << local_ms
                << " ms\n";

            std::cout
                << "Final local vertices: "
                << local_result.output_vertices
                << '\n';

            std::cout
                << "Final local triangles: "
                << local_result.output_triangles
                << '\n';

            pmp::SurfaceMesh
                local_mesh =
                    to_pmp_mesh(
                        local_result.mesh);

            std::vector<std::size_t>
                local_components;

            compute_connected_components(
                local_mesh,
                local_components);

            std::cout
                << "\nLocal refined topology:\n";

            print_topology_report(
                local_result.mesh,
                local_mesh,
                local_components);

            const std::string
                c_string =
                    format_parameter(c_val);

            const std::filesystem::path
                local_ply =
                    std::filesystem::path(
                        "outputs") /
                    ("local_unfolded_c_" +
                     c_string +
                     "_l" +
                     std::to_string(
                         local_levels) +
                     ".ply");

            const std::filesystem::path
                local_obj =
                    std::filesystem::path(
                        "outputs") /
                    ("local_unfolded_c_" +
                     c_string +
                     "_l" +
                     std::to_string(
                         local_levels) +
                     ".obj");

            const std::vector<pmp::Color>
                local_colors =
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
                << "\nSaved local PLY:\n  "
                << local_ply.string()
                << '\n';

            std::cout
                << "Saved local OBJ:\n  "
                << local_obj.string()
                << '\n';

            colors.clear();

            // Open the locally refined mesh for the singular experiment.
            iso::MyViewer window(
                "Local Unfolding / Refinement",
                1024,
                768);

            window.load_mesh(
                local_obj.string().c_str());

            return window.run();
        }

        // ---------------------------------------------------------------------
        // No singularity: save and view only the baseline mesh.
        // ---------------------------------------------------------------------
        std::cout
            << "\nClassification Results:\n";

        std::cout
            << "  Regular (Gray): "
            << regular_count
            << '\n';

        std::cout
            << "  Morse/NonDegenerate (Yellow): "
            << nondegenerate_count
            << '\n';

        std::cout
            << "  Degenerate (Red): "
            << degenerate_count
            << '\n';

        const std::filesystem::path
            output_dir =
                "outputs";

        std::filesystem::create_directories(
            output_dir);

        const std::string
            c_string =
                format_parameter(c_val);

        const std::filesystem::path
            out_ply =
                output_dir /
                ("colored_result_c_" +
                 c_string +
                 ".ply");

        const std::filesystem::path
            out_obj =
                output_dir /
                ("viewer_mesh_c_" +
                 c_string +
                 ".obj");

        write_colored_ply(
            mesh,
            colors,
            out_ply.string());

        std::cout
            << "\nSaved colored PLY:\n  "
            << out_ply.string()
            << '\n';

        pmp::write(
            mesh,
            out_obj.string());

        std::cout
            << "Saved viewer OBJ:\n  "
            << out_obj.string()
            << '\n';

        std::cout
            << "\nOpening Viewer...\n";

        iso::MyViewer window(
            "Marching Cubes / Singularity Viewer",
            1024,
            768);

        window.load_mesh(
            out_obj.string().c_str());

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