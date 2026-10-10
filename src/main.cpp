#include "ImplicitSurface.hpp"
#include "LocalMarchingCubes.hpp"
#include "LocalUnfolder.hpp"
#include "LocalVoxelStitcher.hpp"
#include "MarchingCubes.hpp"
#include "SingularityDetector.hpp"
#include "app/MyViewer.h"
#include "MeshGeometryCheck.hpp"
#include "ScalarFieldPly.hpp"
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
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <map>
#include <utility>
#include <vector>
#include <array>
#include <queue>
namespace
{
std::unique_ptr<ImplicitSurface> make_surface(
    const std::string& name,
    double c)
{
    if (name == "cone")
        return std::make_unique<ParameterizedConeQuadric>(c);
    if (name == "elliptic-cone")
        return std::make_unique<EllipticConeSurface>(c);
    if (name == "rotated-cone")
        return std::make_unique<RotatedEllipticConeSurface>(c);
    if (name == "x-cone")
        return std::make_unique<XAxisEllipticConeSurface>(c);
    if (name == "y-cone")
        return std::make_unique<YAxisEllipticConeSurface>(c);
    if (name == "tilted-cone")
        return std::make_unique<TiltedEllipticConeSurface>(c);
    if (name == "crossing-planes")
        return std::make_unique<IntersectingPlanesSurface>(c);
    if (name == "quartic")
        return std::make_unique<QuarticSaddleSurface>(c, 0.05);
    if (name == "sphere")
        return std::make_unique<SphereImplicitSurface>(1.0);
    throw std::invalid_argument(
        "Unknown surface. Use cone, elliptic-cone, rotated-cone, x-cone, y-cone, tilted-cone, crossing-planes, quartic or sphere.");
}
pmp::SurfaceMesh to_pmp_mesh(const iso::mc::Mesh& source)
{
    using Index = iso::mc::Index;
    using Edge = std::pair<Index, Index>;
    pmp::SurfaceMesh mesh;
    const std::size_t face_count = source.triangles.size();
    const std::size_t corner_count = face_count * 3;
    if (face_count == 0)
        return mesh;
    /*
     * PMP requires each vertex's incident faces to form a manifold fan.
     *
     * The internal MC mesh may intentionally contain a singular vertex
     * where separate face fans meet at the same geometric position.
     *
     * Use one disjoint-set element per triangle corner. Corners are joined
     * only when their faces share an edge containing the same source vertex.
     * Thus, disconnected fans at a singularity are exported as distinct
     * PMP vertices with identical coordinates.
     *
     * The source mesh itself is not modified.
     */
    std::vector<std::size_t> parent(corner_count);
    std::vector<std::size_t> tree_rank(corner_count, 0);
    for (std::size_t i = 0; i < corner_count; ++i)
        parent[i] = i;
    const auto find_root = [&](std::size_t x)
    {
        while (parent[x] != x)
        {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    };
    const auto unite = [&](std::size_t a, std::size_t b)
    {
        a = find_root(a);
        b = find_root(b);
        if (a == b)
            return;
        if (tree_rank[a] < tree_rank[b])
            std::swap(a, b);
        parent[b] = a;
        if (tree_rank[a] == tree_rank[b])
            ++tree_rank[a];
    };
    const auto make_edge = [](Index a, Index b) -> Edge
    {
        if (a > b)
            std::swap(a, b);
        return {a, b};
    };
    /*
     * Find the triangles incident to every undirected edge.
     */
    std::map<Edge, std::vector<std::size_t>> edge_faces;
    for (std::size_t fi = 0; fi < face_count; ++fi)
    {
        const auto& t = source.triangles[fi];
        edge_faces[make_edge(t[0], t[1])].push_back(fi);
        edge_faces[make_edge(t[1], t[2])].push_back(fi);
        edge_faces[make_edge(t[2], t[0])].push_back(fi);
    }
    /*
     * Return the corner-node index for a source vertex in one triangle.
     */
    const auto corner_for = [&](std::size_t face, Index vertex)
        -> std::size_t
    {
        const auto& t = source.triangles[face];
        for (std::size_t c = 0; c < 3; ++c)
        {
            if (t[c] == vertex)
                return 3 * face + c;
        }
        throw std::runtime_error(
            "PMP conversion: edge endpoint not found in triangle");
    };
    /*
     * Adjacent faces belong to the same vertex fan when they share an edge.
     * Join their corner nodes at both endpoints of each manifold edge.
     */
    for (const auto& [edge, faces] : edge_faces)
    {
        if (faces.size() > 2)
        {
            throw std::runtime_error(
                "PMP conversion cannot represent a non-manifold edge");
        }
        if (faces.size() != 2)
            continue;
        const std::size_t f0 = faces[0];
        const std::size_t f1 = faces[1];
        unite(
            corner_for(f0, edge.first),
            corner_for(f1, edge.first));
        unite(
            corner_for(f0, edge.second),
            corner_for(f1, edge.second));
    }
    /*
     * Create one PMP vertex per connected face fan.
     *
     * Disconnected fans sharing a source vertex receive distinct PMP
     * vertices, but those vertices retain the exact same position.
     */
    std::map<std::size_t, pmp::Vertex> fan_vertices;
    std::vector<pmp::Vertex> corner_vertices(corner_count);
    for (std::size_t fi = 0; fi < face_count; ++fi)
    {
        const auto& t = source.triangles[fi];
        for (std::size_t c = 0; c < 3; ++c)
        {
            const std::size_t corner = 3 * fi + c;
            const std::size_t root = find_root(corner);
            auto it = fan_vertices.find(root);
            if (it == fan_vertices.end())
            {
                const auto& p = source.vertices[t[c]];
                const pmp::Vertex pv = mesh.add_vertex(
                    pmp::Point(
                        static_cast<float>(p[0]),
                        static_cast<float>(p[1]),
                        static_cast<float>(p[2])));
                it = fan_vertices.emplace(root, pv).first;
            }
            corner_vertices[corner] = it->second;
        }
    }
    /*
     * Insert the triangles using the fan-aware vertex mapping.
     */
    for (std::size_t fi = 0; fi < face_count; ++fi)
    {
        mesh.add_triangle(
            corner_vertices[3 * fi],
            corner_vertices[3 * fi + 1],
            corner_vertices[3 * fi + 2]);
    }
    return mesh;
}
bool has_flag(int argc, char** argv, const std::string& flag)
{
    for (int i = 1; i < argc; ++i)
    {
        if (std::string(argv[i]) == flag)
            return true;
    }
    return false;
}
std::string option_value(
    int argc,
    char** argv,
    const std::string& option,
    const std::string& fallback)
{
    for (int i = 1; i + 1 < argc; ++i)
    {
        if (std::string(argv[i]) == option)
            return argv[i + 1];
    }
    return fallback;
}
std::string singularity_type_name(SingularityType type)
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
struct SingularityCandidate
{
    bool found = false;
    iso::mc::Point center{0.0, 0.0, 0.0};
    double gradient_norm =
        std::numeric_limits<double>::infinity();
    double function_value = 0.0;
    SingularityType type = SingularityType::Regular;
};
SingularityCandidate find_surface_singularity(
    const iso::mc::Mesh& mesh,
    const SingularityDetector& detector)
{
    SingularityCandidate candidate;
    for (const auto& p : mesh.vertices)
    {
        const pmp::Point point(
            static_cast<float>(p[0]),
            static_cast<float>(p[1]),
            static_cast<float>(p[2]));
        const auto result = detector.classifyPoint(point);
        if (result.type == SingularityType::Regular)
            continue;
        if (std::abs(result.function_value) > 1e-5)
            continue;
        if (!candidate.found ||
            result.gradient_norm < candidate.gradient_norm)
        {
            candidate.found = true;
            candidate.center = p;
            candidate.gradient_norm = result.gradient_norm;
            candidate.function_value = result.function_value;
            candidate.type = result.type;
        }
    }
    return candidate;
}
// Count edge-connected triangle components that contain at least one vertex
// within tolerance of the detected singular point. This deliberately counts
// separate face fans at the same geometric point as distinct components.
std::size_t count_components_touching_point(
    const iso::mc::Mesh& mesh,
    const iso::mc::Point& point,
    double tolerance)
{
    using Index = iso::mc::Index;
    using Edge = std::pair<Index, Index>;

    const std::size_t face_count = mesh.triangles.size();
    if (face_count == 0)
        return 0;

    const auto make_edge = [](Index a, Index b) -> Edge
    {
        if (a > b)
            std::swap(a, b);
        return {a, b};
    };

    std::map<Edge, std::vector<std::size_t>> edge_faces;
    for (std::size_t fi = 0; fi < face_count; ++fi)
    {
        const auto& t = mesh.triangles[fi];
        edge_faces[make_edge(t[0], t[1])].push_back(fi);
        edge_faces[make_edge(t[1], t[2])].push_back(fi);
        edge_faces[make_edge(t[2], t[0])].push_back(fi);
    }

    std::vector<std::vector<std::size_t>> neighbors(face_count);
    for (const auto& entry : edge_faces)
    {
        const auto& faces = entry.second;
        for (std::size_t i = 0; i < faces.size(); ++i)
        {
            for (std::size_t j = i + 1; j < faces.size(); ++j)
            {
                neighbors[faces[i]].push_back(faces[j]);
                neighbors[faces[j]].push_back(faces[i]);
            }
        }
    }

    const double tolerance_squared = tolerance * tolerance;
    const auto vertex_is_near_point = [&](Index vertex) -> bool
    {
        const auto& p = mesh.vertices[vertex];
        const double dx = p[0] - point[0];
        const double dy = p[1] - point[1];
        const double dz = p[2] - point[2];
        return dx * dx + dy * dy + dz * dz <= tolerance_squared;
    };

    std::vector<char> visited(face_count, 0);
    std::queue<std::size_t> pending;
    std::size_t touching_components = 0;

    for (std::size_t seed = 0; seed < face_count; ++seed)
    {
        if (visited[seed])
            continue;

        bool touches_point = false;
        visited[seed] = 1;
        pending.push(seed);

        while (!pending.empty())
        {
            const std::size_t fi = pending.front();
            pending.pop();

            const auto& triangle = mesh.triangles[fi];
            if (vertex_is_near_point(triangle[0]) ||
                vertex_is_near_point(triangle[1]) ||
                vertex_is_near_point(triangle[2]))
            {
                touches_point = true;
            }

            for (const std::size_t adjacent : neighbors[fi])
            {
                if (visited[adjacent])
                    continue;
                visited[adjacent] = 1;
                pending.push(adjacent);
            }
        }

        if (touches_point)
            ++touching_components;
    }

    return touching_components;
}

void print_topology(
    const std::string& label,
    const iso::stitch::TopologyStatistics& topology)
{
    std::cout
        << label << '\n'
        << "  vertices           : " << topology.vertices << '\n'
        << "  triangles          : " << topology.triangles << '\n'
        << "  edges              : " << topology.edges << '\n'
        << "  boundary edges     : " << topology.boundary_edges << '\n'
        << "  boundary vertices  : " << topology.boundary_vertices << '\n'
        << "  non-manifold edges : " << topology.nonmanifold_edges << '\n'
        << "  components         : "
        << topology.connected_components << '\n';
}
std::filesystem::path output_path(
    const std::filesystem::path& dir,
    const std::string& prefix,
    const std::string& surface,
    double c,
    std::size_t resolution,
    std::size_t level)
{
    std::ostringstream c_string;
    c_string << std::fixed << std::setprecision(6) << c;
    return dir /
        (prefix + "_" +
         surface + "_c_" +
         c_string.str() +
         "_r_" + std::to_string(resolution) +
         "_l_" + std::to_string(level) +
         ".obj");
}
void save_mesh(
    const iso::mc::Mesh& mesh,
    const std::filesystem::path& filename)
{
    const auto pmp_mesh = to_pmp_mesh(mesh);
    pmp::write(pmp_mesh, filename.string());
    std::cout << "Saved: " << filename << '\n';
}
void write_report(
    const std::filesystem::path& filename,
    const std::string& surface_name,
    double c,
    std::size_t resolution,
    double extent,
    double radius,
    std::size_t level,
    const SingularityCandidate& singularity,
    std::size_t global_singular_point_components,
    std::size_t final_singular_point_components,
    const iso::stitch::TopologyStatistics& global_topology,
    const LocalMarchingCubes::Result& patch,
    const iso::stitch::TopologyStatistics& patch_topology,
    const iso::stitch::StitchResult* stitched,
    const std::chrono::duration<double, std::milli>& global_duration)
{
    std::ofstream report(filename);
    if (!report)
        return;
    report << "Implicit isosurface experiment\n";
    report << "==============================\n\n";
    report << "Surface: " << surface_name << '\n';
    report << "c: " << std::setprecision(17) << c << '\n';
    report << "Global resolution: " << resolution << "^3\n";
    report << "Extent: " << extent << '\n';
    report << "Local radius: " << radius << '\n';
    report << "Local level: " << level << "\n\n";
    report << "Global MC\n";
    report << "---------\n";
    report << "Extraction time ms: "
           << global_duration.count() << '\n';
    report << "Vertices: " << global_topology.vertices << '\n';
    report << "Triangles: " << global_topology.triangles << '\n';
    report << "Edges: " << global_topology.edges << '\n';
    report << "Boundary edges: "
           << global_topology.boundary_edges << '\n';
    report << "Non-manifold edges: "
           << global_topology.nonmanifold_edges << '\n';
    report << "Components: "
           << global_topology.connected_components << "\n\n";
    report << "Singularity\n";
    report << "-----------\n";
    report << "Found: " << (singularity.found ? "yes" : "no") << '\n';
    report << "Type: "
           << singularity_type_name(singularity.type) << '\n';
    report << "Position: "
           << singularity.center[0] << ' '
           << singularity.center[1] << ' '
           << singularity.center[2] << '\n';
    report << "Gradient norm: "
           << singularity.gradient_norm << '\n';
    report << "Function value: "
           << singularity.function_value << '\n';
    report << "Global edge-connected components touching singular point: "
           << global_singular_point_components << "\n\n";
    report << "Local voxel MC\n";
    report << "--------------\n";
    report << "Resolution: "
           << patch.resolution_x << 'x'
           << patch.resolution_y << 'x'
           << patch.resolution_z << '\n';
    report << "Vertices: "
           << patch.mesh.vertices.size() << '\n';
    report << "Triangles: "
           << patch.mesh.triangles.size() << '\n';
    report << "Boundary edges: "
           << patch.boundary_edges << '\n';
    report << "Boundary vertices: "
           << patch.boundary_vertices << '\n';
    report << "Non-manifold edges: "
           << patch.nonmanifold_edges << '\n';
    report << "Extraction time ms: "
           << patch.extraction_time_ms << '\n';
    report << "Topology edges: "
           << patch_topology.edges << '\n';
    report << "Topology components: "
           << patch_topology.connected_components << "\n\n";
    if (stitched)
    {
        report << "Stitching\n";
        report << "---------\n";
        report << "Success: "
               << (stitched->stitched ? "yes" : "no") << '\n';
        report << "Global triangles removed: "
               << stitched->global_triangles_removed << '\n';
        report << "Local triangles inserted: "
               << stitched->local_triangles_inserted << '\n';
        report << "Seam loops: "
               << stitched->seam_loops << '\n';
        report << "Seam vertices: "
               << stitched->seam_vertices << '\n';
        report << "Bridge triangles added (direct loop welding uses none): "
               << stitched->seam_triangles << '\n';
        report << "Max seam vertex distance: "
               << stitched->max_seam_vertex_distance << '\n';
        report << "Final vertices: "
               << stitched->topology.vertices << '\n';
        report << "Final triangles: "
               << stitched->topology.triangles << '\n';
        report << "Final edges: "
               << stitched->topology.edges << '\n';
        report << "Final boundary edges: "
               << stitched->topology.boundary_edges << '\n';
        report << "Final non-manifold edges: "
               << stitched->topology.nonmanifold_edges << '\n';
        report << "Final components: "
               << stitched->topology.connected_components << '\n';
        report << "Final edge-connected components touching singular point: "
               << final_singular_point_components << '\n';
    }
}
} // namespace
void orient_components_to_implicit_gradient(
    iso::mc::Mesh& mesh,
    const ImplicitSurface& surface)
{
    using Index = iso::mc::Index;
    using Edge = std::pair<Index, Index>;
    const std::size_t face_count = mesh.triangles.size();
    if (face_count == 0)
        return;
    const auto make_edge = [](Index a, Index b) -> Edge
    {
        if (a > b)
            std::swap(a, b);
        return {a, b};
    };
    // Build face adjacency through shared edges. Using edges rather
    // than vertex adjacency keeps the two cone nappes separate at
    // their singular apex.
    std::map<Edge, std::vector<std::size_t>> edge_faces;
    for (std::size_t fi = 0; fi < face_count; ++fi)
    {
        const auto& t = mesh.triangles[fi];
        edge_faces[make_edge(t[0], t[1])].push_back(fi);
        edge_faces[make_edge(t[1], t[2])].push_back(fi);
        edge_faces[make_edge(t[2], t[0])].push_back(fi);
    }
    std::vector<std::vector<std::size_t>> neighbors(face_count);
    for (const auto& [edge, faces] : edge_faces)
    {
        if (faces.size() != 2)
            continue;
        neighbors[faces[0]].push_back(faces[1]);
        neighbors[faces[1]].push_back(faces[0]);
    }
    // Central finite differences approximate grad(f) at a point.
    const auto gradient = [&](double x, double y, double z)
    {
        const double h = 1e-6 * std::max({
            1.0, std::abs(x), std::abs(y), std::abs(z)
        });
        return std::array<double, 3>{
            (surface.eval(x + h, y, z) -
             surface.eval(x - h, y, z)) / (2.0 * h),
            (surface.eval(x, y + h, z) -
             surface.eval(x, y - h, z)) / (2.0 * h),
            (surface.eval(x, y, z + h) -
             surface.eval(x, y, z - h)) / (2.0 * h)
        };
    };
    std::vector<char> visited(face_count, 0);
    std::queue<std::size_t> queue;
    for (std::size_t seed = 0; seed < face_count; ++seed)
    {
        if (visited[seed])
            continue;
        // First collect one edge-connected component. The current
        // stitcher's orientation pass should already make its faces
        // consistently wound relative to one another.
        std::vector<std::size_t> component;
        visited[seed] = 1;
        queue.push(seed);
        while (!queue.empty())
        {
            const std::size_t fi = queue.front();
            queue.pop();
            component.push_back(fi);
            for (const std::size_t adjacent : neighbors[fi])
            {
                if (visited[adjacent])
                    continue;
                visited[adjacent] = 1;
                queue.push(adjacent);
            }
        }
        // Positive score means the component's normals point toward
        // increasing scalar-field values; negative means the reverse.
        // The cross product weights each face by its area.
        double orientation_score = 0.0;
        for (const std::size_t fi : component)
        {
            const auto& t = mesh.triangles[fi];
            const auto& a = mesh.vertices[t[0]];
            const auto& b = mesh.vertices[t[1]];
            const auto& c = mesh.vertices[t[2]];
            const double ux = b[0] - a[0];
            const double uy = b[1] - a[1];
            const double uz = b[2] - a[2];
            const double vx = c[0] - a[0];
            const double vy = c[1] - a[1];
            const double vz = c[2] - a[2];
            const double nx = uy * vz - uz * vy;
            const double ny = uz * vx - ux * vz;
            const double nz = ux * vy - uy * vx;
            const double x = (a[0] + b[0] + c[0]) / 3.0;
            const double y = (a[1] + b[1] + c[1]) / 3.0;
            const double z = (a[2] + b[2] + c[2]) / 3.0;
            const auto grad = gradient(x, y, z);
            orientation_score +=
                nx * grad[0] +
                ny * grad[1] +
                nz * grad[2];
        }
        // Reverse the entire component, not individual faces, so
        // adjacent faces keep compatible edge directions.
        if (orientation_score < 0.0)
        {
            for (const std::size_t fi : component)
            {
                std::swap(
                    mesh.triangles[fi][1],
                    mesh.triangles[fi][2]);
            }
        }
    }
}
int main(int argc, char** argv)
{
    try
    {
        if (has_flag(argc, argv, "--help"))
        {
            std::cout
                << "Usage:\n"
                << "  Isosurfaces [c] [resolution] [extent] [radius] [level]\n"
                << "              [--surface cone|elliptic-cone|rotated-cone|x-cone|y-cone|tilted-cone|crossing-planes|quartic|sphere]\n"
                << "              [--output DIR]\n"
                << "              [--no-gui]\n"
                << "              [--no-stitch]\n"
                << "              [--selection-sweep]\n"
                << "  --selection-sweep compares sphere and topological selectors\n"
                << "  at six probe centers offset by 0.75 * radius from the\n"
                << "  detected singular point; writes a CSV diagnostic.\n";
            return 0;
        }
        const double c =
            argc > 1 ? std::stod(argv[1]) : 0.0;
        const std::size_t resolution =
            argc > 2
                ? static_cast<std::size_t>(std::stoul(argv[2]))
                : 65;
        const double extent =
            argc > 3 ? std::stod(argv[3]) : 2.0;
        const double radius =
            argc > 4 ? std::stod(argv[4]) : 0.30;
        const std::size_t level =
            argc > 5
                ? static_cast<std::size_t>(std::stoul(argv[5]))
                : 2;
        const std::string surface_name =
            option_value(
                argc,
                argv,
                "--surface",
                "cone");
        const std::filesystem::path output_dir =
            option_value(
                argc,
                argv,
                "--output",
                "outputs");
        const bool no_gui =
            has_flag(argc, argv, "--no-gui");
        const bool no_stitch =
            has_flag(argc, argv, "--no-stitch");
        const bool selection_sweep =
            has_flag(argc, argv, "--selection-sweep");
        if (resolution < 2)
            throw std::invalid_argument(
                "resolution must be >= 2");
        if (!(extent > 0.0))
            throw std::invalid_argument(
                "extent must be positive");
        if (!(radius > 0.0))
            throw std::invalid_argument(
                "radius must be positive");
        if (level == 0)
            throw std::invalid_argument(
                "level must be >= 1");
        if (surface_name != "cone" &&
            surface_name != "elliptic-cone" &&
            surface_name != "rotated-cone" &&
            surface_name != "x-cone" &&
            surface_name != "y-cone" &&
            surface_name != "tilted-cone" &&
            surface_name != "crossing-planes" &&
            surface_name != "quartic" &&
            surface_name != "sphere")
        {
            throw std::invalid_argument(
                "Unknown surface. Use cone, elliptic-cone, rotated-cone, x-cone, y-cone, tilted-cone, crossing-planes, quartic or sphere.");
        }
        std::filesystem::create_directories(output_dir);
        const auto surface =
            make_surface(surface_name, c);
        const iso::mc::Bounds bounds{
            {-extent, -extent, -extent},
            { extent,  extent, extent}};
        iso::mc::Options global_options;
        global_options.nx = resolution;
        global_options.ny = resolution;
        global_options.nz = resolution;
        global_options.isovalue = 0.0;
        const auto global_start =
            std::chrono::steady_clock::now();
        const iso::mc::Mesh global =
            iso::mc::extract(
                [&](double x, double y, double z)
                {
                    return surface->eval(x, y, z);
                },
                bounds,
                global_options);
        const auto global_end =
            std::chrono::steady_clock::now();
        const auto global_duration =
            std::chrono::duration<double, std::milli>(
                global_end - global_start);
        const auto global_topology =
            iso::stitch::compute_topology(global);
        std::cout
            << "========================================================\n"
            << "  IMPLICIT ISOSURFACE EXPERIMENT\n"
            << "========================================================\n"
            << "Surface             : " << surface_name << '\n'
            << "c                   : " << c << '\n'
            << "Global resolution   : " << resolution << "^3\n"
            << "Extent              : +/-" << extent << '\n'
            << "Local radius        : " << radius << '\n'
            << "Local level         : " << level << '\n'
            << "Output directory    : " << output_dir << '\n'
            << '\n';
        std::cout
            << "GLOBAL MARCHING CUBES\n"
            << "  extraction time    : "
            << global_duration.count() << " ms\n";
        print_topology("  topology:", global_topology);
        std::cout << '\n';
        const SingularityDetector detector(
            *surface,
            1e-2,
            1e-3);
        const auto singularity =
            find_surface_singularity(
                global,
                detector);
        // Scale the comparison tolerance to the sampled domain. The detected
        // point itself comes from a mesh vertex, while the final mesh can
        // contain independently extracted vertices at the same location.
        const double singular_point_tolerance =
            std::max(1e-8, 1e-6 * std::max(1.0, extent));
        const std::size_t global_singular_point_components =
            singularity.found
                ? count_components_touching_point(
                      global, singularity.center, singular_point_tolerance)
                : 0;
        if (!singularity.found)
        {
            std::cout
                << "No on-surface singularity was detected.\n";
            const auto global_filename =
                output_dir /
                (
                    "global_" +
                    surface_name +
                    "_c_" +
                    std::to_string(c) +
                    "_r_" +
                    std::to_string(resolution) +
                    ".obj");
            save_mesh(global, global_filename);
            if (no_gui)
                return 0;
            iso::MyViewer window(
                "Implicit Isosurface",
                1024,
                768);
            window.load_mesh(
                global_filename.string().c_str());
            window.set_marching_cubes(
                [&](double x, double y, double z)
                {
                    return surface->eval(x, y, z);
                },
                bounds,
                resolution,
                resolution,
                resolution,
                0.0);
            return window.run();
        }
        std::cout
            << "ON-SURFACE SINGULARITY\n"
            << "  position            : ("
            << singularity.center[0] << ", "
            << singularity.center[1] << ", "
            << singularity.center[2] << ")\n"
            << "  type                : "
            << singularity_type_name(
                   singularity.type)
            << '\n'
            << "  gradient norm       : "
            << singularity.gradient_norm << '\n'
            << "  function value      : "
            << singularity.function_value << '\n'
            << "  global components touching singular point: "
            << global_singular_point_components << "\n\n";
        /*
         * Local conforming baseline.
         *
         * We retain the full LevelStatistics reporting here, but do not
         * access a nonexistent LevelStatistics::mesh member. The refined
         * mesh storage remains owned by the LocalUnfolder result.
         */
        LocalUnfolder::Options local_options;
        local_options.radius = radius;
        local_options.levels = level;
        local_options.region_mode =
            LocalUnfolder::RegionMode::TopologicalBfs;
        local_options.sampling_bounds = bounds;
        local_options.base_resolution = resolution;
        local_options.isovalue = 0.0;
        LocalUnfolder unfolder;

        // Optional diagnostic: compare the two selectors in nearby regions
        // whose centers are deliberately offset from the detected singularity.
        // The regular experiment below still refines around the true singular
        // point; these probes are diagnostic only and do not change its mesh.
        if (selection_sweep)
        {
            struct Probe
            {
                const char* label;
                std::array<double, 3> direction;
            };
            const std::array<Probe, 6> probes{{
                {"+X", { 1.0,  0.0,  0.0}},
                {"-X", {-1.0,  0.0,  0.0}},
                {"+Y", { 0.0,  1.0,  0.0}},
                {"-Y", { 0.0, -1.0,  0.0}},
                {"+Z", { 0.0,  0.0,  1.0}},
                {"-Z", { 0.0,  0.0, -1.0}}
            }};
            const double offset_distance = 0.75 * radius;
            const auto sweep_filename = output_dir /
                ("selection_sweep_" + surface_name + "_c_" +
                 std::to_string(c) + "_r_" +
                 std::to_string(resolution) + ".csv");
            std::ofstream sweep_file(sweep_filename);
            if (!sweep_file)
                throw std::runtime_error(
                    "Could not create selection sweep CSV: " +
                    sweep_filename.string());
            sweep_file
                << "probe,center_x,center_y,center_z,radius,"
                << "sphere_faces,topological_faces,selection_changed,"
                << "refinement_level_available\n";

            std::size_t changed_probes = 0;
            std::cout
                << "SELECTOR COMPARISON DIAGNOSTIC\n"
                << "  Probe centers are offset by 0.75 * radius from the\n"
                << "  detected singular point. They do not alter the main\n"
                << "  refinement, which remains centered at the singularity.\n";

            for (const Probe& probe : probes)
            {
                const iso::mc::Point probe_center{
                    singularity.center[0] +
                        offset_distance * probe.direction[0],
                    singularity.center[1] +
                        offset_distance * probe.direction[1],
                    singularity.center[2] +
                        offset_distance * probe.direction[2]};

                LocalUnfolder::Options probe_options = local_options;
                probe_options.radius = radius;
                probe_options.levels = 1;
                // Scalar-grid statistics are unrelated to face-selection
                // comparison. Setting this to zero keeps the sweep cheap.
                probe_options.base_resolution = 0;

                const auto probe_result = unfolder.refine(
                    global, *surface, probe_center, probe_options);

                const bool has_level = !probe_result.levels.empty();
                const bool changed = has_level
                    ? probe_result.levels.front().selection_changed
                    : (probe_result.sphere_faces_in_region !=
                       probe_result.topological_faces_in_region);
                if (changed)
                    ++changed_probes;

                std::cout
                    << "  " << probe.label
                    << " center=(" << probe_center[0] << ", "
                    << probe_center[1] << ", "
                    << probe_center[2] << ")"
                    << " sphere=" << probe_result.sphere_faces_in_region
                    << " topological="
                    << probe_result.topological_faces_in_region
                    << " selection changed=" << (changed ? "YES" : "NO")
                    << (has_level ? "" : " (no refinement level)")
                    << '\n';

                sweep_file << std::setprecision(17)
                    << probe.label << ','
                    << probe_center[0] << ','
                    << probe_center[1] << ','
                    << probe_center[2] << ','
                    << radius << ','
                    << probe_result.sphere_faces_in_region << ','
                    << probe_result.topological_faces_in_region << ','
                    << (changed ? "yes" : "no") << ','
                    << (has_level ? "yes" : "no") << '\n';
            }
            std::cout
                << "  probes with different selected face sets: "
                << changed_probes << "/" << probes.size() << '\n'
                << "  CSV saved: " << sweep_filename << "\n\n";
        }

        const auto baseline_start =
            std::chrono::steady_clock::now();
        const auto baseline =
            unfolder.refine(
                global,
                *surface,
                singularity.center,
                local_options);
        const auto baseline_end =
            std::chrono::steady_clock::now();
        std::cout
            << "LOCAL CONFORMING REFINEMENT BASELINE\n"
            << "  levels returned     : "
            << baseline.levels.size() << "\n\n";
        for (std::size_t i = 0;
             i < baseline.levels.size();
             ++i)
        {
            const auto& result =
                baseline.levels[i];
            std::cout
                << "  Level " << (i + 1) << '\n'
                << "    selected faces      : "
                << result.selected_faces << '\n'
                << "    sphere faces        : "
                << result.selected_faces_sphere << '\n'
                << "    topological faces   : "
                << result.selected_faces_topological << '\n'
                << "    selection changed   : "
                << (result.selection_changed ? "YES" : "NO")
                << '\n'
                << "    output vertices     : "
                << result.output_vertices << '\n'
                << "    output triangles    : "
                << result.output_triangles << '\n'
                << "    interface edges     : "
                << result.interface_edges << '\n'
                << "    interface vertices  : "
                << result.interface_vertices << '\n'
                << "    projection failures : "
                << result.projection_failures << '\n'
                << "    scalar refinement   : "
                << result.scalar.refinement_factor << '\n'
                << "    scalar voxels       : "
                << result.scalar.selected_voxels << '\n'
                << "    scalar samples      : "
                << result.scalar.corner_value_samples << '\n'
                << "    scalar spacing      : "
                << result.scalar.spacing << '\n'
                << "    scalar min          : "
                << result.scalar.min_value << '\n'
                << "    scalar max          : "
                << result.scalar.max_value << '\n'
                << "    scalar mean         : "
                << result.scalar.mean_value << '\n';
        }
        const auto baseline_ms =
            std::chrono::duration<double, std::milli>(
                baseline_end - baseline_start)
                .count();
        std::cout
            << "  total baseline time  : "
            << baseline_ms << " ms\n\n";
        /*
         * Independently extracted local voxel MC patch.
         */
        LocalMarchingCubes local_mc;
        LocalMarchingCubes::Options patch_options;
        patch_options.radius = radius;
        patch_options.level = level;
        patch_options.padding_cells = 1;
        patch_options.sampling_bounds = bounds;
        patch_options.base_resolution = resolution;
        patch_options.isovalue = 0.0;
        const auto patch =
            local_mc.extract(
                *surface,
                singularity.center,
                patch_options);
        const auto patch_topology =
            iso::stitch::compute_topology(
                patch.mesh);
        std::cout
            << "LOCAL VOXEL MARCHING CUBES\n"
            << "  resolution          : "
            << patch.resolution_x << 'x'
            << patch.resolution_y << 'x'
            << patch.resolution_z << '\n'
            << "  vertices             : "
            << patch.mesh.vertices.size() << '\n'
            << "  triangles            : "
            << patch.mesh.triangles.size() << '\n'
            << "  boundary edges       : "
            << patch.boundary_edges << '\n'
            << "  boundary vertices    : "
            << patch.boundary_vertices << '\n'
            << "  non-manifold edges   : "
            << patch.nonmanifold_edges << '\n'
            << "  extraction time      : "
            << patch.extraction_time_ms << " ms\n";
        print_topology(
            "  topology:",
            patch_topology);
        std::cout << '\n';
        const auto local_filename =
            output_path(
                output_dir,
                "local_voxel",
                surface_name,
                c,
                resolution,
                level);
        save_mesh(
            patch.mesh,
            local_filename);
        const auto evaluate_scalar =
            [&](double x, double y, double z)
            {
                return surface->eval(x, y, z);
            };
        auto global_ply_filename =
            output_path(
                output_dir,
                "global_scalar",
                surface_name,
                c,
                resolution,
                level);
        global_ply_filename.replace_extension(".ply");
        iso::field_export::write_scalar_ply(
            global,
            global_ply_filename.string(),
            evaluate_scalar,
            global_options.isovalue,
            singularity.center,
            radius);
        std::cout << "Saved scalar PLY: "
                << global_ply_filename << '\n';
        auto global_vtk_filename = global_ply_filename;
        global_vtk_filename.replace_extension(".vtk");
        iso::field_export::write_scalar_vtk(
            global,
            global_vtk_filename.string(),
            evaluate_scalar,
            global_options.isovalue,
            singularity.center,
            radius);
        std::cout << "Saved scalar VTK: "
                << global_vtk_filename << '\n';
        auto local_ply_filename =
            output_path(
                output_dir,
                "local_voxel_scalar",
                surface_name,
                c,
                resolution,
                level);
        local_ply_filename.replace_extension(".ply");
        iso::field_export::write_scalar_ply(
            patch.mesh,
            local_ply_filename.string(),
            evaluate_scalar,
            patch_options.isovalue,
            singularity.center,
            radius);
        std::cout << "Saved scalar PLY: "
                << local_ply_filename << '\n';
        auto local_vtk_filename = local_ply_filename;
        local_vtk_filename.replace_extension(".vtk");
        iso::field_export::write_scalar_vtk(
            patch.mesh,
            local_vtk_filename.string(),
            evaluate_scalar,
            patch_options.isovalue,
            singularity.center,
            radius);
        std::cout << "Saved scalar VTK: "
                << local_vtk_filename << '\n';
        iso::stitch::StitchResult stitched;
        bool have_stitched_result = false;
        std::size_t final_singular_point_components = 0;
        if (!no_stitch)
        {
            iso::stitch::LocalVoxelStitcher stitcher;
            iso::stitch::LocalVoxelStitcher::Options stitch_options;
            stitch_options.radius = radius;
            stitch_options.use_sphere_seam = true;
            const auto stitch_start =
                std::chrono::steady_clock::now();
            stitched =
                stitcher.stitch(
                    global,
                    patch.mesh,
                    singularity.center,
                    stitch_options);
            orient_components_to_implicit_gradient(
                stitched.mesh,
                *surface);
            stitched.topology =
                iso::stitch::compute_topology(stitched.mesh);
            final_singular_point_components =
                count_components_touching_point(
                    stitched.mesh,
                    singularity.center,
                    singular_point_tolerance);
            const auto stitch_end =
                std::chrono::steady_clock::now();
            const double stitch_ms =
                std::chrono::duration<double, std::milli>(
                    stitch_end - stitch_start)
                    .count();
            have_stitched_result = true;
            if (surface_name == "cone" && std::abs(c) <= 1e-12)
            {
                const auto global_geometry =
                    iso::geometry_check::measure(
                        global,
                        singularity.center,
                        radius,
                        4);
                const auto stitched_geometry =
                    iso::geometry_check::measure(
                        stitched.mesh,
                        singularity.center,
                        radius,
                        4);
                iso::geometry_check::print_report(
                    "GLOBAL MARCHING CUBES",
                    global_geometry);
                iso::geometry_check::print_report(
                    "STITCHED MESH",
                    stitched_geometry);
            }
            std::cout
                << "\nSTITCHED GLOBAL + LOCAL VOXEL MC\n"
                << "  stitched            : "
                << (stitched.stitched ? "YES" : "NO")
                << '\n'
                << "  stitch time         : "
                << stitch_ms << " ms\n"
                << "  global triangles removed: "
                << stitched.global_triangles_removed << '\n'
                << "  local triangles inserted: "
                << stitched.local_triangles_inserted << '\n'
                << "  seam loops           : "
                << stitched.seam_loops << '\n'
                << "  seam vertices        : "
                << stitched.seam_vertices << '\n'
                << "  bridge triangles added (direct loop welding uses none): "
                << stitched.seam_triangles << '\n'
                << "  max seam distance    : "
                << stitched.max_seam_vertex_distance << '\n'
                << "  components touching singular point: "
                << final_singular_point_components << '\n';
            print_topology(
                "  final topology:",
                stitched.topology);
            const auto stitched_filename =
                output_path(
                    output_dir,
                    "stitched_voxel",
                    surface_name,
                    c,
                    resolution,
                    level);
            save_mesh(
                stitched.mesh,
                stitched_filename);
            auto stitched_ply_filename =
            output_path(
                output_dir,
                "stitched_scalar",
                surface_name,
                c,
                resolution,
                level);
        stitched_ply_filename.replace_extension(".ply");
        iso::field_export::write_scalar_ply(
            stitched.mesh,
            stitched_ply_filename.string(),
            evaluate_scalar,
            global_options.isovalue,
            singularity.center,
            radius);
        std::cout << "Saved scalar PLY: "
                << stitched_ply_filename << '\n';
        auto stitched_vtk_filename = stitched_ply_filename;
        stitched_vtk_filename.replace_extension(".vtk");
        iso::field_export::write_scalar_vtk(
            stitched.mesh,
            stitched_vtk_filename.string(),
            evaluate_scalar,
            global_options.isovalue,
            singularity.center,
            radius);
        std::cout << "Saved scalar VTK: "
                << stitched_vtk_filename << '\n';
            const auto report_filename =
                output_dir /
                (
                    "report_" +
                    surface_name +
                    "_c_" +
                    std::to_string(c) +
                    "_r_" +
                    std::to_string(resolution) +
                    "_l_" +
                    std::to_string(level) +
                    ".txt");
            write_report(
                report_filename,
                surface_name,
                c,
                resolution,
                extent,
                radius,
                level,
                singularity,
                global_singular_point_components,
                final_singular_point_components,
                global_topology,
                patch,
                patch_topology,
                &stitched,
                global_duration);
            std::cout
                << "Report saved: "
                << report_filename << '\n';
        }
        if (no_gui)
        {
            if (have_stitched_result)
                return stitched.stitched ? 0 : 2;
            return 0;
        }
        std::filesystem::path viewer_filename;
        if (have_stitched_result &&
            stitched.stitched)
        {
            viewer_filename =
                output_path(
                    output_dir,
                    "stitched_voxel",
                    surface_name,
                    c,
                    resolution,
                    level);
        }
        else
        {
            viewer_filename =
                output_dir /
                (
                    "global_" +
                    surface_name +
                    "_c_" +
                    std::to_string(c) +
                    "_r_" +
                    std::to_string(resolution) +
                    ".obj");
            save_mesh(
                global,
                viewer_filename);
        }
        iso::MyViewer window(
            have_stitched_result &&
                    stitched.stitched
                ? "Stitched Local Voxel MC"
                : "Implicit Isosurface",
            1024,
            768);
        window.load_mesh(
            viewer_filename.string().c_str());
        window.set_marching_cubes(
            [&](double x, double y, double z)
            {
                return surface->eval(x, y, z);
            },
            bounds,
            resolution,
            resolution,
            resolution,
            0.0);
        return window.run();
    }
    catch (const std::exception& e)
    {
        std::cerr
            << "Error: "
            << e.what()
            << '\n';
        return 1;
    }
}
