#include "LocalVoxelStitcher.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace
{
int failures = 0;
#define CHECK(expr) do { if (!(expr)) { std::cerr << "FAIL: " #expr " at " << __FILE__ << ':' << __LINE__ << '\n'; ++failures; } } while (false)

// Creates two separately edge-connected cone nappes. This deterministic mesh
// fixture tests seam stitching without relying on MC's exact triangulation.
iso::mc::Mesh make_double_cone(int angular_samples, int vertical_layers, double z_extent)
{
    iso::mc::Mesh mesh;
    for (int sign : {-1, 1})
    {
        const iso::mc::Index apex = static_cast<iso::mc::Index>(mesh.vertices.size());
        mesh.vertices.push_back({0.0, 0.0, 0.0});
        std::vector<std::vector<iso::mc::Index>> rings;
        rings.reserve(static_cast<std::size_t>(vertical_layers));
        for (int layer = 1; layer <= vertical_layers; ++layer)
        {
            const double z = sign * z_extent * static_cast<double>(layer) / vertical_layers;
            std::vector<iso::mc::Index> ring;
            ring.reserve(static_cast<std::size_t>(angular_samples));
            for (int j = 0; j < angular_samples; ++j)
            {
                const double theta = 2.0 * 3.14159265358979323846 * j / angular_samples;
                ring.push_back(static_cast<iso::mc::Index>(mesh.vertices.size()));
                mesh.vertices.push_back({std::abs(z)*std::cos(theta),
                                         std::abs(z)*std::sin(theta), z});
            }
            rings.push_back(std::move(ring));
        }
        for (int j = 0; j < angular_samples; ++j)
        {
            const int next = (j + 1) % angular_samples;
            mesh.triangles.push_back({apex, rings[0][static_cast<std::size_t>(j)],
                                      rings[0][static_cast<std::size_t>(next)]});
        }
        for (int layer = 0; layer + 1 < vertical_layers; ++layer)
        {
            for (int j = 0; j < angular_samples; ++j)
            {
                const int next = (j + 1) % angular_samples;
                const auto a = rings[static_cast<std::size_t>(layer)][static_cast<std::size_t>(j)];
                const auto b = rings[static_cast<std::size_t>(layer)][static_cast<std::size_t>(next)];
                const auto c = rings[static_cast<std::size_t>(layer + 1)][static_cast<std::size_t>(next)];
                const auto d = rings[static_cast<std::size_t>(layer + 1)][static_cast<std::size_t>(j)];
                mesh.triangles.push_back({a,b,c});
                mesh.triangles.push_back({a,c,d});
            }
        }
    }
    return mesh;
}

void test_boundary_edge_split_preserves_valid_triangle_indices()
{
    iso::mc::Mesh mesh;
    mesh.vertices = {{0,0,0}, {1,0,0}, {0,1,0}};
    mesh.triangles = {{{0,1,2}}};
    std::vector<iso::mc::Index> loop{0,1,2};
    iso::stitch::detail::split_boundary_edge(mesh, loop, 0, 0.25);
    CHECK(mesh.vertices.size() == 4);
    CHECK(mesh.triangles.size() == 2);
    CHECK(loop.size() == 4);
    CHECK(std::abs(mesh.vertices[3][0] - 0.25) < 1e-12);
    CHECK(std::abs(mesh.vertices[3][1]) < 1e-12);
    for (const auto& tri : mesh.triangles)
        for (auto index : tri) CHECK(index < mesh.vertices.size());
    CHECK(iso::stitch::compute_topology(mesh).nonmanifold_edges == 0);
}

void test_different_density_cone_loops_stitch_without_new_boundaries()
{
    const auto global = make_double_cone(32, 32, 1.0);
    const auto local = make_double_cone(96, 96, 0.6);
    iso::stitch::LocalVoxelStitcher stitcher;
    iso::stitch::LocalVoxelStitcher::Options options;
    options.radius = 0.3;
    // Exercise the radius-relative default threshold (0.30 * radius).
    const auto result = stitcher.stitch(global, local, {0,0,0}, options);
    if (!result.stitched)
        std::cerr << "Stitch diagnostic: " << result.failure_reason << '\n';
    CHECK(result.stitched);
    CHECK(result.local_triangles_inserted > 0);
    CHECK(result.seam_loops == 2);
    CHECK(result.seam_vertices > 0);
    CHECK(result.seam_triangles == 0);
    const double allowed_gap = options.max_stitch_distance > 0.0
        ? options.max_stitch_distance : 0.30 * options.radius;
    CHECK(result.max_seam_vertex_distance <= allowed_gap);
    CHECK(result.topology.nonmanifold_edges == 0);
    CHECK(result.topology.boundary_edges == iso::stitch::compute_topology(global).boundary_edges);
    CHECK(result.topology.connected_components == 2);
}

void test_too_small_gap_threshold_fails_safely()
{
    const auto global = make_double_cone(32, 32, 1.0);
    const auto local = make_double_cone(96, 96, 0.6);
    iso::stitch::LocalVoxelStitcher stitcher;
    iso::stitch::LocalVoxelStitcher::Options options;
    options.radius = 0.3;
    options.max_stitch_distance = 1e-9;
    const auto result = stitcher.stitch(global, local, {0,0,0}, options);
    CHECK(!result.stitched);
    CHECK(!result.failure_reason.empty());
    CHECK(result.mesh.vertices.size() == global.vertices.size());
    CHECK(result.mesh.triangles.size() == global.triangles.size());
}
}

int main()
{
    test_boundary_edge_split_preserves_valid_triangle_indices();
    test_different_density_cone_loops_stitch_without_new_boundaries();
    test_too_small_gap_threshold_fails_safely();
    if (failures != 0)
    {
        std::cerr << "Stitching regression tests failed: " << failures << '\n';
        return EXIT_FAILURE;
    }
    std::cout << "All stitching regression tests passed.\n";
    return EXIT_SUCCESS;
}
