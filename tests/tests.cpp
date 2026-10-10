#include "GradientHessian.hpp"
#include "ImplicitSurface.hpp"
#include "LocalUnfolder.hpp"
#include "LocalVoxelStitcher.hpp"
#include "MarchingCubes.hpp"
#include "SingularityDetector.hpp"

#include <Eigen/Eigenvalues>

#include <cmath>
#include <cstddef>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace
{

int failures = 0;

// Fail fast within each test function. This prevents a failed precondition
// (for example, an empty levels vector) from causing a later crash and
// hiding the original test failure.
#define CHECK(condition)                                                       \
    do                                                                         \
    {                                                                          \
        if (!(condition))                                                      \
        {                                                                      \
            ++failures;                                                        \
            std::cerr << __FILE__ << ':' << __LINE__                           \
                      << ": CHECK failed: " << #condition << '\n';            \
            return;                                                            \
        }                                                                      \
    } while (false)

bool near(double a, double b, double eps = 1e-6)
{
    return std::abs(a - b) <= eps;
}

class SphereSurface final : public ImplicitSurface
{
public:
    double eval(double x, double y, double z) const override
    {
        return x * x + y * y + z * z - 1.0;
    }
};

using EdgeKey = std::pair<iso::mc::Index, iso::mc::Index>;

EdgeKey edge(iso::mc::Index a, iso::mc::Index b)
{
    if (a > b)
        std::swap(a, b);
    return {a, b};
}

using std_size_t = std::size_t;

long long euler_characteristic(const iso::mc::Mesh& mesh)
{
    std::map<EdgeKey, std_size_t> edges;
    for (const auto& t : mesh.triangles)
    {
        ++edges[edge(t[0], t[1])];
        ++edges[edge(t[1], t[2])];
        ++edges[edge(t[2], t[0])];
    }

    return static_cast<long long>(mesh.vertices.size()) -
           static_cast<long long>(edges.size()) +
           static_cast<long long>(mesh.triangles.size());
}

void test_cone_classification()
{
    ParameterizedConeQuadric cone(0.0);
    SingularityDetector detector(cone, 1e-8, 1e-8);

    const auto result = detector.classifyPoint(
        pmp::Point(0.0f, 0.0f, 0.0f));

    CHECK(result.type == SingularityType::NonDegenerateSingular);
    CHECK(near(result.function_value, 0.0));
    CHECK(near(result.gradient_norm, 0.0));
    CHECK(near(result.hessian_determinant, -8.0));
    CHECK(near(result.eigenvalues[0], -2.0));
    CHECK(near(result.eigenvalues[1], 2.0));
    CHECK(near(result.eigenvalues[2], 2.0));
}

void test_hyperboloids_have_regular_surface_points()
{
    const double c_positive = 0.25;
    ParameterizedConeQuadric one_sheet(c_positive);
    SingularityDetector detector_positive(one_sheet, 1e-8, 1e-8);

    const auto positive = detector_positive.classifyPoint(
        pmp::Point(0.5f, 0.0f, 0.0f));

    CHECK(positive.type == SingularityType::Regular);
    CHECK(near(positive.function_value, 0.0));
    CHECK(positive.gradient_norm > 0.99);

    const double c_negative = -0.25;
    ParameterizedConeQuadric two_sheet(c_negative);
    SingularityDetector detector_negative(two_sheet, 1e-8, 1e-8);

    const auto negative = detector_negative.classifyPoint(
        pmp::Point(0.0f, 0.0f, 0.5f));

    CHECK(negative.type == SingularityType::Regular);
    CHECK(near(negative.function_value, 0.0));
    CHECK(negative.gradient_norm > 0.99);

    const Eigen::Matrix3d H = two_sheet.hessian(0.0, 0.0, 0.5);
    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(H);

    CHECK(solver.info() == Eigen::Success);
    const Eigen::Vector3d evals = solver.eigenvalues();
    CHECK(near(evals[0], -2.0));
    CHECK(near(evals[1], 2.0));
    CHECK(near(evals[2], 2.0));
}

void test_quartic_saddle_uses_generic_derivatives()
{
    QuarticSaddleSurface surface(0.0, 0.05);
    SingularityDetector detector(surface, 1e-8, 1e-3);

    const auto result = detector.classifyPoint(
        pmp::Point(0.0f, 0.0f, 0.0f));

    CHECK(result.type == SingularityType::NonDegenerateSingular);
    CHECK(near(result.function_value, 0.0, 1e-10));
    CHECK(result.gradient_norm < 1e-7);
    CHECK(near(result.hessian_determinant, -8.0, 1e-3));
    CHECK(near(result.eigenvalues[0], -2.0, 1e-3));
    CHECK(near(result.eigenvalues[1], 2.0, 1e-3));
    CHECK(near(result.eigenvalues[2], 2.0, 1e-3));
}

void test_gradient_hessian_no_silent_zero_fallback()
{
    SphereSurface sphere;
    GradientHessianCalculator calculator(sphere);
    const pmp::Point p(1.0f, 0.0f, 0.0f);

    const pmp::Point g = calculator.computeGradient(p);
    CHECK(std::abs(g[0] - 2.0f) < 1e-3f);
    CHECK(std::abs(g[1]) < 1e-3f);
    CHECK(std::abs(g[2]) < 1e-3f);

    const Eigen::Matrix3d H = calculator.computeHessian(p);
    CHECK(std::abs(H(0, 0) - 2.0) < 1e-3);
    CHECK(std::abs(H(1, 1) - 2.0) < 1e-3);
    CHECK(std::abs(H(2, 2) - 2.0) < 1e-3);
}

void test_closed_mc_sphere_euler_characteristic()
{
    iso::mc::Options options;
    options.nx = 25;
    options.ny = 25;
    options.nz = 25;
    options.isovalue = 0.0;

    const iso::mc::Bounds bounds{
        {-1.5, -1.5, -1.5},
        { 1.5,  1.5,  1.5}
    };

    SphereSurface sphere;
    const auto mesh = iso::mc::extract(
        [&](double x, double y, double z)
        {
            return sphere.eval(x, y, z);
        },
        bounds,
        options);

    CHECK(!mesh.vertices.empty());
    CHECK(!mesh.triangles.empty());
    CHECK(euler_characteristic(mesh) == 2);
}

void test_local_unfolder_region_and_boundaries()
{
    iso::mc::Options options;
    options.nx = 33;
    options.ny = 33;
    options.nz = 33;
    options.isovalue = 0.0;

    const iso::mc::Bounds bounds{
        {-1.5, -1.5, -1.5},
        { 1.5,  1.5,  1.5}
    };

    SphereSurface sphere;
    const auto input = iso::mc::extract(
        [&](double x, double y, double z)
        {
            return sphere.eval(x, y, z);
        },
        bounds,
        options);

    LocalUnfolder::Options local_options;
    local_options.radius = 0.55;
    local_options.levels = 1;
    local_options.region_mode =
        LocalUnfolder::RegionMode::TopologicalBfs;

    LocalUnfolder unfolder;
    // The sphere is centered at the origin but its surface is one unit away.
    // Selecting around the origin with radius 0.55 selects no surface faces.
    // Use a point on the surface so the local region is non-empty.
    const auto result = unfolder.refine(
        input,
        sphere,
        iso::mc::Point{0.0, 0.0, 1.0},
        local_options);

    CHECK(result.sphere_faces_in_region > 0);
    CHECK(result.topological_faces_in_region > 0);
    CHECK(result.faces_in_region == result.topological_faces_in_region);
    CHECK(result.vertices_in_region > 0);
    CHECK(result.interface_edges > 0);
    CHECK(result.interface_vertices > 0);
    CHECK(result.levels.size() == 1);
    // CHECK returns from this test immediately on failure, so the indexing
    // below cannot run when levels is empty.
    CHECK(result.levels[0].split_edges > 0);
    CHECK(result.levels[0].global_boundary_edges ==
          result.global_boundary_edges);
    CHECK(result.output_vertices > result.input_vertices);
    CHECK(result.output_triangles > result.input_triangles);
    CHECK(result.levels[0].projection_failures == 0);
}

void test_selection_methods_can_differ_near_cone_apex()
{
    ParameterizedConeQuadric cone(0.0);

    iso::mc::Options options;
    options.nx = 33;
    options.ny = 33;
    options.nz = 33;
    options.isovalue = 0.0;

    const iso::mc::Bounds bounds{
        {-1.0, -1.0, -1.0},
        { 1.0,  1.0,  1.0}
    };

    const auto input = iso::mc::extract(
        [&](double x, double y, double z)
        {
            return cone.eval(x, y, z);
        },
        bounds,
        options);

    LocalUnfolder::Options local_options;
    local_options.radius = 0.30;
    local_options.levels = 1;
    local_options.region_mode =
        LocalUnfolder::RegionMode::TopologicalBfs;
    // This test checks face selection, not voxel scalar statistics.
    local_options.base_resolution = 0;

    // Offset the selection center toward the upper cone by 0.75 * radius.
    // The sphere selector can select nearby faces on both disconnected
    // nappes, while TopologicalBfs starts from the nearest mesh vertex and
    // follows only the connected component containing that seed.
    const iso::mc::Point probe_center{0.0, 0.0, 0.225};

    LocalUnfolder unfolder;
    const auto result = unfolder.refine(
        input, cone, probe_center, local_options);

    CHECK(result.levels.size() == 1);
    CHECK(result.sphere_faces_in_region > 0);
    CHECK(result.topological_faces_in_region > 0);
    CHECK(result.levels[0].selection_changed);
    CHECK(result.sphere_faces_in_region >
          result.topological_faces_in_region);
}

void test_local_voxel_patch_stitches_into_global_mesh()
{
    ParameterizedConeQuadric cone(0.0);

    const iso::mc::Bounds global_bounds{
        {-1.0, -1.0, -1.0},
        { 1.0,  1.0,  1.0}
    };

    iso::mc::Options global_options;
    global_options.nx = 33;
    global_options.ny = 33;
    global_options.nz = 33;

    const auto global = iso::mc::extract(
        [&](double x, double y, double z)
        {
            return cone.eval(x, y, z);
        },
        global_bounds,
        global_options);

    // Independently extracted finer local voxel patch covering the seam radius.
    const iso::mc::Bounds local_bounds{
        {-0.6, -0.6, -0.6},
        { 0.6,  0.6,  0.6}
    };

    iso::mc::Options local_options;
    local_options.nx = 97;
    local_options.ny = 97;
    local_options.nz = 97;

    const auto local_patch = iso::mc::extract(
        [&](double x, double y, double z)
        {
            return cone.eval(x, y, z);
        },
        local_bounds,
        local_options);

    iso::stitch::LocalVoxelStitcher stitcher;
    iso::stitch::LocalVoxelStitcher::Options stitch_options;
    stitch_options.radius = 0.3;

    const auto result = stitcher.stitch(
        global,
        local_patch,
        iso::mc::Point{0.0, 0.0, 0.0},
        stitch_options);

    CHECK(result.stitched);
    CHECK(result.global_triangles_removed > 0);
    CHECK(result.local_triangles_inserted > 0);
    CHECK(result.seam_loops == 2);
    CHECK(result.seam_vertices > 0);

    // LocalVoxelStitcher welds the matched boundary loops by remapping their
    // vertices to shared indices. It deliberately creates no bridge strip,
    // so zero seam triangles is the expected result for this implementation.
    CHECK(result.seam_triangles == 0);

    CHECK(result.max_seam_vertex_distance >= 0.0);
    CHECK(result.topology.nonmanifold_edges == 0);
    CHECK(result.topology.connected_components == 2);

    // Stitching must not introduce new boundary edges. The global mesh already
    // has open boundaries at the edge of the extraction domain.
    CHECK(result.topology.boundary_edges ==
          iso::stitch::compute_topology(global).boundary_edges);
}

} // namespace

int main()
{
    test_cone_classification();
    test_hyperboloids_have_regular_surface_points();
    test_quartic_saddle_uses_generic_derivatives();
    test_gradient_hessian_no_silent_zero_fallback();
    test_closed_mc_sphere_euler_characteristic();
    test_local_unfolder_region_and_boundaries();
    test_selection_methods_can_differ_near_cone_apex();
    test_local_voxel_patch_stitches_into_global_mesh();

    if (failures != 0)
    {
        std::cerr << "\nTests failed: " << failures << '\n';
        return 1;
    }

    std::cout << "All tests passed.\n";
    return 0;
}
