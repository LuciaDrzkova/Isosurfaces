#include "GradientHessian.hpp"
#include "ImplicitSurface.hpp"
#include "LocalMarchingCubes.hpp"
#include "LocalUnfolder.hpp"
#include "MarchingCubes.hpp"
#include "SingularityDetector.hpp"

#include <Eigen/Eigenvalues>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <map>
#include <utility>

namespace
{

int failures = 0;

#define CHECK(condition)                                                     \
    do                                                                       \
    {                                                                        \
        if (!(condition))                                                    \
        {                                                                    \
            ++failures;                                                      \
            std::cerr << __FILE__ << ':' << __LINE__                         \
                      << ": CHECK failed: " << #condition << '\n';           \
        }                                                                    \
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

long long euler_characteristic(const iso::mc::Mesh& mesh)
{
    std::map<EdgeKey, std::size_t> edges;

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
    ParameterizedConeQuadric one_sheet(0.25);
    SingularityDetector detector_positive(one_sheet, 1e-8, 1e-8);

    const auto positive = detector_positive.classifyPoint(
        pmp::Point(0.5f, 0.0f, 0.0f));

    CHECK(positive.type == SingularityType::Regular);
    CHECK(near(positive.function_value, 0.0));
    CHECK(positive.gradient_norm > 0.99);

    ParameterizedConeQuadric two_sheet(-0.25);
    SingularityDetector detector_negative(two_sheet, 1e-8, 1e-8);

    const auto negative = detector_negative.classifyPoint(
        pmp::Point(0.0f, 0.0f, 0.5f));

    CHECK(negative.type == SingularityType::Regular);
    CHECK(near(negative.function_value, 0.0));
    CHECK(negative.gradient_norm > 0.99);

    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(
        two_sheet.hessian(0.0, 0.0, 0.5));

    CHECK(solver.info() == Eigen::Success);

    const Eigen::Vector3d evals = solver.eigenvalues();

    CHECK(near(evals[0], -2.0));
    CHECK(near(evals[1], 2.0));
    CHECK(near(evals[2], 2.0));
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

    CHECK(std::abs(H(0, 1)) < 1e-3);
    CHECK(std::abs(H(0, 2)) < 1e-3);
    CHECK(std::abs(H(1, 0)) < 1e-3);
    CHECK(std::abs(H(1, 2)) < 1e-3);
    CHECK(std::abs(H(2, 0)) < 1e-3);
    CHECK(std::abs(H(2, 1)) < 1e-3);
}

void test_closed_mc_sphere_euler_characteristic()
{
    iso::mc::Options options;
    options.nx = 25;
    options.ny = 25;
    options.nz = 25;

    const iso::mc::Bounds bounds{
        {-1.5, -1.5, -1.5},
        {1.5, 1.5, 1.5}
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

    const iso::mc::Bounds bounds{
        {-1.5, -1.5, -1.5},
        {1.5, 1.5, 1.5}
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
    local_options.sampling_bounds = bounds;
    local_options.base_resolution = 33;
    local_options.isovalue = 0.0;

    LocalUnfolder unfolder;

    // IMPORTANT:
    // (0,0,0) is the center of the sphere and is not on the surface.
    // LocalUnfolder expects a point on/near the surface.
    const iso::mc::Point surface_point{1.0, 0.0, 0.0};

    const auto result = unfolder.refine(
        input,
        sphere,
        surface_point,
        local_options);

    CHECK(result.vertices_in_region > 0);
    CHECK(result.faces_in_region > 0);
    CHECK(result.levels.size() == 1);
    CHECK(result.boundary_edges > 0);
    CHECK(result.boundary_vertices > 0);

    // Do not access result.levels[0] if region selection failed.
    if (!result.levels.empty())
    {
        CHECK(result.levels[0].selected_faces > 0);
        CHECK(result.levels[0].new_vertices > 0);
        CHECK(result.levels[0].new_triangles > 0);
        CHECK(result.levels[0].projection_failures == 0);
    }
}

void test_local_voxel_marching_cubes_patch()
{
    SphereSurface sphere;

    const iso::mc::Bounds bounds{
        {-1.5, -1.5, -1.5},
        {1.5, 1.5, 1.5}
    };

    LocalMarchingCubes::Options options;
    options.radius = 0.55;
    options.level = 1;
    options.padding_cells = 1;
    options.sampling_bounds = bounds;
    options.base_resolution = 33;
    options.isovalue = 0.0;

    LocalMarchingCubes extractor;

    // IMPORTANT:
    // The local MC patch is centered around a point on the surface,
    // not at the center of the sphere.
    const iso::mc::Point surface_point{1.0, 0.0, 0.0};

    const auto result = extractor.extract(
        sphere,
        surface_point,
        options);

    CHECK(!result.mesh.vertices.empty());
    CHECK(!result.mesh.triangles.empty());

    CHECK(result.refinement_factor == 2);
    CHECK(result.refined_spacing < result.base_spacing);

    CHECK(result.resolution_x > result.base_cells_x);
    CHECK(result.resolution_y > result.base_cells_y);
    CHECK(result.resolution_z > result.base_cells_z);

    CHECK(result.boundary_edges > 0);
    CHECK(result.boundary_vertices > 0);

    CHECK(result.nonmanifold_edges == 0);
    CHECK(result.selected_cells > 0);
    CHECK(result.extraction_time_ms >= 0.0);
}

} // namespace

int main()
{
    test_cone_classification();
    test_hyperboloids_have_regular_surface_points();
    test_gradient_hessian_no_silent_zero_fallback();
    test_closed_mc_sphere_euler_characteristic();
    test_local_unfolder_region_and_boundaries();
    test_local_voxel_marching_cubes_patch();

    if (failures == 0)
    {
        std::cout << "All tests passed.\n";
        return 0;
    }

    std::cerr << failures << " test(s) failed.\n";
    return 1;
}