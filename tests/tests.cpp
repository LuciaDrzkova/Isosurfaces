// Minimal dependency-free test runner: each CHECK failure is reported and the
// process exits non-zero at the end.

#include "method/ImplicitSurface.h"
#include "method/Job.h"
#include "method/MeshMetrics.h"
#include "method/Remesher.h"
#include "app/Settings.h"

#include <pmp/algorithms/shapes.h>
#include <pmp/io/io.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

using namespace std;

namespace {

int failures = 0;

#define CHECK(cond)                                                          \
    do                                                                       \
    {                                                                        \
        if (!(cond))                                                         \
        {                                                                    \
            ++failures;                                                      \
            cerr << __FILE__ << ':' << __LINE__ << ": CHECK failed: "   \
                      << #cond << '\n';                                      \
        }                                                                    \
    } while (0)

using P = pmp::Point;

bool near(double a, double b, double eps = 1e-4)
{
    return abs(a - b) <= eps;
}

void test_median()
{
    CHECK(iso::median({}) == 0.0f);
    CHECK(iso::median({5.0f}) == 5.0f);
    CHECK(iso::median({3.0f, 1.0f, 2.0f}) == 2.0f);
    CHECK(iso::median({4.0f, 1.0f, 3.0f, 2.0f}) == 2.5f);
}

void test_triangle_angles()
{
    const auto eq = iso::triangle_angles(P(0, 0, 0), P(1, 0, 0),
                                         P(0.5f, sqrt(3.0f) / 2, 0));
    for (float a : eq)
        CHECK(near(a, 60.0, 1e-2));

    const auto right = iso::triangle_angles(P(0, 0, 0), P(1, 0, 0), P(0, 1, 0));
    CHECK(near(right[0], 90.0, 1e-2));

    // Coincident points must not produce NaN
    for (float a : iso::triangle_angles(P(1, 1, 1), P(1, 1, 1), P(2, 1, 1)))
        CHECK(!isnan(a));
}

void test_gradient()
{
    // ghost: f = x^3 + y^2 + z^2 - 0.5  =>  grad f = (3x^2, 2y, 2z)
    const auto* ghost = iso::find_surface("ghost");
    CHECK(ghost != nullptr);
    const pmp::Point p(0.7f, -1.2f, 0.4f);
    const auto s = ghost->evaluate(p);
    CHECK(near(s.value, 0.343 + 1.44 + 0.16 - 0.5));
    CHECK(near(s.gradient[0], 3 * 0.49));
    CHECK(near(s.gradient[1], -2.4));
    CHECK(near(s.gradient[2], 0.8));

    CHECK(iso::find_surface("no-such-surface") == nullptr);
}

void test_projection()
{
    const auto* sphere = iso::find_surface("sphere"); // radius sqrt(5)
    CHECK(sphere != nullptr);
    const pmp::Point q = sphere->project(pmp::Point(2.0f, 0.5f, 1.0f));
    CHECK(abs(sphere->value(q)) < 1e-4);

    // Already on the surface: unchanged
    const pmp::Point on(sqrt(5.0f), 0.0f, 0.0f);
    CHECK(pmp::distance(sphere->project(on), on) < 1e-4);
}

void test_remesh_sphere()
{
    const auto* sphere = iso::find_surface("sphere");
    pmp::SurfaceMesh mesh = pmp::icosphere(4);
    iso::project_onto_surface(mesh, *sphere);

    const size_t before = mesh.n_vertices();
    const auto m0 = iso::compute_mesh_metrics(mesh, *sphere, 25.0f);
    CHECK(m0.vertices_off_surface == 0);

    // A divisor below 1 makes (almost) every triangle a collapse candidate
    iso::RemeshOptions options;
    options.area_divisor = 0.5f;
    const auto result = iso::remesh(mesh, *sphere, options);

    CHECK(result.collapses_per_iteration.size() == 6);
    CHECK(result.collapses_per_iteration[0] > 0);
    CHECK(mesh.n_vertices() < before);
    CHECK(mesh.n_faces() > 0);
    CHECK(mesh.is_triangle_mesh());

    // Collapsed vertices are re-projected, so the result is still on the surface
    const auto m1 = iso::compute_mesh_metrics(mesh, *sphere, 25.0f);
    CHECK(m1.max_distance < 1e-3);
    CHECK(m1.n_vertices == mesh.n_vertices());
}

void test_job_report()
{
    const auto dir = filesystem::temp_directory_path() / "iso_test_job";
    filesystem::create_directories(dir);

    pmp::write(pmp::icosphere(3), dir / "in.obj"); // radius 1, not on f = 0

    iso::JobSettings settings;
    settings.input = (dir / "in.obj").string();
    settings.output = (dir / "out.obj").string();
    settings.function = "sphere";
    settings.remesh.area_divisor = 0.8f;
    const auto report = iso::run_job(settings);

    CHECK(report.input.max_distance > 1.0); // radius 1 vs sqrt(5)
    CHECK(report.output.n_vertices < report.input.n_vertices);
    CHECK(filesystem::exists(settings.output));

    // "Projected" is the input mesh after the first projection: the same mesh
    // as the input, but on the surface. The output stays on the surface too.
    CHECK(report.projected.n_vertices == report.input.n_vertices);
    CHECK(report.projected.n_faces == report.input.n_faces);
    CHECK(report.projected.max_distance < 1e-5);
    CHECK(report.output.max_distance < 1e-5);

    filesystem::remove_all(dir);
}

void test_remesh_rejects_bad_input()
{
    const auto* sphere = iso::find_surface("sphere");

    pmp::SurfaceMesh quads = pmp::quad_sphere(1);
    bool threw = false;
    try
    {
        iso::remesh(quads, *sphere);
    }
    catch (const invalid_argument&)
    {
        threw = true;
    }
    CHECK(threw);

    // Empty mesh is fine
    pmp::SurfaceMesh empty;
    iso::remesh(empty, *sphere);
    const auto m = iso::compute_mesh_metrics(empty, *sphere, 25.0f);
    CHECK(m.n_faces == 0 && m.avg_distance == 0.0);
}

void test_settings_round_trip()
{
    const auto file = filesystem::temp_directory_path() / "iso_test_settings" /
                      "settings.ini";
    filesystem::remove_all(file.parent_path());

    // No file yet: defaults
    iso::JobSettings defaults = iso::load_settings(file);
    CHECK(defaults.input.empty() && defaults.function.empty());
    CHECK(defaults.remesh.iterations == iso::RemeshOptions{}.iterations);

    iso::JobSettings s;
    s.input = "/some dir/mesh with spaces.obj";
    s.output = "result.obj";
    s.function = "gyroid";
    s.experiment = "singularities";
    s.remesh.iterations = 3;
    s.remesh.min_angle = 30.5f;
    s.remesh.area_divisor = 2.25f;
    CHECK(iso::save_settings(file, s)); // also creates the directory

    const auto loaded = iso::load_settings(file);
    CHECK(loaded.input == s.input);
    CHECK(loaded.output == s.output);
    CHECK(loaded.function == "gyroid");
    CHECK(loaded.experiment == "singularities");
    CHECK(loaded.remesh.iterations == 3);
    CHECK(near(loaded.remesh.min_angle, 30.5));
    CHECK(near(loaded.remesh.area_divisor, 2.25));

    // Garbage, unknown names and out-of-range values are ignored
    {
        ofstream out(file);
        out << "function=no-such-surface\nexperiment=nope\niterations=abc\n"
               "area_divisor=-4\nmin_angle=999\nnot a line\ninput=kept.obj\n";
    }
    const auto bad = iso::load_settings(file);
    CHECK(bad.function.empty() && bad.experiment.empty());
    CHECK(bad.remesh.iterations == iso::RemeshOptions{}.iterations);
    CHECK(near(bad.remesh.area_divisor, iso::RemeshOptions{}.area_divisor));
    CHECK(near(bad.remesh.min_angle, iso::RemeshOptions{}.min_angle));
    CHECK(bad.input == "kept.obj");

    // Empty path (no home directory): saving is a harmless no-op
    CHECK(!iso::save_settings({}, s));

    filesystem::remove_all(file.parent_path());
}

} // namespace

int main()
{
    test_median();
    test_triangle_angles();
    test_gradient();
    test_projection();
    test_remesh_sphere();
    test_job_report();
    test_remesh_rejects_bad_input();
    test_settings_round_trip();

    if (failures)
    {
        cerr << failures << " check(s) failed\n";
        return 1;
    }
    cout << "All tests passed\n";
    return 0;
}
