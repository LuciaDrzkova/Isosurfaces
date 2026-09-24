#include "app/MyViewer.h"

// Needed for Hessian-based singularity classification. Deliberately using
// FORWARD mode's second-order type (dual2nd) rather than reverse mode
// (var): mixing forward's wrt/at/derivative with reverse's own wrt/at/
// derivative overloads in the same translation unit causes ambiguous
// overload resolution. The Hessian entries are built from six calls to
// the SAME derivative()/wrt()/at() combo already used for the gradient
// below (dual2nd supports mixed second partials via wrt(x, y)), rather
// than a hessian() convenience wrapper, to avoid depending on an autodiff
// Eigen-support API whose exact signature varies between versions.
#include <autodiff/forward/dual.hpp>
#include <Eigen/Dense>
#include <fstream>
#include <pmp/io/io.h>
#include <pmp/surface_mesh.h>
#include <pmp/algorithms/differential_geometry.h>

#include <autodiff/forward/dual.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <numeric>
#include <string>
#include <unordered_set>
#include <vector>

using namespace pmp;
using namespace autodiff;

// Allow user to specify
int num_iterations = 6; // You can modify this number or pass it as a command-line argument
float angle_treshold = 25.0f; // Angle in which triangles it should collapse
float divisor_median = 1.5f;  // divides the median_volume

// --- Singularity classification parameters ---
// Percentile of the per-mesh gradient-norm distribution used to flag
// "near-singular" vertices. 0.01 = bottom 1% of |grad f| on this mesh.
float singularity_grad_percentile = 0.01f;
// Relative eigenvalue tolerance used to decide Hessian rank-deficiency.
// An eigenvalue is treated as "zero" if |eig| < rel_hessian_tol * max|eig|.
double rel_hessian_tol = 1e-3;

// Required for unordered_set to work with pmp::Edge and pmp::Vertex
namespace std {
template <>
struct hash<pmp::Edge>
{
    size_t operator()(const pmp::Edge &e) const
    {
        return std::hash<int>()(e.idx()); // Hash by Edge index
    }
};

template <>
struct hash<pmp::Vertex>
{
    size_t operator()(const pmp::Vertex &v) const
    {
        return std::hash<int>()(v.idx()); // Hash by Vertex index
    }
};
} // namespace std

// Function to compute the angles of a triangle
std::vector<float> compute_triangle_angles(const Point &p0, const Point &p1,
                                           const Point &p2)
{
    std::vector<float> angles;
    float a = distance(p1, p2);
    float b = distance(p0, p2);
    float c = distance(p0, p1);

    float alpha = acos(std::clamp((b * b + c * c - a * a) / (2.0f * b * c),
                                  -1.0f, 1.0f)) *
                  180.0f / M_PI;
    float beta = acos(std::clamp((a * a + c * c - b * b) / (2.0f * a * c),
                                 -1.0f, 1.0f)) *
                 180.0f / M_PI;
    float gamma = 180.0f - alpha - beta;

    angles.push_back(alpha);
    angles.push_back(beta);
    angles.push_back(gamma);

    return angles;
}

std::vector<Vertex> getRingVertices(const SurfaceMesh &mesh, Vertex old_v)
{
    std::vector<Vertex> vertices;
    for (auto he : mesh.halfedges(old_v))
    {
        Vertex v = mesh.to_vertex(he);
        if (v != old_v)
            vertices.push_back(v);
    }
    return vertices;
}

Point getCentroid(const SurfaceMesh &mesh, const std::vector<Vertex> &vertices)
{
    Point centroid(0.0f, 0.0f, 0.0f);
    for (const auto &v : vertices)
    {
        centroid += mesh.position(v);
    }
    centroid /= static_cast<float>(vertices.size());
    return centroid;
}

void print_vector(const std::vector<Vertex> &vec)
{
    std::cout << "[ ";
    for (const auto &element : vec)
    {
        std::cout << element << " ";
    }
    std::cout << "]\n";
}

// Function to compute median of a vector of floats
float compute_median(std::vector<float> &values)
{
    std::nth_element(values.begin(), values.begin() + values.size() / 2,
                     values.end());
    if (values.size() % 2 == 0)
    {
        return (values[values.size() / 2 - 1] + values[values.size() / 2]) /
               2.0f;
    }
    else
    {
        return values[values.size() / 2];
    }
}

// Implicit function in the style of autodiff example.
// Templated so the SAME formula can be evaluated with autodiff::dual
// (first derivatives, used by the existing gradient/projection code)
// and with autodiff::dual2nd (used below for the Hessian / singularity
// classification). Keep exactly one active `return` line, as before.
template <typename T>
T implicitFunctionImpl(T x, T y, T z)
{
    //return x * x * x + y * y + z * z - 0.5; //ghost
    //return x * x + y * y + z * z - 5.0; //sphere
    //return sin(x) * cos(y) + sin(y) * cos(z) + sin(z) * cos(x) + 1; // Gyroid_surf.obj
    //return x * x * x - 3 * x * y * y + z * z - 1; // Julia_surface.obj
    //return sin(3 * x) + sin(3 * y) - z - 1; // sinusoidal_surf
    //return x * x + y * y + (z - 1) * (z - 1) -pow((0.7 + 0.3 * sin(5 * atan2(y, x))),2.0); // turbine.obj
    //return x * x + y * y + z * z + 0.1 * cos(4 * x) + 0.1 * cos(4 * y) + 0.1 * cos(4 * z) - 5; // weird_blobby_shape.obj
    return -(cos(x + y) * cos(y + z) * cos(z + x) - 0.5); // wavy_network
    //return abs(sin(x) + sin(y) + sin(z)) - 1.2 + 0.5 * cos(3 * x * y * z); // Frozen_Lattice_Shard
    //return pow((x *x + y *y + z *z - 1) , 2) - 0.3 * cos(6 * atan2(y, x)) * exp(-z *z); // Lotus_Form
    //return x *x + y *y + z *z + 0.81 * sin(6 * (x + y + z)) - 1; // wavz_shell
    //return x * x + y * y - z * z; // cone: singular ONLY at origin (0,0,0) - use this to sanity-check the singularity detector
    //moebius_3D_surface
}

// Thin wrapper preserving the original signature/call sites (used by
// computeGradient / projectOntoImplicitSurface below, unchanged).
autodiff::dual implicitFunction(autodiff::dual x, autodiff::dual y,
                                autodiff::dual z)
{
    return implicitFunctionImpl(x, y, z);
}

// Second wrapper, same formula, for forward-mode second derivatives
// (used only by the Hessian / singularity classification code below).
autodiff::dual2nd implicitFunctionHessian(autodiff::dual2nd x,
                                          autodiff::dual2nd y,
                                          autodiff::dual2nd z)
{
    return implicitFunctionImpl(x, y, z);
}

// Compute gradient using autodiff's derivative tools
Point computeGradient(const Point &p)
{
    using namespace autodiff;
    autodiff::dual x = p[0], y = p[1], z = p[2];

    double dx = derivative(implicitFunction, wrt(x), at(x, y, z));
    double dy = derivative(implicitFunction, wrt(y), at(x, y, z));
    double dz = derivative(implicitFunction, wrt(z), at(x, y, z));

    return Point(static_cast<float>(dx), static_cast<float>(dy),
                 static_cast<float>(dz));
}

// Newton-Raphson projection using the autodiff gradient
Point projectOntoImplicitSurface(Point p, int max_iterations = 10,
                                 float tolerance = 1e-6f)
{
    for (int i = 0; i < max_iterations; ++i)
    {
        autodiff::dual x = p[0], y = p[1], z = p[2];
        autodiff::dual fx = implicitFunction(x, y, z);
        float fx_val = val(fx);

        Point grad = computeGradient(p);
        float grad_norm2 = dot(grad, grad);

        if (grad_norm2 < 1e-10f)
            break;

        // Newton-Raphson update
        p = p - (fx_val / grad_norm2) * grad;

        if (std::abs(fx_val) < tolerance)
            break;
    }

    return p;
}

// Classification labels for a vertex, based on local behaviour of f.
enum class SingularityClass
{
    Regular = 0,           // |grad f| above the per-mesh threshold
    NonDegenerateSingular, // |grad f| ~ 0, Hessian full rank (Morse-type: cone/node-like)
    DegenerateSingular     // |grad f| ~ 0, Hessian rank-deficient (needs blowup/unfolding treatment)
};

// Returns the value at the given percentile (0.01 = bottom 1%) of `values`.
// Generic helper: does NOT recompute anything, just sorts and indexes into
// whatever distribution the caller already has (used below on a
// once-computed gradient-norm distribution, to avoid a second pass over
// every vertex).
float compute_percentile(std::vector<float> values, float percentile)
{
    if (values.empty())
        return 0.0f;
    std::sort(values.begin(), values.end());
    size_t idx = static_cast<size_t>(percentile * values.size());
    idx = std::min(idx, values.size() - 1);
    return values[idx];
}

// Classifies a single point once it has already been flagged as
// near-singular (|grad f| below the per-mesh threshold): builds the 3x3
// Hessian from six second-order derivative() calls (dual2nd supports pure
// second partials via wrt(x, x) and mixed partials via wrt(x, y)), then
// checks its rank relative to its own largest eigenvalue magnitude, for
// the same scale-independence reason as the gradient threshold above.
SingularityClass classifyFlaggedPoint(const Point &p)
{
    autodiff::dual2nd x = p[0], y = p[1], z = p[2];

    double Hxx = derivative(implicitFunctionHessian, wrt(x, x), at(x, y, z));
    double Hyy = derivative(implicitFunctionHessian, wrt(y, y), at(x, y, z));
    double Hzz = derivative(implicitFunctionHessian, wrt(z, z), at(x, y, z));
    double Hxy = derivative(implicitFunctionHessian, wrt(x, y), at(x, y, z));
    double Hxz = derivative(implicitFunctionHessian, wrt(x, z), at(x, y, z));
    double Hyz = derivative(implicitFunctionHessian, wrt(y, z), at(x, y, z));

    Eigen::Matrix3d H;
    H << Hxx, Hxy, Hxz,
         Hxy, Hyy, Hyz,
         Hxz, Hyz, Hzz;

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(H);
    Eigen::Vector3d eig = solver.eigenvalues(); // ascending order

    double max_abs_eig = eig.cwiseAbs().maxCoeff();
    if (max_abs_eig <= 0.0)
        return SingularityClass::DegenerateSingular; // Hessian itself is ~0

    for (int i = 0; i < 3; ++i)
    {
        if (std::abs(eig[i]) < rel_hessian_tol * max_abs_eig)
            return SingularityClass::DegenerateSingular;
    }

    return SingularityClass::NonDegenerateSingular;
}

// Writes a categorically colored mesh as ASCII PLY directly, WITHOUT going
// through pmp::write()'s property-based IO: in this environment, pmp's
// writer hits an internal assertion/debug-break when the mesh carries a
// "v:color" property. This is a small, standard ASCII PLY writer with no
// pmp involvement in the color path, and colors are passed in as a plain
// std::vector (NOT stored as a pmp mesh property), so nothing here can
// affect the later pmp::write() of the actual output .obj mesh.
void write_colored_ply(const SurfaceMesh &mesh,
                       const std::vector<Color> &vertex_colors,
                       const std::string &filename)
{
    auto points = mesh.get_vertex_property<Point>("v:point");

    std::ofstream out(filename);
    if (!out.is_open())
    {
        std::cerr << "Failed to open " << filename << " for writing."
                  << std::endl;
        return;
    }

    out << "ply\n";
    out << "format ascii 1.0\n";
    out << "element vertex " << mesh.n_vertices() << "\n";
    out << "property float x\n";
    out << "property float y\n";
    out << "property float z\n";
    out << "property uchar red\n";
    out << "property uchar green\n";
    out << "property uchar blue\n";
    out << "element face " << mesh.n_faces() << "\n";
    out << "property list uchar int vertex_indices\n";
    out << "end_header\n";

    for (auto v : mesh.vertices())
    {
        Point p = points[v];
        Color c = vertex_colors[v.idx()];
        int r = static_cast<int>(std::clamp(c[0], 0.0f, 1.0f) * 255.0f);
        int g = static_cast<int>(std::clamp(c[1], 0.0f, 1.0f) * 255.0f);
        int b = static_cast<int>(std::clamp(c[2], 0.0f, 1.0f) * 255.0f);
        out << p[0] << " " << p[1] << " " << p[2] << " " << r << " " << g
            << " " << b << "\n";
    }

    for (auto f : mesh.faces())
    {
        std::vector<int> idx;
        for (auto v : mesh.vertices(f))
            idx.push_back(v.idx());

        out << idx.size();
        for (int i : idx)
            out << " " << i;
        out << "\n";
    }

    out.close();
}

// Runs the classification over every vertex of `mesh`, colors the vertices
// categorically (gray = regular, yellow = non-degenerate singular,
// red = degenerate singular), writes a colored PLY for visualization, and
// prints the statistics block. Diagnostic only: does not modify geometry
// or connectivity. Colors are kept in a plain std::vector, NOT added as a
// pmp mesh property, so the mesh object itself is untouched (see
// write_colored_ply for why).
void classify_and_export_singularities(SurfaceMesh &mesh,
                                       const std::string &label,
                                       const std::string &ply_filename)
{
    std::cout << "Classifying singularities (" << label << ")..." << std::flush;
    auto t0 = std::chrono::high_resolution_clock::now();

    auto points = mesh.get_vertex_property<Point>("v:point");

    // Single pass over all vertices: compute each vertex's gradient once
    // and cache it, instead of computing it once for the threshold and
    // again for the classification loop (that redundancy roughly doubled
    // the cost of this step on large meshes).
    std::vector<Vertex> verts;
    std::vector<float> grad_norms;
    verts.reserve(mesh.n_vertices());
    grad_norms.reserve(mesh.n_vertices());

    for (auto v : mesh.vertices())
    {
        Point grad = computeGradient(points[v]);
        verts.push_back(v);
        grad_norms.push_back(std::sqrt(dot(grad, grad)));
    }

    float grad_threshold =
        compute_percentile(grad_norms, singularity_grad_percentile);

    std::vector<Color> vertex_colors(mesh.n_vertices());

    const Color color_regular(0.6f, 0.6f, 0.6f);
    const Color color_nondegenerate(1.0f, 0.85f, 0.0f);
    const Color color_degenerate(1.0f, 0.0f, 0.0f);

    int count_regular = 0, count_nondegenerate = 0, count_degenerate = 0;

    for (size_t i = 0; i < verts.size(); ++i)
    {
        Vertex v = verts[i];

        if (grad_norms[i] >= grad_threshold)
        {
            vertex_colors[v.idx()] = color_regular;
            count_regular++;
            continue;
        }

        // Only the flagged (near-singular) vertices pay for a Hessian.
        SingularityClass cls = classifyFlaggedPoint(points[v]);
        if (cls == SingularityClass::NonDegenerateSingular)
        {
            vertex_colors[v.idx()] = color_nondegenerate;
            count_nondegenerate++;
        }
        else
        {
            vertex_colors[v.idx()] = color_degenerate;
            count_degenerate++;
        }
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << " done (" << std::chrono::duration<double>(t1 - t0).count()
              << "s)" << std::endl;

    std::cout << "=== SINGULARITY CLASSIFICATION (" << label << ") ==="
              << std::endl;
    std::cout << "Gradient threshold (per-mesh, "
              << (singularity_grad_percentile * 100.0f)
              << "th percentile): " << grad_threshold << std::endl;
    std::cout << "Regular vertices: " << count_regular << std::endl;
    std::cout << "Non-degenerate singular vertices (Hessian full rank): "
              << count_nondegenerate << std::endl;
    std::cout << "Degenerate singular vertices (Hessian rank-deficient): "
              << count_degenerate << std::endl;

    std::cout << "Writing " << ply_filename << " ..." << std::flush;
    auto t2 = std::chrono::high_resolution_clock::now();
    write_colored_ply(mesh, vertex_colors, ply_filename);
    auto t3 = std::chrono::high_resolution_clock::now();
    std::cout << " done (" << std::chrono::duration<double>(t3 - t2).count()
              << "s)" << std::endl;
}

void compute_mesh_metrics(const SurfaceMesh &mesh, const std::string &label)
{
    auto t0 = std::chrono::high_resolution_clock::now();
    int triangle_count = 0;
    int triangles_less_than_median = 0;
    int vertex_count = 0;
    std::vector<float> angles_less_than_threshold;
    std::vector<float> volumes;
    float sum_squared_deviations = 0.0f;

    std::vector<float> distances_from_surface;
    int vertices_off_surface = 0;
    const float surface_tolerance = 1e-4f;

    // Precompute valences
    std::vector<int> vertex_valences(mesh.n_vertices(), 0);

    vertex_count = mesh.n_vertices(); // Count total number of vertices

    // Compute distances from surface and valences
    auto points = mesh.get_vertex_property<Point>("v:point");
    for (auto v : mesh.vertices())
    {
        Point p = points[v];
        autodiff::dual x = p[0], y = p[1], z = p[2];
        autodiff::dual fx = implicitFunction(x, y, z);
        float dist = std::abs(val(fx)); // Unsigned distance
        distances_from_surface.push_back(dist);

        if (dist > surface_tolerance)
        {
            vertices_off_surface++;
        }

        // Count valence (number of adjacent vertices)
        int idx = v.idx(); // get integer index
        for (auto vv : mesh.vertices(v))
        {
            vertex_valences[idx]++;
        }
    }

    // Regular vertex stats
    int regular_vertices = 0;
    int interior_vertices = 0;
    for (auto v : mesh.vertices())
    {
        if (!mesh.is_boundary(v))
        {
            interior_vertices++;
            int idx = v.idx();
            if (vertex_valences[idx] == 6)
            {
                regular_vertices++;
            }
        }
    }

    // Main loop over faces
    for (auto f : mesh.faces())
    {
        // Compute triangle area (volume in this context)
        float vol = face_area(mesh, f);
        volumes.push_back(vol);

        // Get triangle angles
        std::vector<Point> points_triangle;
        for (auto v : mesh.vertices(f))
            points_triangle.push_back(mesh.position(v));

        auto angles = compute_triangle_angles(
            points_triangle[0], points_triangle[1], points_triangle[2]);

        triangle_count++;

        // Angle and area thresholds
        if (vol < compute_median(volumes))
            triangles_less_than_median++;

        for (float angle : angles)
        {
            if (angle < angle_treshold)
            {
                angles_less_than_threshold.push_back(angle);
                break;
            }
        }
    }

    // Compute the median volume
    float median_volume = compute_median(volumes);

    // Calculate the squared deviations from the median volume
    for (float vol : volumes)
    {
        float deviation = vol - median_volume;
        sum_squared_deviations += deviation * deviation;
    }

    // Compute the quadric deviation (sum of squared deviations)
    float quadric_deviation = std::sqrt(sum_squared_deviations);

    // Compute distance stats
    float max_dist = *std::max_element(distances_from_surface.begin(),
                                       distances_from_surface.end());
    float avg_dist = std::accumulate(distances_from_surface.begin(),
                                     distances_from_surface.end(), 0.0f) /
                     distances_from_surface.size();

    // Compute regular vertex percentage
    float regular_vertex_percentage = 0.0f;
    if (interior_vertices > 0)
        regular_vertex_percentage =
            (static_cast<float>(regular_vertices) / interior_vertices) * 100.0f;

    // Output results
    std::cout << "=== " << label << " ===" << std::endl;
    std::cout << "Total number of triangles: " << triangle_count << std::endl;
    std::cout << "Total number of vertices: " << vertex_count << std::endl;
    std::cout << "Median Volume: " << median_volume << std::endl;
    std::cout << "Triangles with volume less than median: "
              << triangles_less_than_median << std::endl;
    std::cout << "Triangles with at least one angle less than "
              << angle_treshold
              << " degrees: " << angles_less_than_threshold.size() << std::endl;
    std::cout << "Quadric Deviation from Median Volume: " << quadric_deviation
              << std::endl;

    std::cout << "Average unsigned distance from implicit surface: " << avg_dist
              << std::endl;
    std::cout << "Max unsigned distance from implicit surface: " << max_dist
              << std::endl;
    std::cout << "Number of vertices NOT on the surface (tol > "
              << surface_tolerance << "): " << vertices_off_surface
              << std::endl;
    std::cout << "Percentage of regular (valence-6) interior vertices: "
              << regular_vertex_percentage << "%" << std::endl;

    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Metric computation time: "
              << std::chrono::duration<double>(t1 - t0).count() << "s"
              << std::endl;
}

int main(int argc, char **argv)
{
    iso::MyViewer window("View", 1024, 768);

    const char *output_name = "output.obj";
    SurfaceMesh mymesh;

    // Accept input filename from argv, fallback to bundled test meshes
    const char *input_name = "assets/mesh_catalog/wavy_network.obj";
    if (argc > 1)
        input_name = argv[1];

    try
    {
        pmp::read(mymesh, input_name);
    }
    catch (const pmp::IOException &)
    {
        std::cerr << "Warning: failed to open '" << input_name << "'. Trying 'assets/sphere.obj' fallback." << std::endl;
        try
        {
            pmp::read(mymesh, "assets/sphere.obj");
        }
        catch (const pmp::IOException &e)
        {
            std::cerr << "Error: failed to open fallback 'assets/sphere.obj': " << e.what() << std::endl;
            return 1;
        }
    }

    // Compute and print input metrics
    compute_mesh_metrics(mymesh, "INPUT MESH");
    classify_and_export_singularities(mymesh, "INPUT MESH", "input_singularities.ply");

    //  NEW: Project all vertices onto the implicit surface
    auto points = mymesh.get_vertex_property<Point>("v:point");
    for (auto v : mymesh.vertices())
    {
        Point p = points[v];
        points[v] = projectOntoImplicitSurface(p);
    }

    // Compute and print input metrics
    compute_mesh_metrics(mymesh, "PROJECTED MESH");
    classify_and_export_singularities(mymesh, "PROJECTED MESH", "projected_singularities.ply");

    auto start_time = std::chrono::high_resolution_clock::now(); // START

    for (int iter = 0; iter < num_iterations; iter++)
    {
        printf("Iter: %d\n", iter + 1);

        /////////////////////////////////////////////////////////////////////////////// - MEDIAN
        std::vector<float> volumes;

        for (auto f : mymesh.faces())
        {
            float vol = face_area(mymesh, f);
            volumes.push_back(vol);
        }

        // Calculate the median volume instead of the mean volume
        std::sort(volumes.begin(), volumes.end());

        float median_volume;
        size_t n = volumes.size();

        if (n % 2 == 0) // If even number of elements
        {
            median_volume = (volumes[n / 2 - 1] + volumes[n / 2]) / 2.0f;
        }
        else // If odd number of elements
        {
            median_volume = volumes[n / 2];
        }

        float threshold = median_volume / divisor_median;

        /////////////////////////////////////////////////////////////////////////////// - END MEDIAN

        int count_triangles = 0, small_iter = 0;
        bool do_next_cycle;
        do
        {
            do_next_cycle = false;
            small_iter++;
            //printf("Small Iter: %d\n", small_iter);

            //printf("Faces: %d\n", std::distance(mymesh.faces().begin(), mymesh.faces().end()));

            // Step 1: Precompute boundary edges and vertices
            std::unordered_set<Vertex> disabled_vertices;

            for (auto e : mymesh.edges())
            {
                if (mymesh.is_boundary(e))
                {
                    disabled_vertices.insert(mymesh.vertex(e, 0));
                    disabled_vertices.insert(mymesh.vertex(e, 1));
                }
            }

            auto points = mymesh.get_vertex_property<Point>("v:point");

            for (Face f : mymesh.faces())
            {
                //////////////////////////////////////////////////////////
                // Pocita min angle

                std::vector<Point> vertices;
                for (auto v : mymesh.vertices(f))
                {
                    vertices.push_back(points[v]);
                }

                std::vector<float> angles = compute_triangle_angles(
                    vertices[0], vertices[1], vertices[2]);

                float min_angle =
                    *std::min_element(angles.begin(), angles.end());
                //////////////////////////////////////////////////////////

                float vol = face_area(mymesh, f);

                if (vol < threshold || min_angle < angle_treshold)
                {
                    Halfedge shortest_he;

                    float shortest = std::numeric_limits<float>::max();

                    for (auto he : mymesh.halfedges(f))
                    {
                        float length = distance(points[mymesh.from_vertex(he)],
                                                points[mymesh.to_vertex(he)]);
                        if (length < shortest)
                        {
                            shortest = length;
                            shortest_he = he;
                        }
                    }

                    Vertex v0 = mymesh.from_vertex(shortest_he);
                    Vertex v1 = mymesh.to_vertex(shortest_he);

                    if (disabled_vertices.count(v0) ||
                        disabled_vertices.count(v1))
                        continue;

                    if (mymesh.is_collapse_ok(shortest_he))
                    {
                        // STEP 1
                        Vertex old_v = mymesh.to_vertex(shortest_he);

                        mymesh.collapse(shortest_he);
                        do_next_cycle = true;

                        // STEP 2
                        std::vector<Vertex> ring_vertices =
                            getRingVertices(mymesh, old_v);

                        // STEP 3
                        for (Vertex v : ring_vertices)
                        {
                            disabled_vertices.insert(v);
                        }
                        disabled_vertices.insert(old_v);

                        // STEP 4
                        Point centroid = getCentroid(mymesh, ring_vertices);

                        // STEP 5 Apply Newton-Raphson projection onto the implicit surface
                        Point new_p = projectOntoImplicitSurface(centroid);

                        // STEP 6
                        mymesh.position(old_v) = new_p;
                    }
                }
            }
        } while (do_next_cycle);
    }

    mymesh.garbage_collection();
    pmp::write(mymesh, output_name);

    auto end_time = std::chrono::high_resolution_clock::now(); // END
    std::chrono::duration<double> runtime = end_time - start_time;
    std::cout << "=== TOTAL RUNTIME: " << runtime.count()
              << " seconds ===" << std::endl;

    // Compute and print output metrics
    compute_mesh_metrics(mymesh, "OUTPUT MESH");
    classify_and_export_singularities(mymesh, "OUTPUT MESH", "output_singularities.ply");

    window.load_mesh(output_name);
    return window.run();
}