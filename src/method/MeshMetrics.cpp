// Measurements of a mesh (see MeshMetrics.h).

#include "method/MeshMetrics.h"

#include <pmp/algorithms/differential_geometry.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <ostream>

using namespace std;

namespace iso {

namespace {
constexpr float kRadToDeg = 180.0f / 3.14159265358979323846f;

// Angle at the vertex opposite to side `a`, given the three side lengths.
// Law of cosines; returns 0 if either adjacent side has zero length.
float angle_opposite(float a, float b, float c)
{
    const float denom = 2.0f * b * c;
    if (denom <= 0.0f)
        return 0.0f;
    return acos(clamp((b * b + c * c - a * a) / denom, -1.0f, 1.0f)) *
           kRadToDeg;
}
} // namespace

float median(vector<float> values)
{
    if (values.empty())
        return 0.0f;

    const size_t mid = values.size() / 2;
    nth_element(values.begin(), values.begin() + mid, values.end());
    if (values.size() % 2 == 1)
        return values[mid];

    // Even count: the lower middle value is the max of the lower partition.
    const float lower = *max_element(values.begin(), values.begin() + mid);
    return 0.5f * (lower + values[mid]);
}

array<float, 3> triangle_angles(const pmp::Point& p0, const pmp::Point& p1,
                                     const pmp::Point& p2)
{
    const float a = pmp::distance(p1, p2);
    const float b = pmp::distance(p0, p2);
    const float c = pmp::distance(p0, p1);

    const float alpha = angle_opposite(a, b, c);
    const float beta = angle_opposite(b, a, c);
    const float gamma = max(0.0f, 180.0f - alpha - beta);
    return {alpha, beta, gamma};
}

MeshMetrics compute_mesh_metrics(const pmp::SurfaceMesh& mesh,
                                 const ImplicitSurface& surface,
                                 float min_angle_threshold,
                                 float surface_tolerance)
{
    const auto t0 = chrono::high_resolution_clock::now();
    MeshMetrics m;
    m.n_vertices = mesh.n_vertices();
    m.n_faces = mesh.n_faces();
    m.min_angle_threshold = min_angle_threshold;
    m.surface_tolerance = surface_tolerance;

    // --- Vertices: distance to the surface and valence regularity ----------
    double distance_sum = 0.0;
    size_t interior = 0, regular = 0;
    for (auto v : mesh.vertices())
    {
        const double d = abs(surface.value(mesh.position(v)));
        distance_sum += d;
        m.max_distance = max(m.max_distance, d);
        if (d > surface_tolerance)
            ++m.vertices_off_surface;

        if (!mesh.is_boundary(v))
        {
            ++interior;
            if (mesh.valence(v) == 6)
                ++regular;
        }
    }
    if (m.n_vertices > 0)
        m.avg_distance = distance_sum / static_cast<double>(m.n_vertices);
    if (interior > 0)
        m.regular_vertex_percentage =
            100.0f * static_cast<float>(regular) / static_cast<float>(interior);

    // --- Faces: area distribution and small angles -------------------------
    vector<float> areas;
    areas.reserve(m.n_faces);
    for (auto f : mesh.faces())
    {
        areas.push_back(pmp::face_area(mesh, f));

        array<pmp::Point, 3> p;
        int i = 0;
        for (auto v : mesh.vertices(f))
            if (i < 3)
                p[i++] = mesh.position(v);

        const auto angles = triangle_angles(p[0], p[1], p[2]);
        const float min_angle = *min_element(angles.begin(), angles.end());
        
        if (min_angle < m.min_angle_threshold && min_angle > 0.0f) 
            ++m.faces_with_small_angle;
    }

    m.median_area = median(areas);
    double sum_sq = 0.0;
    for (float a : areas)
    {
        if (a < m.median_area)
            ++m.faces_below_median_area;
        const double dev = a - m.median_area;
        sum_sq += dev * dev;
    }
    m.area_deviation = static_cast<float>(sqrt(sum_sq));
    m.computation_time =
        chrono::duration<double>(chrono::high_resolution_clock::now() - t0)
            .count();

    return m;
}

void print_metrics(ostream& out, const string& label,
                   const MeshMetrics& m)
{
    out << "=== " << label << " ===\n"
        << "Total number of triangles: " << m.n_faces << '\n'
        << "Total number of vertices: " << m.n_vertices << '\n'
        << "Median Volume: " << m.median_area << '\n'
        << "Triangles with volume less than median: "
        << m.faces_below_median_area << '\n'
        << "Triangles with at least one angle less than "
        << m.min_angle_threshold << " degrees: " << m.faces_with_small_angle
        << '\n'
        << "Quadric Deviation from Median Volume: " << m.area_deviation
        << '\n'
        << "Average unsigned distance from implicit surface: "
        << m.avg_distance << '\n'
        << "Max unsigned distance from implicit surface: " << m.max_distance
        << '\n'
        << "Number of vertices NOT on the surface (tol > "
        << m.surface_tolerance << "): " << m.vertices_off_surface << '\n'
        << "Percentage of regular (valence-6) interior vertices: "
        << m.regular_vertex_percentage << "%\n"
        << "Metric computation time: " << m.computation_time << "s\n";
}

} // namespace iso
