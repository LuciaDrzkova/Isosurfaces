// The remeshing algorithm (see Remesher.h). Read remesh() first, from the
// bottom of the file upwards.

#include "method/Remesher.h"

#include "method/MeshMetrics.h"

#include <pmp/algorithms/differential_geometry.h>

#include <algorithm>
#include <array>
#include <limits>
#include <stdexcept>
#include <unordered_set>
#include <numeric>

using namespace std;

// Provide std::hash specializations for pmp::Vertex and pmp::Edge so the
// unordered_set<Vertex> approach in kod.txt works identically.
namespace std {
template <>
struct hash<pmp::Vertex>
{
    size_t operator()(const pmp::Vertex &v) const noexcept
    {
        return std::hash<int>()(v.idx());
    }
};

template <>
struct hash<pmp::Edge>
{
    size_t operator()(const pmp::Edge &e) const noexcept
    {
        return std::hash<int>()(e.idx());
    }
};
} // namespace std

namespace iso {

using namespace pmp;

namespace {

float median_face_area(const SurfaceMesh& mesh)
{
    // Use sorting to compute the median (matches the kod.txt approach).
    vector<float> areas;
    areas.reserve(mesh.n_faces());
    for (auto f : mesh.faces())
        areas.push_back(face_area(mesh, f));
    if (areas.empty())
        return 0.0f;
    sort(areas.begin(), areas.end());
    const size_t n = areas.size();
    if (n % 2 == 0)
        return 0.5f * (areas[n / 2 - 1] + areas[n / 2]);
    return areas[n / 2];
}

Halfedge shortest_halfedge(const SurfaceMesh& mesh, Face f)
{
    Halfedge shortest;
    float shortest_length = numeric_limits<float>::max();
    for (auto h : mesh.halfedges(f))
    {
        const float length = distance(mesh.position(mesh.from_vertex(h)),
                                      mesh.position(mesh.to_vertex(h)));
        if (length < shortest_length)
        {
            shortest_length = length;
            shortest = h;
        }
    }
    return shortest;
}

// One sweep over all faces. Returns the number of collapses performed.
// Forward declarations for helpers defined below (used by collapse_sweep)
static std::vector<Vertex> getRingVertices(const SurfaceMesh& mesh, Vertex v);
static pmp::Point getCentroid(const SurfaceMesh& mesh, const std::vector<Vertex>& verts);

int collapse_sweep(SurfaceMesh& mesh, const ImplicitSurface& surface,
                   const RemeshOptions& options, float area_threshold)
{
    // Match the kod.txt approach: maintain a set of disabled vertex indices
    // (starts with boundary vertices) and reset it for each sweep.
    std::unordered_set<int> disabled;
    for (auto e : mesh.edges())
    {
        if (mesh.is_boundary(e))
        {
            disabled.insert(mesh.vertex(e, 0).idx());
            disabled.insert(mesh.vertex(e, 1).idx());
        }
    }

    int collapses = 0;
    for (Face f : mesh.faces())
    {
        array<Point, 3> p;
        int i = 0;
        for (auto v : mesh.vertices(f))
            p[i++] = mesh.position(v);

        const auto angles = triangle_angles(p[0], p[1], p[2]);
        const float min_angle = *min_element(angles.begin(), angles.end());

        const float area = face_area(mesh, f);
        if (!(area < area_threshold || min_angle < options.min_angle))
            continue;

        const Halfedge h = shortest_halfedge(mesh, f);
        const Vertex v0 = mesh.from_vertex(h);
        const Vertex v1 = mesh.to_vertex(h);

        if (disabled.count(v0.idx()) || disabled.count(v1.idx()))
            continue;
        if (!mesh.is_collapse_ok(h))
            continue;

        // As in kod.txt, keep the 'old' vertex (the survivor) and collapse
        Vertex old_v = mesh.to_vertex(h);
        mesh.collapse(h);
        ++collapses;

        // Collect ring around the survivor and disable those vertices
        std::vector<Vertex> ring = getRingVertices(mesh, old_v);
        for (const Vertex& rv : ring)
            disabled.insert(rv.idx());
        disabled.insert(old_v.idx());

        // Move survivor to centroid of the ring and project onto surface
        Point centroid = getCentroid(mesh, ring);
        Point new_p = surface.project(centroid);
        mesh.position(old_v) = new_p;
    }
    return collapses;
}

// Helpers matching kod.txt behavior: ring vertices and centroid
static std::vector<Vertex> getRingVertices(const SurfaceMesh& mesh, Vertex v)
{
    std::vector<Vertex> vertices;
    for (auto he : mesh.halfedges(v))
    {
        Vertex vv = mesh.to_vertex(he);
        if (vv != v)
            vertices.push_back(vv);
    }
    return vertices;
}

static pmp::Point getCentroid(const SurfaceMesh& mesh, const std::vector<Vertex>& verts)
{
    pmp::Point centroid(0.0f, 0.0f, 0.0f);
    for (const auto& v : verts)
        centroid += mesh.position(v);
    if (!verts.empty())
        centroid /= static_cast<pmp::Scalar>(verts.size());
    return centroid;
}

} // namespace

void project_onto_surface(SurfaceMesh& mesh, const ImplicitSurface& surface)
{
    for (auto v : mesh.vertices())
        mesh.position(v) = surface.project(mesh.position(v));
}

RemeshResult remesh(SurfaceMesh& mesh, const ImplicitSurface& surface,
                    const RemeshOptions& options)
{
    if (!mesh.is_triangle_mesh())
        throw invalid_argument("remesh: input must be a triangle mesh");
    if (options.area_divisor <= 0.0f)
        throw invalid_argument("remesh: area_divisor must be positive");

    RemeshResult result;
    for (int iter = 0; iter < options.iterations; ++iter)
    {
        const float threshold =
            median_face_area(mesh) / options.area_divisor;

        // Repeat until a sweep finds nothing left to collapse
        int total = 0, collapses;
        do
        {
            collapses = collapse_sweep(mesh, surface, options, threshold);
            total += collapses;
        } while (collapses > 0);

        result.collapses_per_iteration.push_back(total);
    }

    mesh.garbage_collection();
    return result;
}

} // namespace iso
