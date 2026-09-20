// The remeshing algorithm (see Remesher.h). Read remesh() first, from the
// bottom of the file upwards.

#include "method/Remesher.h"

#include "method/MeshMetrics.h"

#include <pmp/algorithms/differential_geometry.h>

#include <algorithm>
#include <array>
#include <limits>
#include <stdexcept>

using namespace std;

namespace iso {

using namespace pmp;

namespace {

// Vertices on the mesh boundary, indexed by Vertex::idx(). Sized with
// vertices_size(), not n_vertices(): the latter excludes deleted vertices, whose
// slots still exist until garbage collection.
vector<char> boundary_vertices(const SurfaceMesh& mesh)
{
    vector<char> locked(mesh.vertices_size(), 0);
    for (auto e : mesh.edges())
    {
        if (mesh.is_boundary(e))
        {
            locked[mesh.vertex(e, 0).idx()] = 1;
            locked[mesh.vertex(e, 1).idx()] = 1;
        }
    }
    return locked;
}

float median_face_area(const SurfaceMesh& mesh)
{
    vector<float> areas;
    areas.reserve(mesh.n_faces());
    for (auto f : mesh.faces())
        areas.push_back(face_area(mesh, f));
    return median(std::move(areas));
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
int collapse_sweep(SurfaceMesh& mesh, const ImplicitSurface& surface,
                   const RemeshOptions& options, float area_threshold)
{
    // Boundary vertices stay fixed for the whole run; vertices touched by a
    // collapse are additionally locked until the next sweep.
    vector<char> locked = boundary_vertices(mesh);
    int collapses = 0;

    for (Face f : mesh.faces())
    {
        array<Point, 3> p;
        int i = 0;
        for (auto v : mesh.vertices(f))
            p[i++] = mesh.position(v);

        const auto angles = triangle_angles(p[0], p[1], p[2]);
        const float min_angle = *min_element(angles.begin(), angles.end());

        if (face_area(mesh, f) >= area_threshold &&
            min_angle >= options.min_angle)
            continue;

        const Halfedge h = shortest_halfedge(mesh, f);
        const Vertex from = mesh.from_vertex(h);
        const Vertex kept = mesh.to_vertex(h); // survives the collapse

        if (locked[from.idx()] || locked[kept.idx()])
            continue;
        if (!mesh.is_collapse_ok(h))
            continue;

        mesh.collapse(h);
        ++collapses;

        // Move the survivor to the centroid of its new neighbourhood, back
        // onto the surface, and freeze that neighbourhood for this sweep.
        Point centroid(0, 0, 0);
        int n_ring = 0;
        for (auto v : mesh.vertices(kept))
        {
            centroid += mesh.position(v);
            locked[v.idx()] = 1;
            ++n_ring;
        }
        locked[kept.idx()] = 1;

        if (n_ring > 0)
            mesh.position(kept) =
                surface.project(centroid / static_cast<Scalar>(n_ring));
    }
    return collapses;
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
