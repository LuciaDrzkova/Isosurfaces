#pragma once

#include "MarchingCubes.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace iso::stitch
{

struct TopologyStatistics
{
    std::size_t vertices = 0;
    std::size_t triangles = 0;
    std::size_t edges = 0;
    std::size_t boundary_edges = 0;
    std::size_t boundary_vertices = 0;
    std::size_t nonmanifold_edges = 0;
    std::size_t connected_components = 0;
};

struct StitchResult
{
    iso::mc::Mesh mesh;
    bool stitched = false;
    std::size_t global_triangles_removed = 0;
    std::size_t local_triangles_inserted = 0;
    std::size_t seam_loops = 0;
    std::size_t seam_vertices = 0;
    std::size_t seam_triangles = 0;
    double max_seam_vertex_distance = 0.0;
    TopologyStatistics topology;
};

namespace detail
{

using Index = iso::mc::Index;
using Point = iso::mc::Point;
using Triangle = iso::mc::Triangle;
using Edge = std::pair<Index, Index>;

inline constexpr Index invalid_index = std::numeric_limits<Index>::max();

inline Edge make_edge(Index a, Index b)
{
    if (a > b)
        std::swap(a, b);
    return {a, b};
}

inline double squared_distance(const Point& a, const Point& b)
{
    const double dx = a[0] - b[0];
    const double dy = a[1] - b[1];
    const double dz = a[2] - b[2];
    return dx * dx + dy * dy + dz * dz;
}

inline double distance(const Point& a, const Point& b)
{
    return std::sqrt(squared_distance(a, b));
}

inline Point midpoint(const Point& a, const Point& b)
{
    return {0.5 * (a[0] + b[0]),
            0.5 * (a[1] + b[1]),
            0.5 * (a[2] + b[2])};
}

inline Point interpolate(const Point& a, const Point& b, double t)
{
    return {a[0] + t * (b[0] - a[0]),
            a[1] + t * (b[1] - a[1]),
            a[2] + t * (b[2] - a[2])};
}

inline Point centroid(const Point& a, const Point& b, const Point& c)
{
    return {(a[0] + b[0] + c[0]) / 3.0,
            (a[1] + b[1] + c[1]) / 3.0,
            (a[2] + b[2] + c[2]) / 3.0};
}

inline double radius_squared(const Point& p, const Point& center)
{
    return squared_distance(p, center);
}

inline TopologyStatistics topology(const iso::mc::Mesh& mesh)
{
    TopologyStatistics result;
    result.vertices = mesh.vertices.size();
    result.triangles = mesh.triangles.size();

    std::map<Edge, std::size_t> uses;
    for (const Triangle& t : mesh.triangles)
    {
        if (t[0] >= mesh.vertices.size() ||
            t[1] >= mesh.vertices.size() ||
            t[2] >= mesh.vertices.size())
            throw std::runtime_error("Topology: triangle index is out of range");

        ++uses[make_edge(t[0], t[1])];
        ++uses[make_edge(t[1], t[2])];
        ++uses[make_edge(t[2], t[0])];
    }

    result.edges = uses.size();
    std::set<Index> boundary_vertices;
    std::vector<std::vector<Index>> neighbors(mesh.vertices.size());

    for (const auto& [edge, count] : uses)
    {
        neighbors[edge.first].push_back(edge.second);
        neighbors[edge.second].push_back(edge.first);

        if (count == 1)
        {
            ++result.boundary_edges;
            boundary_vertices.insert(edge.first);
            boundary_vertices.insert(edge.second);
        }
        else if (count > 2)
        {
            ++result.nonmanifold_edges;
        }
    }
    result.boundary_vertices = boundary_vertices.size();

    std::vector<char> visited(mesh.vertices.size(), 0);
    for (Index start = 0;
         start < static_cast<Index>(mesh.vertices.size());
         ++start)
    {
        if (visited[start])
            continue;

        ++result.connected_components;
        std::queue<Index> q;
        q.push(start);
        visited[start] = 1;

        while (!q.empty())
        {
            const Index v = q.front();
            q.pop();
            for (const Index n : neighbors[v])
            {
                if (!visited[n])
                {
                    visited[n] = 1;
                    q.push(n);
                }
            }
        }
    }

    return result;
}

/*
 * Clip every triangle against a sphere using the signed value
 * |p-center|^2-radius^2 linearly interpolated over each source edge.
 * Edge intersections are cached by the source edge, so adjacent triangles
 * reuse the exact same seam vertex. The returned seam is therefore a true
 * triangle-edge cut, rather than a jagged selection by triangle centroids.
 */
inline iso::mc::Mesh clip_by_sphere(
    const iso::mc::Mesh& source,
    const Point& center,
    double r2,
    bool keep_inside)
{
    iso::mc::Mesh out;
    if (source.triangles.empty())
        return out;

    std::vector<double> signed_value(source.vertices.size(), 0.0);
    for (std::size_t i = 0; i < source.vertices.size(); ++i)
        signed_value[i] = radius_squared(source.vertices[i], center) - r2;

    std::vector<Index> source_to_out(source.vertices.size(), invalid_index);
    std::map<Edge, Index> crossing_vertices;
    const double zero_tol = 1e-12 * std::max(1.0, r2);

    const auto copy_original = [&](Index source_index) -> Index
    {
        Index& mapped = source_to_out[source_index];
        if (mapped != invalid_index)
            return mapped;
        mapped = static_cast<Index>(out.vertices.size());
        out.vertices.push_back(source.vertices[source_index]);
        return mapped;
    };

    const auto intersection = [&](Index a, Index b) -> Index
    {
        if (std::abs(signed_value[a]) <= zero_tol)
            return copy_original(a);
        if (std::abs(signed_value[b]) <= zero_tol)
            return copy_original(b);

        const Edge key = make_edge(a, b);
        const auto existing = crossing_vertices.find(key);
        if (existing != crossing_vertices.end())
            return existing->second;

        const double denominator = signed_value[a] - signed_value[b];
        if (std::abs(denominator) <= std::numeric_limits<double>::epsilon())
            return copy_original(a);

        const double t = std::clamp(signed_value[a] / denominator, 0.0, 1.0);
        const Index new_index = static_cast<Index>(out.vertices.size());
        out.vertices.push_back(interpolate(
            source.vertices[a], source.vertices[b], t));
        crossing_vertices.emplace(key, new_index);
        return new_index;
    };

    const auto is_kept = [&](double value)
    {
        return keep_inside ? value <= 0.0 : value >= 0.0;
    };

    for (const Triangle& triangle : source.triangles)
    {
        std::vector<Index> polygon;
        polygon.reserve(4);

        for (std::size_t edge_i = 0; edge_i < 3; ++edge_i)
        {
            const Index a = triangle[edge_i];
            const Index b = triangle[(edge_i + 1) % 3];
            const bool a_in = is_kept(signed_value[a]);
            const bool b_in = is_kept(signed_value[b]);

            if (a_in && b_in)
            {
                polygon.push_back(copy_original(b));
            }
            else if (a_in && !b_in)
            {
                polygon.push_back(intersection(a, b));
            }
            else if (!a_in && b_in)
            {
                polygon.push_back(intersection(a, b));
                polygon.push_back(copy_original(b));
            }
        }

        // Remove duplicate consecutive indices (possible at an exact seam vertex).
        std::vector<Index> cleaned;
        cleaned.reserve(polygon.size());
        for (const Index v : polygon)
        {
            if (cleaned.empty() || cleaned.back() != v)
                cleaned.push_back(v);
        }
        if (cleaned.size() > 1 && cleaned.front() == cleaned.back())
            cleaned.pop_back();

        if (cleaned.size() < 3)
            continue;

        // Convex clipping of one triangle yields at most a quadrilateral.
        for (std::size_t i = 1; i + 1 < cleaned.size(); ++i)
        {
            const Index a = cleaned[0];
            const Index b = cleaned[i];
            const Index c = cleaned[i + 1];
            if (a == b || b == c || c == a)
                continue;
            out.triangles.push_back({a, b, c});
        }
    }

    // Compact away any vertices that belonged only to clipped-away polygons.
    std::vector<char> used(out.vertices.size(), 0);
    for (const Triangle& t : out.triangles)
    {
        used[t[0]] = used[t[1]] = used[t[2]] = 1;
    }
    std::vector<Index> compact_map(out.vertices.size(), invalid_index);
    std::vector<Point> compact_vertices;
    compact_vertices.reserve(out.vertices.size());
    for (Index i = 0; i < static_cast<Index>(out.vertices.size()); ++i)
    {
        if (!used[i])
            continue;
        compact_map[i] = static_cast<Index>(compact_vertices.size());
        compact_vertices.push_back(out.vertices[i]);
    }
    for (Triangle& t : out.triangles)
    {
        t[0] = compact_map[t[0]];
        t[1] = compact_map[t[1]];
        t[2] = compact_map[t[2]];
    }
    out.vertices = std::move(compact_vertices);
    return out;
}

inline std::vector<std::vector<Index>> boundary_loops(const iso::mc::Mesh& mesh)
{
    std::map<Edge, std::size_t> uses;
    for (const Triangle& t : mesh.triangles)
    {
        ++uses[make_edge(t[0], t[1])];
        ++uses[make_edge(t[1], t[2])];
        ++uses[make_edge(t[2], t[0])];
    }

    std::vector<std::vector<Index>> adjacency(mesh.vertices.size());
    std::set<Edge> boundary_edges;
    for (const auto& [edge, count] : uses)
    {
        if (count != 1)
            continue;
        boundary_edges.insert(edge);
        adjacency[edge.first].push_back(edge.second);
        adjacency[edge.second].push_back(edge.first);
    }

    for (const Edge& e : boundary_edges)
    {
        if (adjacency[e.first].size() != 2 || adjacency[e.second].size() != 2)
            throw std::runtime_error(
                "Clipped mesh boundary is not a collection of simple loops");
    }

    std::set<Edge> visited_edges;
    std::vector<std::vector<Index>> loops;
    for (const Edge& start_edge : boundary_edges)
    {
        if (visited_edges.count(start_edge))
            continue;

        const Index start = start_edge.first;
        Index previous = invalid_index;
        Index current = start;
        std::vector<Index> loop;

        for (std::size_t guard = 0; guard <= boundary_edges.size() + 1; ++guard)
        {
            loop.push_back(current);
            Index next = invalid_index;
            for (const Index candidate : adjacency[current])
            {
                const Edge e = make_edge(current, candidate);
                if (visited_edges.count(e))
                    continue;
                if (candidate != previous || adjacency[current].size() <= 1)
                {
                    next = candidate;
                    break;
                }
            }

            if (next == invalid_index)
                break;

            visited_edges.insert(make_edge(current, next));
            previous = current;
            current = next;
            if (current == start)
                break;
        }

        if (loop.size() >= 3 && current == start)
            loops.push_back(std::move(loop));
    }
    return loops;
}

inline Point loop_centroid(
    const iso::mc::Mesh& mesh,
    const std::vector<Index>& loop)
{
    Point c{0.0, 0.0, 0.0};
    if (loop.empty())
        return c;
    for (const Index v : loop)
    {
        c[0] += mesh.vertices[v][0];
        c[1] += mesh.vertices[v][1];
        c[2] += mesh.vertices[v][2];
    }
    const double inv = 1.0 / static_cast<double>(loop.size());
    c[0] *= inv;
    c[1] *= inv;
    c[2] *= inv;
    return c;
}

inline double loop_mean_radius(
    const iso::mc::Mesh& mesh,
    const std::vector<Index>& loop,
    const Point& center)
{
    if (loop.empty())
        return 0.0;
    double sum = 0.0;
    for (const Index v : loop)
        sum += distance(mesh.vertices[v], center);
    return sum / static_cast<double>(loop.size());
}

inline void split_boundary_edge(
    iso::mc::Mesh& mesh,
    std::vector<Index>& loop,
    std::size_t edge_position,
    double fraction = 0.5)
{
    if (loop.size() < 3 || edge_position >= loop.size())
        throw std::runtime_error("Cannot split an invalid boundary loop edge");

    fraction = std::clamp(fraction, 1e-8, 1.0 - 1e-8);
    const std::size_t next_position = (edge_position + 1) % loop.size();
    const Index a = loop[edge_position];
    const Index b = loop[next_position];
    const Index new_vertex = static_cast<Index>(mesh.vertices.size());
    mesh.vertices.push_back(interpolate(mesh.vertices[a], mesh.vertices[b], fraction));

    bool found = false;
    std::vector<Triangle> rebuilt;
    rebuilt.reserve(mesh.triangles.size() + 1);
    for (const Triangle& t : mesh.triangles)
    {
        bool matched = false;
        for (std::size_t k = 0; k < 3; ++k)
        {
            const Index first = t[k];
            const Index second = t[(k + 1) % 3];
            const Index opposite = t[(k + 2) % 3];
            if ((first == a && second == b) || (first == b && second == a))
            {
                // Keep the original directed edge order and triangle winding.
                rebuilt.push_back({first, new_vertex, opposite});
                rebuilt.push_back({new_vertex, second, opposite});
                matched = true;
                found = true;
                break;
            }
        }
        if (!matched)
            rebuilt.push_back(t);
    }

    if (!found)
    {
        mesh.vertices.pop_back();
        throw std::runtime_error("Boundary edge could not be split");
    }

    mesh.triangles = std::move(rebuilt);
    loop.insert(loop.begin() + static_cast<std::ptrdiff_t>(edge_position + 1),
                new_vertex);
}

inline void split_boundary_edge_near_other_loop(
    iso::mc::Mesh& target_mesh,
    std::vector<Index>& target_loop,
    const iso::mc::Mesh& guide_mesh,
    const std::vector<Index>& guide_loop)
{
    if (target_loop.size() < 3 || guide_loop.empty())
        throw std::runtime_error("Cannot equalize invalid boundary loops");

    std::size_t best_edge = 0;
    double best_fraction = 0.5;
    double best_d2 = std::numeric_limits<double>::infinity();
    bool found_candidate = false;

    for (const Index gv : guide_loop)
    {
        const Point& p = guide_mesh.vertices[gv];
        for (std::size_t i = 0; i < target_loop.size(); ++i)
        {
            const Point& a = target_mesh.vertices[target_loop[i]];
            const Point& b = target_mesh.vertices[target_loop[(i + 1) % target_loop.size()]];
            const double dx = b[0] - a[0];
            const double dy = b[1] - a[1];
            const double dz = b[2] - a[2];
            const double length2 = dx * dx + dy * dy + dz * dz;
            if (length2 <= 1e-24)
                continue;

            double t = ((p[0] - a[0]) * dx +
                        (p[1] - a[1]) * dy +
                        (p[2] - a[2]) * dz) / length2;
            t = std::clamp(t, 0.0, 1.0);
            if (t <= 0.10 || t >= 0.90)
                continue;

            const Point projected = interpolate(a, b, t);
            const double d2 = squared_distance(p, projected);
            if (d2 < best_d2)
            {
                best_d2 = d2;
                best_edge = i;
                best_fraction = t;
                found_candidate = true;
            }
        }
    }

    if (!found_candidate)
    {
        double longest2 = -1.0;
        for (std::size_t i = 0; i < target_loop.size(); ++i)
        {
            const double d2 = squared_distance(
                target_mesh.vertices[target_loop[i]],
                target_mesh.vertices[target_loop[(i + 1) % target_loop.size()]]);
            if (d2 > longest2)
            {
                longest2 = d2;
                best_edge = i;
            }
        }
        best_fraction = 0.5;
    }

    split_boundary_edge(target_mesh, target_loop, best_edge, best_fraction);
}

inline void equalize_loop_sizes(
    iso::mc::Mesh& mesh_a,
    std::vector<Index>& loop_a,
    iso::mc::Mesh& mesh_b,
    std::vector<Index>& loop_b)
{
    const std::size_t target = std::max(loop_a.size(), loop_b.size());
    while (loop_a.size() < target)
        split_boundary_edge_near_other_loop(mesh_a, loop_a, mesh_b, loop_b);
    while (loop_b.size() < target)
        split_boundary_edge_near_other_loop(mesh_b, loop_b, mesh_a, loop_a);
}

struct LoopMatch
{
    std::size_t shift = 0;
    bool reversed = false;
    double cost = std::numeric_limits<double>::infinity();
};

inline LoopMatch best_loop_match(
    const iso::mc::Mesh& a_mesh,
    const std::vector<Index>& a,
    const iso::mc::Mesh& b_mesh,
    const std::vector<Index>& b)
{
    if (a.size() != b.size() || a.empty())
        throw std::runtime_error("Boundary loops must have equal non-zero size");

    LoopMatch best;
    for (const bool reversed : {false, true})
    {
        for (std::size_t shift = 0; shift < b.size(); ++shift)
        {
            double cost = 0.0;
            for (std::size_t i = 0; i < a.size(); ++i)
            {
                const std::size_t j = reversed
                    ? (shift + b.size() - i) % b.size()
                    : (shift + i) % b.size();
                cost += squared_distance(a_mesh.vertices[a[i]], b_mesh.vertices[b[j]]);
            }
            if (cost < best.cost)
            {
                best.cost = cost;
                best.shift = shift;
                best.reversed = reversed;
            }
        }
    }
    return best;
}

inline void orient_consistently_from_global(
    iso::mc::Mesh& mesh,
    std::size_t global_triangle_count)
{
    struct Incidence
    {
        Index triangle = 0;
        bool forward = false;
    };

    std::map<Edge, std::vector<Incidence>> incidence;
    for (Index ti = 0; ti < static_cast<Index>(mesh.triangles.size()); ++ti)
    {
        const Triangle& t = mesh.triangles[ti];
        for (const auto& [a, b] : std::array<std::pair<Index, Index>, 3>{{
                 {t[0], t[1]}, {t[1], t[2]}, {t[2], t[0]}}})
        {
            const Edge e = make_edge(a, b);
            incidence[e].push_back({ti, a == e.first && b == e.second});
        }
    }
    for (const auto& [edge, uses] : incidence)
    {
        (void)edge;
        if (uses.size() > 2)
            throw std::runtime_error("Cannot orient stitched mesh: non-manifold edge");
    }

    std::vector<int> flip(mesh.triangles.size(), -1);
    std::queue<Index> q;
    const Index global_count = static_cast<Index>(
        std::min(global_triangle_count, mesh.triangles.size()));

    const auto propagate = [&](Index seed)
    {
        while (!q.empty())
        {
            const Index current = q.front();
            q.pop();
            const Triangle& t = mesh.triangles[current];
            for (const auto& [a, b] : std::array<std::pair<Index, Index>, 3>{{
                     {t[0], t[1]}, {t[1], t[2]}, {t[2], t[0]}}})
            {
                const Edge e = make_edge(a, b);
                const auto it = incidence.find(e);
                if (it == incidence.end() || it->second.size() != 2)
                    continue;

                const auto& uses = it->second;
                const Incidence* other = uses[0].triangle == current
                    ? &uses[1] : &uses[0];
                const bool current_forward = a == e.first && b == e.second;
                const int required = (current_forward == other->forward)
                    ? 1 - flip[current] : flip[current];

                if (flip[other->triangle] == -1)
                {
                    flip[other->triangle] = required;
                    q.push(other->triangle);
                }
                else if (flip[other->triangle] != required)
                {
                    throw std::runtime_error("Cannot consistently orient stitched mesh");
                }
            }
        }
        (void)seed;
    };

    // Preserve global-face winding and propagate it into the local patch.
    for (Index ti = 0; ti < global_count; ++ti)
    {
        if (flip[ti] != -1)
            continue;
        flip[ti] = 0;
        q.push(ti);
        propagate(ti);
    }
    // Orient disconnected local components consistently as well.
    for (Index ti = 0; ti < static_cast<Index>(mesh.triangles.size()); ++ti)
    {
        if (flip[ti] != -1)
            continue;
        flip[ti] = 0;
        q.push(ti);
        propagate(ti);
    }

    for (Index ti = 0; ti < static_cast<Index>(mesh.triangles.size()); ++ti)
    {
        if (flip[ti] == 1)
            std::swap(mesh.triangles[ti][1], mesh.triangles[ti][2]);
    }
}

} // namespace detail

inline TopologyStatistics compute_topology(const iso::mc::Mesh& mesh)
{
    return detail::topology(mesh);
}

class LocalVoxelStitcher
{
public:
    using Index = detail::Index;
    using Triangle = detail::Triangle;

    struct Options
    {
        double radius = 0.3;
        // Kept for API compatibility. Stitching uses a sphere-cut seam.
        bool use_sphere_seam = true;
    };

    StitchResult stitch(
        const iso::mc::Mesh& global,
        const iso::mc::Mesh& local_patch,
        const iso::mc::Point& center) const
    {
        return stitch(global, local_patch, center, Options{});
    }

    StitchResult stitch(
        const iso::mc::Mesh& global,
        const iso::mc::Mesh& local_patch,
        const iso::mc::Point& center,
        const Options& options) const
    {
        if (!(options.radius > 0.0))
            throw std::invalid_argument("LocalVoxelStitcher radius must be positive");
        if (!options.use_sphere_seam)
            throw std::invalid_argument(
                "LocalVoxelStitcher currently requires use_sphere_seam=true");

        StitchResult result;
        result.mesh = global;
        result.topology = detail::topology(global);
        if (global.triangles.empty() || local_patch.triangles.empty())
            return result;

        const double r2 = options.radius * options.radius;
        const double loop_radius_tolerance = std::max(0.5 * options.radius, 1e-8);

        // Keep the outside of the global mesh and the inside of the local patch.
        // Clipping each face produces a conforming seam on the source triangle edges.
        iso::mc::Mesh outside;
        iso::mc::Mesh inside;
        try
        {
            outside = detail::clip_by_sphere(global, center, r2, false);
            inside = detail::clip_by_sphere(local_patch, center, r2, true);
        }
        catch (const std::exception&)
        {
            return result;
        }

        if (outside.triangles.empty() || inside.triangles.empty())
            return result;

        // Keep this statistic based on original faces (not clipped triangle counts,
        // which may increase when a triangle is split by the seam).
        for (const Triangle& t : global.triangles)
        {
            const detail::Point c = detail::centroid(
                global.vertices[t[0]], global.vertices[t[1]], global.vertices[t[2]]);
            if (detail::squared_distance(c, center) < r2)
                ++result.global_triangles_removed;
        }

        std::vector<std::vector<Index>> global_loops;
        std::vector<std::vector<Index>> local_loops;
        try
        {
            global_loops = detail::boundary_loops(outside);
            local_loops = detail::boundary_loops(inside);
        }
        catch (const std::exception&)
        {
            return result;
        }

        // The global outside mesh also has its original domain-boundary loops.
        // Only pair loops close to the requested sphere radius; do not require
        // the unrelated outer-boundary loops to have local counterparts.
        std::vector<std::size_t> global_seams;
        std::vector<std::size_t> local_seams;
        for (std::size_t i = 0; i < global_loops.size(); ++i)
        {
            const double mean_r = detail::loop_mean_radius(outside, global_loops[i], center);
            if (std::abs(mean_r - options.radius) <= loop_radius_tolerance)
                global_seams.push_back(i);
        }
        for (std::size_t i = 0; i < local_loops.size(); ++i)
        {
            const double mean_r = detail::loop_mean_radius(inside, local_loops[i], center);
            if (std::abs(mean_r - options.radius) <= loop_radius_tolerance)
                local_seams.push_back(i);
        }

        if (global_seams.empty() || global_seams.size() != local_seams.size())
            return result;

        struct MatchedLoops
        {
            std::vector<Index> global_loop;
            std::vector<Index> local_loop;
            detail::LoopMatch match;
        };
        std::vector<MatchedLoops> matched;
        matched.reserve(local_seams.size());
        std::set<std::size_t> used_global;

        try
        {
            for (const std::size_t local_index : local_seams)
            {
                const auto& local_loop = local_loops[local_index];
                const detail::Point lc = detail::loop_centroid(inside, local_loop);
                const double lr = detail::loop_mean_radius(inside, local_loop, center);
                double best_score = std::numeric_limits<double>::infinity();
                std::size_t best_global = std::numeric_limits<std::size_t>::max();

                for (const std::size_t global_index : global_seams)
                {
                    if (used_global.count(global_index))
                        continue;
                    const auto& global_loop = global_loops[global_index];
                    const detail::Point gc = detail::loop_centroid(outside, global_loop);
                    const double gr = detail::loop_mean_radius(outside, global_loop, center);
                    const double score = detail::distance(lc, gc) + 0.5 * std::abs(lr - gr);
                    if (score < best_score)
                    {
                        best_score = score;
                        best_global = global_index;
                    }
                }

                if (best_global == std::numeric_limits<std::size_t>::max())
                    return result;
                used_global.insert(best_global);

                auto gloop = global_loops[best_global];
                auto lloop = local_loops[local_index];
                detail::equalize_loop_sizes(outside, gloop, inside, lloop);
                if (gloop.size() != lloop.size() || gloop.size() < 3)
                    throw std::runtime_error("Unable to equalize seam loop sizes");

                const detail::LoopMatch match = detail::best_loop_match(
                    outside, gloop, inside, lloop);
                matched.push_back({std::move(gloop), std::move(lloop), match});
            }
        }
        catch (const std::exception&)
        {
            return result;
        }

        // At this point the cut loops are compatible. The output shares one
        // vertex index for every pair on the seam; no bridge strip is necessary.
        result.mesh = std::move(outside);
        const std::size_t global_triangle_count = result.mesh.triangles.size();
        std::vector<Index> local_remap(inside.vertices.size(), detail::invalid_index);

        for (const MatchedLoops& loops : matched)
        {
            const std::size_t n = loops.global_loop.size();
            ++result.seam_loops;
            result.seam_vertices += n;

            for (std::size_t i = 0; i < n; ++i)
            {
                const std::size_t j = loops.match.reversed
                    ? (loops.match.shift + n - i) % n
                    : (loops.match.shift + i) % n;
                const Index gv = loops.global_loop[i];
                const Index lv = loops.local_loop[j];

                if (local_remap[lv] != detail::invalid_index && local_remap[lv] != gv)
                {
                    result.mesh = global;
                    result.topology = detail::topology(global);
                    result.stitched = false;
                    return result;
                }

                const detail::Point gp = result.mesh.vertices[gv];
                const detail::Point lp = inside.vertices[lv];
                result.max_seam_vertex_distance = std::max(
                    result.max_seam_vertex_distance, detail::distance(gp, lp));
                result.mesh.vertices[gv] = detail::midpoint(gp, lp);
                local_remap[lv] = gv;
            }
        }

        for (Index lv = 0; lv < static_cast<Index>(inside.vertices.size()); ++lv)
        {
            if (local_remap[lv] != detail::invalid_index)
                continue;
            local_remap[lv] = static_cast<Index>(result.mesh.vertices.size());
            result.mesh.vertices.push_back(inside.vertices[lv]);
        }

        result.local_triangles_inserted = 0;
        for (const Triangle& t : inside.triangles)
        {
            const Index a = local_remap[t[0]];
            const Index b = local_remap[t[1]];
            const Index c = local_remap[t[2]];
            if (a == b || b == c || c == a)
                continue;
            result.mesh.triangles.push_back({a, b, c});
            ++result.local_triangles_inserted;
        }
        result.seam_triangles = 0;

        try
        {
            detail::orient_consistently_from_global(
                result.mesh, global_triangle_count);
        }
        catch (const std::exception&)
        {
            result.stitched = false;
            result.topology = detail::topology(result.mesh);
            return result;
        }

        result.topology = detail::topology(result.mesh);
        const TopologyStatistics original_topology = detail::topology(global);
        result.stitched = result.seam_loops == matched.size() &&
                          result.seam_loops == local_seams.size() &&
                          result.topology.nonmanifold_edges == 0 &&
                          result.topology.boundary_edges == original_topology.boundary_edges;

        if (!result.stitched)
        {
            // Do not silently present the clipped partial result as a valid stitch.
            // Preserve the failed output for diagnosis in the saved OBJ/report.
        }
        return result;
    }
};

} // namespace iso::stitch
