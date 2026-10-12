#pragma once

#include "MarchingCubes.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <stdexcept>
#include <string>
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
    std::string failure_reason;
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
    if (a > b) std::swap(a, b);
    return {a, b};
}

inline double squared_distance(const Point& a, const Point& b)
{
    const double x = a[0] - b[0];
    const double y = a[1] - b[1];
    const double z = a[2] - b[2];
    return x*x + y*y + z*z;
}
inline double distance(const Point& a, const Point& b)
{
    return std::sqrt(squared_distance(a, b));
}
inline Point interpolate(const Point& a, const Point& b, double t)
{
    return {a[0] + t*(b[0]-a[0]),
            a[1] + t*(b[1]-a[1]),
            a[2] + t*(b[2]-a[2])};
}
inline Point midpoint(const Point& a, const Point& b)
{
    return interpolate(a, b, 0.5);
}
inline Point centroid(const Point& a, const Point& b, const Point& c)
{
    return {(a[0]+b[0]+c[0])/3.0,
            (a[1]+b[1]+c[1])/3.0,
            (a[2]+b[2]+c[2])/3.0};
}
inline double dot(const Point& a, const Point& b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline TopologyStatistics topology(const iso::mc::Mesh& mesh)
{
    TopologyStatistics out;
    out.vertices = mesh.vertices.size();
    out.triangles = mesh.triangles.size();

    std::map<Edge, std::size_t> uses;
    for (const Triangle& t : mesh.triangles)
    {
        if (t[0] >= mesh.vertices.size() ||
            t[1] >= mesh.vertices.size() ||
            t[2] >= mesh.vertices.size())
            throw std::runtime_error("Topology check: triangle index out of range");
        ++uses[make_edge(t[0], t[1])];
        ++uses[make_edge(t[1], t[2])];
        ++uses[make_edge(t[2], t[0])];
    }

    out.edges = uses.size();
    std::set<Index> boundary_vertices;
    std::vector<std::vector<Index>> neighbors(mesh.vertices.size());
    for (const auto& [edge, count] : uses)
    {
        neighbors[edge.first].push_back(edge.second);
        neighbors[edge.second].push_back(edge.first);
        if (count == 1)
        {
            ++out.boundary_edges;
            boundary_vertices.insert(edge.first);
            boundary_vertices.insert(edge.second);
        }
        else if (count > 2)
        {
            ++out.nonmanifold_edges;
        }
    }
    out.boundary_vertices = boundary_vertices.size();

    std::vector<char> visited(mesh.vertices.size(), 0);
    for (Index start = 0; start < static_cast<Index>(mesh.vertices.size()); ++start)
    {
        if (visited[start]) continue;
        ++out.connected_components;
        std::queue<Index> q;
        q.push(start);
        visited[start] = 1;
        while (!q.empty())
        {
            const Index v = q.front(); q.pop();
            for (Index n : neighbors[v])
            {
                if (!visited[n])
                {
                    visited[n] = 1;
                    q.push(n);
                }
            }
        }
    }
    return out;
}

/*
 * Clip triangles against a sphere using a linearly interpolated signed
 * radius value. Intersections are cached by source edge so adjacent faces
 * share exactly the same seam vertex. This produces a piecewise-linear seam.
 */
inline iso::mc::Mesh clip_by_sphere(
    const iso::mc::Mesh& source,
    const Point& center,
    double radius_squared,
    bool keep_inside)
{
    iso::mc::Mesh out;
    if (source.triangles.empty()) return out;

    std::vector<double> value(source.vertices.size(), 0.0);
    for (std::size_t i = 0; i < source.vertices.size(); ++i)
        value[i] = squared_distance(source.vertices[i], center) - radius_squared;

    std::vector<Index> original_map(source.vertices.size(), invalid_index);
    std::map<Edge, Index> cut_map;
    const double zero_tol = 1e-12 * std::max(1.0, radius_squared);

    const auto copy_original = [&](Index i) -> Index
    {
        Index& mapped = original_map[i];
        if (mapped == invalid_index)
        {
            mapped = static_cast<Index>(out.vertices.size());
            out.vertices.push_back(source.vertices[i]);
        }
        return mapped;
    };
    const auto intersection = [&](Index a, Index b) -> Index
    {
        if (std::abs(value[a]) <= zero_tol) return copy_original(a);
        if (std::abs(value[b]) <= zero_tol) return copy_original(b);
        const Edge key = make_edge(a, b);
        auto found = cut_map.find(key);
        if (found != cut_map.end()) return found->second;
        const double denominator = value[a] - value[b];
        if (std::abs(denominator) <= std::numeric_limits<double>::epsilon())
            return copy_original(a);
        const double t = std::clamp(value[a] / denominator, 0.0, 1.0);
        const Index id = static_cast<Index>(out.vertices.size());
        out.vertices.push_back(interpolate(source.vertices[a], source.vertices[b], t));
        cut_map.emplace(key, id);
        return id;
    };
    const auto kept = [&](double v) { return keep_inside ? v <= 0.0 : v >= 0.0; };

    for (const Triangle& tri : source.triangles)
    {
        std::vector<Index> polygon;
        polygon.reserve(4);
        for (std::size_t e = 0; e < 3; ++e)
        {
            const Index a = tri[e];
            const Index b = tri[(e + 1) % 3];
            const bool ai = kept(value[a]);
            const bool bi = kept(value[b]);
            if (ai && bi) polygon.push_back(copy_original(b));
            else if (ai && !bi) polygon.push_back(intersection(a, b));
            else if (!ai && bi)
            {
                polygon.push_back(intersection(a, b));
                polygon.push_back(copy_original(b));
            }
        }

        std::vector<Index> clean;
        clean.reserve(polygon.size());
        for (Index v : polygon)
            if (clean.empty() || clean.back() != v) clean.push_back(v);
        if (clean.size() > 1 && clean.front() == clean.back()) clean.pop_back();
        if (clean.size() < 3) continue;
        for (std::size_t i = 1; i + 1 < clean.size(); ++i)
        {
            const Index a = clean[0], b = clean[i], c = clean[i+1];
            if (a != b && b != c && c != a) out.triangles.push_back({a,b,c});
        }
    }

    // Remove vertices that are referenced by no retained triangle.
    std::vector<char> used(out.vertices.size(), 0);
    for (const Triangle& t : out.triangles) used[t[0]] = used[t[1]] = used[t[2]] = 1;
    std::vector<Index> remap(out.vertices.size(), invalid_index);
    std::vector<Point> compact;
    compact.reserve(out.vertices.size());
    for (Index i = 0; i < static_cast<Index>(out.vertices.size()); ++i)
    {
        if (!used[i]) continue;
        remap[i] = static_cast<Index>(compact.size());
        compact.push_back(out.vertices[i]);
    }
    for (Triangle& t : out.triangles)
        t = {remap[t[0]], remap[t[1]], remap[t[2]]};
    out.vertices = std::move(compact);
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
    for (const auto& [e, count] : uses)
    {
        if (count != 1) continue;
        boundary_edges.insert(e);
        adjacency[e.first].push_back(e.second);
        adjacency[e.second].push_back(e.first);
    }
    for (const Edge& e : boundary_edges)
    {
        if (adjacency[e.first].size() != 2 || adjacency[e.second].size() != 2)
            throw std::runtime_error("Clipped mesh boundary is not a collection of simple loops");
    }

    std::set<Edge> visited;
    std::vector<std::vector<Index>> loops;
    for (const Edge& start_edge : boundary_edges)
    {
        if (visited.count(start_edge)) continue;
        const Index start = start_edge.first;
        Index previous = invalid_index;
        Index current = start;
        std::vector<Index> loop;
        for (std::size_t guard = 0; guard <= boundary_edges.size() + 1; ++guard)
        {
            loop.push_back(current);
            Index next = invalid_index;
            for (Index candidate : adjacency[current])
            {
                const Edge e = make_edge(current, candidate);
                if (visited.count(e)) continue;
                if (candidate != previous || adjacency[current].size() <= 1)
                {
                    next = candidate;
                    break;
                }
            }
            if (next == invalid_index) break;
            visited.insert(make_edge(current, next));
            previous = current;
            current = next;
            if (current == start) break;
        }
        if (loop.size() >= 3 && current == start) loops.push_back(std::move(loop));
    }
    return loops;
}

inline Point loop_centroid(const iso::mc::Mesh& mesh, const std::vector<Index>& loop)
{
    Point c{0.0,0.0,0.0};
    for (Index v : loop)
        for (std::size_t k = 0; k < 3; ++k) c[k] += mesh.vertices[v][k];
    if (!loop.empty())
        for (double& x : c) x /= static_cast<double>(loop.size());
    return c;
}
inline double loop_mean_radius(
    const iso::mc::Mesh& mesh, const std::vector<Index>& loop, const Point& center)
{
    if (loop.empty()) return 0.0;
    double total = 0.0;
    for (Index v : loop) total += distance(mesh.vertices[v], center);
    return total / static_cast<double>(loop.size());
}

/* Split a boundary edge at a point on it; this is the edge-splitting step
 * from the article, performed in the triangle mesh as well as in the loop. */
inline void split_boundary_edge(
    iso::mc::Mesh& mesh,
    std::vector<Index>& loop,
    std::size_t edge_position,
    double fraction)
{
    if (loop.size() < 3 || edge_position >= loop.size())
        throw std::runtime_error("Cannot split invalid boundary loop edge");
    fraction = std::clamp(fraction, 1e-7, 1.0 - 1e-7);
    const Index a = loop[edge_position];
    const Index b = loop[(edge_position + 1) % loop.size()];
    const Index inserted = static_cast<Index>(mesh.vertices.size());
    mesh.vertices.push_back(interpolate(mesh.vertices[a], mesh.vertices[b], fraction));

    bool found = false;
    std::vector<Triangle> rebuilt;
    rebuilt.reserve(mesh.triangles.size() + 2);
    for (const Triangle& t : mesh.triangles)
    {
        bool split = false;
        for (std::size_t k = 0; k < 3; ++k)
        {
            const Index first = t[k];
            const Index second = t[(k+1)%3];
            const Index opposite = t[(k+2)%3];
            if ((first == a && second == b) || (first == b && second == a))
            {
                rebuilt.push_back({first, inserted, opposite});
                rebuilt.push_back({inserted, second, opposite});
                split = true;
                found = true;
                break;
            }
        }
        if (!split) rebuilt.push_back(t);
    }
    if (!found)
    {
        mesh.vertices.pop_back();
        throw std::runtime_error("Boundary edge could not be split in the mesh");
    }
    mesh.triangles = std::move(rebuilt);
    loop.insert(loop.begin() + static_cast<std::ptrdiff_t>(edge_position + 1), inserted);
}

/*
 * Article-inspired feature search: for a vertex on the guide loop, project to
 * every eligible edge on the target loop and select the closest vertex-edge
 * feature. Unlike the old code, we do not fall back to splitting an arbitrary
 * longest edge: no acceptable near feature means stitching must fail safely.
 */
inline bool split_target_loop_near_guide_vertex(
    iso::mc::Mesh& target_mesh,
    std::vector<Index>& target_loop,
    const iso::mc::Mesh& guide_mesh,
    const std::vector<Index>& guide_loop,
    double max_distance)
{
    if (target_loop.size() < 3 || guide_loop.empty()) return false;
    constexpr double endpoint_margin = 0.04;
    double best_d2 = max_distance * max_distance;
    std::size_t best_edge = 0;
    double best_t = 0.5;
    bool found = false;

    for (Index gv : guide_loop)
    {
        const Point& p = guide_mesh.vertices[gv];
        for (std::size_t i = 0; i < target_loop.size(); ++i)
        {
            const Point& a = target_mesh.vertices[target_loop[i]];
            const Point& b = target_mesh.vertices[target_loop[(i+1)%target_loop.size()]];
            const Point ab{b[0]-a[0], b[1]-a[1], b[2]-a[2]};
            const double length2 = dot(ab, ab);
            if (length2 <= 1e-24) continue;
            const Point ap{p[0]-a[0], p[1]-a[1], p[2]-a[2]};
            const double t = dot(ap, ab) / length2;
            // Endpoint matches are handled as vertex-vertex candidates. Do not
            // create tiny edges (the T-junction/degenerate case in the article).
            if (t <= endpoint_margin || t >= 1.0 - endpoint_margin) continue;
            const Point projected = interpolate(a, b, t);
            const double d2 = squared_distance(p, projected);
            if (d2 <= best_d2)
            {
                best_d2 = d2;
                best_edge = i;
                best_t = t;
                found = true;
            }
        }
    }
    if (found)
    {
        split_boundary_edge(target_mesh, target_loop, best_edge, best_t);
        return true;
    }

    // A very fine loop may have no vertex whose projection falls safely inside
    // a coarse edge (the article's vertex-edge candidate). In that case look
    // for the closest edge-to-edge feature: split the coarse edge midpoint only
    // when it is close to an edge of the guide loop. This is deliberately not
    // the old unconditional "split the longest edge" fallback.
    best_d2 = max_distance * max_distance;
    found = false;
    best_t = 0.5;
    for (std::size_t i = 0; i < target_loop.size(); ++i)
    {
        const Point& a = target_mesh.vertices[target_loop[i]];
        const Point& b = target_mesh.vertices[target_loop[(i+1)%target_loop.size()]];
        const Point midpoint_ab = midpoint(a, b);
        for (std::size_t j = 0; j < guide_loop.size(); ++j)
        {
            const Point& c = guide_mesh.vertices[guide_loop[j]];
            const Point& d = guide_mesh.vertices[guide_loop[(j+1)%guide_loop.size()]];
            const Point cd{d[0]-c[0], d[1]-c[1], d[2]-c[2]};
            const double length2 = dot(cd, cd);
            if (length2 <= 1e-24) continue;
            const Point cm{midpoint_ab[0]-c[0], midpoint_ab[1]-c[1], midpoint_ab[2]-c[2]};
            const double u = std::clamp(dot(cm, cd) / length2, 0.0, 1.0);
            const Point nearest = interpolate(c, d, u);
            const double d2 = squared_distance(midpoint_ab, nearest);
            if (d2 <= best_d2)
            {
                best_d2 = d2;
                best_edge = i;
                found = true;
            }
        }
    }
    if (!found) return false;
    split_boundary_edge(target_mesh, target_loop, best_edge, 0.5);
    return true;
}

inline bool equalize_loop_sizes(
    iso::mc::Mesh& mesh_a, std::vector<Index>& loop_a,
    iso::mc::Mesh& mesh_b, std::vector<Index>& loop_b,
    double max_distance)
{
    const std::size_t target_size = std::max(loop_a.size(), loop_b.size());
    // Only split the lower-resolution loop. Every split must be supported by
    // a close projection from the opposite loop; arbitrary fallback splits are forbidden.
    while (loop_a.size() < target_size)
    {
        if (!split_target_loop_near_guide_vertex(
                mesh_a, loop_a, mesh_b, loop_b, max_distance)) return false;
    }
    while (loop_b.size() < target_size)
    {
        if (!split_target_loop_near_guide_vertex(
                mesh_b, loop_b, mesh_a, loop_a, max_distance)) return false;
    }
    return loop_a.size() == loop_b.size();
}

struct LoopMatch
{
    std::size_t shift = 0;
    bool reversed = false;
    double total_squared_distance = std::numeric_limits<double>::infinity();
    double max_distance = std::numeric_limits<double>::infinity();
};

inline bool best_loop_match(
    const iso::mc::Mesh& a_mesh, const std::vector<Index>& a,
    const iso::mc::Mesh& b_mesh, const std::vector<Index>& b,
    double max_distance, LoopMatch& result)
{
    if (a.size() != b.size() || a.empty()) return false;
    const double max_d2 = max_distance * max_distance;
    LoopMatch best;
    for (const bool reversed : {false, true})
    {
        for (std::size_t shift = 0; shift < b.size(); ++shift)
        {
            double sum = 0.0;
            double max_d2_this = 0.0;
            bool acceptable = true;
            for (std::size_t i = 0; i < a.size(); ++i)
            {
                const std::size_t j = reversed
                    ? (shift + b.size() - i) % b.size()
                    : (shift + i) % b.size();
                const double d2 = squared_distance(a_mesh.vertices[a[i]], b_mesh.vertices[b[j]]);
                sum += d2;
                max_d2_this = std::max(max_d2_this, d2);
                if (d2 > max_d2) { acceptable = false; break; }
            }
            if (acceptable && sum < best.total_squared_distance)
            {
                best.shift = shift;
                best.reversed = reversed;
                best.total_squared_distance = sum;
                best.max_distance = std::sqrt(max_d2_this);
            }
        }
    }
    if (!std::isfinite(best.total_squared_distance)) return false;
    result = best;
    return true;
}

inline void orient_consistently_from_global(
    iso::mc::Mesh& mesh, std::size_t global_triangle_count)
{
    struct Incidence { Index triangle = 0; bool forward = false; };
    std::map<Edge, std::vector<Incidence>> incidence;
    for (Index ti = 0; ti < static_cast<Index>(mesh.triangles.size()); ++ti)
    {
        const Triangle& t = mesh.triangles[ti];
        for (const auto& [a,b] : std::array<std::pair<Index,Index>,3>{{
                 {t[0],t[1]}, {t[1],t[2]}, {t[2],t[0]}}})
        {
            const Edge e = make_edge(a,b);
            incidence[e].push_back({ti, a == e.first && b == e.second});
        }
    }
    for (const auto& [e, uses] : incidence)
    {
        (void)e;
        if (uses.size() > 2) throw std::runtime_error("Cannot orient stitched mesh: non-manifold edge");
    }

    std::vector<int> flip(mesh.triangles.size(), -1);
    std::queue<Index> q;
    auto flood = [&]()
    {
        while (!q.empty())
        {
            const Index current = q.front(); q.pop();
            const Triangle& t = mesh.triangles[current];
            for (const auto& [a,b] : std::array<std::pair<Index,Index>,3>{{
                     {t[0],t[1]}, {t[1],t[2]}, {t[2],t[0]}}})
            {
                auto it = incidence.find(make_edge(a,b));
                if (it == incidence.end() || it->second.size() != 2) continue;
                const Incidence* other = it->second[0].triangle == current
                    ? &it->second[1] : &it->second[0];
                const bool current_forward = (a < b);
                const int required = (current_forward == other->forward)
                    ? 1 - flip[current] : flip[current];
                if (flip[other->triangle] == -1)
                {
                    flip[other->triangle] = required;
                    q.push(other->triangle);
                }
                else if (flip[other->triangle] != required)
                    throw std::runtime_error("Cannot consistently orient stitched mesh");
            }
        }
    };

    const Index global_count = static_cast<Index>(std::min(global_triangle_count, mesh.triangles.size()));
    for (Index i = 0; i < global_count; ++i)
    {
        if (flip[i] != -1) continue;
        flip[i] = 0; q.push(i); flood();
    }
    for (Index i = 0; i < static_cast<Index>(mesh.triangles.size()); ++i)
    {
        if (flip[i] != -1) continue;
        flip[i] = 0; q.push(i); flood();
    }
    for (Index i = 0; i < static_cast<Index>(mesh.triangles.size()); ++i)
        if (flip[i] == 1) std::swap(mesh.triangles[i][1], mesh.triangles[i][2]);
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
        bool use_sphere_seam = true;
        // 0 selects a radius-relative default of 0.30 * radius.
        // Set explicitly when running controlled experiments.
        double max_stitch_distance = 0.0;
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
        StitchResult result;
        result.mesh = global;
        result.topology = detail::topology(global);

        if (!(options.radius > 0.0) || !std::isfinite(options.radius))
            throw std::invalid_argument("LocalVoxelStitcher radius must be finite and positive");
        if (!options.use_sphere_seam)
            throw std::invalid_argument("LocalVoxelStitcher currently requires a sphere seam");
        if (options.max_stitch_distance < 0.0 || !std::isfinite(options.max_stitch_distance))
            throw std::invalid_argument("max_stitch_distance must be finite and non-negative");
        if (global.triangles.empty() || local_patch.triangles.empty())
        {
            result.failure_reason = "Global mesh or local patch is empty";
            return result;
        }

        const double max_distance = options.max_stitch_distance > 0.0
            ? options.max_stitch_distance : 0.30 * options.radius;
        const double radius2 = options.radius * options.radius;
        // Much narrower than the old 0.5 * radius band: outer domain boundaries
        // must not accidentally be treated as the spherical stitching seam.
        const double radial_band = std::max(0.20 * options.radius, 1e-8);

        iso::mc::Mesh outside;
        iso::mc::Mesh inside;
        try
        {
            outside = detail::clip_by_sphere(global, center, radius2, false);
            inside = detail::clip_by_sphere(local_patch, center, radius2, true);
        }
        catch (const std::exception& e)
        {
            result.failure_reason = std::string("Sphere clipping failed: ") + e.what();
            return result;
        }
        if (outside.triangles.empty() || inside.triangles.empty())
        {
            result.failure_reason = "Sphere clipping produced an empty mesh";
            return result;
        }

        for (const Triangle& t : global.triangles)
        {
            const detail::Point c = detail::centroid(
                global.vertices[t[0]], global.vertices[t[1]], global.vertices[t[2]]);
            if (detail::squared_distance(c, center) < radius2)
                ++result.global_triangles_removed;
        }

        std::vector<std::vector<Index>> global_loops, local_loops;
        try
        {
            global_loops = detail::boundary_loops(outside);
            local_loops = detail::boundary_loops(inside);
        }
        catch (const std::exception& e)
        {
            result.failure_reason = std::string("Boundary extraction failed: ") + e.what();
            return result;
        }

        const auto find_seams = [&](const iso::mc::Mesh& mesh,
                                    const std::vector<std::vector<Index>>& loops)
        {
            std::vector<std::size_t> seam_ids;
            for (std::size_t i = 0; i < loops.size(); ++i)
            {
                const double mean_r = detail::loop_mean_radius(mesh, loops[i], center);
                if (std::abs(mean_r - options.radius) <= radial_band)
                    seam_ids.push_back(i);
            }
            return seam_ids;
        };
        const std::vector<std::size_t> global_seams = find_seams(outside, global_loops);
        const std::vector<std::size_t> local_seams = find_seams(inside, local_loops);
        if (global_seams.empty() || global_seams.size() != local_seams.size())
        {
            result.failure_reason = "Global/local seam-loop counts differ or no spherical seam was found";
            return result;
        }

        struct MatchedLoops
        {
            std::vector<Index> global_loop;
            std::vector<Index> local_loop;
            detail::LoopMatch match;
        };
        std::vector<MatchedLoops> matched;
        std::set<std::size_t> used_global;
        try
        {
            for (std::size_t local_id : local_seams)
            {
                const auto& local_loop = local_loops[local_id];
                const detail::Point lc = detail::loop_centroid(inside, local_loop);
                const double lr = detail::loop_mean_radius(inside, local_loop, center);
                double best_score = std::numeric_limits<double>::infinity();
                std::size_t best_global_id = std::numeric_limits<std::size_t>::max();
                for (std::size_t global_id : global_seams)
                {
                    if (used_global.count(global_id)) continue;
                    const auto& global_loop = global_loops[global_id];
                    const detail::Point gc = detail::loop_centroid(outside, global_loop);
                    const double gr = detail::loop_mean_radius(outside, global_loop, center);
                    const double score = detail::distance(lc, gc) + 0.5 * std::abs(lr-gr);
                    if (score < best_score) { best_score = score; best_global_id = global_id; }
                }
                if (best_global_id == std::numeric_limits<std::size_t>::max())
                {
                    result.failure_reason = "Could not pair a local seam loop with a global seam loop";
                    return result;
                }
                used_global.insert(best_global_id);

                auto gloop = global_loops[best_global_id];
                auto lloop = local_loop;
                // The article's feature-based idea adapted to two cut loops:
                // split only at the closest admissible vertex-to-edge projection,
                // and refuse to create a split when the gap exceeds the threshold.
                if (!detail::equalize_loop_sizes(outside, gloop, inside, lloop, max_distance))
                {
                    result.failure_reason = "No nearby vertex-edge feature could equalize seam loops within max_stitch_distance";
                    return result;
                }
                detail::LoopMatch loop_match;
                if (!detail::best_loop_match(outside, gloop, inside, lloop,
                                             max_distance, loop_match))
                {
                    result.failure_reason = "No cyclic seam correspondence satisfies max_stitch_distance";
                    return result;
                }
                matched.push_back({std::move(gloop), std::move(lloop), loop_match});
            }
        }
        catch (const std::exception& e)
        {
            result.failure_reason = std::string("Feature-based seam matching failed: ") + e.what();
            return result;
        }

        // Build a candidate result in isolation. Until topology checks pass, the
        // returned mesh remains the original global mesh (no partial patch leak).
        iso::mc::Mesh candidate = std::move(outside);
        const std::size_t global_triangle_count = candidate.triangles.size();
        std::vector<Index> local_remap(inside.vertices.size(), detail::invalid_index);
        std::size_t seam_loop_count = 0;
        std::size_t seam_vertex_count = 0;
        double maximum_pair_distance = 0.0;

        for (const MatchedLoops& loops : matched)
        {
            ++seam_loop_count;
            const std::size_t n = loops.global_loop.size();
            seam_vertex_count += n;
            for (std::size_t i = 0; i < n; ++i)
            {
                const std::size_t j = loops.match.reversed
                    ? (loops.match.shift + n - i) % n
                    : (loops.match.shift + i) % n;
                const Index gv = loops.global_loop[i];
                const Index lv = loops.local_loop[j];
                if (local_remap[lv] != detail::invalid_index && local_remap[lv] != gv)
                {
                    result.failure_reason = "A local seam vertex maps to more than one global vertex";
                    return result;
                }
                const detail::Point gp = candidate.vertices[gv];
                const detail::Point lp = inside.vertices[lv];
                const double gap = detail::distance(gp, lp);
                if (gap > max_distance)
                {
                    result.failure_reason = "Seam vertex pair exceeded max_stitch_distance";
                    return result;
                }
                maximum_pair_distance = std::max(maximum_pair_distance, gap);
                // Symmetric merge as in the article's average-position collapse.
                candidate.vertices[gv] = detail::midpoint(gp, lp);
                local_remap[lv] = gv;
            }
        }

        for (Index lv = 0; lv < static_cast<Index>(inside.vertices.size()); ++lv)
        {
            if (local_remap[lv] != detail::invalid_index) continue;
            local_remap[lv] = static_cast<Index>(candidate.vertices.size());
            candidate.vertices.push_back(inside.vertices[lv]);
        }
        std::size_t inserted_triangles = 0;
        for (const Triangle& t : inside.triangles)
        {
            const Index a = local_remap[t[0]], b = local_remap[t[1]], c = local_remap[t[2]];
            if (a == b || b == c || c == a) continue;
            candidate.triangles.push_back({a,b,c});
            ++inserted_triangles;
        }

        try
        {
            detail::orient_consistently_from_global(candidate, global_triangle_count);
            const TopologyStatistics candidate_topology = detail::topology(candidate);
            const TopologyStatistics original_topology = detail::topology(global);
            if (candidate_topology.nonmanifold_edges != 0)
            {
                result.failure_reason = "Candidate stitch creates non-manifold edges";
                return result;
            }
            if (candidate_topology.boundary_edges != original_topology.boundary_edges)
            {
                result.failure_reason = "Candidate stitch changes the original domain-boundary edge count";
                return result;
            }
            if (inserted_triangles == 0 || seam_loop_count != local_seams.size())
            {
                result.failure_reason = "Candidate stitch did not insert every local seam patch";
                return result;
            }
            result.mesh = std::move(candidate);
            result.stitched = true;
            result.local_triangles_inserted = inserted_triangles;
            result.seam_loops = seam_loop_count;
            result.seam_vertices = seam_vertex_count;
            result.seam_triangles = 0; // A welded seam needs no bridge strip.
            result.max_seam_vertex_distance = maximum_pair_distance;
            result.topology = candidate_topology;
            result.failure_reason.clear();
            return result;
        }
        catch (const std::exception& e)
        {
            result.failure_reason = std::string("Candidate topology validation failed: ") + e.what();
            return result;
        }
    }
};

} // namespace iso::stitch
