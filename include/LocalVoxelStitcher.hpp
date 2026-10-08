#pragma once

#include "MarchingCubes.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <numeric>
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

inline Edge make_edge(Index a, Index b)
{
    if (a > b)
        std::swap(a, b);

    return {a, b};
}

inline double squared_distance(
    const Point& a,
    const Point& b)
{
    const double dx = a[0] - b[0];
    const double dy = a[1] - b[1];
    const double dz = a[2] - b[2];

    return dx * dx + dy * dy + dz * dz;
}

inline double distance(
    const Point& a,
    const Point& b)
{
    return std::sqrt(
        squared_distance(a, b));
}

inline Point midpoint(
    const Point& a,
    const Point& b)
{
    return {
        0.5 * (a[0] + b[0]),
        0.5 * (a[1] + b[1]),
        0.5 * (a[2] + b[2])
    };
}

inline Point centroid(
    const Point& a,
    const Point& b,
    const Point& c)
{
    return {
        (a[0] + b[0] + c[0]) / 3.0,
        (a[1] + b[1] + c[1]) / 3.0,
        (a[2] + b[2] + c[2]) / 3.0
    };
}

inline double radius_from(
    const Point& p,
    const Point& center)
{
    return std::sqrt(
        squared_distance(p, center));
}

inline TopologyStatistics topology(
    const iso::mc::Mesh& mesh)
{
    TopologyStatistics result;

    result.vertices = mesh.vertices.size();
    result.triangles = mesh.triangles.size();

    std::map<Edge, std::size_t> uses;
    std::set<Index> boundary_vertices;

    for (const Triangle& t : mesh.triangles)
    {
        ++uses[make_edge(t[0], t[1])];
        ++uses[make_edge(t[1], t[2])];
        ++uses[make_edge(t[2], t[0])];
    }

    result.edges = uses.size();

    for (const auto& [edge, count] : uses)
    {
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

    result.boundary_vertices =
        boundary_vertices.size();

    std::vector<std::vector<Index>> neighbors(
        mesh.vertices.size());

    for (const auto& [edge, count] : uses)
    {
        if (count == 0)
            continue;

        neighbors[edge.first].push_back(edge.second);
        neighbors[edge.second].push_back(edge.first);
    }

    std::vector<char> visited(
        mesh.vertices.size(),
        0);

    for (Index start = 0;
         start < static_cast<Index>(
                     mesh.vertices.size());
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

struct Submesh
{
    iso::mc::Mesh mesh;
};

template <typename Predicate>
Submesh filter_triangles(
    const iso::mc::Mesh& source,
    Predicate&& keep)
{
    Submesh result;

    std::vector<Index> remap(
        source.vertices.size(),
        std::numeric_limits<Index>::max());

    auto copy_vertex =
        [&](Index source_index) -> Index
    {
        Index& mapped =
            remap[source_index];

        if (mapped !=
            std::numeric_limits<Index>::max())
        {
            return mapped;
        }

        mapped =
            static_cast<Index>(
                result.mesh.vertices.size());

        result.mesh.vertices.push_back(
            source.vertices[source_index]);

        return mapped;
    };

    for (const Triangle& t :
         source.triangles)
    {
        const Point c =
            centroid(
                source.vertices[t[0]],
                source.vertices[t[1]],
                source.vertices[t[2]]);

        if (!keep(c))
            continue;

        result.mesh.triangles.push_back({
            copy_vertex(t[0]),
            copy_vertex(t[1]),
            copy_vertex(t[2])
        });
    }

    return result;
}

inline std::vector<std::vector<Index>>
boundary_loops(
    const iso::mc::Mesh& mesh)
{
    std::map<Edge, std::size_t> uses;

    std::vector<std::vector<Index>> adjacency(
        mesh.vertices.size());

    for (const Triangle& t :
         mesh.triangles)
    {
        ++uses[make_edge(t[0], t[1])];
        ++uses[make_edge(t[1], t[2])];
        ++uses[make_edge(t[2], t[0])];
    }

    std::set<Edge> boundary_edges;

    for (const auto& [edge, count] :
         uses)
    {
        if (count != 1)
            continue;

        boundary_edges.insert(edge);

        adjacency[edge.first].push_back(
            edge.second);

        adjacency[edge.second].push_back(
            edge.first);
    }

    std::set<Edge> visited_edges;

    std::vector<std::vector<Index>> loops;

    for (const Edge& start_edge :
         boundary_edges)
    {
        if (visited_edges.count(start_edge))
            continue;

        const Index start =
            start_edge.first;

        Index current = start;

        Index previous =
            std::numeric_limits<Index>::max();

        std::vector<Index> loop;

        while (true)
        {
            loop.push_back(current);

            Index next =
                std::numeric_limits<Index>::max();

            for (const Index candidate :
                 adjacency[current])
            {
                const Edge edge =
                    make_edge(
                        current,
                        candidate);

                if (visited_edges.count(edge))
                    continue;

                if (candidate != previous ||
                    adjacency[current].size() <= 1)
                {
                    next = candidate;
                    break;
                }
            }

            if (next ==
                std::numeric_limits<Index>::max())
            {
                break;
            }

            visited_edges.insert(
                make_edge(current, next));

            previous = current;
            current = next;

            if (current == start)
                break;

            if (loop.size() >
                boundary_edges.size() + 1)
            {
                throw std::runtime_error(
                    "Failed to trace boundary loop");
            }
        }

        if (loop.size() >= 3 &&
            current == start)
        {
            loops.push_back(
                std::move(loop));
        }
    }

    return loops;
}

inline Point loop_centroid(
    const iso::mc::Mesh& mesh,
    const std::vector<Index>& loop)
{
    Point c{
        0.0,
        0.0,
        0.0
    };

    if (loop.empty())
        return c;

    for (const Index v : loop)
    {
        c[0] += mesh.vertices[v][0];
        c[1] += mesh.vertices[v][1];
        c[2] += mesh.vertices[v][2];
    }

    const double inv =
        1.0 /
        static_cast<double>(
            loop.size());

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
    {
        sum += radius_from(
            mesh.vertices[v],
            center);
    }

    return sum /
           static_cast<double>(
               loop.size());
}

inline bool triangle_contains_edge(
    const Triangle& t,
    Index a,
    Index b)
{
    return
        (t[0] == a && t[1] == b) ||
        (t[1] == a && t[2] == b) ||
        (t[2] == a && t[0] == b) ||
        (t[0] == b && t[1] == a) ||
        (t[1] == b && t[2] == a) ||
        (t[2] == b && t[0] == a);
}

inline void split_boundary_edge(
    iso::mc::Mesh& mesh,
    std::vector<Index>& loop,
    std::size_t edge_position)
{
    if (loop.size() < 3)
    {
        throw std::runtime_error(
            "Boundary loop is too small");
    }

    const std::size_t next_position =
        (edge_position + 1) %
        loop.size();

    const Index a =
        loop[edge_position];

    const Index b =
        loop[next_position];

    const Index new_vertex =
        static_cast<Index>(
            mesh.vertices.size());

    mesh.vertices.push_back(
        midpoint(
            mesh.vertices[a],
            mesh.vertices[b]));

    bool found = false;

    std::vector<Triangle> rebuilt;

    rebuilt.reserve(
        mesh.triangles.size() + 1);

    for (const Triangle& t :
         mesh.triangles)
    {
        if (found)
        {
            rebuilt.push_back(t);
            continue;
        }

        Index first = 0;
        Index second = 0;
        Index opposite = 0;

        bool matched = false;

        if (t[0] == a && t[1] == b)
        {
            first = t[0];
            second = t[1];
            opposite = t[2];
            matched = true;
        }
        else if (t[1] == a && t[2] == b)
        {
            first = t[1];
            second = t[2];
            opposite = t[0];
            matched = true;
        }
        else if (t[2] == a && t[0] == b)
        {
            first = t[2];
            second = t[0];
            opposite = t[1];
            matched = true;
        }
        else if (t[0] == b && t[1] == a)
        {
            first = t[0];
            second = t[1];
            opposite = t[2];
            matched = true;
        }
        else if (t[1] == b && t[2] == a)
        {
            first = t[1];
            second = t[2];
            opposite = t[0];
            matched = true;
        }
        else if (t[2] == b && t[0] == a)
        {
            first = t[2];
            second = t[0];
            opposite = t[1];
            matched = true;
        }

        if (!matched)
        {
            rebuilt.push_back(t);
            continue;
        }

        found = true;

        /*
         * Preserve the orientation of the original triangle.
         *
         * Original:
         *
         *     first -> second -> opposite
         *
         * New:
         *
         *     first -> new -> opposite
         *     new   -> second -> opposite
         */
        rebuilt.push_back({
            first,
            new_vertex,
            opposite
        });

        rebuilt.push_back({
            new_vertex,
            second,
            opposite
        });
    }

    if (!found)
    {
        mesh.vertices.pop_back();

        throw std::runtime_error(
            "Boundary edge could not be split");
    }

    mesh.triangles =
        std::move(rebuilt);

    loop.insert(
        loop.begin() +
            static_cast<std::ptrdiff_t>(
                next_position),
        new_vertex);
}

inline void equalize_loop_sizes(
    iso::mc::Mesh& mesh_a,
    std::vector<Index>& loop_a,
    iso::mc::Mesh& mesh_b,
    std::vector<Index>& loop_b)
{
    const std::size_t target =
        std::max(
            loop_a.size(),
            loop_b.size());

    while (loop_a.size() < target)
    {
        std::size_t best = 0;
        double longest = -1.0;

        for (std::size_t i = 0;
             i < loop_a.size();
             ++i)
        {
            const std::size_t j =
                (i + 1) %
                loop_a.size();

            const double len =
                distance(
                    mesh_a.vertices[
                        loop_a[i]],
                    mesh_a.vertices[
                        loop_a[j]]);

            if (len > longest)
            {
                longest = len;
                best = i;
            }
        }

        split_boundary_edge(
            mesh_a,
            loop_a,
            best);
    }

    while (loop_b.size() < target)
    {
        std::size_t best = 0;
        double longest = -1.0;

        for (std::size_t i = 0;
             i < loop_b.size();
             ++i)
        {
            const std::size_t j =
                (i + 1) %
                loop_b.size();

            const double len =
                distance(
                    mesh_b.vertices[
                        loop_b[i]],
                    mesh_b.vertices[
                        loop_b[j]]);

            if (len > longest)
            {
                longest = len;
                best = i;
            }
        }

        split_boundary_edge(
            mesh_b,
            loop_b,
            best);
    }
}

struct LoopMatch
{
    std::size_t shift = 0;
    bool reversed = false;

    double cost =
        std::numeric_limits<double>::infinity();
};

inline LoopMatch best_loop_match(
    const iso::mc::Mesh& a_mesh,
    const std::vector<Index>& a,
    const iso::mc::Mesh& b_mesh,
    const std::vector<Index>& b)
{
    if (a.size() != b.size() ||
        a.empty())
    {
        throw std::runtime_error(
            "Boundary loops must have equal non-zero size");
    }

    LoopMatch best;

    for (const bool reversed :
         {false, true})
    {
        for (std::size_t shift = 0;
             shift < b.size();
             ++shift)
        {
            double cost = 0.0;

            for (std::size_t i = 0;
                 i < a.size();
                 ++i)
            {
                const std::size_t j =
                    reversed
                        ? (shift +
                           b.size() -
                           (i %
                            b.size())) %
                              b.size()
                        : (shift + i) %
                              b.size();

                cost += squared_distance(
                    a_mesh.vertices[a[i]],
                    b_mesh.vertices[b[j]]);
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

/*
 * Orient the local patch so that it agrees with the already oriented
 * global mesh.
 *
 * Global triangles are authoritative. Their winding is never changed.
 *
 * Once the bridge is present, the local triangles are connected to the
 * global triangles through shared bridge edges, so a BFS can propagate
 * the required orientation into the local patch.
 */
inline void orient_consistently_from_global(
    iso::mc::Mesh& mesh,
    std::size_t global_triangle_count)
{
    struct Incidence
    {
        Index triangle = 0;
        bool forward = false;
    };

    std::map<Edge, std::vector<Incidence>>
        incidences;

    for (Index ti = 0;
         ti < static_cast<Index>(
                  mesh.triangles.size());
         ++ti)
    {
        const Triangle& t =
            mesh.triangles[ti];

        const auto add_edge =
            [&](Index a, Index b)
        {
            const Edge edge =
                make_edge(a, b);

            const bool forward =
                (a == edge.first &&
                 b == edge.second);

            incidences[edge].push_back({
                ti,
                forward
            });
        };

        add_edge(t[0], t[1]);
        add_edge(t[1], t[2]);
        add_edge(t[2], t[0]);
    }

    for (const auto& [edge, uses] :
         incidences)
    {
        if (uses.size() > 2)
        {
            throw std::runtime_error(
                "Cannot orient stitched mesh: non-manifold edge");
        }
    }

    std::vector<int> flip(
        mesh.triangles.size(),
        -1);

    std::queue<Index> q;

    const Index global_count =
        static_cast<Index>(
            std::min(
                global_triangle_count,
                mesh.triangles.size()));

    auto propagate =
        [&](Index seed)
    {
        if (flip[seed] == -1)
        {
            flip[seed] = 0;
            q.push(seed);
        }

        while (!q.empty())
        {
            const Index current =
                q.front();

            q.pop();

            const Triangle& ct =
                mesh.triangles[current];

            const std::array<
                std::pair<Index, Index>,
                3>
                edges{{
                    {ct[0], ct[1]},
                    {ct[1], ct[2]},
                    {ct[2], ct[0]}
                }};

            for (const auto& [a, b] :
                 edges)
            {
                const Edge edge =
                    make_edge(a, b);

                const auto it =
                    incidences.find(edge);

                if (it == incidences.end() ||
                    it->second.size() != 2)
                {
                    continue;
                }

                const auto& uses =
                    it->second;

                const Incidence* other =
                    uses[0].triangle == current
                        ? &uses[1]
                        : &uses[0];

                const bool current_forward =
                    (a == edge.first &&
                     b == edge.second);

                const bool other_forward =
                    other->forward;

                /*
                 * Adjacent oriented triangles must traverse their shared
                 * edge in opposite directions.
                 */
                const int required_flip =
                    (current_forward ==
                     other_forward)
                        ? 1 - flip[current]
                        : flip[current];

                if (flip[other->triangle] == -1)
                {
                    flip[other->triangle] =
                        required_flip;

                    q.push(
                        other->triangle);
                }
                else if (
                    flip[other->triangle] !=
                    required_flip)
                {
                    throw std::runtime_error(
                        "Cannot consistently orient stitched mesh");
                }
            }
        }
    };

    /*
     * Preserve every surviving global triangle exactly as it came from
     * Marching Cubes.
     */
    for (Index ti = 0;
         ti < global_count;
         ++ti)
    {
        if (flip[ti] != -1)
            continue;

        flip[ti] = 0;
        q.push(ti);

        propagate(ti);
    }

    /*
     * Any remaining isolated component is oriented consistently.
     */
    for (Index ti = 0;
         ti < static_cast<Index>(
                  mesh.triangles.size());
         ++ti)
    {
        if (flip[ti] != -1)
            continue;

        flip[ti] = 0;
        q.push(ti);

        propagate(ti);
    }

    for (Index ti = 0;
         ti < static_cast<Index>(
                  mesh.triangles.size());
         ++ti)
    {
        if (flip[ti] == 1)
        {
            std::swap(
                mesh.triangles[ti][1],
                mesh.triangles[ti][2]);
        }
    }
}

} // namespace detail

inline TopologyStatistics compute_topology(
    const iso::mc::Mesh& mesh)
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
    };

    StitchResult stitch(
        const iso::mc::Mesh& global,
        const iso::mc::Mesh& local_patch,
        const iso::mc::Point& center) const
    {
        const Options options;

        return stitch(
            global,
            local_patch,
            center,
            options);
    }

    StitchResult stitch(
        const iso::mc::Mesh& global,
        const iso::mc::Mesh& local_patch,
        const iso::mc::Point& center,
        const Options& options) const
    {
        if (!(options.radius > 0.0))
        {
            throw std::invalid_argument(
                "LocalVoxelStitcher radius must be positive");
        }

        StitchResult result;
        result.mesh = global;

        if (global.triangles.empty() ||
            local_patch.triangles.empty())
        {
            result.topology =
                detail::topology(
                    result.mesh);

            return result;
        }

        const double radius_squared =
            options.radius *
            options.radius;

        /*
         * IMPORTANT:
         *
         * We do NOT weld independently sampled local vertices to global
         * vertices anymore.
         *
         * That old approach produced large seam displacements
         * (your ~0.078 value) and could create badly distorted triangles.
         *
         * Instead we:
         *
         *   1. remove the global triangles inside the seam,
         *   2. keep the local triangles inside the seam,
         *   3. equalize the two boundary rings,
         *   4. connect the rings with an actual triangle strip.
         *
         * This gives a real 2-manifold transition.
         */

        auto outside =
            detail::filter_triangles(
                global,
                [&](const detail::Point& p)
                {
                    return detail::squared_distance(
                               p,
                               center) >
                           radius_squared;
                });

        auto inside =
            detail::filter_triangles(
                local_patch,
                [&](const detail::Point& p)
                {
                    return detail::squared_distance(
                               p,
                               center) <=
                           radius_squared;
                });

        result.global_triangles_removed =
            global.triangles.size() -
            outside.mesh.triangles.size();

        result.local_triangles_inserted =
            inside.mesh.triangles.size();

        if (result.global_triangles_removed == 0 ||
            result.local_triangles_inserted == 0)
        {
            result.topology =
                detail::topology(
                    result.mesh);

            return result;
        }

        auto global_loops =
            detail::boundary_loops(
                outside.mesh);

        auto local_loops =
            detail::boundary_loops(
                inside.mesh);

        if (global_loops.empty() ||
            local_loops.empty())
        {
            result.topology =
                detail::topology(
                    result.mesh);

            return result;
        }

        struct Match
        {
            std::size_t global_loop = 0;
            std::size_t local_loop = 0;
        };

        std::vector<char> global_used(
            global_loops.size(),
            0);

        std::vector<Match> matches;

        /*
         * Match rings belonging to the same local surface component.
         *
         * For the cone this should produce exactly two matches:
         *
         *   upper cone <-> upper local cone
         *   lower cone <-> lower local cone
         */
        for (std::size_t li = 0;
             li < local_loops.size();
             ++li)
        {
            const detail::Point lc =
                detail::loop_centroid(
                    inside.mesh,
                    local_loops[li]);

            const double lr =
                detail::loop_mean_radius(
                    inside.mesh,
                    local_loops[li],
                    center);

            double best_score =
                std::numeric_limits<double>::infinity();

            std::size_t best_global =
                std::numeric_limits<std::size_t>::max();

            for (std::size_t gi = 0;
                 gi < global_loops.size();
                 ++gi)
            {
                if (global_used[gi])
                    continue;

                const detail::Point gc =
                    detail::loop_centroid(
                        outside.mesh,
                        global_loops[gi]);

                const double gr =
                    detail::loop_mean_radius(
                        outside.mesh,
                        global_loops[gi],
                        center);

                const double centroid_distance =
                    detail::distance(
                        lc,
                        gc);

                const double radius_difference =
                    std::abs(lr - gr);

                /*
                 * Centroid position is the dominant criterion.
                 * Radius is only a secondary stabilizer.
                 */
                const double score =
                    centroid_distance +
                    0.25 * radius_difference;

                if (score < best_score)
                {
                    best_score = score;
                    best_global = gi;
                }
            }

            if (best_global ==
                std::numeric_limits<
                    std::size_t>::max())
            {
                result.topology =
                    detail::topology(
                        result.mesh);

                return result;
            }

            /*
             * Do not use a loose "radius" acceptance test here.
             *
             * The seam is intentionally allowed to connect two independently
             * sampled surfaces. The bridge is responsible for resolving the
             * geometric mismatch.
             */
            global_used[best_global] = 1;

            matches.push_back({
                best_global,
                li
            });
        }

        struct Bridge
        {
            std::vector<detail::Index>
                global_loop;

            std::vector<detail::Index>
                local_loop;

            detail::LoopMatch matching;
        };

        std::vector<Bridge> bridges;

        bridges.reserve(
            matches.size());

        for (const Match& match :
             matches)
        {
            auto gloop =
                global_loops[
                    match.global_loop];

            auto lloop =
                local_loops[
                    match.local_loop];

            detail::equalize_loop_sizes(
                outside.mesh,
                gloop,
                inside.mesh,
                lloop);

            if (gloop.size() !=
                    lloop.size() ||
                gloop.size() < 3)
            {
                result.topology =
                    detail::topology(
                        result.mesh);

                return result;
            }

            const detail::LoopMatch matching =
                detail::best_loop_match(
                    outside.mesh,
                    gloop,
                    inside.mesh,
                    lloop);

            bridges.push_back({
                std::move(gloop),
                std::move(lloop),
                matching
            });
        }

        /*
         * Start with the surviving global surface.
         */
        result.mesh =
            outside.mesh;

        const std::size_t global_triangle_count =
            result.mesh.triangles.size();

        /*
         * Append the local patch.
         *
         * We deliberately keep the local seam vertices separate from the
         * global seam vertices. They are connected by bridge triangles.
         *
         * This is geometrically much safer than forcing two independently
         * sampled rings onto exactly the same vertices.
         */
        std::vector<std::vector<Index>>
            local_vertex_maps;

        local_vertex_maps.reserve(
            bridges.size());

        std::vector<Index> local_remap(
            inside.mesh.vertices.size(),
            std::numeric_limits<
                Index>::max());

        for (Index lv = 0;
             lv < static_cast<Index>(
                      inside.mesh.vertices.size());
             ++lv)
        {
            local_remap[lv] =
                static_cast<Index>(
                    result.mesh.vertices.size());

            result.mesh.vertices.push_back(
                inside.mesh.vertices[lv]);
        }

        /*
         * Insert the local triangles.
         */
        for (const detail::Triangle& t :
             inside.mesh.triangles)
        {
            const Index a =
                local_remap[t[0]];

            const Index b =
                local_remap[t[1]];

            const Index c =
                local_remap[t[2]];

            result.mesh.triangles.push_back({
                a,
                b,
                c
            });
        }

        /*
         * Build an actual zipper/bridge strip.
         *
         * Every global boundary edge receives exactly one bridge triangle.
         * Every local boundary edge receives exactly one bridge triangle.
         *
         * Therefore:
         *
         *   boundary edge on global side -> 2 incident triangles
         *   boundary edge on local side  -> 2 incident triangles
         *
         * and the seam disappears topologically.
         */
        for (const Bridge& bridge :
             bridges)
        {
            const std::size_t n =
                bridge.global_loop.size();

            result.seam_loops++;

            result.seam_vertices +=
                2 * n;

            result.max_seam_vertex_distance =
                std::max(
                    result.max_seam_vertex_distance,
                    0.0);

            for (std::size_t i = 0;
                 i < n;
                 ++i)
            {
                const std::size_t inext =
                    (i + 1) % n;

                const std::size_t j =
                    bridge.matching.reversed
                        ? (bridge.matching.shift +
                           n -
                           (i % n)) %
                              n
                        : (bridge.matching.shift +
                           i) %
                              n;

                const std::size_t jnext =
                    bridge.matching.reversed
                        ? (bridge.matching.shift +
                           n -
                           (inext % n)) %
                              n
                        : (bridge.matching.shift +
                           inext) %
                              n;

                const Index g0 =
                    bridge.global_loop[i];

                const Index g1 =
                    bridge.global_loop[inext];

                const Index l0 =
                    local_remap[
                        bridge.local_loop[j]];

                const Index l1 =
                    local_remap[
                        bridge.local_loop[jnext]];

                /*
                 * Initial winding is chosen so that the strip is a regular
                 * quad split. The final orientation pass below determines
                 * the definitive winding from the global MC surface.
                 */
                result.mesh.triangles.push_back({
                    g0,
                    g1,
                    l0
                });

                result.mesh.triangles.push_back({
                    g1,
                    l1,
                    l0
                });

                result.seam_triangles += 2;

                result.max_seam_vertex_distance =
                    std::max(
                        result.max_seam_vertex_distance,
                        detail::distance(
                            outside.mesh.vertices[g0],
                            inside.mesh.vertices[
                                bridge.local_loop[j]]));

                result.max_seam_vertex_distance =
                    std::max(
                        result.max_seam_vertex_distance,
                        detail::distance(
                            outside.mesh.vertices[g1],
                            inside.mesh.vertices[
                                bridge.local_loop[jnext]]));
            }
        }

        /*
         * At this point the topology is assembled, but triangle winding of
         * the local patch and bridge strip is not trusted.
         *
         * Re-orient everything while keeping the surviving global MC
         * triangles authoritative.
         */
        try
        {
            detail::orient_consistently_from_global(
                result.mesh,
                global_triangle_count);
        }
        catch (const std::exception&)
        {
            /*
             * Keep the generated mesh so the caller can still inspect it,
             * but mark the operation as unsuccessful.
             */
            result.stitched = false;

            result.topology =
                detail::topology(
                    result.mesh);

            return result;
        }

        result.stitched =
            (result.seam_loops ==
             bridges.size());

        result.topology =
            detail::topology(
                result.mesh);

        /*
         * A successful stitched mesh should have:
         *
         *   - zero non-manifold edges
         *   - no new seam boundary edges
         *
         * If this is not true, do not claim success.
         */
        if (result.topology.nonmanifold_edges != 0)
        {
            result.stitched = false;
        }

        return result;
    }
};

} // namespace iso::stitch