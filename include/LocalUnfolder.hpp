#pragma once

#include "ImplicitSurface.hpp"
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

class LocalUnfolder
{
public:
    enum class RegionMode
    {
        Sphere,
        TopologicalBfs
    };

    struct ScalarStatistics
    {
        std::size_t refinement_level = 0;
        std::size_t refinement_factor = 1;
        std::size_t selected_voxels = 0;
        std::size_t corner_value_samples = 0;
        double spacing = 0.0;
        double min_value = 0.0;
        double max_value = 0.0;
        double mean_value = 0.0;
        bool has_representative_voxel = false;
        std::array<double, 8> representative_values{};
    };

    struct Options
    {
        double radius = 0.3;
        std::size_t levels = 2;
        double projection_tolerance = 1e-10;
        std::size_t projection_iterations = 12;
        RegionMode region_mode = RegionMode::TopologicalBfs;
        iso::mc::Bounds sampling_bounds{{0.0, 0.0, 0.0},
                                         {0.0, 0.0, 0.0}};
        std::size_t base_resolution = 0;
        double isovalue = 0.0;
    };

    struct LevelStatistics
    {
        std::size_t level = 0;
        std::size_t selected_faces = 0;
        std::size_t selected_faces_sphere = 0;
        std::size_t selected_faces_topological = 0;
        bool selection_changed = false;
        std::size_t output_vertices = 0;
        std::size_t output_triangles = 0;
        std::size_t split_edges = 0;

        // Interface between refined and untouched surface regions.
        std::size_t interface_edges = 0;
        std::size_t interface_vertices = 0;

        // Boundary produced by clipping the complete surface to the finite
        // global Marching Cubes sampling box.
        std::size_t global_boundary_edges = 0;
        std::size_t global_boundary_vertices = 0;

        std::size_t new_vertices = 0;
        std::size_t new_triangles = 0;
        std::size_t projection_failures = 0;
        ScalarStatistics scalar;

        // Compatibility aliases for older callers. These mean refinement
        // interface counts, not the global surface boundary.
        std::size_t boundary_edges = 0;
        std::size_t boundary_vertices = 0;
    };

    struct Result
    {
        iso::mc::Mesh mesh;
        std::size_t input_vertices = 0;
        std::size_t input_triangles = 0;
        std::size_t output_vertices = 0;
        std::size_t output_triangles = 0;
        std::size_t vertices_in_region = 0;
        std::size_t faces_in_region = 0;
        std::size_t sphere_faces_in_region = 0;
        std::size_t topological_faces_in_region = 0;

        std::size_t interface_edges = 0;
        std::size_t interface_vertices = 0;
        std::size_t global_boundary_edges = 0;
        std::size_t global_boundary_vertices = 0;

        // Compatibility aliases. These mean refinement interface counts.
        std::size_t boundary_edges = 0;
        std::size_t boundary_vertices = 0;

        ScalarStatistics initial_scalar;
        std::vector<LevelStatistics> levels;
    };

private:
    using Index = iso::mc::Index;
    using Triangle = iso::mc::Triangle;
    using Point = iso::mc::Point;
    using EdgeKey = std::pair<Index, Index>;

    struct EdgeAdjacency
    {
        std::vector<std::size_t> faces;
    };

    struct RegionSelection
    {
        std::vector<bool> sphere;
        std::vector<bool> topological;
        std::size_t sphere_count = 0;
        std::size_t topological_count = 0;
    };

    static EdgeKey make_edge(Index a, Index b)
    {
        if (a > b)
            std::swap(a, b);
        return {a, b};
    }

    static double squared_distance(const Point& a, const Point& b)
    {
        const double dx = a[0] - b[0];
        const double dy = a[1] - b[1];
        const double dz = a[2] - b[2];
        return dx * dx + dy * dy + dz * dz;
    }

    static double distance(const Point& a)
    {
        return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    }

    static Point centroid(const Point& a, const Point& b, const Point& c)
    {
        return {(a[0] + b[0] + c[0]) / 3.0,
                (a[1] + b[1] + c[1]) / 3.0,
                (a[2] + b[2] + c[2]) / 3.0};
    }

    static Point midpoint(const Point& a, const Point& b)
    {
        return {0.5 * (a[0] + b[0]),
                0.5 * (a[1] + b[1]),
                0.5 * (a[2] + b[2])};
    }

    static Eigen::Vector3d numerical_gradient(
        const Point& p,
        const ImplicitSurface& surface)
    {
        const double h = 1e-6 * std::max(1.0, distance(p));
        Point pxp = p, pxm = p, pyp = p, pym = p, pzp = p, pzm = p;
        pxp[0] += h; pxm[0] -= h;
        pyp[1] += h; pym[1] -= h;
        pzp[2] += h; pzm[2] -= h;

        return Eigen::Vector3d(
            (surface.eval(pxp[0], pxp[1], pxp[2]) -
             surface.eval(pxm[0], pxm[1], pxm[2])) / (2.0 * h),
            (surface.eval(pyp[0], pyp[1], pyp[2]) -
             surface.eval(pym[0], pym[1], pym[2])) / (2.0 * h),
            (surface.eval(pzp[0], pzp[1], pzp[2]) -
             surface.eval(pzm[0], pzm[1], pzm[2])) / (2.0 * h));
    }

    static Eigen::Vector3d gradient(
        const Point& p,
        const ImplicitSurface& surface)
    {
        if (const auto* quadric =
                dynamic_cast<const ParameterizedConeQuadric*>(&surface))
            return quadric->gradient(p[0], p[1], p[2]);
        return numerical_gradient(p, surface);
    }

    static bool project_to_surface(
        Point& p,
        const ImplicitSurface& surface,
        const Options& options)
    {
        for (std::size_t i = 0; i < options.projection_iterations; ++i)
        {
            const double value =
                surface.eval(p[0], p[1], p[2]) - options.isovalue;
            if (std::abs(value) <= options.projection_tolerance)
                return true;

            const Eigen::Vector3d g = gradient(p, surface);
            const double g2 = g.squaredNorm();
            if (g2 <= std::numeric_limits<double>::epsilon())
                return false;

            p[0] -= value * g[0] / g2;
            p[1] -= value * g[1] / g2;
            p[2] -= value * g[2] / g2;
        }

        return std::abs(
                   surface.eval(p[0], p[1], p[2]) - options.isovalue) <=
               options.projection_tolerance;
    }

    static std::map<EdgeKey, EdgeAdjacency>
    build_edge_adjacency(const iso::mc::Mesh& mesh)
    {
        std::map<EdgeKey, EdgeAdjacency> adjacency;
        for (std::size_t fi = 0; fi < mesh.triangles.size(); ++fi)
        {
            const Triangle& t = mesh.triangles[fi];
            adjacency[make_edge(t[0], t[1])].faces.push_back(fi);
            adjacency[make_edge(t[1], t[2])].faces.push_back(fi);
            adjacency[make_edge(t[2], t[0])].faces.push_back(fi);
        }
        return adjacency;
    }

    static std::vector<bool> select_faces_sphere(
        const iso::mc::Mesh& mesh,
        const Point& center,
        double radius)
    {
        std::vector<bool> selected(mesh.triangles.size(), false);
        const double radius2 = radius * radius;
        for (std::size_t i = 0; i < mesh.triangles.size(); ++i)
        {
            const Triangle& t = mesh.triangles[i];
            selected[i] = squared_distance(
                centroid(mesh.vertices[t[0]], mesh.vertices[t[1]], mesh.vertices[t[2]]),
                center) <= radius2;
        }
        return selected;
    }

    static std::vector<bool> select_faces_topological(
        const iso::mc::Mesh& mesh,
        const Point& center,
        double radius)
    {
        std::vector<bool> selected(mesh.triangles.size(), false);
        if (mesh.vertices.empty() || mesh.triangles.empty())
            return selected;

        constexpr double seed_tolerance = 1e-8;
        const double seed2 = seed_tolerance * seed_tolerance;
        std::vector<Index> seeds;

        for (Index vi = 0; vi < static_cast<Index>(mesh.vertices.size()); ++vi)
        {
            if (squared_distance(mesh.vertices[vi], center) <= seed2)
                seeds.push_back(vi);
        }

        if (seeds.empty())
        {
            Index nearest = 0;
            double best = std::numeric_limits<double>::infinity();
            for (Index vi = 0; vi < static_cast<Index>(mesh.vertices.size()); ++vi)
            {
                const double d2 = squared_distance(mesh.vertices[vi], center);
                if (d2 < best)
                {
                    best = d2;
                    nearest = vi;
                }
            }
            seeds.push_back(nearest);
        }

        const auto adjacency = build_edge_adjacency(mesh);
        std::vector<std::vector<std::size_t>> vertex_faces(mesh.vertices.size());
        for (std::size_t fi = 0; fi < mesh.triangles.size(); ++fi)
        {
            const Triangle& t = mesh.triangles[fi];
            vertex_faces[t[0]].push_back(fi);
            vertex_faces[t[1]].push_back(fi);
            vertex_faces[t[2]].push_back(fi);
        }

        std::queue<std::size_t> queue;
        std::vector<bool> visited(mesh.triangles.size(), false);
        for (Index seed : seeds)
        {
            for (std::size_t fi : vertex_faces[seed])
            {
                if (!visited[fi])
                {
                    visited[fi] = true;
                    queue.push(fi);
                }
            }
        }

        const double radius2 = radius * radius;
        while (!queue.empty())
        {
            const std::size_t fi = queue.front();
            queue.pop();

            const Triangle& t = mesh.triangles[fi];
            const Point c = centroid(mesh.vertices[t[0]],
                                     mesh.vertices[t[1]],
                                     mesh.vertices[t[2]]);
            if (squared_distance(c, center) > radius2)
                continue;

            selected[fi] = true;
            for (const EdgeKey e : {
                     make_edge(t[0], t[1]),
                     make_edge(t[1], t[2]),
                     make_edge(t[2], t[0])})
            {
                const auto it = adjacency.find(e);
                if (it == adjacency.end())
                    continue;
                for (std::size_t n : it->second.faces)
                {
                    if (!visited[n])
                    {
                        visited[n] = true;
                        queue.push(n);
                    }
                }
            }
        }

        return selected;
    }

    static RegionSelection select_region(
        const iso::mc::Mesh& mesh,
        const Point& center,
        double radius)
    {
        RegionSelection region;
        region.sphere = select_faces_sphere(mesh, center, radius);
        region.topological = select_faces_topological(mesh, center, radius);
        region.sphere_count = static_cast<std::size_t>(
            std::count(region.sphere.begin(), region.sphere.end(), true));
        region.topological_count = static_cast<std::size_t>(
            std::count(region.topological.begin(), region.topological.end(), true));
        return region;
    }

    static std::size_t count_region_vertices(
        const iso::mc::Mesh& mesh,
        const std::vector<bool>& selected)
    {
        std::set<Index> vertices;
        for (std::size_t fi = 0; fi < mesh.triangles.size(); ++fi)
        {
            if (!selected[fi])
                continue;
            const Triangle& t = mesh.triangles[fi];
            vertices.insert(t[0]);
            vertices.insert(t[1]);
            vertices.insert(t[2]);
        }
        return vertices.size();
    }

    static ScalarStatistics collect_scalar_statistics(
        const ImplicitSurface& surface,
        const Point& center,
        double radius,
        const iso::mc::Bounds& bounds,
        std::size_t base_resolution,
        double isovalue,
        std::size_t refinement_level)
    {
        ScalarStatistics s;
        s.refinement_level = refinement_level;
        if (base_resolution < 2)
            return s;
        const std::size_t cells_per_axis =
            (base_resolution - 1) * (std::size_t(1) << refinement_level);
        const double hx = (bounds.max[0] - bounds.min[0]) / cells_per_axis;
        const double hy = (bounds.max[1] - bounds.min[1]) / cells_per_axis;
        const double hz = (bounds.max[2] - bounds.min[2]) / cells_per_axis;
        s.spacing = std::min({hx, hy, hz});
        s.refinement_factor = std::size_t(1) << refinement_level;

        const double radius2 = radius * radius;
        double min_v = std::numeric_limits<double>::infinity();
        double max_v = -std::numeric_limits<double>::infinity();
        long double sum = 0.0L;
        double best_d2 = std::numeric_limits<double>::infinity();

        constexpr int off[8][3] = {
            {0,0,0}, {1,0,0}, {1,0,1}, {0,0,1},
            {0,1,0}, {1,1,0}, {1,1,1}, {0,1,1}};

        for (std::size_t z = 0; z < cells_per_axis; ++z)
        {
            for (std::size_t y = 0; y < cells_per_axis; ++y)
            {
                for (std::size_t x = 0; x < cells_per_axis; ++x)
                {
                    const Point c{
                        bounds.min[0] + (x + 0.5) * hx,
                        bounds.min[1] + (y + 0.5) * hy,
                        bounds.min[2] + (z + 0.5) * hz};
                    if (squared_distance(c, center) > radius2)
                        continue;
                    ++s.selected_voxels;

                    std::array<double, 8> values{};
                    for (int i = 0; i < 8; ++i)
                    {
                        const double px = bounds.min[0] + (x + off[i][0]) * hx;
                        const double py = bounds.min[1] + (y + off[i][1]) * hy;
                        const double pz = bounds.min[2] + (z + off[i][2]) * hz;
                        values[i] = surface.eval(px, py, pz) - isovalue;
                        min_v = std::min(min_v, values[i]);
                        max_v = std::max(max_v, values[i]);
                        sum += static_cast<long double>(values[i]);
                        ++s.corner_value_samples;
                    }

                    const double d2 = squared_distance(c, center);
                    if (d2 < best_d2)
                    {
                        best_d2 = d2;
                        s.representative_values = values;
                        s.has_representative_voxel = true;
                    }
                }
            }
        }

        if (s.corner_value_samples > 0)
        {
            s.min_value = min_v;
            s.max_value = max_v;
            s.mean_value = static_cast<double>(
                sum / static_cast<long double>(s.corner_value_samples));
        }
        return s;
    }

    static Index midpoint_vertex(
        Index a,
        Index b,
        const iso::mc::Mesh& current,
        iso::mc::Mesh& next,
        const ImplicitSurface& surface,
        const Options& options,
        std::map<EdgeKey, Index>& cache,
        std::size_t& projection_failures)
    {
        const EdgeKey e = make_edge(a, b);
        const auto found = cache.find(e);
        if (found != cache.end())
            return found->second;

        Point p = midpoint(current.vertices[a], current.vertices[b]);
        if (!project_to_surface(p, surface, options))
            ++projection_failures;

        const Index index = static_cast<Index>(next.vertices.size());
        next.vertices.push_back(p);
        cache[e] = index;
        return index;
    }

    static void add_triangle(
        std::vector<Triangle>& triangles,
        Index a,
        Index b,
        Index c)
    {
        triangles.push_back(Triangle{a, b, c});
    }

    static void subdivide_face(
        const Triangle& t,
        unsigned char mask,
        Index m0,
        Index m1,
        Index m2,
        std::vector<Triangle>& output)
    {
        const Index a = t[0], b = t[1], c = t[2];
        switch (mask)
        {
        case 0:
            add_triangle(output, a, b, c); break;
        case 1:
            add_triangle(output, a, m0, c);
            add_triangle(output, m0, b, c); break;
        case 2:
            add_triangle(output, b, m1, a);
            add_triangle(output, m1, c, a); break;
        case 3:
            add_triangle(output, b, m1, m0);
            add_triangle(output, a, m0, m1);
            add_triangle(output, a, m1, c); break;
        case 4:
            add_triangle(output, c, m2, b);
            add_triangle(output, m2, a, b); break;
        case 5:
            add_triangle(output, a, m0, m2);
            add_triangle(output, m0, b, c);
            add_triangle(output, m0, c, m2); break;
        case 6:
            add_triangle(output, c, m2, m1);
            add_triangle(output, a, b, m1);
            add_triangle(output, a, m1, m2); break;
        case 7:
            add_triangle(output, a, m0, m2);
            add_triangle(output, m0, b, m1);
            add_triangle(output, m2, m1, c);
            add_triangle(output, m0, m1, m2); break;
        default:
            throw std::runtime_error("Invalid local refinement mask");
        }
    }

    static void global_boundary_statistics(
        const iso::mc::Mesh& mesh,
        std::size_t& edges,
        std::size_t& vertices)
    {
        std::map<EdgeKey, std::size_t> uses;
        for (const Triangle& t : mesh.triangles)
        {
            ++uses[make_edge(t[0], t[1])];
            ++uses[make_edge(t[1], t[2])];
            ++uses[make_edge(t[2], t[0])];
        }

        std::set<Index> boundary_vertices;
        edges = 0;
        for (const auto& [e, count] : uses)
        {
            if (count == 1)
            {
                ++edges;
                boundary_vertices.insert(e.first);
                boundary_vertices.insert(e.second);
            }
        }
        vertices = boundary_vertices.size();
    }

public:
    Result refine(
        const iso::mc::Mesh& input,
        const ImplicitSurface& surface,
        const Point& center) const
    {
        return refine(input, surface, center, Options{});
    }

    Result refine(
        const iso::mc::Mesh& input,
        const ImplicitSurface& surface,
        const Point& center,
        const Options& options) const
    {
        if (!(options.radius > 0.0))
            throw std::invalid_argument("LocalUnfolder radius must be positive");

        Result result;
        result.input_vertices = input.vertices.size();
        result.input_triangles = input.triangles.size();
        result.mesh = input;

        const RegionSelection initial =
            select_region(input, center, options.radius);
        result.sphere_faces_in_region = initial.sphere_count;
        result.topological_faces_in_region = initial.topological_count;

        const std::vector<bool>& chosen =
            options.region_mode == RegionMode::Sphere
                ? initial.sphere
                : initial.topological;

        result.faces_in_region = static_cast<std::size_t>(
            std::count(chosen.begin(), chosen.end(), true));
        result.vertices_in_region = count_region_vertices(input, chosen);
        result.initial_scalar = collect_scalar_statistics(
            surface, center, options.radius, options.sampling_bounds,
            options.base_resolution, options.isovalue, 0);

        for (std::size_t level = 1; level <= options.levels; ++level)
        {
            const iso::mc::Mesh current = result.mesh;
            const RegionSelection region =
                select_region(current, center, options.radius);
            const std::vector<bool>& selected =
                options.region_mode == RegionMode::Sphere
                    ? region.sphere
                    : region.topological;

            const std::size_t selected_count = static_cast<std::size_t>(
                std::count(selected.begin(), selected.end(), true));
            if (selected_count == 0)
                break;

            const auto adjacency = build_edge_adjacency(current);
            std::map<EdgeKey, bool> refine_edges;
            std::set<EdgeKey> interface_edges;
            std::set<Index> interface_vertices;

            for (const auto& [edge, a] : adjacency)
            {
                bool has_selected = false;
                bool has_unselected = false;
                for (std::size_t fi : a.faces)
                {
                    has_selected |= selected[fi];
                    has_unselected |= !selected[fi];
                }
                if (has_selected)
                    refine_edges[edge] = true;
                if (has_selected && has_unselected)
                {
                    interface_edges.insert(edge);
                    interface_vertices.insert(edge.first);
                    interface_vertices.insert(edge.second);
                }
            }

            iso::mc::Mesh next;
            next.vertices = current.vertices;
            next.triangles.reserve(current.triangles.size() * 2);

            std::map<EdgeKey, Index> cache;
            std::size_t projection_failures = 0;

            for (std::size_t fi = 0; fi < current.triangles.size(); ++fi)
            {
                const Triangle& t = current.triangles[fi];
                const EdgeKey e0 = make_edge(t[0], t[1]);
                const EdgeKey e1 = make_edge(t[1], t[2]);
                const EdgeKey e2 = make_edge(t[2], t[0]);
                const bool s0 = refine_edges[e0];
                const bool s1 = refine_edges[e1];
                const bool s2 = refine_edges[e2];

                unsigned char mask = 0;
                if (s0) mask |= 1;
                if (s1) mask |= 2;
                if (s2) mask |= 4;

                Index m0 = 0, m1 = 0, m2 = 0;
                if (s0)
                    m0 = midpoint_vertex(t[0], t[1], current, next,
                                         surface, options, cache, projection_failures);
                if (s1)
                    m1 = midpoint_vertex(t[1], t[2], current, next,
                                         surface, options, cache, projection_failures);
                if (s2)
                    m2 = midpoint_vertex(t[2], t[0], current, next,
                                         surface, options, cache, projection_failures);

                subdivide_face(t, mask, m0, m1, m2, next.triangles);
            }

            LevelStatistics stats;
            stats.level = level;
            stats.selected_faces = selected_count;
            stats.selected_faces_sphere = region.sphere_count;
            stats.selected_faces_topological = region.topological_count;
            stats.selection_changed = region.sphere != region.topological;
            stats.output_vertices = next.vertices.size();
            stats.output_triangles = next.triangles.size();
            stats.split_edges = cache.size();
            stats.interface_edges = interface_edges.size();
            stats.interface_vertices = interface_vertices.size();
            stats.boundary_edges = stats.interface_edges;
            stats.boundary_vertices = stats.interface_vertices;
            stats.new_vertices = next.vertices.size() - current.vertices.size();
            stats.new_triangles = next.triangles.size() - current.triangles.size();
            stats.projection_failures = projection_failures;
            global_boundary_statistics(
                next, stats.global_boundary_edges, stats.global_boundary_vertices);
            stats.scalar = collect_scalar_statistics(
                surface, center, options.radius, options.sampling_bounds,
                options.base_resolution, options.isovalue, level);

            result.interface_edges = stats.interface_edges;
            result.interface_vertices = stats.interface_vertices;
            result.boundary_edges = stats.interface_edges;
            result.boundary_vertices = stats.interface_vertices;
            result.global_boundary_edges = stats.global_boundary_edges;
            result.global_boundary_vertices = stats.global_boundary_vertices;
            result.levels.push_back(stats);
            result.mesh = std::move(next);
        }

        result.output_vertices = result.mesh.vertices.size();
        result.output_triangles = result.mesh.triangles.size();
        global_boundary_statistics(
            result.mesh,
            result.global_boundary_edges,
            result.global_boundary_vertices);
        return result;
    }
};
