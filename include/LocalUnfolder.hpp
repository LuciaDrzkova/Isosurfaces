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

        iso::mc::Bounds sampling_bounds{
            {0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0}
        };

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

        // Interface between refined and non-refined regions.
        std::size_t interface_edges = 0;
        std::size_t interface_vertices = 0;

        // Boundary of the complete surface caused by clipping
        // the implicit surface to the finite MC sampling box.
        std::size_t global_boundary_edges = 0;
        std::size_t global_boundary_vertices = 0;

        std::size_t new_vertices = 0;
        std::size_t new_triangles = 0;
        std::size_t projection_failures = 0;

        ScalarStatistics scalar;
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

        // Interface between refined and non-refined regions.
        std::size_t interface_edges = 0;
        std::size_t interface_vertices = 0;

        // Boundary of the complete surface caused by clipping
        // to the finite MC sampling box.
        std::size_t global_boundary_edges = 0;
        std::size_t global_boundary_vertices = 0;

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

        std::vector<bool> selected;

        std::size_t sphere_count = 0;
        std::size_t topological_count = 0;
        std::size_t selected_count = 0;
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

    static double distance(const Point& a, const Point& b)
    {
        return std::sqrt(squared_distance(a, b));
    }

    static Point midpoint(const Point& a, const Point& b)
    {
        return {
            0.5 * (a[0] + b[0]),
            0.5 * (a[1] + b[1]),
            0.5 * (a[2] + b[2])
        };
    }

    static Point centroid(
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

    static Eigen::Vector3d numerical_gradient(
        const Point& p,
        const ImplicitSurface& surface)
    {
        const double scale =
            std::max(1.0, distance(p, Point{0.0, 0.0, 0.0}));

        const double h = 1e-6 * scale;

        Point pxp = p;
        Point pxm = p;
        Point pyp = p;
        Point pym = p;
        Point pzp = p;
        Point pzm = p;

        pxp[0] += h;
        pxm[0] -= h;

        pyp[1] += h;
        pym[1] -= h;

        pzp[2] += h;
        pzm[2] -= h;

        const double dfdx =
            (surface.eval(pxp[0], pxp[1], pxp[2]) -
             surface.eval(pxm[0], pxm[1], pxm[2])) /
            (2.0 * h);

        const double dfdy =
            (surface.eval(pyp[0], pyp[1], pyp[2]) -
             surface.eval(pym[0], pym[1], pym[2])) /
            (2.0 * h);

        const double dfdz =
            (surface.eval(pzp[0], pzp[1], pzp[2]) -
             surface.eval(pzm[0], pzm[1], pzm[2])) /
            (2.0 * h);

        return Eigen::Vector3d(dfdx, dfdy, dfdz);
    }

    static Eigen::Vector3d gradient(
        const Point& p,
        const ImplicitSurface& surface)
    {
        const auto* quadric =
            dynamic_cast<const ParameterizedConeQuadric*>(&surface);

        if (quadric)
        {
            return quadric->gradient(
                p[0],
                p[1],
                p[2]);
        }

        return numerical_gradient(p, surface);
    }

    static bool project_to_surface(
        Point& p,
        const ImplicitSurface& surface,
        const Options& options)
    {
        for (std::size_t iteration = 0;
             iteration < options.projection_iterations;
             ++iteration)
        {
            const double value =
                surface.eval(p[0], p[1], p[2]) -
                options.isovalue;

            if (std::abs(value) <= options.projection_tolerance)
                return true;

            const Eigen::Vector3d grad =
                gradient(p, surface);

            const double grad_squared =
                grad.squaredNorm();

            if (grad_squared <=
                std::numeric_limits<double>::epsilon())
            {
                return false;
            }

            p[0] -= value * grad[0] / grad_squared;
            p[1] -= value * grad[1] / grad_squared;
            p[2] -= value * grad[2] / grad_squared;
        }

        return std::abs(
                   surface.eval(p[0], p[1], p[2]) -
                   options.isovalue) <=
               options.projection_tolerance;
    }

    static std::map<EdgeKey, EdgeAdjacency>
    build_edge_adjacency(const iso::mc::Mesh& mesh)
    {
        std::map<EdgeKey, EdgeAdjacency> adjacency;

        for (std::size_t fi = 0;
             fi < mesh.triangles.size();
             ++fi)
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
        std::vector<bool> selected(
            mesh.triangles.size(),
            false);

        const double radius_squared =
            radius * radius;

        for (std::size_t i = 0;
             i < mesh.triangles.size();
             ++i)
        {
            const Triangle& t =
                mesh.triangles[i];

            const Point c =
                centroid(
                    mesh.vertices[t[0]],
                    mesh.vertices[t[1]],
                    mesh.vertices[t[2]]);

            selected[i] =
                squared_distance(c, center) <=
                radius_squared;
        }

        return selected;
    }

    /*
     * Multi-source topological BFS.
     *
     * The singularity is represented by one or more coincident mesh
     * vertices. For the cone at c = 0 these are two different vertices
     * belonging to two different connected components.
     *
     * Starting BFS from only one vertex therefore refines only one cone
     * sheet. We instead find ALL vertices geometrically coincident with
     * the singularity center and use all of them as BFS seeds.
     */
    static std::vector<bool> select_faces_topological(
        const iso::mc::Mesh& mesh,
        const Point& center,
        double radius)
    {
        const std::size_t face_count =
            mesh.triangles.size();

        std::vector<bool> selected(
            face_count,
            false);

        if (mesh.vertices.empty() ||
            mesh.triangles.empty())
        {
            return selected;
        }

        const double radius_squared =
            radius * radius;

        /*
         * First find every mesh vertex geometrically coincident with
         * the singularity center.
         *
         * A small tolerance is necessary because MC interpolation can
         * produce values such as -2.7e-17 instead of exact zero.
         */
        constexpr double seed_tolerance = 1e-8;
        const double seed_tolerance_squared =
            seed_tolerance * seed_tolerance;

        std::vector<Index> seeds;

        for (Index vi = 0;
             vi < static_cast<Index>(mesh.vertices.size());
             ++vi)
        {
            if (squared_distance(
                    mesh.vertices[vi],
                    center) <= seed_tolerance_squared)
            {
                seeds.push_back(vi);
            }
        }

        /*
         * If no exact/coincident singular mesh vertex exists, fall back
         * to the nearest vertex. This keeps the method useful for
         * resolutions where the singularity is not sampled exactly.
         */
        if (seeds.empty())
        {
            Index nearest = 0;
            double nearest_distance =
                std::numeric_limits<double>::infinity();

            for (Index vi = 0;
                 vi < static_cast<Index>(mesh.vertices.size());
                 ++vi)
            {
                const double d2 =
                    squared_distance(
                        mesh.vertices[vi],
                        center);

                if (d2 < nearest_distance)
                {
                    nearest_distance = d2;
                    nearest = vi;
                }
            }

            seeds.push_back(nearest);
        }

        const auto adjacency =
            build_edge_adjacency(mesh);

        /*
         * Build vertex -> incident-face adjacency.
         */
        std::vector<std::vector<std::size_t>>
            vertex_faces(mesh.vertices.size());

        for (std::size_t fi = 0;
             fi < face_count;
             ++fi)
        {
            const Triangle& t =
                mesh.triangles[fi];

            vertex_faces[t[0]].push_back(fi);
            vertex_faces[t[1]].push_back(fi);
            vertex_faces[t[2]].push_back(fi);
        }

        /*
         * Multi-source BFS over the triangle adjacency graph.
         *
         * The radius is used as a geometric stopping condition on the
         * face centroid. Importantly, BFS is performed independently
         * from ALL singularity representatives.
         */
        std::queue<std::size_t> queue;
        std::vector<bool> visited(
            face_count,
            false);

        for (const Index seed : seeds)
        {
            for (const std::size_t face_index :
                 vertex_faces[seed])
            {
                if (!visited[face_index])
                {
                    visited[face_index] = true;
                    queue.push(face_index);
                }
            }
        }

        while (!queue.empty())
        {
            const std::size_t face_index =
                queue.front();

            queue.pop();

            const Triangle& face =
                mesh.triangles[face_index];

            const Point c =
                centroid(
                    mesh.vertices[face[0]],
                    mesh.vertices[face[1]],
                    mesh.vertices[face[2]]);

            /*
             * Do not expand outside the requested geometric region.
             *
             * A face outside the radius is not selected and is not
             * traversed further. This produces a connected topological
             * neighbourhood around every seed while preserving the
             * requested locality.
             */
            if (squared_distance(c, center) >
                radius_squared)
            {
                continue;
            }

            selected[face_index] = true;

            /*
             * Traverse all adjacent faces across triangle edges.
             */
            for (const EdgeKey edge : {
                     make_edge(face[0], face[1]),
                     make_edge(face[1], face[2]),
                     make_edge(face[2], face[0])})
            {
                const auto it =
                    adjacency.find(edge);

                if (it == adjacency.end())
                    continue;

                for (const std::size_t neighbour :
                     it->second.faces)
                {
                    if (!visited[neighbour])
                    {
                        visited[neighbour] = true;
                        queue.push(neighbour);
                    }
                }
            }
        }

        return selected;
    }

    static std::size_t count_vertices_in_region(
        const iso::mc::Mesh& mesh,
        const Point& center,
        double radius)
    {
        std::size_t count = 0;

        const double radius_squared =
            radius * radius;

        for (const Point& p : mesh.vertices)
        {
            if (squared_distance(p, center) <=
                radius_squared)
            {
                ++count;
            }
        }

        return count;
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
        ScalarStatistics statistics;
        statistics.refinement_level =
            refinement_level;

        if (base_resolution < 2)
            return statistics;

        if (refinement_level >=
            sizeof(std::size_t) * 8 - 1)
        {
            throw std::invalid_argument(
                "Local scalar refinement level is too large");
        }

        const std::size_t factor =
            std::size_t(1) << refinement_level;

        statistics.refinement_factor =
            factor;

        const std::size_t base_cells =
            base_resolution - 1;

        const std::size_t cells =
            base_cells * factor;

        const double hx =
            (bounds.max[0] - bounds.min[0]) /
            static_cast<double>(cells);

        const double hy =
            (bounds.max[1] - bounds.min[1]) /
            static_cast<double>(cells);

        const double hz =
            (bounds.max[2] - bounds.min[2]) /
            static_cast<double>(cells);

        statistics.spacing =
            std::min({hx, hy, hz});

        const double radius_squared =
            radius * radius;

        double minimum =
            std::numeric_limits<double>::infinity();

        double maximum =
            -std::numeric_limits<double>::infinity();

        long double sum = 0.0L;

        double best_d2 =
            std::numeric_limits<double>::infinity();

        /*
         * Same corner convention as MarchingCubes.hpp.
         */
        constexpr int corner_offset[8][3] = {
            {0, 0, 0},
            {1, 0, 0},
            {1, 0, 1},
            {0, 0, 1},
            {0, 1, 0},
            {1, 1, 0},
            {1, 1, 1},
            {0, 1, 1}
        };

        for (std::size_t z = 0;
             z < cells;
             ++z)
        {
            for (std::size_t y = 0;
                 y < cells;
                 ++y)
            {
                for (std::size_t x = 0;
                     x < cells;
                     ++x)
                {
                    const Point cell_center{
                        bounds.min[0] +
                            (static_cast<double>(x) + 0.5) * hx,
                        bounds.min[1] +
                            (static_cast<double>(y) + 0.5) * hy,
                        bounds.min[2] +
                            (static_cast<double>(z) + 0.5) * hz
                    };

                    if (squared_distance(
                            cell_center,
                            center) > radius_squared)
                    {
                        continue;
                    }

                    ++statistics.selected_voxels;

                    std::array<double, 8> values{};

                    for (int c = 0; c < 8; ++c)
                    {
                        const double px =
                            bounds.min[0] +
                            static_cast<double>(
                                x + corner_offset[c][0]) *
                                hx;

                        const double py =
                            bounds.min[1] +
                            static_cast<double>(
                                y + corner_offset[c][1]) *
                                hy;

                        const double pz =
                            bounds.min[2] +
                            static_cast<double>(
                                z + corner_offset[c][2]) *
                                hz;

                        values[c] =
                            surface.eval(
                                px,
                                py,
                                pz) -
                            isovalue;

                        minimum =
                            std::min(
                                minimum,
                                values[c]);

                        maximum =
                            std::max(
                                maximum,
                                values[c]);

                        sum +=
                            static_cast<long double>(
                                values[c]);

                        ++statistics.corner_value_samples;
                    }

                    const double d2 =
                        squared_distance(
                            cell_center,
                            center);

                    if (d2 < best_d2)
                    {
                        best_d2 = d2;
                        statistics.representative_values =
                            values;
                        statistics.has_representative_voxel =
                            true;
                    }
                }
            }
        }

        if (statistics.corner_value_samples > 0)
        {
            statistics.min_value = minimum;
            statistics.max_value = maximum;

            statistics.mean_value =
                static_cast<double>(
                    sum /
                    static_cast<long double>(
                        statistics.corner_value_samples));
        }

        return statistics;
    }

    static Index get_midpoint(
        Index a,
        Index b,
        const iso::mc::Mesh& current,
        iso::mc::Mesh& next,
        const ImplicitSurface& surface,
        const Options& options,
        std::map<EdgeKey, Index>& midpoint_cache,
        std::size_t& projection_failures)
    {
        const EdgeKey key =
            make_edge(a, b);

        const auto it =
            midpoint_cache.find(key);

        if (it != midpoint_cache.end())
            return it->second;

        Point p =
            midpoint(
                current.vertices[a],
                current.vertices[b]);

        if (!project_to_surface(
                p,
                surface,
                options))
        {
            ++projection_failures;
        }

        const Index index =
            static_cast<Index>(
                next.vertices.size());

        next.vertices.push_back(p);
        midpoint_cache[key] = index;

        return index;
    }

    static void append_triangle(
        std::vector<Triangle>& triangles,
        Index a,
        Index b,
        Index c)
    {
        triangles.push_back(
            Triangle{a, b, c});
    }

    static void split_face(
        const Triangle& face,
        unsigned char mask,
        Index m0,
        Index m1,
        Index m2,
        std::vector<Triangle>& output)
    {
        const Index v0 = face[0];
        const Index v1 = face[1];
        const Index v2 = face[2];

        switch (mask)
        {
        case 0:
            append_triangle(
                output,
                v0,
                v1,
                v2);
            break;

        case 1:
            append_triangle(
                output,
                v0,
                m0,
                v2);

            append_triangle(
                output,
                m0,
                v1,
                v2);
            break;

        case 2:
            append_triangle(
                output,
                v1,
                m1,
                v0);

            append_triangle(
                output,
                m1,
                v2,
                v0);
            break;

        case 3:
            append_triangle(
                output,
                v1,
                m1,
                m0);

            append_triangle(
                output,
                v0,
                m0,
                m1);

            append_triangle(
                output,
                v0,
                m1,
                v2);
            break;

        case 4:
            append_triangle(
                output,
                v2,
                m2,
                v1);

            append_triangle(
                output,
                m2,
                v0,
                v1);
            break;

        case 5:
            append_triangle(
                output,
                v0,
                m0,
                m2);

            append_triangle(
                output,
                m0,
                v1,
                v2);

            append_triangle(
                output,
                m0,
                v2,
                m2);
            break;

        case 6:
            append_triangle(
                output,
                v2,
                m2,
                m1);

            append_triangle(
                output,
                v0,
                v1,
                m1);

            append_triangle(
                output,
                v0,
                m1,
                m2);
            break;

        case 7:
            append_triangle(
                output,
                v0,
                m0,
                m2);

            append_triangle(
                output,
                m0,
                v1,
                m1);

            append_triangle(
                output,
                m2,
                m1,
                v2);

            append_triangle(
                output,
                m0,
                m1,
                m2);
            break;

        default:
            throw std::runtime_error(
                "Invalid local refinement split mask");
        }
    }

    static std::set<EdgeKey>
    compute_global_boundary_edges(
        const iso::mc::Mesh& mesh)
    {
        std::map<EdgeKey, std::size_t> uses;

        for (const Triangle& t :
             mesh.triangles)
        {
            ++uses[make_edge(t[0], t[1])];
            ++uses[make_edge(t[1], t[2])];
            ++uses[make_edge(t[2], t[0])];
        }

        std::set<EdgeKey> boundary;

        for (const auto& entry : uses)
        {
            if (entry.second == 1)
                boundary.insert(entry.first);
        }

        return boundary;
    }

    static std::set<Index>
    compute_global_boundary_vertices(
        const iso::mc::Mesh& mesh)
    {
        const auto boundary_edges =
            compute_global_boundary_edges(mesh);

        std::set<Index> vertices;

        for (const EdgeKey& edge :
             boundary_edges)
        {
            vertices.insert(edge.first);
            vertices.insert(edge.second);
        }

        return vertices;
    }

    static RegionSelection select_region(
        const iso::mc::Mesh& mesh,
        const Point& center,
        double radius)
    {
        RegionSelection result;

        result.sphere =
            select_faces_sphere(
                mesh,
                center,
                radius);

        result.topological =
            select_faces_topological(
                mesh,
                center,
                radius);

        result.sphere_count =
            static_cast<std::size_t>(
                std::count(
                    result.sphere.begin(),
                    result.sphere.end(),
                    true));

        result.topological_count =
            static_cast<std::size_t>(
                std::count(
                    result.topological.begin(),
                    result.topological.end(),
                    true));

        result.selected =
            result.topological;

        result.selected_count =
            result.topological_count;

        return result;
    }

public:
    Result refine(
        const iso::mc::Mesh& input,
        const ImplicitSurface& surface,
        const Point& center) const
    {
        return refine(
            input,
            surface,
            center,
            Options{});
    }

    Result refine(
        const iso::mc::Mesh& input,
        const ImplicitSurface& surface,
        const Point& center,
        const Options& options) const
    {
        if (options.radius <= 0.0)
        {
            throw std::invalid_argument(
                "LocalUnfolder radius must be positive");
        }

        Result result;

        result.input_vertices =
            input.vertices.size();

        result.input_triangles =
            input.triangles.size();

        result.mesh = input;

        result.initial_scalar =
            collect_scalar_statistics(
                surface,
                center,
                options.radius,
                options.sampling_bounds,
                options.base_resolution,
                options.isovalue,
                0);

        /*
         * The initial region is selected using the SAME multi-source
         * topological BFS used during refinement.
         */
        {
            const RegionSelection initial =
                select_region(
                    input,
                    center,
                    options.radius);

            result.sphere_faces_in_region =
                initial.sphere_count;

            result.topological_faces_in_region =
                initial.topological_count;

            if (options.region_mode ==
                RegionMode::Sphere)
            {
                result.faces_in_region =
                    initial.sphere_count;
            }
            else
            {
                result.faces_in_region =
                    initial.topological_count;
            }

            std::set<Index> region_vertices;

            for (std::size_t fi = 0;
                 fi < input.triangles.size();
                 ++fi)
            {
                if (!initial.selected[fi])
                    continue;

                const Triangle& t =
                    input.triangles[fi];

                region_vertices.insert(t[0]);
                region_vertices.insert(t[1]);
                region_vertices.insert(t[2]);
            }

            result.vertices_in_region =
                region_vertices.size();

            /*
             * Compute initial refinement interface.
             */
            const auto adjacency =
                build_edge_adjacency(input);

            std::set<EdgeKey> interface_edges;
            std::set<Index> interface_vertices;

            for (const auto& entry :
                 adjacency)
            {
                bool has_selected = false;
                bool has_unselected = false;

                for (const std::size_t fi :
                     entry.second.faces)
                {
                    if (initial.selected[fi])
                        has_selected = true;
                    else
                        has_unselected = true;
                }

                if (has_selected &&
                    has_unselected)
                {
                    interface_edges.insert(
                        entry.first);

                    interface_vertices.insert(
                        entry.first.first);

                    interface_vertices.insert(
                        entry.first.second);
                }
            }

            result.interface_edges =
                interface_edges.size();

            result.interface_vertices =
                interface_vertices.size();
        }

        if (options.levels == 0 ||
            input.triangles.empty())
        {
            result.output_vertices =
                result.mesh.vertices.size();

            result.output_triangles =
                result.mesh.triangles.size();

            const auto global_boundary_edges =
                compute_global_boundary_edges(
                    result.mesh);

            const auto global_boundary_vertices =
                compute_global_boundary_vertices(
                    result.mesh);

            result.global_boundary_edges =
                global_boundary_edges.size();

            result.global_boundary_vertices =
                global_boundary_vertices.size();

            return result;
        }

        for (std::size_t level = 1;
             level <= options.levels;
             ++level)
        {
            const iso::mc::Mesh current =
                result.mesh;

            const RegionSelection region =
                select_region(
                    current,
                    center,
                    options.radius);

            const std::vector<bool>& sphere_selected =
                region.sphere;

            const std::vector<bool>& bfs_selected =
                region.topological;

            std::vector<bool> selected;

            if (options.region_mode ==
                RegionMode::Sphere)
            {
                selected = sphere_selected;
            }
            else
            {
                /*
                 * IMPORTANT:
                 *
                 * bfs_selected is a multi-source BFS. At c=0 the two
                 * cone sheets have separate mesh components, but both
                 * contain a vertex at the same geometric singularity.
                 *
                 * Consequently both components are seeded and both
                 * sheets are selected here.
                 */
                selected = bfs_selected;
            }

            const std::size_t selected_count =
                static_cast<std::size_t>(
                    std::count(
                        selected.begin(),
                        selected.end(),
                        true));

            if (selected_count == 0)
                break;

            const auto adjacency =
                build_edge_adjacency(current);

            std::map<EdgeKey, bool> refine_edges;

            std::set<EdgeKey>
                interface_edges;

            std::set<Index>
                interface_vertices;

            for (const auto& entry :
                 adjacency)
            {
                const EdgeKey& edge =
                    entry.first;

                const auto& faces =
                    entry.second.faces;

                bool has_selected = false;
                bool has_unselected = false;

                for (const std::size_t face_index :
                     faces)
                {
                    if (selected[face_index])
                        has_selected = true;
                    else
                        has_unselected = true;
                }

                /*
                 * Every edge of a selected triangle is split.
                 */
                if (has_selected)
                {
                    refine_edges[edge] = true;
                }

                /*
                 * Interface means one selected face and one
                 * unselected face sharing the edge.
                 */
                if (has_selected &&
                    has_unselected)
                {
                    interface_edges.insert(edge);

                    interface_vertices.insert(
                        edge.first);

                    interface_vertices.insert(
                        edge.second);
                }
            }

            iso::mc::Mesh next;

            next.vertices =
                current.vertices;

            next.triangles.reserve(
                current.triangles.size() * 2);

            std::map<EdgeKey, Index>
                midpoint_cache;

            std::size_t projection_failures = 0;

            for (std::size_t fi = 0;
                 fi < current.triangles.size();
                 ++fi)
            {
                const Triangle& face =
                    current.triangles[fi];

                const EdgeKey e0 =
                    make_edge(
                        face[0],
                        face[1]);

                const EdgeKey e1 =
                    make_edge(
                        face[1],
                        face[2]);

                const EdgeKey e2 =
                    make_edge(
                        face[2],
                        face[0]);

                const bool split0 =
                    refine_edges[e0];

                const bool split1 =
                    refine_edges[e1];

                const bool split2 =
                    refine_edges[e2];

                unsigned char mask = 0;

                if (split0)
                    mask |= 1;

                if (split1)
                    mask |= 2;

                if (split2)
                    mask |= 4;

                Index m0 = 0;
                Index m1 = 0;
                Index m2 = 0;

                if (split0)
                {
                    m0 =
                        get_midpoint(
                            face[0],
                            face[1],
                            current,
                            next,
                            surface,
                            options,
                            midpoint_cache,
                            projection_failures);
                }

                if (split1)
                {
                    m1 =
                        get_midpoint(
                            face[1],
                            face[2],
                            current,
                            next,
                            surface,
                            options,
                            midpoint_cache,
                            projection_failures);
                }

                if (split2)
                {
                    m2 =
                        get_midpoint(
                            face[2],
                            face[0],
                            current,
                            next,
                            surface,
                            options,
                            midpoint_cache,
                            projection_failures);
                }

                split_face(
                    face,
                    mask,
                    m0,
                    m1,
                    m2,
                    next.triangles);
            }

            const auto global_boundary_edges =
                compute_global_boundary_edges(
                    next);

            const auto global_boundary_vertices =
                compute_global_boundary_vertices(
                    next);

            LevelStatistics statistics;

            statistics.level = level;

            statistics.selected_faces =
                selected_count;

            statistics.selected_faces_sphere =
                region.sphere_count;

            statistics.selected_faces_topological =
                region.topological_count;

            statistics.selection_changed =
                region.sphere !=
                region.topological;

            statistics.output_vertices =
                next.vertices.size();

            statistics.output_triangles =
                next.triangles.size();

            statistics.split_edges =
                midpoint_cache.size();

            statistics.interface_edges =
                interface_edges.size();

            statistics.interface_vertices =
                interface_vertices.size();

            statistics.global_boundary_edges =
                global_boundary_edges.size();

            statistics.global_boundary_vertices =
                global_boundary_vertices.size();

            statistics.new_vertices =
                next.vertices.size() -
                current.vertices.size();

            statistics.new_triangles =
                next.triangles.size() -
                current.triangles.size();

            statistics.projection_failures =
                projection_failures;

            statistics.scalar =
                collect_scalar_statistics(
                    surface,
                    center,
                    options.radius,
                    options.sampling_bounds,
                    options.base_resolution,
                    options.isovalue,
                    level);

            result.interface_edges =
                statistics.interface_edges;

            result.interface_vertices =
                statistics.interface_vertices;

            result.global_boundary_edges =
                statistics.global_boundary_edges;

            result.global_boundary_vertices =
                statistics.global_boundary_vertices;

            result.levels.push_back(
                statistics);

            result.mesh =
                std::move(next);
        }

        result.output_vertices =
            result.mesh.vertices.size();

        result.output_triangles =
            result.mesh.triangles.size();

        const auto final_global_boundary_edges =
            compute_global_boundary_edges(
                result.mesh);

        const auto final_global_boundary_vertices =
            compute_global_boundary_vertices(
                result.mesh);

        result.global_boundary_edges =
            final_global_boundary_edges.size();

        result.global_boundary_vertices =
            final_global_boundary_vertices.size();

        return result;
    }
};