#pragma once

#include "ImplicitSurface.hpp"
#include "MarchingCubes.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

class LocalUnfolder
{
public:
    struct Options
    {
        double radius = 0.3;
        std::size_t levels = 2;

        double projection_tolerance = 1e-10;
        std::size_t projection_iterations = 12;
    };

    struct LevelStatistics
    {
        std::size_t level = 0;
        std::size_t selected_faces = 0;
        std::size_t output_vertices = 0;
        std::size_t output_triangles = 0;
        std::size_t new_vertices = 0;
        std::size_t new_triangles = 0;
        std::size_t projection_failures = 0;
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

    struct EdgeSplitMap
    {
        std::map<EdgeKey, Index> midpoint;
    };

    static EdgeKey make_edge(
        Index a,
        Index b)
    {
        if (a > b)
            std::swap(a, b);

        return {a, b};
    }

    static double squared_distance(
        const Point& a,
        const Point& b)
    {
        const double dx = a[0] - b[0];
        const double dy = a[1] - b[1];
        const double dz = a[2] - b[2];

        return dx * dx + dy * dy + dz * dz;
    }

    static double distance(
        const Point& a,
        const Point& b)
    {
        return std::sqrt(
            squared_distance(a, b));
    }

    static Point midpoint(
        const Point& a,
        const Point& b)
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
            std::max(
                1.0,
                distance(
                    p,
                    Point{0.0, 0.0, 0.0}));

        const double h =
            1e-6 * scale;

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
            (surface.eval(
                 pxp[0],
                 pxp[1],
                 pxp[2]) -
             surface.eval(
                 pxm[0],
                 pxm[1],
                 pxm[2])) /
            (2.0 * h);

        const double dfdy =
            (surface.eval(
                 pyp[0],
                 pyp[1],
                 pyp[2]) -
             surface.eval(
                 pym[0],
                 pym[1],
                 pym[2])) /
            (2.0 * h);

        const double dfdz =
            (surface.eval(
                 pzp[0],
                 pzp[1],
                 pzp[2]) -
             surface.eval(
                 pzm[0],
                 pzm[1],
                 pzm[2])) /
            (2.0 * h);

        return Eigen::Vector3d(
            dfdx,
            dfdy,
            dfdz);
    }

    static Eigen::Vector3d gradient(
        const Point& p,
        const ImplicitSurface& surface)
    {
        const auto* quadric =
            dynamic_cast<
                const ParameterizedConeQuadric*>(
                &surface);

        if (quadric)
        {
            return quadric->gradient(
                p[0],
                p[1],
                p[2]);
        }

        return numerical_gradient(
            p,
            surface);
    }

    static bool project_to_surface(
        Point& p,
        const ImplicitSurface& surface,
        const Options& options)
    {
        for (std::size_t iteration = 0;
             iteration <
             options.projection_iterations;
             ++iteration)
        {
            const double value =
                surface.eval(
                    p[0],
                    p[1],
                    p[2]);

            if (std::abs(value) <=
                options.projection_tolerance)
            {
                return true;
            }

            const Eigen::Vector3d grad =
                gradient(
                    p,
                    surface);

            const double grad_squared =
                grad.squaredNorm();

            if (grad_squared <=
                std::numeric_limits<double>::epsilon())
            {
                return false;
            }

            p[0] -=
                value * grad[0] /
                grad_squared;

            p[1] -=
                value * grad[1] /
                grad_squared;

            p[2] -=
                value * grad[2] /
                grad_squared;
        }

        return std::abs(
                   surface.eval(
                       p[0],
                       p[1],
                       p[2])) <=
               options.projection_tolerance;
    }

    static Index get_midpoint(
        Index a,
        Index b,
        const iso::mc::Mesh& current,
        iso::mc::Mesh& next,
        const ImplicitSurface& surface,
        const Options& options,
        EdgeSplitMap& midpoint_cache,
        std::size_t& projection_failures)
    {
        const EdgeKey key =
            make_edge(a, b);

        const auto it =
            midpoint_cache.midpoint.find(key);

        if (it !=
            midpoint_cache.midpoint.end())
        {
            return it->second;
        }

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

        midpoint_cache.midpoint[key] =
            index;

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
        unsigned char split_mask,
        Index m0,
        Index m1,
        Index m2,
        std::vector<Triangle>& output)
    {
        const Index v0 = face[0];
        const Index v1 = face[1];
        const Index v2 = face[2];

        switch (split_mask)
        {
        case 0:
            append_triangle(
                output,
                v0,
                v1,
                v2);
            break;

        // Edge 0: v0-v1.
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

        // Edge 1: v1-v2.
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

        // Edges 0 and 1.
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

        // Edge 2: v2-v0.
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

        // Edges 0 and 2.
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

        // Edges 1 and 2.
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

        // All three edges.
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

    static std::map<EdgeKey, EdgeAdjacency>
    build_edge_adjacency(
        const iso::mc::Mesh& mesh)
    {
        std::map<
            EdgeKey,
            EdgeAdjacency> adjacency;

        for (std::size_t fi = 0;
             fi < mesh.triangles.size();
             ++fi)
        {
            const Triangle& t =
                mesh.triangles[fi];

            adjacency[
                make_edge(
                    t[0],
                    t[1])]
                .faces.push_back(fi);

            adjacency[
                make_edge(
                    t[1],
                    t[2])]
                .faces.push_back(fi);

            adjacency[
                make_edge(
                    t[2],
                    t[0])]
                .faces.push_back(fi);
        }

        return adjacency;
    }

    static std::vector<bool>
    select_faces(
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
                squared_distance(
                    c,
                    center) <=
                radius_squared;
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

        for (const Point& p :
             mesh.vertices)
        {
            if (squared_distance(
                    p,
                    center) <=
                radius_squared)
            {
                ++count;
            }
        }

        return count;
    }

    static std::size_t count_faces_in_region(
        const iso::mc::Mesh& mesh,
        const Point& center,
        double radius)
    {
        const std::vector<bool> selected =
            select_faces(
                mesh,
                center,
                radius);

        return static_cast<std::size_t>(
            std::count(
                selected.begin(),
                selected.end(),
                true));
    }

public:
    Result refine(
        const iso::mc::Mesh& input,
        const ImplicitSurface& surface,
        const Point& center) const
    {
        Options options;

        return refine(
            input,
            surface,
            center,
            options);
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

        result.mesh =
            input;

        result.vertices_in_region =
            count_vertices_in_region(
                input,
                center,
                options.radius);

        result.faces_in_region =
            count_faces_in_region(
                input,
                center,
                options.radius);

        if (options.levels == 0 ||
            input.triangles.empty())
        {
            result.output_vertices =
                result.mesh.vertices.size();

            result.output_triangles =
                result.mesh.triangles.size();

            return result;
        }

        for (std::size_t level = 0;
             level < options.levels;
             ++level)
        {
            const iso::mc::Mesh current =
                result.mesh;

            const std::vector<bool> selected =
                select_faces(
                    current,
                    center,
                    options.radius);

            const std::size_t selected_count =
                static_cast<std::size_t>(
                    std::count(
                        selected.begin(),
                        selected.end(),
                        true));

            if (selected_count == 0)
            {
                break;
            }

            const auto adjacency =
                build_edge_adjacency(
                    current);

            std::map<EdgeKey, bool>
                split_edges;

            for (const auto& entry :
                 adjacency)
            {
                const auto& faces =
                    entry.second.faces;

                bool split =
                    false;

                for (std::size_t face_index :
                     faces)
                {
                    if (selected[face_index])
                    {
                        split = true;
                        break;
                    }
                }

                split_edges[entry.first] =
                    split;
            }

            iso::mc::Mesh next;

            next.vertices =
                current.vertices;

            next.triangles.reserve(
                current.triangles.size() * 2);

            EdgeSplitMap midpoint_cache;

            std::size_t projection_failures =
                0;

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
                    split_edges[e0];

                const bool split1 =
                    split_edges[e1];

                const bool split2 =
                    split_edges[e2];

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

            LevelStatistics statistics;

            statistics.level =
                level + 1;

            statistics.selected_faces =
                selected_count;

            statistics.output_vertices =
                next.vertices.size();

            statistics.output_triangles =
                next.triangles.size();

            statistics.new_vertices =
                next.vertices.size() -
                current.vertices.size();

            statistics.new_triangles =
                next.triangles.size() -
                current.triangles.size();

            statistics.projection_failures =
                projection_failures;

            result.levels.push_back(
                statistics);

            result.mesh =
                std::move(next);
        }

        result.output_vertices =
            result.mesh.vertices.size();

        result.output_triangles =
            result.mesh.triangles.size();

        return result;
    }
};