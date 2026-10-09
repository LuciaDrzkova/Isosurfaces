#pragma once

#include "ImplicitSurface.hpp"
#include "MarchingCubes.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>

// Independent local voxel extraction used as the next-stage research path.
// The local box is aligned to the original MC lattice, so level k uses the
// same coarse box but 2^k times more samples per original grid cell.
// Stitching into the global surface is deliberately kept separate from this
// extractor; LocalUnfolder remains the conforming triangle-transition baseline.
class LocalMarchingCubes
{
public:
    struct Options
    {
        double radius = 0.3;
        std::size_t level = 1;
        std::size_t padding_cells = 1;
        iso::mc::Bounds sampling_bounds{{0.0, 0.0, 0.0},
                                         {0.0, 0.0, 0.0}};
        std::size_t base_resolution = 0;
        double isovalue = 0.0;
    };

    struct Result
    {
        iso::mc::Mesh mesh;
        iso::mc::Bounds local_bounds{{0.0, 0.0, 0.0},
                                      {0.0, 0.0, 0.0}};

        std::size_t base_cells_x = 0;
        std::size_t base_cells_y = 0;
        std::size_t base_cells_z = 0;
        std::size_t resolution_x = 0;
        std::size_t resolution_y = 0;
        std::size_t resolution_z = 0;
        std::size_t refinement_factor = 1;
        double base_spacing = 0.0;
        double refined_spacing = 0.0;

        std::size_t boundary_edges = 0;
        std::size_t boundary_vertices = 0;
        std::size_t nonmanifold_edges = 0;
        std::size_t selected_cells = 0;
        std::size_t active_cells = 0;
        double extraction_time_ms = 0.0;

        // Kept for source compatibility with earlier experiments. This
        // extractor deliberately does not weld coincident vertices; value is 0.
        std::size_t welded_vertices = 0;
    };

private:
    using Index = iso::mc::Index;
    using EdgeKey = std::pair<Index, Index>;

    static EdgeKey make_edge(Index a, Index b)
    {
        if (a > b)
            std::swap(a, b);
        return {a, b};
    }

    static std::size_t clamp_index(long long value, std::size_t maximum)
    {
        if (value < 0)
            return 0;
        return std::min(static_cast<std::size_t>(value), maximum);
    }

    static std::size_t floor_index(
        double coordinate,
        double global_min,
        double spacing,
        std::size_t maximum)
    {
        const long long value = static_cast<long long>(
            std::floor((coordinate - global_min) / spacing + 1e-12));
        return clamp_index(value, maximum);
    }

    static std::size_t ceil_index(
        double coordinate,
        double global_min,
        double spacing,
        std::size_t maximum)
    {
        const long long value = static_cast<long long>(
            std::ceil((coordinate - global_min) / spacing - 1e-12));
        return clamp_index(value, maximum);
    }

    static void boundary_statistics(
        const iso::mc::Mesh& mesh,
        std::size_t& boundary_edges,
        std::size_t& boundary_vertices,
        std::size_t& nonmanifold_edges)
    {
        std::map<EdgeKey, std::size_t> uses;
        for (const auto& t : mesh.triangles)
        {
            ++uses[make_edge(t[0], t[1])];
            ++uses[make_edge(t[1], t[2])];
            ++uses[make_edge(t[2], t[0])];
        }

        std::set<Index> boundary_vertex_set;
        boundary_edges = 0;
        nonmanifold_edges = 0;

        for (const auto& [edge, count] : uses)
        {
            if (count == 1)
            {
                ++boundary_edges;
                boundary_vertex_set.insert(edge.first);
                boundary_vertex_set.insert(edge.second);
            }
            else if (count > 2)
            {
                ++nonmanifold_edges;
            }
        }

        boundary_vertices = boundary_vertex_set.size();
    }

public:
    Result extract(
        const ImplicitSurface& surface,
        const iso::mc::Point& center,
        const Options& options) const
    {
        if (!(options.radius > 0.0))
            throw std::invalid_argument(
                "LocalMarchingCubes radius must be positive");

        if (options.base_resolution < 2)
            throw std::invalid_argument(
                "LocalMarchingCubes requires base_resolution >= 2");

        if (!(options.sampling_bounds.min[0] < options.sampling_bounds.max[0] &&
              options.sampling_bounds.min[1] < options.sampling_bounds.max[1] &&
              options.sampling_bounds.min[2] < options.sampling_bounds.max[2]))
        {
            throw std::invalid_argument(
                "LocalMarchingCubes received invalid sampling bounds");
        }

        if (options.level >= sizeof(std::size_t) * 8 - 2)
            throw std::invalid_argument(
                "LocalMarchingCubes refinement level is too large");

        const std::size_t base_cells = options.base_resolution - 1;
        const double hx =
            (options.sampling_bounds.max[0] - options.sampling_bounds.min[0]) /
            static_cast<double>(base_cells);
        const double hy =
            (options.sampling_bounds.max[1] - options.sampling_bounds.min[1]) /
            static_cast<double>(base_cells);
        const double hz =
            (options.sampling_bounds.max[2] - options.sampling_bounds.min[2]) /
            static_cast<double>(base_cells);

        const std::size_t rx0 = floor_index(
            center[0] - options.radius,
            options.sampling_bounds.min[0], hx, base_cells);
        const std::size_t ry0 = floor_index(
            center[1] - options.radius,
            options.sampling_bounds.min[1], hy, base_cells);
        const std::size_t rz0 = floor_index(
            center[2] - options.radius,
            options.sampling_bounds.min[2], hz, base_cells);
        const std::size_t rx1 = ceil_index(
            center[0] + options.radius,
            options.sampling_bounds.min[0], hx, base_cells);
        const std::size_t ry1 = ceil_index(
            center[1] + options.radius,
            options.sampling_bounds.min[1], hy, base_cells);
        const std::size_t rz1 = ceil_index(
            center[2] + options.radius,
            options.sampling_bounds.min[2], hz, base_cells);

        const std::size_t x0i = rx0 > options.padding_cells
            ? rx0 - options.padding_cells : 0;
        const std::size_t y0i = ry0 > options.padding_cells
            ? ry0 - options.padding_cells : 0;
        const std::size_t z0i = rz0 > options.padding_cells
            ? rz0 - options.padding_cells : 0;
        const std::size_t x1i = std::min(
            base_cells, rx1 + options.padding_cells);
        const std::size_t y1i = std::min(
            base_cells, ry1 + options.padding_cells);
        const std::size_t z1i = std::min(
            base_cells, rz1 + options.padding_cells);

        const std::size_t coarse_x = std::max<std::size_t>(1, x1i - x0i);
        const std::size_t coarse_y = std::max<std::size_t>(1, y1i - y0i);
        const std::size_t coarse_z = std::max<std::size_t>(1, z1i - z0i);

        const std::size_t factor = std::size_t(1) << options.level;
        const std::size_t nx = coarse_x * factor + 1;
        const std::size_t ny = coarse_y * factor + 1;
        const std::size_t nz = coarse_z * factor + 1;

        Result result;
        result.local_bounds = iso::mc::Bounds{
            {options.sampling_bounds.min[0] + x0i * hx,
             options.sampling_bounds.min[1] + y0i * hy,
             options.sampling_bounds.min[2] + z0i * hz},
            {options.sampling_bounds.min[0] + x1i * hx,
             options.sampling_bounds.min[1] + y1i * hy,
             options.sampling_bounds.min[2] + z1i * hz}};
        result.base_cells_x = coarse_x;
        result.base_cells_y = coarse_y;
        result.base_cells_z = coarse_z;
        result.resolution_x = nx;
        result.resolution_y = ny;
        result.resolution_z = nz;
        result.refinement_factor = factor;
        result.base_spacing = std::min({hx, hy, hz});
        result.refined_spacing = result.base_spacing /
                                 static_cast<double>(factor);

        const double radius2 = options.radius * options.radius;
        for (std::size_t z = 0; z < coarse_z; ++z)
        {
            for (std::size_t y = 0; y < coarse_y; ++y)
            {
                for (std::size_t x = 0; x < coarse_x; ++x)
                {
                    const iso::mc::Point p{
                        result.local_bounds.min[0] + (x + 0.5) * hx,
                        result.local_bounds.min[1] + (y + 0.5) * hy,
                        result.local_bounds.min[2] + (z + 0.5) * hz};
                    const double dx = p[0] - center[0];
                    const double dy = p[1] - center[1];
                    const double dz = p[2] - center[2];
                    if (dx * dx + dy * dy + dz * dz <= radius2)
                        ++result.selected_cells;
                }
            }
        }

        iso::mc::Options mc_options;
        mc_options.nx = nx;
        mc_options.ny = ny;
        mc_options.nz = nz;
        mc_options.isovalue = options.isovalue;

        const auto start = std::chrono::steady_clock::now();
        result.mesh = iso::mc::extract(
            [&](double x, double y, double z)
            {
                return surface.eval(x, y, z);
            },
            result.local_bounds,
            mc_options);
        const auto end = std::chrono::steady_clock::now();

        result.extraction_time_ms =
            std::chrono::duration<double, std::milli>(end - start).count();

        // Important: do not weld equal-position vertices here. At the cone
        // apex and other grid-aligned zeroes, that merges otherwise distinct
        // Marching-Cubes fans and drops degenerate triangles. The resulting
        // pinched/non-manifold topology prevents the seam loops from matching.
        // Mesh export should split disconnected fans if the viewer's mesh
        // representation cannot encode the singular vertex directly.
        result.active_cells = result.mesh.triangles.size();

        boundary_statistics(
            result.mesh,
            result.boundary_edges,
            result.boundary_vertices,
            result.nonmanifold_edges);

        return result;
    }
};
