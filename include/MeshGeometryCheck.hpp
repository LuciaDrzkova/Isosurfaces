#pragma once

#include "MarchingCubes.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

// Exact geometric-error check for the double cone x^2 + y^2 - z^2 = 0.
// This is an analytic distance, not the first-order |f|/|grad f| estimate,
// so it remains well-defined at the singular apex.
namespace iso::geometry_check
{

using Point = iso::mc::Point;
using Mesh = iso::mc::Mesh;
using Triangle = iso::mc::Triangle;

struct Statistics
{
    std::size_t count = 0;
    double mean = 0.0;
    double rms = 0.0;
    double p95 = 0.0;
    double p99 = 0.0;
    double maximum = 0.0;
};

struct Report
{
    Statistics all_vertices;
    Statistics all_triangle_samples;
    Statistics local_vertices;
    Statistics local_triangle_samples;
};

inline double distance_to_double_cone(const Point& p)
{
    const double radial = std::hypot(p[0], p[1]);
    return std::abs(radial - std::abs(p[2])) / std::sqrt(2.0);
}

inline double squared_distance(const Point& a, const Point& b)
{
    const double dx = a[0] - b[0];
    const double dy = a[1] - b[1];
    const double dz = a[2] - b[2];
    return dx * dx + dy * dy + dz * dz;
}

inline Statistics summarize(std::vector<double> values)
{
    Statistics result;
    result.count = values.size();
    if (values.empty())
        return result;

    double sum = 0.0;
    double sum_squared = 0.0;
    for (const double value : values)
    {
        sum += value;
        sum_squared += value * value;
    }

    result.mean = sum / static_cast<double>(values.size());
    result.rms = std::sqrt(sum_squared / static_cast<double>(values.size()));
    result.maximum = *std::max_element(values.begin(), values.end());

    std::sort(values.begin(), values.end());
    const auto percentile = [&values](double p)
    {
        const std::size_t index = static_cast<std::size_t>(
            std::ceil(p * static_cast<double>(values.size()))) - 1;
        return values[std::min(index, values.size() - 1)];
    };
    result.p95 = percentile(0.95);
    result.p99 = percentile(0.99);
    return result;
}

inline Report measure(const Mesh& mesh,
                      const Point& center,
                      double local_radius,
                      std::size_t triangle_subdivisions = 4)
{
    Report report;
    std::vector<double> all_vertex_errors;
    std::vector<double> local_vertex_errors;
    std::vector<double> all_face_errors;
    std::vector<double> local_face_errors;

    const double radius_squared = local_radius * local_radius;
    const auto is_local = [&](const Point& p)
    {
        return squared_distance(p, center) <= radius_squared;
    };

    // Unique mesh vertices: a direct check that catches vertex displacement.
    all_vertex_errors.reserve(mesh.vertices.size());
    local_vertex_errors.reserve(mesh.vertices.size());
    for (const Point& p : mesh.vertices)
    {
        const double error = distance_to_double_cone(p);
        all_vertex_errors.push_back(error);
        if (is_local(p))
            local_vertex_errors.push_back(error);
    }

    // Sample the interior of every planar triangle on a barycentric grid.
    // Subdivision 4 gives 15 samples per triangle, including its corners.
    // This detects planar-facet deviation that vertex-only checks miss.
    const std::size_t n = std::max<std::size_t>(1, triangle_subdivisions);
    all_face_errors.reserve(mesh.triangles.size() *
                            ((n + 1) * (n + 2) / 2));

    for (const Triangle& t : mesh.triangles)
    {
        const Point& a = mesh.vertices[t[0]];
        const Point& b = mesh.vertices[t[1]];
        const Point& c = mesh.vertices[t[2]];

        for (std::size_t i = 0; i <= n; ++i)
        {
            for (std::size_t j = 0; j <= n - i; ++j)
            {
                const std::size_t k = n - i - j;
                const double inv_n = 1.0 / static_cast<double>(n);
                const Point p{
                    (static_cast<double>(i) * a[0] +
                     static_cast<double>(j) * b[0] +
                     static_cast<double>(k) * c[0]) * inv_n,
                    (static_cast<double>(i) * a[1] +
                     static_cast<double>(j) * b[1] +
                     static_cast<double>(k) * c[1]) * inv_n,
                    (static_cast<double>(i) * a[2] +
                     static_cast<double>(j) * b[2] +
                     static_cast<double>(k) * c[2]) * inv_n};

                const double error = distance_to_double_cone(p);
                all_face_errors.push_back(error);
                if (is_local(p))
                    local_face_errors.push_back(error);
            }
        }
    }

    report.all_vertices = summarize(std::move(all_vertex_errors));
    report.all_triangle_samples = summarize(std::move(all_face_errors));
    report.local_vertices = summarize(std::move(local_vertex_errors));
    report.local_triangle_samples = summarize(std::move(local_face_errors));
    return report;
}

inline void print_statistics(const std::string& label, const Statistics& s)
{
    std::cout << "  " << label << " (" << s.count << " samples)\n"
              << "    mean distance : " << s.mean << '\n'
              << "    RMS distance  : " << s.rms << '\n'
              << "    95th percentile: " << s.p95 << '\n'
              << "    99th percentile: " << s.p99 << '\n'
              << "    maximum       : " << s.maximum << '\n';
}

inline void print_report(const std::string& label, const Report& r)
{
    std::cout << "\nGEOMETRIC ERROR TO EXACT DOUBLE CONE: " << label << '\n';
    std::cout << "  Distances are Euclidean and reported in model units.\n";
    print_statistics("all mesh vertices", r.all_vertices);
    print_statistics("all triangle samples", r.all_triangle_samples);
    print_statistics("vertices inside local radius", r.local_vertices);
    print_statistics("triangle samples inside local radius",
                     r.local_triangle_samples);
}

} // namespace iso::geometry_check
