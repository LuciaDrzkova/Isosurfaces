#pragma once

#include "MarchingCubes.hpp"

#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace iso::field_export
{

/*
 * Write an ASCII PLY mesh with per-vertex scalar-field metadata.
 *
 * Properties:
 *   scalar_value       = f(x,y,z), evaluated at the exported vertex
 *   iso_residual       = f(x,y,z) - isovalue
 *   distance_to_center = Euclidean distance from the requested local center
 *   in_local_region    = 1 when the vertex lies within local_radius, else 0
 *   red/green/blue     = visualization colors for the spatial local-region mask
 *
 * The region flag describes spatial membership, not vertex provenance. In a
 * stitched mesh, seam vertices are shared by both patches and therefore do
 * not have a unique global/local source.
 */
template <class ScalarFunction>
void write_scalar_ply(
    const iso::mc::Mesh& mesh,
    const std::string& filename,
    ScalarFunction&& scalar_function,
    double isovalue,
    const iso::mc::Point& center,
    double local_radius)
{
    if (!(local_radius > 0.0) || !std::isfinite(local_radius))
        throw std::invalid_argument(
            "write_scalar_ply: local_radius must be finite and positive");

    std::ofstream out(filename);
    if (!out)
        throw std::runtime_error(
            "write_scalar_ply: could not open output file: " + filename);

    out << "ply\n"
        << "format ascii 1.0\n"
        << "comment Per-vertex implicit scalar-field measurements\n"
        << "comment scalar_value is the field evaluated at each mesh vertex\n"
        << "comment iso_residual is scalar_value minus the target isovalue\n"
        << "comment in_local_region is a spatial mask based on distance_to_center\n"
        << "element vertex " << mesh.vertices.size() << "\n"
        << "property double x\n"
        << "property double y\n"
        << "property double z\n"
        << "property uchar red\n"
        << "property uchar green\n"
        << "property uchar blue\n"
        << "property double scalar_value\n"
        << "property double iso_residual\n"
        << "property double distance_to_center\n"
        << "property uchar in_local_region\n"
        << "element face " << mesh.triangles.size() << "\n"
        << "property list uchar uint vertex_indices\n"
        << "end_header\n";

    out << std::setprecision(17);
    const double radius_squared = local_radius * local_radius;

    for (const auto& p : mesh.vertices)
    {
        const double dx = p[0] - center[0];
        const double dy = p[1] - center[1];
        const double dz = p[2] - center[2];
        const double distance = std::sqrt(dx * dx + dy * dy + dz * dz);
        const bool in_local_region = distance * distance <= radius_squared;
        const double value = scalar_function(p[0], p[1], p[2]);
        const double residual = value - isovalue;

        // Blue outside and orange inside provide a useful default view in
        // common PLY viewers. The numeric mask remains available for analysis.
        const unsigned int red = in_local_region ? 255u : 55u;
        const unsigned int green = in_local_region ? 135u : 145u;
        const unsigned int blue = in_local_region ? 35u : 235u;

        out << p[0] << ' ' << p[1] << ' ' << p[2] << ' '
            << red << ' ' << green << ' ' << blue << ' '
            << value << ' ' << residual << ' ' << distance << ' '
            << (in_local_region ? 1 : 0) << '\n';
    }

    for (const auto& t : mesh.triangles)
    {
        out << "3 "
            << static_cast<unsigned long long>(t[0]) << ' '
            << static_cast<unsigned long long>(t[1]) << ' '
            << static_cast<unsigned long long>(t[2]) << '\n';
    }

    if (!out)
        throw std::runtime_error(
            "write_scalar_ply: failed while writing output file: " + filename);
}



/*
 * Write a VTK legacy POLYDATA mesh with point-data arrays. Unlike arbitrary
 * extra PLY properties, these arrays are directly exposed by ParaView's VTK
 * reader and can be selected in the Coloring menu.
 */
template <class ScalarFunction>
void write_scalar_vtk(
    const iso::mc::Mesh& mesh,
    const std::string& filename,
    ScalarFunction&& scalar_function,
    double isovalue,
    const iso::mc::Point& center,
    double local_radius)
{
    if (!(local_radius > 0.0) || !std::isfinite(local_radius))
        throw std::invalid_argument(
            "write_scalar_vtk: local_radius must be finite and positive");

    std::ofstream out(filename);
    if (!out)
        throw std::runtime_error(
            "write_scalar_vtk: could not open output file: " + filename);

    out << "# vtk DataFile Version 3.0\n"
        << "Implicit scalar field on isosurface vertices\n"
        << "ASCII\n"
        << "DATASET POLYDATA\n"
        << "POINTS " << mesh.vertices.size() << " double\n"
        << std::setprecision(17);

    for (const auto& p : mesh.vertices)
        out << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';

    out << "POLYGONS " << mesh.triangles.size() << ' '
        << mesh.triangles.size() * 4 << '\n';
    for (const auto& t : mesh.triangles)
    {
        out << "3 "
            << static_cast<unsigned long long>(t[0]) << ' '
            << static_cast<unsigned long long>(t[1]) << ' '
            << static_cast<unsigned long long>(t[2]) << '\n';
    }

    out << "POINT_DATA " << mesh.vertices.size() << '\n';

    const auto write_scalar_array = [&](const char* name, auto value_at)
    {
        out << "SCALARS " << name << " double 1\n"
            << "LOOKUP_TABLE default\n";
        for (const auto& p : mesh.vertices)
            out << value_at(p) << '\n';
    };

    write_scalar_array("scalar_value", [&](const iso::mc::Point& p)
    {
        return scalar_function(p[0], p[1], p[2]);
    });
    write_scalar_array("iso_residual", [&](const iso::mc::Point& p)
    {
        return scalar_function(p[0], p[1], p[2]) - isovalue;
    });
    write_scalar_array("distance_to_center", [&](const iso::mc::Point& p)
    {
        const double dx = p[0] - center[0];
        const double dy = p[1] - center[1];
        const double dz = p[2] - center[2];
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    });

    out << "SCALARS in_local_region int 1\n"
        << "LOOKUP_TABLE default\n";
    const double radius_squared = local_radius * local_radius;
    for (const auto& p : mesh.vertices)
    {
        const double dx = p[0] - center[0];
        const double dy = p[1] - center[1];
        const double dz = p[2] - center[2];
        out << ((dx * dx + dy * dy + dz * dz <= radius_squared) ? 1 : 0)
            << '\n';
    }

    if (!out)
        throw std::runtime_error(
            "write_scalar_vtk: failed while writing output file: " + filename);
}

} // namespace iso::field_export
