// Measurements of a mesh: triangle angles, median area, distance to the surface,
// number of skinny triangles, ...

#pragma once

#include "method/ImplicitSurface.h"

#include <pmp/surface_mesh.h>

#include <array>
#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

namespace iso {

//! Median of \p values (mean of the two middle values for even counts);
//! 0 for an empty vector.
float median(std::vector<float> values);

//! Interior angles in degrees at p0, p1, p2. Degenerate triangles (a repeated
//! point) yield zeros for the undefined angles instead of NaN.
std::array<float, 3> triangle_angles(const pmp::Point& p0, const pmp::Point& p1,
                                     const pmp::Point& p2);

//! Quality summary of a triangle mesh with respect to an implicit surface.
struct MeshMetrics
{
    std::size_t n_vertices = 0;
    std::size_t n_faces = 0;

    float median_area = 0;
    std::size_t faces_below_median_area = 0;
    //! L2 norm of the face areas' deviation from the median area
    float area_deviation = 0;

    float min_angle_threshold = 0; //!< degrees, echoed from the request
    std::size_t faces_with_small_angle = 0;

    double avg_distance = 0; //!< mean |f(v)| over all vertices
    double max_distance = 0; //!< max  |f(v)|
    float surface_tolerance = 0;
    std::size_t vertices_off_surface = 0; //!< |f(v)| > surface_tolerance

    //! Share of interior vertices with valence 6, in percent
    float regular_vertex_percentage = 0;

    double computation_time = 0; //!< seconds spent computing these metrics
};

MeshMetrics compute_mesh_metrics(const pmp::SurfaceMesh& mesh,
                                 const ImplicitSurface& surface,
                                 float min_angle_threshold,
                                 float surface_tolerance = 1e-4f);

void print_metrics(std::ostream& out, const std::string& label,
                   const MeshMetrics& m);

} // namespace iso
