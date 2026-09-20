// The remeshing algorithm: project the mesh onto the surface, then remove small
// and skinny triangles by collapsing their shortest edge.

#pragma once

#include "method/ImplicitSurface.h"

#include <pmp/surface_mesh.h>

#include <vector>

namespace iso {

struct RemeshOptions
{
    //! Number of outer passes. The target area is re-measured on each one.
    int iterations = 6;

    //! A triangle is a collapse candidate if any of its angles is below this
    //! (degrees) ...
    float min_angle = 25.0f;

    //! ... or if its area is below (median area / area_divisor).
    float area_divisor = 1.5f;
};

struct RemeshResult
{
    //! Edge collapses performed in each outer pass
    std::vector<int> collapses_per_iteration;
};

//! Coarsens and regularises \p mesh in place. Small or skinny triangles are
//! removed by collapsing their shortest edge; the surviving vertex is moved to
//! the centroid of its neighbours and projected back onto \p surface.
//! Boundary vertices are never touched.
//!
//! \p mesh must be a triangle mesh. Deleted elements are removed
//! (garbage collection) before returning.
RemeshResult remesh(pmp::SurfaceMesh& mesh, const ImplicitSurface& surface,
                    const RemeshOptions& options = {});

//! Projects every vertex of \p mesh onto \p surface.
void project_onto_surface(pmp::SurfaceMesh& mesh,
                          const ImplicitSurface& surface);

} // namespace iso
