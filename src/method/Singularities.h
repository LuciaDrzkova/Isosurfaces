// Research code for the singularities work (implementation in Singularities.cpp).

#pragma once

#include "method/ImplicitSurface.h"

#include <pmp/surface_mesh.h>

#include <ostream>

namespace iso {

//! Starting point for the singularities research. Registered as the
//! experiment "singularities" (src/Experiments.cpp).
void analyze_singularities(pmp::SurfaceMesh& mesh,
                           const ImplicitSurface& surface, std::ostream& log);

} // namespace iso
