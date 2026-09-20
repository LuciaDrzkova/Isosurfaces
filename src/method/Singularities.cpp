#include "method/Singularities.h"

#include <cmath>

using namespace std;

namespace iso {

// ---------------------------------------------------------------------------
// Your playground for the singularities research.
//
// What you get:
//   mesh     the remeshed triangle mesh (pmp::SurfaceMesh). You may modify it;
//            the modified mesh is what gets written to the output file.
//   surface  the implicit surface f, with
//              surface.value(p)     -> f(p)
//              surface.evaluate(p)  -> { f(p), grad f(p) } via autodiff
//              surface.project(p)   -> Newton projection onto f = 0
//   log      anything written here appears in the GUI and on the console.
//
// The body below is only a small example of using that API (counting vertices
// where the gradient nearly vanishes, i.e. candidate singular points of f = 0).
// Replace it with your own work.
// ---------------------------------------------------------------------------
void analyze_singularities(pmp::SurfaceMesh& mesh,
                           const ImplicitSurface& surface, ostream& log)
{
    constexpr double kGradientEpsilon = 1e-2;

    size_t candidates = 0;
    for (auto v : mesh.vertices())
    {
        const auto sample = surface.evaluate(mesh.position(v));
        if (norm(sample.gradient) < kGradientEpsilon)
        {
            ++candidates;
            log << "  candidate at vertex " << v.idx() << ", |grad f| = "
                << norm(sample.gradient) << '\n';
        }
    }
    log << candidates << " of " << mesh.n_vertices()
        << " vertices with |grad f| < " << kGradientEpsilon << '\n';
}

} // namespace iso
