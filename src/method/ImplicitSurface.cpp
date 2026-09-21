// How the value, gradient and projection are computed.
// To add a surface you don't need to change anything here: see Surfaces.cpp.

#include "method/ImplicitSurface.h"

#include <cmath>

using namespace std;

namespace iso {

using autodiff::dual;
using autodiff::at;
using autodiff::derivatives;
using autodiff::wrt;

double ImplicitSurface::value(const pmp::Point& p) const
{
    return static_cast<double>(f_(p[0], p[1], p[2]));
}

ImplicitSurface::Sample ImplicitSurface::evaluate(const pmp::Point& p) const
{
    dual x = p[0], y = p[1], z = p[2];

    // derivatives() returns {f(p), df/dvar}; seed one variable per pass
    const auto fx = derivatives(f_, wrt(x), at(x, y, z));
    const auto fy = derivatives(f_, wrt(y), at(x, y, z));
    const auto fz = derivatives(f_, wrt(z), at(x, y, z));

    return {fx[0], pmp::Point(static_cast<pmp::Scalar>(fx[1]),
                              static_cast<pmp::Scalar>(fy[1]),
                              static_cast<pmp::Scalar>(fz[1]))};
}

pmp::Point ImplicitSurface::project(pmp::Point p, int max_iterations,
                                    float tolerance) const
{
    for (int i = 0; i < max_iterations; ++i)
    {
        const Sample s = evaluate(p);
        const float f = static_cast<float>(s.value);
        const float grad_norm2 = pmp::dot(s.gradient, s.gradient);

        if (grad_norm2 < 1e-10f)
            break; // critical point: no direction to move in

        // Newton step towards f = 0 along the gradient
        p = p - (f / grad_norm2) * s.gradient;

        // Checked after the step (with the value from before it): once the
        // point is within tolerance, one more step is still taken. This is what
        // the original program did, so results stay identical.
        if (abs(f) < tolerance)
            break;
    }
    return p;
}

} // namespace iso
