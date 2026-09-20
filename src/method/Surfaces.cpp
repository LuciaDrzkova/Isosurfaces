// ===========================================================================
// THE IMPLICIT SURFACES  -  add your own formulas here.
//
// A surface is the set of points (x, y, z) where a function f(x, y, z) is zero.
// To add one:
//
//   1. Write the function below. Use  sin cos tan exp sqrt abs atan2  and
//      sq(v) for "v squared". Write numbers normally (0.5, 2, ...).
//   2. Add one line to the list at the bottom:
//          {"name", "text shown in the program", function},
//
// That's all: the new surface shows up in the window and on the command line,
// and its gradient is computed automatically.
//
// Example (a torus):
//
//     dual torus(dual x, dual y, dual z)
//     {
//         return sq(sqrt(x * x + y * y) - 2) + z * z - 0.5;
//     }
//     ...
//     {"torus", "(sqrt(x^2 + y^2) - 2)^2 + z^2 = 0.5", torus},
// ===========================================================================

#include "method/ImplicitSurface.h"

using namespace std;
using autodiff::dual;

namespace iso {

namespace {

// v squared
dual sq(dual v)
{
    return v * v;
}

// --- The formulas ---------------------------------------------------------

dual ghost(dual x, dual y, dual z)
{
    return x * x * x + y * y + z * z - 0.5;
}

dual sphere(dual x, dual y, dual z)
{
    return x * x + y * y + z * z - 5.0;
}

dual gyroid(dual x, dual y, dual z)
{
    return sin(x) * cos(y) + sin(y) * cos(z) + sin(z) * cos(x) + 1;
}

dual julia(dual x, dual y, dual z)
{
    return x * x * x - 3 * x * y * y + z * z - 1;
}

dual sinusoidal(dual x, dual y, dual z)
{
    return sin(3 * x) + sin(3 * y) - z - 1;
}

dual turbine(dual x, dual y, dual z)
{
    return x * x + y * y + (z - 1) * (z - 1) -
           sq(0.7 + 0.3 * sin(5 * atan2(y, x)));
}

dual blobby(dual x, dual y, dual z)
{
    return x * x + y * y + z * z + 0.1 * cos(4 * x) + 0.1 * cos(4 * y) +
           0.1 * cos(4 * z) - 5;
}

dual wavy_network(dual x, dual y, dual z)
{
    return -(cos(x + y) * cos(y + z) * cos(z + x) - 0.5);
}

dual frozen_lattice(dual x, dual y, dual z)
{
    return abs(sin(x) + sin(y) + sin(z)) - 1.2 + 0.5 * cos(3 * x * y * z);
}

dual lotus(dual x, dual y, dual z)
{
    return sq(x * x + y * y + z * z - 1) -
           0.3 * cos(6 * atan2(y, x)) * exp(-z * z);
}

dual wavy_shell(dual x, dual y, dual z)
{
    return x * x + y * y + z * z + 0.81 * sin(6 * (x + y + z)) - 1;
}

} // namespace

// --- The list of all surfaces (name, description, formula) ----------------

const vector<ImplicitSurface>& builtin_surfaces()
{
    static const vector<ImplicitSurface> surfaces = {
        {"ghost", "x^3 + y^2 + z^2 = 0.5", ghost},
        {"sphere", "x^2 + y^2 + z^2 = 5", sphere},
        {"gyroid", "sin x cos y + sin y cos z + sin z cos x = -1", gyroid},
        {"julia", "x^3 - 3xy^2 + z^2 = 1", julia},
        {"sinusoidal", "sin 3x + sin 3y - z = 1", sinusoidal},
        {"turbine", "x^2 + y^2 + (z-1)^2 = (0.7 + 0.3 sin(5 atan2(y,x)))^2",
         turbine},
        {"blobby", "x^2 + y^2 + z^2 + 0.1(cos 4x + cos 4y + cos 4z) = 5",
         blobby},
        {"wavy_network", "cos(x+y) cos(y+z) cos(z+x) = 0.5", wavy_network},
        {"frozen_lattice", "|sin x + sin y + sin z| + 0.5 cos(3xyz) = 1.2",
         frozen_lattice},
        {"lotus", "(x^2 + y^2 + z^2 - 1)^2 = 0.3 cos(6 atan2(y,x)) exp(-z^2)",
         lotus},
        {"wavy_shell", "x^2 + y^2 + z^2 + 0.81 sin(6(x+y+z)) = 1", wavy_shell},
    };
    return surfaces;
}

// Look a surface up by its name; nullptr if there is none.
const ImplicitSurface* find_surface(const string& name)
{
    for (const auto& s : builtin_surfaces())
        if (s.name() == name)
            return &s;
    return nullptr;
}

} // namespace iso
