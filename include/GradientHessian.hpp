#pragma once

#include "ImplicitSurface.hpp"

#include <Eigen/Dense>
#include <pmp/types.h>

#include <algorithm>
#include <cmath>
#include <limits>

// Numerical derivative helper used by code paths that operate on arbitrary
// ImplicitSurface implementations. The specialized cone/quadric derivative
// routines live on ParameterizedConeQuadric itself.
class GradientHessianCalculator
{
private:
    const ImplicitSurface& surface_;

    double step(double x, double y, double z) const
    {
        const double scale = std::max(
            1.0,
            std::sqrt(x * x + y * y + z * z));
        return 1e-6 * scale;
    }

public:
    explicit GradientHessianCalculator(const ImplicitSurface& surface)
        : surface_(surface)
    {
    }

    pmp::Point computeGradient(const pmp::Point& p) const
    {
        const double x = static_cast<double>(p[0]);
        const double y = static_cast<double>(p[1]);
        const double z = static_cast<double>(p[2]);
        const double h = step(x, y, z);

        const double dx =
            (surface_.eval(x + h, y, z) -
             surface_.eval(x - h, y, z)) /
            (2.0 * h);
        const double dy =
            (surface_.eval(x, y + h, z) -
             surface_.eval(x, y - h, z)) /
            (2.0 * h);
        const double dz =
            (surface_.eval(x, y, z + h) -
             surface_.eval(x, y, z - h)) /
            (2.0 * h);

        return pmp::Point(
            static_cast<float>(dx),
            static_cast<float>(dy),
            static_cast<float>(dz));
    }

    Eigen::Matrix3d computeHessian(const pmp::Point& p) const
    {
        const double x = static_cast<double>(p[0]);
        const double y = static_cast<double>(p[1]);
        const double z = static_cast<double>(p[2]);
        const double h = 1e-4 *
            std::max(1.0, std::sqrt(x * x + y * y + z * z));
        const double h2 = h * h;

        const double f000 = surface_.eval(x, y, z);
        Eigen::Matrix3d H = Eigen::Matrix3d::Zero();

        H(0, 0) =
            (surface_.eval(x + h, y, z) - 2.0 * f000 +
             surface_.eval(x - h, y, z)) /
            h2;
        H(1, 1) =
            (surface_.eval(x, y + h, z) - 2.0 * f000 +
             surface_.eval(x, y - h, z)) /
            h2;
        H(2, 2) =
            (surface_.eval(x, y, z + h) - 2.0 * f000 +
             surface_.eval(x, y, z - h)) /
            h2;

        H(0, 1) = H(1, 0) =
            (surface_.eval(x + h, y + h, z) -
             surface_.eval(x + h, y - h, z) -
             surface_.eval(x - h, y + h, z) +
             surface_.eval(x - h, y - h, z)) /
            (4.0 * h2);

        H(0, 2) = H(2, 0) =
            (surface_.eval(x + h, y, z + h) -
             surface_.eval(x + h, y, z - h) -
             surface_.eval(x - h, y, z + h) +
             surface_.eval(x - h, y, z - h)) /
            (4.0 * h2);

        H(1, 2) = H(2, 1) =
            (surface_.eval(x, y + h, z + h) -
             surface_.eval(x, y + h, z - h) -
             surface_.eval(x, y - h, z + h) +
             surface_.eval(x, y - h, z - h)) /
            (4.0 * h2);

        return H;
    }

    // Newton projection onto f(x,y,z)=0.
    bool projectOntoSurface(
        pmp::Point& p,
        double tolerance = 1e-10,
        std::size_t iterations = 12) const
    {
        for (std::size_t i = 0; i < iterations; ++i)
        {
            const double x = static_cast<double>(p[0]);
            const double y = static_cast<double>(p[1]);
            const double z = static_cast<double>(p[2]);
            const double value = surface_.eval(x, y, z);

            if (std::abs(value) <= tolerance)
                return true;

            const pmp::Point g = computeGradient(p);
            const double gx = static_cast<double>(g[0]);
            const double gy = static_cast<double>(g[1]);
            const double gz = static_cast<double>(g[2]);
            const double norm2 = gx * gx + gy * gy + gz * gz;

            if (norm2 <= std::numeric_limits<double>::epsilon())
                return false;

            p[0] = static_cast<float>(x - value * gx / norm2);
            p[1] = static_cast<float>(y - value * gy / norm2);
            p[2] = static_cast<float>(z - value * gz / norm2);
        }

        return std::abs(surface_.eval(
            static_cast<double>(p[0]),
            static_cast<double>(p[1]),
            static_cast<double>(p[2]))) <= tolerance;
    }
};
