#pragma once

#include <Eigen/Dense>
#include <autodiff/forward/dual.hpp>

// -----------------------------------------------------------------------------
// Base interface for implicit surfaces
// -----------------------------------------------------------------------------

class ImplicitSurface
{
public:
    virtual ~ImplicitSurface() = default;

    virtual double eval(
        double x,
        double y,
        double z) const = 0;
};

// -----------------------------------------------------------------------------
// Parameterized cone / quadric
//
//     f(x,y,z) = x^2 + y^2 - z^2 - c
//
// For c = 0:
//
//     f(x,y,z) = x^2 + y^2 - z^2
//
// The singular point is:
//
//     p = (0,0,0)
//
// Gradient:
//
//     grad f = (2x, 2y, -2z)
//
// Hessian:
//
//     H = [ 2  0  0 ]
//         [ 0  2  0 ]
//         [ 0  0 -2 ]
//
// Hessian determinant:
//
//     det(H) = -8
//
// Therefore the origin is a NON-DEGENERATE critical point.
// -----------------------------------------------------------------------------

class ParameterizedConeQuadric : public ImplicitSurface
{
private:
    double c_;

public:
    explicit ParameterizedConeQuadric(double c = 0.0)
        : c_(c)
    {
    }

    void setParameter(double c)
    {
        c_ = c;
    }

    double getParameter() const
    {
        return c_;
    }

    // -------------------------------------------------------------------------
    // Implicit function
    // -------------------------------------------------------------------------

    double eval(
        double x,
        double y,
        double z) const override
    {
        return x * x + y * y - z * z - c_;
    }

    // -------------------------------------------------------------------------
    // Exact analytic gradient
    // -------------------------------------------------------------------------

    Eigen::Vector3d gradient(
        double x,
        double y,
        double z) const
    {
        return Eigen::Vector3d(
            2.0 * x,
            2.0 * y,
            -2.0 * z);
    }

    // -------------------------------------------------------------------------
    // Exact analytic Hessian
    //
    // It does not depend on x, y or z for this quadratic.
    // -------------------------------------------------------------------------

    Eigen::Matrix3d hessian(
        double,
        double,
        double) const
    {
        Eigen::Matrix3d H;

        H << 2.0, 0.0, 0.0,
             0.0, 2.0, 0.0,
             0.0, 0.0, -2.0;

        return H;
    }

    // -------------------------------------------------------------------------
    // Exact Hessian determinant
    // -------------------------------------------------------------------------

    double hessianDeterminant(
        double x,
        double y,
        double z) const
    {
        return hessian(x, y, z).determinant();
    }

    // -------------------------------------------------------------------------
    // Template version retained for the existing autodiff infrastructure.
    // -------------------------------------------------------------------------

    template <typename T>
    T evalTemplate(
        T x,
        T y,
        T z) const
    {
        return x * x + y * y - z * z - c_;
    }
};