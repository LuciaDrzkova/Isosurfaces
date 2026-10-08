#pragma once

#include <Eigen/Dense>

#include <cmath>

// Base interface for scalar implicit surfaces f(x,y,z)=0.
class ImplicitSurface
{
public:
    virtual ~ImplicitSurface() = default;

    virtual double eval(double x, double y, double z) const = 0;
};

// -----------------------------------------------------------------------------
// Parameterized quadric used throughout the original experiments:
//     f(x,y,z) = x^2 + y^2 - z^2 - c
//
// c < 0 : two-sheet hyperboloid
// c = 0 : double cone with a non-degenerate Morse saddle at the origin
// c > 0 : one-sheet hyperboloid
// -----------------------------------------------------------------------------
class ParameterizedConeQuadric : public ImplicitSurface
{
private:
    double c_ = 0.0;

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

    double eval(double x, double y, double z) const override
    {
        return x * x + y * y - z * z - c_;
    }

    Eigen::Vector3d gradient(double x, double y, double z) const
    {
        return Eigen::Vector3d(
            2.0 * x,
            2.0 * y,
            -2.0 * z);
    }

    Eigen::Matrix3d hessian(double, double, double) const
    {
        Eigen::Matrix3d H;
        H << 2.0, 0.0, 0.0,
             0.0, 2.0, 0.0,
             0.0, 0.0, -2.0;
        return H;
    }

    double hessianDeterminant(double x, double y, double z) const
    {
        return hessian(x, y, z).determinant();
    }

    template <typename T>
    T evalTemplate(T x, T y, T z) const
    {
        return x * x + y * y - z * z - c_;
    }
};

// -----------------------------------------------------------------------------
// A genuinely non-quadric family used to test that singularity detection and
// local refinement do not depend on the ParameterizedConeQuadric type.
//
//     f(x,y,z) = x^2 + y^2 - z^2
//                + alpha * (x^4 + y^4 + z^4) - c
//
// At c = 0 the origin is still a non-degenerate Morse saddle, but away from
// the origin the surface is no longer quadratic. SingularityDetector and
// LocalUnfolder therefore exercise their generic numerical derivative paths.
// -----------------------------------------------------------------------------
class QuarticSaddleSurface : public ImplicitSurface
{
private:
    double c_ = 0.0;
    double alpha_ = 0.05;

public:
    explicit QuarticSaddleSurface(double c = 0.0,
                                  double alpha = 0.05)
        : c_(c),
          alpha_(alpha)
    {
    }

    double getParameter() const
    {
        return c_;
    }

    double getAlpha() const
    {
        return alpha_;
    }

    double eval(double x, double y, double z) const override
    {
        const double x2 = x * x;
        const double y2 = y * y;
        const double z2 = z * z;

        return x2 + y2 - z2 +
               alpha_ * (x2 * x2 + y2 * y2 + z2 * z2) - c_;
    }
};

// -----------------------------------------------------------------------------
// Closed regular surface used as a topology/control family.
// There is no critical point on the level set for a non-zero radius, so the
// local-singularity stage should correctly skip it.
// -----------------------------------------------------------------------------
class SphereImplicitSurface : public ImplicitSurface
{
private:
    double radius_ = 1.0;

public:
    explicit SphereImplicitSurface(double radius = 1.0)
        : radius_(radius)
    {
        if (!(radius_ > 0.0))
            radius_ = 1.0;
    }

    double eval(double x, double y, double z) const override
    {
        return x * x + y * y + z * z - radius_ * radius_;
    }

    double radius() const
    {
        return radius_;
    }
};
