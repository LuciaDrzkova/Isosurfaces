#pragma once

#include <Eigen/Dense>

// Base interface for scalar implicit surfaces f(x,y,z)=0.
class ImplicitSurface
{
public:
    virtual ~ImplicitSurface() = default;

    virtual double eval(
        double x,
        double y,
        double z) const = 0;
};

// Parameterized quadric used throughout the experiments:
//     f(x,y,z) = x^2 + y^2 - z^2 - c
//
// c < 0 : two-sheet hyperboloid
// c = 0 : double cone with a non-degenerate Morse saddle at the origin
// c > 0 : one-sheet hyperboloid
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

    double eval(
        double x,
        double y,
        double z) const override
    {
        return x * x + y * y - z * z - c_;
    }

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

    double hessianDeterminant(
        double x,
        double y,
        double z) const
    {
        return hessian(x, y, z).determinant();
    }

    template <typename T>
    T evalTemplate(T x, T y, T z) const
    {
        return x * x + y * y - z * z - c_;
    }
};
