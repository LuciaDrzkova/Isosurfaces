#pragma once

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <cmath>
#include <stdexcept>

// Base interface for scalar implicit surfaces f(x,y,z)=0.
class ImplicitSurface
{
public:
    virtual ~ImplicitSurface() = default;
    virtual double eval(double x, double y, double z) const = 0;
};

// Parameterized quadric used throughout the original experiments:
// f(x,y,z) = x^2 + y^2 - z^2 - c.
// c < 0: two-sheet hyperboloid; c = 0: double cone; c > 0: one-sheet hyperboloid.
class ParameterizedConeQuadric : public ImplicitSurface
{
private:
    double c_ = 0.0;
public:
    explicit ParameterizedConeQuadric(double c = 0.0) : c_(c) {}

    void setParameter(double c) { c_ = c; }
    double getParameter() const { return c_; }

    double eval(double x, double y, double z) const override
    {
        return x * x + y * y - z * z - c_;
    }

    Eigen::Vector3d gradient(double x, double y, double z) const
    {
        return Eigen::Vector3d(2.0 * x, 2.0 * y, -2.0 * z);
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

// Degree-2 elliptic cone with an isolated, non-degenerate singularity at
// the origin when c=0:
//     x^2/a^2 + y^2/b^2 - z^2 - c = 0.
// The unequal transverse scales make this a useful non-circular cone test.
class EllipticConeSurface final : public ImplicitSurface
{
private:
    double c_;
    double a_;
    double b_;
public:
    explicit EllipticConeSurface(double c = 0.0,
                                 double semi_axis_x = 1.5,
                                 double semi_axis_y = 0.75)
        : c_(c), a_(semi_axis_x), b_(semi_axis_y)
    {
        if (!(a_ > 0.0) || !(b_ > 0.0))
            throw std::invalid_argument("Elliptic cone semi-axes must be positive");
    }

    double eval(double x, double y, double z) const override
    {
        return (x * x) / (a_ * a_) + (y * y) / (b_ * b_) - z * z - c_;
    }
};

// Degree-2 elliptic cone rotated in the XY plane. It retains an isolated,
// non-degenerate singularity at the origin for c=0, but its principal axes
// are not aligned with the grid. This tests the algorithms on cross terms
// that arise when the equation is expanded in x and y.
class RotatedEllipticConeSurface final : public ImplicitSurface
{
private:
    double c_;
    double a_;
    double b_;
    double theta_;
public:
    explicit RotatedEllipticConeSurface(
        double c = 0.0,
        double semi_axis_x = 1.5,
        double semi_axis_y = 0.75,
        double angle_radians = 0.52359877559829887308)
        : c_(c), a_(semi_axis_x), b_(semi_axis_y), theta_(angle_radians)
    {
        if (!(a_ > 0.0) || !(b_ > 0.0))
            throw std::invalid_argument("Rotated cone semi-axes must be positive");
    }

    double eval(double x, double y, double z) const override
    {
        const double ct = std::cos(theta_);
        const double st = std::sin(theta_);
        const double xr = ct * x + st * y;
        const double yr = -st * x + ct * y;
        return (xr * xr) / (a_ * a_) +
               (yr * yr) / (b_ * b_) - z * z - c_;
    }
};

// Elliptic cone whose singular axis is the X axis:
//     y^2/a^2 + z^2/b^2 - x^2 - c = 0.
class XAxisEllipticConeSurface final : public ImplicitSurface
{
private:
    double c_;
    double a_;
    double b_;
public:
    explicit XAxisEllipticConeSurface(double c = 0.0,
                                      double semi_axis_y = 1.5,
                                      double semi_axis_z = 0.75)
        : c_(c), a_(semi_axis_y), b_(semi_axis_z)
    {
        if (!(a_ > 0.0) || !(b_ > 0.0))
            throw std::invalid_argument("X-axis cone semi-axes must be positive");
    }

    double eval(double x, double y, double z) const override
    {
        return (y * y) / (a_ * a_) + (z * z) / (b_ * b_) - x * x - c_;
    }
};

// Elliptic cone whose singular axis is the Y axis:
//     x^2/a^2 + z^2/b^2 - y^2 - c = 0.
class YAxisEllipticConeSurface final : public ImplicitSurface
{
private:
    double c_;
    double a_;
    double b_;
public:
    explicit YAxisEllipticConeSurface(double c = 0.0,
                                      double semi_axis_x = 1.5,
                                      double semi_axis_z = 0.75)
        : c_(c), a_(semi_axis_x), b_(semi_axis_z)
    {
        if (!(a_ > 0.0) || !(b_ > 0.0))
            throw std::invalid_argument("Y-axis cone semi-axes must be positive");
    }

    double eval(double x, double y, double z) const override
    {
        return (x * x) / (a_ * a_) + (z * z) / (b_ * b_) - y * y - c_;
    }
};

// Elliptic cone rotated about all three world axes. In local coordinates it
// has the usual isolated cone singularity at the origin; in world coordinates
// the equation has mixed terms and the singular axis is not grid-aligned.
class TiltedEllipticConeSurface final : public ImplicitSurface
{
private:
    double c_;
    double a_;
    double b_;
    Eigen::Matrix3d rotation_;
public:
    explicit TiltedEllipticConeSurface(double c = 0.0,
                                       double semi_axis_x = 1.5,
                                       double semi_axis_y = 0.75)
        : c_(c), a_(semi_axis_x), b_(semi_axis_y)
    {
        if (!(a_ > 0.0) || !(b_ > 0.0))
            throw std::invalid_argument("Tilted cone semi-axes must be positive");

        const Eigen::AngleAxisd yaw(0.47, Eigen::Vector3d::UnitZ());
        const Eigen::AngleAxisd pitch(0.61, Eigen::Vector3d::UnitY());
        const Eigen::AngleAxisd roll(0.28, Eigen::Vector3d::UnitX());
        rotation_ = (yaw * pitch * roll).toRotationMatrix();
    }

    double eval(double x, double y, double z) const override
    {
        const Eigen::Vector3d world(x, y, z);
        const Eigen::Vector3d local = rotation_.transpose() * world;
        return (local.x() * local.x()) / (a_ * a_) +
               (local.y() * local.y()) / (b_ * b_) -
               local.z() * local.z() - c_;
    }
};

// Pair of intersecting planes x=y and x=-y when c=0:
//     x^2 - y^2 - c = 0.
// Unlike the cone examples, its singular set is the entire z axis and the
// Hessian is rank deficient, so the origin is a degenerate singularity.
class IntersectingPlanesSurface final : public ImplicitSurface
{
private:
    double c_;
public:
    explicit IntersectingPlanesSurface(double c = 0.0) : c_(c) {}

    double eval(double x, double y, double) const override
    {
        return x * x - y * y - c_;
    }
};

// A genuinely non-quadric family used to test generic numerical derivatives:
// f = x^2 + y^2 - z^2 + alpha*(x^4+y^4+z^4) - c.
class QuarticSaddleSurface : public ImplicitSurface
{
private:
    double c_ = 0.0;
    double alpha_ = 0.05;
public:
    explicit QuarticSaddleSurface(double c = 0.0, double alpha = 0.05)
        : c_(c), alpha_(alpha) {}

    double getParameter() const { return c_; }
    double getAlpha() const { return alpha_; }

    double eval(double x, double y, double z) const override
    {
        const double x2 = x * x;
        const double y2 = y * y;
        const double z2 = z * z;
        return x2 + y2 - z2 + alpha_ * (x2 * x2 + y2 * y2 + z2 * z2) - c_;
    }
};

// Closed regular control surface. For positive radius its level set has no
// critical point on the surface, so the local-singularity stage should skip it.
class SphereImplicitSurface : public ImplicitSurface
{
private:
    double radius_ = 1.0;
public:
    explicit SphereImplicitSurface(double radius = 1.0) : radius_(radius)
    {
        if (!(radius_ > 0.0))
            radius_ = 1.0;
    }

    double eval(double x, double y, double z) const override
    {
        return x * x + y * y + z * z - radius_ * radius_;
    }

    double radius() const { return radius_; }
};
