#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <stdexcept>

// Base interface for scalar implicit surfaces f(x,y,z)=0.
class ImplicitSurface
{
public:
    virtual ~ImplicitSurface() = default;
    virtual double eval(double x, double y, double z) const = 0;
};

// Original benchmark: f=x^2+y^2-z^2-c. The singularity is at the origin only
// when c=0; for c!=0 the level set is regular.
class ParameterizedConeQuadric : public ImplicitSurface
{
private:
    double c_ = 0.0;
public:
    explicit ParameterizedConeQuadric(double c = 0.0) : c_(c) {}
    void setParameter(double c) { c_ = c; }
    double getParameter() const { return c_; }
    double eval(double x, double y, double z) const override
    { return x*x + y*y - z*z - c_; }
    Eigen::Vector3d gradient(double x, double y, double z) const
    { return Eigen::Vector3d(2.0*x, 2.0*y, -2.0*z); }
    Eigen::Matrix3d hessian(double, double, double) const
    {
        Eigen::Matrix3d H;
        H << 2.0, 0.0, 0.0,
             0.0, 2.0, 0.0,
             0.0, 0.0, -2.0;
        return H;
    }
    double hessianDeterminant(double x, double y, double z) const
    { return hessian(x,y,z).determinant(); }
    template<class T> T evalTemplate(T x, T y, T z) const
    { return x*x + y*y - z*z - c_; }
};

// Quartic control family. For c=0 the origin is a non-degenerate Morse saddle.
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
        const double x2=x*x, y2=y*y, z2=z*z;
        return x2+y2-z2 + alpha_*(x2*x2+y2*y2+z2*z2)-c_;
    }
};

// Closed regular control surface: f=x^2+y^2+z^2-r^2.
class SphereImplicitSurface : public ImplicitSurface
{
private:
    double radius_ = 1.0;
public:
    explicit SphereImplicitSurface(double radius = 1.0) : radius_(radius)
    { if (!(radius_ > 0.0) || !std::isfinite(radius_)) radius_=1.0; }
    double eval(double x, double y, double z) const override
    { return x*x+y*y+z*z-radius_*radius_; }
    double radius() const { return radius_; }
};

// Isolated non-degenerate elliptic cone, z-axis aligned.
class EllipticConeSurface : public ImplicitSurface
{
private:
    double c_=0.0, a_=1.5, b_=0.75;
public:
    explicit EllipticConeSurface(double c=0.0, double a=1.5, double b=0.75)
        : c_(c), a_(a), b_(b) {}
    double getParameter() const { return c_; }
    double eval(double x,double y,double z) const override
    { return x*x/(a_*a_)+y*y/(b_*b_)-z*z-c_; }
};

// Same elliptic cone rotated in the XY plane (rotation does not change its apex).
class RotatedEllipticConeSurface : public ImplicitSurface
{
private:
    double c_=0.0, a_=1.5, b_=0.75, angle_=0.52359877559829887308;
public:
    explicit RotatedEllipticConeSurface(double c=0.0, double angle=0.52359877559829887308,
                                        double a=1.5, double b=0.75)
        : c_(c), a_(a), b_(b), angle_(angle) {}
    double eval(double x,double y,double z) const override
    {
        const double ca=std::cos(angle_), sa=std::sin(angle_);
        const double u=ca*x+sa*y, v=-sa*x+ca*y;
        return u*u/(a_*a_)+v*v/(b_*b_)-z*z-c_;
    }
};

// Elliptic cone with axis X: -x^2 + y^2/a^2 + z^2/b^2 - c = 0.
class XAxisEllipticConeSurface : public ImplicitSurface
{
private:
    double c_=0.0, a_=1.5, b_=0.75;
public:
    explicit XAxisEllipticConeSurface(double c=0.0,double a=1.5,double b=0.75)
        : c_(c), a_(a), b_(b) {}
    double eval(double x,double y,double z) const override
    { return -x*x+y*y/(a_*a_)+z*z/(b_*b_)-c_; }
};
using XConeSurface = XAxisEllipticConeSurface;
using XAxisConeSurface = XAxisEllipticConeSurface;

// Elliptic cone with axis Y: x^2/a^2 - y^2 + z^2/b^2 - c = 0.
class YAxisEllipticConeSurface : public ImplicitSurface
{
private:
    double c_=0.0, a_=1.5, b_=0.75;
public:
    explicit YAxisEllipticConeSurface(double c=0.0,double a=1.5,double b=0.75)
        : c_(c), a_(a), b_(b) {}
    double eval(double x,double y,double z) const override
    { return x*x/(a_*a_)-y*y+z*z/(b_*b_)-c_; }
};
using YConeSurface = YAxisEllipticConeSurface;
using YAxisConeSurface = YAxisEllipticConeSurface;

// Fully tilted elliptic cone. The Euler rotation maps world coordinates into a
// local orthonormal frame; the Hessian at the apex is full-rank and indefinite.
class TiltedEllipticConeSurface : public ImplicitSurface
{
private:
    double c_=0.0, a_=1.5, b_=0.75;
    double rx_=0.47, ry_=-0.61, rz_=0.38;
public:
    explicit TiltedEllipticConeSurface(double c=0.0,double rx=0.47,
                                       double ry=-0.61,double rz=0.38,
                                       double a=1.5,double b=0.75)
        : c_(c),a_(a),b_(b),rx_(rx),ry_(ry),rz_(rz) {}
    double eval(double x,double y,double z) const override
    {
        const double cx=std::cos(rx_), sx=std::sin(rx_);
        const double cy=std::cos(ry_), sy=std::sin(ry_);
        const double cz=std::cos(rz_), sz=std::sin(rz_);
        // R^T p for R=Rz*Ry*Rx.
        const double x1=cz*x+sz*y;
        const double y1=-sz*x+cz*y;
        const double u=cy*x1-sy*z;
        const double v=cx*y1-sx*(sy*x1+cy*z);
        const double w=sx*y1+cx*(sy*x1+cy*z);
        return u*u/(a_*a_)+v*v/(b_*b_)-w*w-c_;
    }
};
using TiltedConeSurface = TiltedEllipticConeSurface;

// A true degenerate case: f=x^2-y^2-c. At c=0, the zero set is the union
// x=y and x=-y; their intersection (the singular set) is the entire Z-axis.
class IntersectingPlanesSurface : public ImplicitSurface
{
private:
    double c_=0.0;
public:
    explicit IntersectingPlanesSurface(double c=0.0):c_(c){}
    double eval(double x,double y,double) const override { return x*x-y*y-c_; }
};
using CrossingPlanesSurface = IntersectingPlanesSurface;
