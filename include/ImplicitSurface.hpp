#pragma once

#include <autodiff/forward/dual.hpp>

// Base interface for generic implicit surfaces
class ImplicitSurface {
public:
    virtual ~ImplicitSurface() = default;
    virtual double eval(double x, double y, double z) const = 0;
};

// Target family: f_c(x,y,z) = x^2 + y^2 - z^2 - c
class ParameterizedConeQuadric : public ImplicitSurface {
private:
    double c_;

public:
    explicit ParameterizedConeQuadric(double c = 0.0) : c_(c) {}

    void setParameter(double c) { c_ = c; }
    double getParameter() const { return c_; }

    double eval(double x, double y, double z) const override {
        return x * x + y * y - z * z - c_;
    }

    // Template methods for autodiff
    template <typename T>
    T evalTemplate(T x, T y, T z) const {
        return x * x + y * y - z * z - c_;
    }
};

