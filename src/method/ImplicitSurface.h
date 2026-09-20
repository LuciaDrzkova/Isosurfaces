// One implicit surface f(x, y, z) = 0: its value, its gradient (computed exactly
// by automatic differentiation) and the projection of a point onto it.
// The formulas themselves are in Surfaces.cpp.

#pragma once

#include <pmp/mat_vec.h>
#include <pmp/types.h>

#include <autodiff/forward/dual.hpp>

#include <string>
#include <utility>
#include <vector>

namespace iso {

//! A named implicit surface { p : f(p) = 0 }, differentiated with autodiff.
class ImplicitSurface
{
public:
    using Function = autodiff::dual (*)(autodiff::dual, autodiff::dual,
                                        autodiff::dual);

    //! Value of f and its gradient at one point.
    struct Sample
    {
        double value;
        pmp::Point gradient;
    };

    ImplicitSurface(std::string name, std::string description, Function f)
        : name_(std::move(name)),
          description_(std::move(description)),
          f_(f)
    {
    }

    const std::string& name() const { return name_; }
    const std::string& description() const { return description_; }

    //! f(p)
    double value(const pmp::Point& p) const;

    //! f(p) and grad f(p), from a single forward-mode pass per axis.
    Sample evaluate(const pmp::Point& p) const;

    //! Newton-Raphson projection of \p p onto the surface along the gradient.
    //! Stops when |f| < \p tolerance, when the gradient vanishes, or after
    //! \p max_iterations steps.
    pmp::Point project(pmp::Point p, int max_iterations = 10,
                       double tolerance = 1e-6) const;

private:
    std::string name_;
    std::string description_;
    Function f_;
};

//! All built-in surfaces, in a stable order.
const std::vector<ImplicitSurface>& builtin_surfaces();

//! Look up a built-in surface by name; nullptr if there is none.
const ImplicitSurface* find_surface(const std::string& name);

} // namespace iso
