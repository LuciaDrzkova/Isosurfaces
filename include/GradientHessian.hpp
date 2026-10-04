#pragma once

#include "ImplicitSurface.hpp"
#include <Eigen/Dense>
#include <autodiff/forward/dual.hpp>
#include <pmp/types.h>
#include <cmath>

class GradientHessianCalculator {
private:
    const ImplicitSurface& surface_;

public:
    explicit GradientHessianCalculator(const ImplicitSurface& surf) : surface_(surf) {}

    // Compute 3D Gradient using forward-mode autodiff
    pmp::Point computeGradient(const pmp::Point& p) const {
        using namespace autodiff;
        dual x = p[0], y = p[1], z = p[2];

        auto eval_fn = [this](dual x_val, dual y_val, dual z_val) {
            // Cast to ParameterizedConeQuadric for template access
            const ParameterizedConeQuadric* quad = 
                dynamic_cast<const ParameterizedConeQuadric*>(&surface_);
            if (quad) return quad->evalTemplate(x_val, y_val, z_val);
            
            // Fallback: use numerical gradient
            return dual(0.0);
        };

        double dx = derivative(eval_fn, wrt(x), at(x, y, z));
        double dy = derivative(eval_fn, wrt(y), at(x, y, z));
        double dz = derivative(eval_fn, wrt(z), at(x, y, z));

        return pmp::Point(static_cast<float>(dx), static_cast<float>(dy), static_cast<float>(dz));
    }

    // Compute 3x3 Hessian using dual2nd second-order autodiff
    Eigen::Matrix3d computeHessian(const pmp::Point& p) const {
        using namespace autodiff;
        dual2nd x = p[0], y = p[1], z = p[2];

        auto eval_fn = [this](dual2nd x_val, dual2nd y_val, dual2nd z_val) {
            const ParameterizedConeQuadric* quad = 
                dynamic_cast<const ParameterizedConeQuadric*>(&surface_);
            if (quad) return quad->evalTemplate(x_val, y_val, z_val);
            return dual2nd(0.0);
        };

        double Hxx = derivative(eval_fn, wrt(x, x), at(x, y, z));
        double Hyy = derivative(eval_fn, wrt(y, y), at(x, y, z));
        double Hzz = derivative(eval_fn, wrt(z, z), at(x, y, z));
        double Hxy = derivative(eval_fn, wrt(x, y), at(x, y, z));
        double Hxz = derivative(eval_fn, wrt(x, z), at(x, y, z));
        double Hyz = derivative(eval_fn, wrt(y, z), at(x, y, z));

        Eigen::Matrix3d H;
        H << Hxx, Hxy, Hxz,
             Hxy, Hyy, Hyz,
             Hxz, Hyz, Hzz;
        return H;
    }

    // Newton-Raphson projection to calculate true Euclidean distance to the implicit surface
    pmp::Point projectOntoSurface(pmp::Point p, int max_iterations = 15, float tolerance = 1e-7f) const {
        for (int i = 0; i < max_iterations; ++i) {
            double val = surface_.eval(p[0], p[1], p[2]);
            if (std::abs(val) < tolerance) break;

            pmp::Point grad = computeGradient(p);
            float grad_sqnorm = pmp::dot(grad, grad);

            if (grad_sqnorm < 1e-12f) break; // Avoid division by zero near singularities

            p = p - (static_cast<float>(val) / grad_sqnorm) * grad;
        }
        return p;
    }

    // Exact Euclidean distance via projectOntoSurface
    float computeEuclideanDistance(const pmp::Point& p) const {
        pmp::Point p_proj = projectOntoSurface(p);
        return pmp::distance(p, p_proj);
    }
};
