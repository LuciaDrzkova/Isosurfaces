#pragma once

#include "ImplicitSurface.hpp"

#include <Eigen/Dense>
#include <pmp/types.h>

#include <cmath>
#include <limits>

// -----------------------------------------------------------------------------
// Singularity classification
// -----------------------------------------------------------------------------

enum class SingularityType
{
    Regular,

    // Critical point whose Hessian is full rank:
    // det(H) != 0
    NonDegenerateSingular,

    // Critical point whose Hessian is rank deficient:
    // det(H) ~= 0
    DegenerateSingular
};

// -----------------------------------------------------------------------------
// Complete diagnostic information for one mesh vertex
// -----------------------------------------------------------------------------

struct SingularityClassificationResult
{
    SingularityType type =
        SingularityType::Regular;

    // Vertex position p
    pmp::Point position;

    // f(p)
    double function_value = 0.0;

    // grad f(p)
    Eigen::Vector3d gradient =
        Eigen::Vector3d::Zero();

    // ||grad f(p)||
    double gradient_norm = 0.0;

    // Hessian H(p)
    Eigen::Matrix3d hessian =
        Eigen::Matrix3d::Zero();

    // det(H(p))
    double hessian_determinant = 0.0;

    // Eigenvalues of H(p)
    Eigen::Vector3d eigenvalues =
        Eigen::Vector3d::Zero();
};

// -----------------------------------------------------------------------------
// Singularity detector
// -----------------------------------------------------------------------------

class SingularityDetector
{
private:
    const ImplicitSurface& surface_;

    // A vertex is considered a critical point when
    //
    //     ||grad f(p)|| <= grad_threshold_
    //
    double grad_threshold_;

    // A critical point is considered non-degenerate when
    //
    //     |det(H)| > hessian_threshold_
    //
    double hessian_threshold_;

public:
    SingularityDetector(
        const ImplicitSurface& surface,
        double grad_threshold = 1e-2,
        double hessian_threshold = 1e-3)
        : surface_(surface),
          grad_threshold_(grad_threshold),
          hessian_threshold_(hessian_threshold)
    {
    }

    // -------------------------------------------------------------------------
    // Analyze one mesh vertex.
    //
    // IMPORTANT:
    // We use the analytic Hessian of ParameterizedConeQuadric here instead of
    // GradientHessian.hpp's dual2nd Hessian path.
    // -------------------------------------------------------------------------

    SingularityClassificationResult classifyPoint(
        const pmp::Point& p) const
    {
        SingularityClassificationResult result;

        result.position = p;

        const double x =
            static_cast<double>(p[0]);

        const double y =
            static_cast<double>(p[1]);

        const double z =
            static_cast<double>(p[2]);

        // ---------------------------------------------------------------------
        // 1. Function value f(p)
        // ---------------------------------------------------------------------

        result.function_value =
            surface_.eval(x, y, z);

        // ---------------------------------------------------------------------
        // 2. Gradient grad f(p)
        // ---------------------------------------------------------------------

        if (const auto* quadric =
                dynamic_cast<const ParameterizedConeQuadric*>(
                    &surface_))
        {
            result.gradient =
                quadric->gradient(x, y, z);
        }
        else
        {
            // Generic fallback: numerical central difference.
            result.gradient =
                numericalGradient(x, y, z);
        }

        result.gradient_norm =
            result.gradient.norm();

        // ---------------------------------------------------------------------
        // Regular point:
        //
        // gradient is not approximately zero.
        // ---------------------------------------------------------------------

        if (result.gradient_norm >
            grad_threshold_)
        {
            result.type =
                SingularityType::Regular;

            return result;
        }

        // ---------------------------------------------------------------------
        // 3. Hessian H(p)
        // ---------------------------------------------------------------------

        if (const auto* quadric =
                dynamic_cast<const ParameterizedConeQuadric*>(
                    &surface_))
        {
            result.hessian =
                quadric->hessian(x, y, z);
        }
        else
        {
            // Generic fallback: numerical Hessian.
            result.hessian =
                numericalHessian(x, y, z);
        }

        // ---------------------------------------------------------------------
        // 4. Hessian determinant
        // ---------------------------------------------------------------------

        result.hessian_determinant =
            result.hessian.determinant();

        // ---------------------------------------------------------------------
        // 5. Hessian eigenvalues
        // ---------------------------------------------------------------------

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>
            solver(result.hessian);

        if (solver.info() == Eigen::Success)
        {
            result.eigenvalues =
                solver.eigenvalues();
        }

        // ---------------------------------------------------------------------
        // 6. Final singularity classification
        //
        // Critical point + det(H) != 0
        //     -> NonDegenerateSingular
        //
        // Critical point + det(H) ~= 0
        //     -> DegenerateSingular
        // ---------------------------------------------------------------------

        if (std::abs(result.hessian_determinant) >
            hessian_threshold_)
        {
            result.type =
                SingularityType::NonDegenerateSingular;
        }
        else
        {
            result.type =
                SingularityType::DegenerateSingular;
        }

        return result;
    }

private:

    // -------------------------------------------------------------------------
    // Numerical gradient fallback
    // Used only if another ImplicitSurface implementation is added later.
    // -------------------------------------------------------------------------

    Eigen::Vector3d numericalGradient(
        double x,
        double y,
        double z) const
    {
        const double h = 1e-6;

        const double fx1 =
            surface_.eval(x + h, y, z);

        const double fx0 =
            surface_.eval(x - h, y, z);

        const double fy1 =
            surface_.eval(x, y + h, z);

        const double fy0 =
            surface_.eval(x, y - h, z);

        const double fz1 =
            surface_.eval(x, y, z + h);

        const double fz0 =
            surface_.eval(x, y, z - h);

        return Eigen::Vector3d(
            (fx1 - fx0) / (2.0 * h),
            (fy1 - fy0) / (2.0 * h),
            (fz1 - fz0) / (2.0 * h));
    }

    // -------------------------------------------------------------------------
    // Numerical Hessian fallback
    // -------------------------------------------------------------------------

    Eigen::Matrix3d numericalHessian(
        double x,
        double y,
        double z) const
    {
        const double h = 1e-4;

        Eigen::Matrix3d H =
            Eigen::Matrix3d::Zero();

        const double f000 =
            surface_.eval(x, y, z);

        // Diagonal entries

        H(0, 0) =
            (surface_.eval(x + h, y, z)
             - 2.0 * f000
             + surface_.eval(x - h, y, z))
            / (h * h);

        H(1, 1) =
            (surface_.eval(x, y + h, z)
             - 2.0 * f000
             + surface_.eval(x, y - h, z))
            / (h * h);

        H(2, 2) =
            (surface_.eval(x, y, z + h)
             - 2.0 * f000
             + surface_.eval(x, y, z - h))
            / (h * h);

        // Mixed xy

        H(0, 1) =
            (surface_.eval(x + h, y + h, z)
             - surface_.eval(x + h, y - h, z)
             - surface_.eval(x - h, y + h, z)
             + surface_.eval(x - h, y - h, z))
            / (4.0 * h * h);

        // Mixed xz

        H(0, 2) =
            (surface_.eval(x + h, y, z + h)
             - surface_.eval(x + h, y, z - h)
             - surface_.eval(x - h, y, z + h)
             + surface_.eval(x - h, y, z - h))
            / (4.0 * h * h);

        // Mixed yz

        H(1, 2) =
            (surface_.eval(x, y + h, z + h)
             - surface_.eval(x, y + h, z - h)
             - surface_.eval(x, y - h, z + h)
             + surface_.eval(x, y - h, z - h))
            / (4.0 * h * h);

        // Symmetry

        H(1, 0) = H(0, 1);
        H(2, 0) = H(0, 2);
        H(2, 1) = H(1, 2);

        return H;
    }
};