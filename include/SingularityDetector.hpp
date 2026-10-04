#pragma once

#include "GradientHessian.hpp"
#include <Eigen/Dense>

enum class SingularityType {
    Regular,
    NonDegenerateSingular, // Morse point: Hessian full-rank
    DegenerateSingular     // Higher-order: Hessian rank-deficient
};

struct SingularityClassificationResult {
    SingularityType type;
    Eigen::Vector3d eigenvalues;
    double max_grad_norm;
};

class SingularityDetector {
private:
    GradientHessianCalculator calc_;
    double grad_threshold_;
    double rel_hessian_tol_;

public:
    SingularityDetector(const ImplicitSurface& surf, double grad_thresh = 1e-3, double hess_tol = 1e-3)
        : calc_(surf), grad_threshold_(grad_thresh), rel_hessian_tol_(hess_tol) {}

    SingularityClassificationResult classifyPoint(const pmp::Point& p) const {
        pmp::Point grad = calc_.computeGradient(p);
        double grad_norm = std::sqrt(pmp::dot(grad, grad));

        SingularityClassificationResult res;
        res.max_grad_norm = grad_norm;

        if (grad_norm > grad_threshold_) {
            res.type = SingularityType::Regular;
            res.eigenvalues.setZero();
            return res;
        }

        // Candidate singular point: check Hessian rank via eigenvalues
        Eigen::Matrix3d H = calc_.computeHessian(p);
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(H);
        res.eigenvalues = solver.eigenvalues();

        double max_abs_eig = res.eigenvalues.cwiseAbs().maxCoeff();
        if (max_abs_eig <= 1e-12) {
            res.type = SingularityType::DegenerateSingular;
            return res;
        }

        bool rank_deficient = false;
        for (int i = 0; i < 3; ++i) {
            if (std::abs(res.eigenvalues[i]) < rel_hessian_tol_ * max_abs_eig) {
                rank_deficient = true;
                break;
            }
        }

        res.type = rank_deficient ? SingularityType::DegenerateSingular : SingularityType::NonDegenerateSingular;
        return res;
    }
};
