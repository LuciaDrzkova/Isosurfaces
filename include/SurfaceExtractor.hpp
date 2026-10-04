#pragma once

#include "ImplicitSurface.hpp"
#include <pmp/surface_mesh.h>
#include <vector>

// Baseline grid sampler and iso-surface generator
class VoxelGridExtractor {
private:
    double min_bound_;
    double max_bound_;
    int resolution_;

public:
    VoxelGridExtractor(double min_b = -2.0, double max_b = 2.0, int res = 32)
        : min_bound_(min_b), max_bound_(max_b), resolution_(res) {}

    pmp::SurfaceMesh extractMesh(const ImplicitSurface& surface) const {
        pmp::SurfaceMesh mesh;
        double step = (max_bound_ - min_bound_) / resolution_;

        // Uniform voxel sampling grid generating quad-patches across sign-crossings
        for (int i = 0; i < resolution_; ++i) {
            double x = min_bound_ + i * step;
            for (int j = 0; j < resolution_; ++j) {
                double y = min_bound_ + j * step;
                for (int k = 0; k < resolution_; ++k) {
                    double z = min_bound_ + k * step;

                    double v0 = surface.eval(x, y, z);
                    double v1 = surface.eval(x + step, y, z);
                    double v2 = surface.eval(x, y + step, z);
                    double v3 = surface.eval(x, y, z + step);

                    // Check for sign change along voxel edges
                    if ((v0 * v1 < 0) || (v0 * v2 < 0) || (v0 * v3 < 0)) {
                        pmp::Point p(static_cast<float>(x), static_cast<float>(y), static_cast<float>(z));
                        
                        // Construct a local quad cell approximating local surface
                        pmp::Vertex v_a = mesh.add_vertex(p);
                        pmp::Vertex v_b = mesh.add_vertex(p + pmp::Point(step, 0, 0));
                        pmp::Vertex v_c = mesh.add_vertex(p + pmp::Point(step, step, 0));
                        pmp::Vertex v_d = mesh.add_vertex(p + pmp::Point(0, step, 0));

                        try {
                            mesh.add_quad(v_a, v_b, v_c, v_d);
                        } catch (...) {
                            // Skip invalid quads
                        }
                    }
                }
            }
        }
        return mesh;
    }
};
