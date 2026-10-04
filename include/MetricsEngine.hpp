#pragma once

#include "GradientHessian.hpp"
#include <pmp/surface_mesh.h>
#include <pmp/algorithms/differential_geometry.h>
#include <algorithm>
#include <vector>
#include <numeric>
#include <cmath>

struct MeshMetrics {
    int num_vertices = 0;
    int num_edges = 0;
    int num_faces = 0;
    
    // Geometry
    float min_area = 0.0f;
    float max_area = 0.0f;
    float mean_area = 0.0f;
    float median_area = 0.0f;
    int skinny_triangles_count = 0; // min interior angle < threshold
    
    // Accuracy
    float max_implicit_error = 0.0f;   // |f(p)|
    float mean_implicit_error = 0.0f;
    float max_euclidean_dist = 0.0f;   // ||p - p_proj||
    float mean_euclidean_dist = 0.0f;
    
    // Topology
    int num_connected_components = 0;
    int num_boundary_edges = 0;
    int euler_characteristic = 0;
};

class MetricsEngine {
public:
    static std::vector<float> computeTriangleAngles(const pmp::Point& p0, const pmp::Point& p1, const pmp::Point& p2) {
        float a = pmp::distance(p1, p2);
        float b = pmp::distance(p0, p2);
        float c = pmp::distance(p0, p1);

        float alpha = std::acos(std::clamp((b * b + c * c - a * a) / (2.0f * b * c), -1.0f, 1.0f)) * 180.0f / M_PI;
        float beta  = std::acos(std::clamp((a * a + c * c - b * b) / (2.0f * a * c), -1.0f, 1.0f)) * 180.0f / M_PI;
        float gamma = 180.0f - alpha - beta;

        return {alpha, beta, gamma};
    }

    static int countConnectedComponents(const pmp::SurfaceMesh& mesh) {
        std::vector<bool> visited(mesh.n_vertices(), false);
        int components = 0;

        for (auto v : mesh.vertices()) {
            if (visited[v.idx()]) continue;
            components++;
            
            std::vector<pmp::Vertex> queue = {v};
            visited[v.idx()] = true;
            size_t head = 0;

            while (head < queue.size()) {
                pmp::Vertex curr = queue[head++];
                for (auto neighbor : mesh.vertices(curr)) {
                    if (!visited[neighbor.idx()]) {
                        visited[neighbor.idx()] = true;
                        queue.push_back(neighbor);
                    }
                }
            }
        }
        return components;
    }

    static MeshMetrics evaluateMesh(const pmp::SurfaceMesh& mesh, const ImplicitSurface& surface, float angle_threshold = 25.0f) {
        MeshMetrics m;
        m.num_vertices = mesh.n_vertices();
        m.num_edges = mesh.n_edges();
        m.num_faces = mesh.n_faces();
        m.euler_characteristic = m.num_vertices - m.num_edges + m.num_faces;

        if (m.num_faces == 0) return m;

        // 1. Boundary Edges
        for (auto e : mesh.edges()) {
            if (mesh.is_boundary(e)) m.num_boundary_edges++;
        }

        // 2. Connected Components
        m.num_connected_components = countConnectedComponents(mesh);

        // 3. Triangle Statistics
        std::vector<float> areas;
        areas.reserve(mesh.n_faces());

        for (auto f : mesh.faces()) {
            float area = pmp::face_area(mesh, f);
            areas.push_back(area);

            std::vector<pmp::Point> pts;
            for (auto v : mesh.vertices(f)) pts.push_back(mesh.position(v));

            if (pts.size() == 3) {
                auto angles = computeTriangleAngles(pts[0], pts[1], pts[2]);
                if (*std::min_element(angles.begin(), angles.end()) < angle_threshold) {
                    m.skinny_triangles_count++;
                }
            }
        }

        std::sort(areas.begin(), areas.end());
        m.min_area = areas.front();
        m.max_area = areas.back();
        m.mean_area = std::accumulate(areas.begin(), areas.end(), 0.0f) / areas.size();
        m.median_area = areas[areas.size() / 2];

        // 4. Accuracy Metrics (|f(p)| vs True Distance)
        GradientHessianCalculator calc(surface);
        float sum_imp = 0.0f, sum_euc = 0.0f;

        for (auto v : mesh.vertices()) {
            pmp::Point p = mesh.position(v);
            float imp_err = std::abs(static_cast<float>(surface.eval(p[0], p[1], p[2])));
            float euc_err = calc.computeEuclideanDistance(p);

            m.max_implicit_error = std::max(m.max_implicit_error, imp_err);
            m.max_euclidean_dist = std::max(m.max_euclidean_dist, euc_err);
            sum_imp += imp_err;
            sum_euc += euc_err;
        }

        m.mean_implicit_error = sum_imp / m.num_vertices;
        m.mean_euclidean_dist = sum_euc / m.num_vertices;

        return m;
    }
};
