#include "ImplicitSurface.hpp"
#include "SingularityDetector.hpp"
#include "app/MyViewer.h"

#include <pmp/surface_mesh.h>
#include <pmp/io/io.h>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>

// Your custom PLY writer from bachelor's thesis[cite: 1]
void write_colored_ply(const pmp::SurfaceMesh& mesh,
                       const std::vector<pmp::Color>& vertex_colors,
                       const std::string& filename)
{
    auto points = mesh.get_vertex_property<pmp::Point>("v:point");
    std::ofstream out(filename);
    if (!out.is_open()) return;

    out << "ply\nformat ascii 1.0\n";
    out << "element vertex " << mesh.n_vertices() << "\n";
    out << "property float x\nproperty float y\nproperty float z\n";
    out << "property uchar red\nproperty uchar green\nproperty uchar blue\n";
    out << "element face " << mesh.n_faces() << "\n";
    out << "property list uchar int vertex_indices\nend_header\n";

    for (auto v : mesh.vertices()) {
        pmp::Point p = points[v];
        pmp::Color c = vertex_colors[v.idx()];
        int r = static_cast<int>(std::clamp(c[0], 0.0f, 1.0f) * 255.0f);
        int g = static_cast<int>(std::clamp(c[1], 0.0f, 1.0f) * 255.0f);
        int b = static_cast<int>(std::clamp(c[2], 0.0f, 1.0f) * 255.0f);
        out << p[0] << " " << p[1] << " " << p[2] << " " << r << " " << g << " " << b << "\n";
    }

    for (auto f : mesh.faces()) {
        std::vector<int> idx;
        for (auto v : mesh.vertices(f)) idx.push_back(v.idx());
        out << idx.size();
        for (int i : idx) out << " " << i;
        out << "\n";
    }
    out.close();
}

int main(int argc, char** argv) {
    // 1. Define the parameter matching your MeshLab export
    double c_val = 0.0; // Change this to match the mesh you loaded
    ParameterizedConeQuadric quadric(c_val);
    SingularityDetector detector(quadric, 1e-2, 1e-3); // threshold, hessian_tol

    // 2. Load the MeshLab-generated mesh
    std::string input_mesh = "cone.obj"; 
    pmp::SurfaceMesh mesh;
    try {
        pmp::read(mesh, input_mesh);
        std::cout << "Successfully loaded " << input_mesh << " (" << mesh.n_vertices() << " vertices)\n";
    } catch (const std::exception& e) {
        std::cerr << "Failed to load mesh: " << e.what() << "\n";
        return 1;
    }

    // 3. Classify and Color every vertex
    std::vector<pmp::Color> colors(mesh.n_vertices());
    int reg_count = 0, nondeg_count = 0, deg_count = 0;

    for (auto v : mesh.vertices()) {
        pmp::Point p = mesh.position(v);
        SingularityClassificationResult res = detector.classifyPoint(p);

        if (res.type == SingularityType::Regular) {
            colors[v.idx()] = pmp::Color(0.6f, 0.6f, 0.6f); // Gray
            reg_count++;
        } else if (res.type == SingularityType::NonDegenerateSingular) {
            colors[v.idx()] = pmp::Color(1.0f, 0.85f, 0.0f); // Yellow (Morse)
            nondeg_count++;
        } else {
            colors[v.idx()] = pmp::Color(1.0f, 0.0f, 0.0f); // Red (Degenerate)
            deg_count++;
        }
    }

    std::cout << "Classification Results:\n"
              << "  Regular (Gray): " << reg_count << "\n"
              << "  Morse/NonDegenerate (Yellow): " << nondeg_count << "\n"
              << "  Degenerate (Red): " << deg_count << "\n";

    // 4. Export colored PLY and visualize
    std::string out_ply = "colored_result.ply";
    write_colored_ply(mesh, colors, out_ply);
    
    std::cout << "Opening Viewer...\n";
    iso::MyViewer window("Singularity Viewer", 1024, 768);
    window.load_mesh(out_ply.c_str());
    return window.run();
}