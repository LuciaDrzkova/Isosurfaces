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
#include "MarchingCubes.hpp"

// Store vertex colors directly to a PLY file bypassing PMP's property IO
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

// Helper to find mesh file in multiple locations
std::string find_mesh(const std::string& filename) {
    std::vector<std::string> search_paths = {
        filename,
        "assets/mesh_catalog/" + filename,
        "../assets/mesh_catalog/" + filename,
        "../../assets/mesh_catalog/" + filename,
    };
    
    for (const auto& path : search_paths) {
        std::ifstream test(path);
        if (test.good()) {
            std::cout << "Found mesh at: " << path << "\n";
            return path;
        }
    }
    return "";
}

int main(int argc, char** argv) {
    // 1. Parse command line or use default
    std::string mesh_filename = (argc > 1) ? argv[1] : "cone.obj";
    double c_val = (argc > 2) ? std::stod(argv[2]) : 0.0;
    
    std::cout << "========================================================\n";
    std::cout << "  PHASE 1: SINGULARITY DETECTION ON MESHLAB MESHES\n";
    std::cout << "========================================================\n\n";

    // 2. Find and load mesh
    std::string mesh_path = find_mesh(mesh_filename);
    if (mesh_path.empty()) {
        std::cerr << "ERROR: Could not find mesh '" << mesh_filename << "'\n";
        std::cerr << "Searched in:\n";
        std::cerr << "  - Current directory\n";
        std::cerr << "  - assets/mesh_catalog/\n";
        std::cerr << "\nUsage: ./Isosurfaces <mesh_file> [c_parameter]\n";
        std::cerr << "Example: ./Isosurfaces cone.obj 0.0\n";
        return 1;
    }

    // 3. Define the parameter matching your MeshLab export
    ParameterizedConeQuadric quadric(c_val);
    SingularityDetector detector(quadric);

    // 4. Load the MeshLab-generated mesh
    pmp::SurfaceMesh mesh;
    try {
        pmp::read(mesh, mesh_path);
        std::cout << "Successfully loaded mesh: " << mesh.n_vertices() << " vertices, " 
                  << mesh.n_faces() << " faces\n";
        std::cout << "Parameter c = " << c_val << "\n\n";
    } catch (const std::exception& e) {
        std::cerr << "Failed to load mesh: " << e.what() << "\n";
        return 1;
    }

    // 5. Classify, Color, and Localize the Singularity
    std::vector<pmp::Color> colors(mesh.n_vertices());
    int reg_count = 0, nondeg_count = 0, deg_count = 0;
    
    pmp::Point singularity_center(0.0f, 0.0f, 0.0f);

    std::cout << "Classifying vertices...\n";
    for (auto v : mesh.vertices()) {
        pmp::Point p = mesh.position(v);
        SingularityClassificationResult res = detector.classifyPoint(p);

        if (res.type == SingularityType::Regular) {
            colors[v.idx()] = pmp::Color(0.6f, 0.6f, 0.6f); // Gray
            reg_count++;
        } else if (res.type == SingularityType::NonDegenerateSingular) {
            colors[v.idx()] = pmp::Color(1.0f, 0.85f, 0.0f); // Yellow (Morse)
            singularity_center += p;
            nondeg_count++;
        } else {
            colors[v.idx()] = pmp::Color(1.0f, 0.0f, 0.0f); // Red (Degenerate)
            singularity_center += p;
            deg_count++;
        }
    }

    int total_singular = nondeg_count + deg_count;
    float local_radius = 0.3f; // The radius of our local unfolding bounding sphere
    int vertices_in_radius = 0;

    if (total_singular > 0) {
        // Calculate the exact centroid of the singularity
        singularity_center /= static_cast<float>(total_singular);
        
        std::cout << "\n--- Localization Results ---\n";
        std::cout << "Singularity Center Estimated at: (" 
                  << singularity_center[0] << ", " 
                  << singularity_center[1] << ", " 
                  << singularity_center[2] << ")\n";

        // Mark the bounding volume for local unfolding
        for (auto v : mesh.vertices()) {
            pmp::Point p = mesh.position(v);
            if (pmp::distance(p, singularity_center) <= local_radius) {
                vertices_in_radius++;
                // If it is a regular vertex, tint it light blue to visualize the unfolding zone
                if (colors[v.idx()][0] == 0.6f) { // Checking if it was gray
                    colors[v.idx()] = pmp::Color(0.2f, 0.6f, 1.0f); 
                }
            }
        }
        std::cout << "Local Unfolding Region (Radius " << local_radius << ") contains " 
                  << vertices_in_radius << " vertices.\n";
    }

    std::cout << "\nClassification Breakdown:\n";
    std::cout << "  🔴 Degenerate (RED):     " << deg_count << " vertices\n";
    std::cout << "  🟡 Morse Saddle (YELLOW): " << nondeg_count << " vertices\n";
    std::cout << "  🔵 Unfolding Zone (BLUE): " << (vertices_in_radius - total_singular) << " vertices\n";
    std::cout << "  ⚪ Regular (GRAY):       " << (reg_count - (vertices_in_radius - total_singular)) << " vertices\n";
    std::cout << "  Total:                   " << mesh.n_vertices() << " vertices\n\n";

    // 6. Export colored PLY for MeshLab visualization (bypassing PMP's IO)
    std::string out_ply = "colored_result_c_" + std::to_string(c_val) + ".ply";
    write_colored_ply(mesh, colors, out_ply);
    std::cout << "Saved colored PLY to: " << out_ply << " (Open this in MeshLab to see colors)\n";
    
    // 7. Export standard OBJ for the C++ viewer
    std::string out_obj = "viewer_mesh_c_" + std::to_string(c_val) + ".obj";
    pmp::write(mesh, out_obj);
    std::cout << "Saved viewer OBJ to: " << out_obj << "\n";
    
    std::cout << "\n>>> Opening Viewer. Close the window to exit.\n";
    iso::MyViewer window("Singularity Detector - Phase 1", 1024, 768);
    window.load_mesh(out_obj.c_str()); // Load the format PMP knows how to read
    return window.run();
}