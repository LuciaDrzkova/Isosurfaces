#include "ImplicitSurface.hpp"
#include "SingularityDetector.hpp"
#include "MarchingCubes.hpp"
#include "app/MyViewer.h"

#include <pmp/io/io.h>
#include <pmp/surface_mesh.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Convert singularity type to readable text.
// -----------------------------------------------------------------------------

const char* singularity_type_name(
    SingularityType type)
{
    switch (type)
    {
    case SingularityType::Regular:
        return "Regular";

    case SingularityType::NonDegenerateSingular:
        return "NonDegenerateSingular";

    case SingularityType::DegenerateSingular:
        return "DegenerateSingular";
    }

    return "Unknown";
}

// -----------------------------------------------------------------------------
// Write colored PLY.
// -----------------------------------------------------------------------------

void write_colored_ply(
    const pmp::SurfaceMesh& mesh,
    const std::vector<pmp::Color>& vertex_colors,
    const std::string& filename)
{
    auto points =
        mesh.get_vertex_property<pmp::Point>("v:point");

    std::ofstream out(filename);

    if (!out.is_open())
    {
        throw std::runtime_error(
            "Could not open PLY output file: " + filename);
    }

    out << "ply\n";
    out << "format ascii 1.0\n";

    out << "element vertex "
        << mesh.n_vertices()
        << "\n";

    out << "property float x\n";
    out << "property float y\n";
    out << "property float z\n";

    out << "property uchar red\n";
    out << "property uchar green\n";
    out << "property uchar blue\n";

    out << "element face "
        << mesh.n_faces()
        << "\n";

    out << "property list uchar int vertex_indices\n";
    out << "end_header\n";

    for (auto v : mesh.vertices())
    {
        const pmp::Point p =
            points[v];

        const pmp::Color c =
            vertex_colors[v.idx()];

        const int r =
            static_cast<int>(
                std::clamp(
                    c[0],
                    0.0f,
                    1.0f) * 255.0f);

        const int g =
            static_cast<int>(
                std::clamp(
                    c[1],
                    0.0f,
                    1.0f) * 255.0f);

        const int b =
            static_cast<int>(
                std::clamp(
                    c[2],
                    0.0f,
                    1.0f) * 255.0f);

        out << p[0] << ' '
            << p[1] << ' '
            << p[2] << ' '
            << r << ' '
            << g << ' '
            << b << '\n';
    }

    for (auto f : mesh.faces())
    {
        std::vector<int> indices;

        indices.reserve(3);

        for (auto v : mesh.vertices(f))
            indices.push_back(v.idx());

        out << indices.size();

        for (int index : indices)
            out << ' ' << index;

        out << '\n';
    }
}

// -----------------------------------------------------------------------------
// Format c for filenames.
// -----------------------------------------------------------------------------

std::string format_parameter(
    double value)
{
    std::ostringstream stream;

    stream
        << std::fixed
        << std::setprecision(6)
        << value;

    return stream.str();
}

// -----------------------------------------------------------------------------
// Main
// -----------------------------------------------------------------------------

int main(
    int argc,
    char** argv)
{
    try
    {
        // ---------------------------------------------------------------------
        // Usage:
        //
        // ./build/Isosurfaces [c] [resolution] [extent]
        //
        // Example:
        //
        // ./build/Isosurfaces 0 129 2
        // ---------------------------------------------------------------------

        const double c_val =
            (argc > 1)
                ? std::stod(argv[1])
                : 0.0;

        const std::size_t resolution =
            (argc > 2)
                ? static_cast<std::size_t>(
                      std::stoul(argv[2]))
                : 129;

        const double extent =
            (argc > 3)
                ? std::stod(argv[3])
                : 2.0;

        if (resolution < 2)
        {
            std::cerr
                << "ERROR: resolution must be at least 2.\n";

            return 1;
        }

        if (extent <= 0.0)
        {
            std::cerr
                << "ERROR: extent must be positive.\n";

            return 1;
        }

        // ---------------------------------------------------------------------
        // 1. Analytic implicit surface
        // ---------------------------------------------------------------------

        ParameterizedConeQuadric quadric(c_val);

        SingularityDetector detector(
            quadric,
            1e-2,   // gradient threshold
            1e-3);  // Hessian determinant threshold

        // ---------------------------------------------------------------------
        // 2. Marching Cubes domain
        // ---------------------------------------------------------------------

        iso::mc::Bounds bounds{
            {-extent, -extent, -extent},
            { extent,  extent,  extent}
        };

        iso::mc::Options mc_options;

        mc_options.nx =
            resolution;

        mc_options.ny =
            resolution;

        mc_options.nz =
            resolution;

        mc_options.isovalue =
            0.0;

        std::cout
            << "========================================================\n"
            << "  MARCHING CUBES ISOSURFACE GENERATION\n"
            << "========================================================\n\n";

        std::cout
            << "Parameter c      : "
            << c_val
            << '\n';

        std::cout
            << "Grid resolution  : "
            << resolution
            << "^3\n";

        std::cout
            << "Sampling bounds  : [-"
            << extent
            << ", "
            << extent
            << "]^3\n\n";

        // ---------------------------------------------------------------------
        // 3. Marching Cubes
        // ---------------------------------------------------------------------

        std::cout
            << "Sampling implicit surface on "
            << resolution
            << "^3 grid...\n";

        const iso::mc::Mesh extracted =
            iso::mc::extract(
                [&](double x,
                    double y,
                    double z)
                {
                    return quadric.eval(
                        x,
                        y,
                        z);
                },
                bounds,
                mc_options);

        std::cout
            << "Marching Cubes produced "
            << extracted.vertices.size()
            << " vertices and "
            << extracted.triangles.size()
            << " triangles.\n";

        // ---------------------------------------------------------------------
        // 4. Convert to PMP SurfaceMesh
        // ---------------------------------------------------------------------

        pmp::SurfaceMesh mesh;

        std::vector<pmp::Vertex> vertices;

        vertices.reserve(
            extracted.vertices.size());

        for (const auto& p :
             extracted.vertices)
        {
            vertices.push_back(
                mesh.add_vertex(
                    pmp::Point(
                        static_cast<float>(p[0]),
                        static_cast<float>(p[1]),
                        static_cast<float>(p[2]))));
        }

        for (const auto& triangle :
             extracted.triangles)
        {
            mesh.add_triangle(
                vertices[triangle[0]],
                vertices[triangle[1]],
                vertices[triangle[2]]);
        }

        // ---------------------------------------------------------------------
        // 5. Classify every mesh vertex
        // ---------------------------------------------------------------------

        for (auto v : mesh.vertices())
        {
            const pmp::Point p = mesh.position(v);

            if (p == pmp::Point(0, 0, 0))
            {
                std::cout << "\nAPEX VERTEX " << v.idx() << "\n";

                std::cout
                    << "Valence: "
                    << mesh.valence(v)
                    << "\n";

                std::cout << "Adjacent vertices:\n";

                for (auto vv : mesh.vertices(v))
                {
                    const pmp::Point q =
                        mesh.position(vv);

                    std::cout
                        << "  "
                        << vv.idx()
                        << ": ("
                        << q[0] << ", "
                        << q[1] << ", "
                        << q[2] << ")\n";
                }
            }
        }

        std::vector<pmp::Color> colors(
            mesh.n_vertices());

        std::vector<SingularityType> classifications(
            mesh.n_vertices(),
            SingularityType::Regular);

        int regular_count = 0;
        int nondegenerate_count = 0;
        int degenerate_count = 0;

        pmp::Point singularity_center(
            0.0f,
            0.0f,
            0.0f);

        std::cout
            << "\n========================================================\n"
            << "  SINGULARITY DIAGNOSTICS\n"
            << "========================================================\n";

        for (auto v :
             mesh.vertices())
        {
            const pmp::Point p =
                mesh.position(v);

            const SingularityClassificationResult result =
                detector.classifyPoint(p);

            classifications[v.idx()] =
                result.type;

            // -----------------------------------------------------------------
            // Regular
            // -----------------------------------------------------------------

            if (result.type ==
                SingularityType::Regular)
            {
                colors[v.idx()] =
                    pmp::Color(
                        0.6f,
                        0.6f,
                        0.6f);

                ++regular_count;

                continue;
            }

            // -----------------------------------------------------------------
            // Singular candidate
            // -----------------------------------------------------------------

            std::cout
                << "\n--- Singular candidate ---\n";

            std::cout
                << std::setprecision(12);

            std::cout
                << "Vertex index: "
                << v.idx()
                << '\n';

            std::cout
                << "Position: ("
                << result.position[0]
                << ", "
                << result.position[1]
                << ", "
                << result.position[2]
                << ")\n";

            std::cout
                << "f(p): "
                << result.function_value
                << '\n';

            std::cout
                << "Gradient: ("
                << result.gradient[0]
                << ", "
                << result.gradient[1]
                << ", "
                << result.gradient[2]
                << ")\n";

            std::cout
                << "|Gradient|: "
                << result.gradient_norm
                << '\n';

            std::cout
                << "Hessian determinant: "
                << result.hessian_determinant
                << '\n';

            std::cout
                << "Hessian eigenvalues: ("
                << result.eigenvalues[0]
                << ", "
                << result.eigenvalues[1]
                << ", "
                << result.eigenvalues[2]
                << ")\n";

            std::cout
                << "Classification: "
                << singularity_type_name(
                       result.type)
                << '\n';

            // -----------------------------------------------------------------
            // Non-degenerate singular point = yellow
            // -----------------------------------------------------------------

            if (result.type ==
                SingularityType::NonDegenerateSingular)
            {
                colors[v.idx()] =
                    pmp::Color(
                        1.0f,
                        0.85f,
                        0.0f);

                singularity_center +=
                    p;

                ++nondegenerate_count;
            }
            // -----------------------------------------------------------------
            // Degenerate singular point = red
            // -----------------------------------------------------------------

            else
            {
                colors[v.idx()] =
                    pmp::Color(
                        1.0f,
                        0.0f,
                        0.0f);

                singularity_center +=
                    p;

                ++degenerate_count;
            }
        }

        // ---------------------------------------------------------------------
        // 6. Localize singular region
        // ---------------------------------------------------------------------

        const int total_singular =
            nondegenerate_count +
            degenerate_count;

        const float local_radius =
            0.3f;

        int vertices_in_radius =
            0;

        if (total_singular > 0)
        {
            singularity_center /=
                static_cast<float>(
                    total_singular);

            for (auto v :
                 mesh.vertices())
            {
                const pmp::Point p =
                    mesh.position(v);

                if (pmp::distance(
                        p,
                        singularity_center)
                    <= local_radius)
                {
                    ++vertices_in_radius;

                    // Only recolor regular vertices blue.
                    if (classifications[v.idx()] ==
                        SingularityType::Regular)
                    {
                        colors[v.idx()] =
                            pmp::Color(
                                0.2f,
                                0.6f,
                                1.0f);
                    }
                }
            }

            std::cout
                << "\n--- Localization Results ---\n";

            std::cout
                << "Singularity Center Estimated at: ("
                << singularity_center[0]
                << ", "
                << singularity_center[1]
                << ", "
                << singularity_center[2]
                << ")\n";

            std::cout
                << "Local Unfolding Region (Radius "
                << local_radius
                << ") contains "
                << vertices_in_radius
                << " vertices.\n";
        }

        // ---------------------------------------------------------------------
        // 7. Summary
        // ---------------------------------------------------------------------

        std::cout
            << "\nClassification Results:\n";

        std::cout
            << "  Regular (Gray): "
            << regular_count
            << '\n';

        std::cout
            << "  Morse/NonDegenerate (Yellow): "
            << nondegenerate_count
            << '\n';

        std::cout
            << "  Degenerate (Red): "
            << degenerate_count
            << '\n';

        // ---------------------------------------------------------------------
        // 8. Create outputs directory
        // ---------------------------------------------------------------------

        const std::filesystem::path output_dir =
            "outputs";

        std::filesystem::create_directories(
            output_dir);

        const std::string c_string =
            format_parameter(c_val);

        const std::filesystem::path out_ply =
            output_dir /
            ("colored_result_c_" +
             c_string +
             ".ply");

        const std::filesystem::path out_obj =
            output_dir /
            ("viewer_mesh_c_" +
             c_string +
             ".obj");

        // ---------------------------------------------------------------------
        // 9. Save colored PLY
        // ---------------------------------------------------------------------

        write_colored_ply(
            mesh,
            colors,
            out_ply.string());

        std::cout
            << "\nSaved colored PLY:\n  "
            << out_ply.string()
            << '\n';

        // ---------------------------------------------------------------------
        // 10. Save OBJ
        // ---------------------------------------------------------------------

        pmp::write(
            mesh,
            out_obj.string());

        std::cout
            << "Saved viewer OBJ:\n  "
            << out_obj.string()
            << '\n';

        // ---------------------------------------------------------------------
        // 11. Open viewer
        // ---------------------------------------------------------------------

        std::cout
            << "\nOpening Viewer...\n";

        iso::MyViewer window(
            "Marching Cubes / Singularity Viewer",
            1024,
            768);

        window.load_mesh(
            out_obj.string().c_str());

        return window.run();
    }
    catch (const std::exception& e)
    {
        std::cerr
            << "\nERROR: "
            << e.what()
            << '\n';

        return 1;
    }
}