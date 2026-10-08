#include "ImplicitSurface.hpp"
#include "LocalMarchingCubes.hpp"
#include "LocalVoxelStitcher.hpp"
#include "MarchingCubes.hpp"
#include "SingularityDetector.hpp"
#include "app/MyViewer.h"

#include <pmp/io/io.h>
#include <pmp/surface_mesh.h>

#include <cmath>
#include <cstddef>
#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

std::unique_ptr<ImplicitSurface> make_surface(
    const std::string& name,
    double c)
{
    if (name == "cone")
        return std::make_unique<ParameterizedConeQuadric>(c);
    if (name == "quartic")
        return std::make_unique<QuarticSaddleSurface>(c, 0.05);
    if (name == "sphere")
        return std::make_unique<SphereImplicitSurface>(1.0);

    throw std::invalid_argument(
        "Unknown surface. Use cone, quartic or sphere.");
}

pmp::SurfaceMesh to_pmp_mesh(const iso::mc::Mesh& source)
{
    pmp::SurfaceMesh mesh;
    std::vector<pmp::Vertex> vertices;
    vertices.reserve(source.vertices.size());

    for (const auto& p : source.vertices)
    {
        vertices.push_back(mesh.add_vertex(
            pmp::Point(
                static_cast<float>(p[0]),
                static_cast<float>(p[1]),
                static_cast<float>(p[2]))));
    }

    for (const auto& triangle : source.triangles)
    {
        mesh.add_triangle(
            vertices[triangle[0]],
            vertices[triangle[1]],
            vertices[triangle[2]]);
    }

    return mesh;
}

bool no_gui(int argc, char** argv)
{
    for (int i = 1; i < argc; ++i)
    {
        if (std::string(argv[i]) == "--no-gui")
            return true;
    }
    return false;
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        const double c =
            argc > 1 ? std::stod(argv[1]) : 0.0;
        const std::size_t resolution =
            argc > 2 ? static_cast<std::size_t>(std::stoul(argv[2])) : 65;
        const double extent =
            argc > 3 ? std::stod(argv[3]) : 2.0;
        const double radius =
            argc > 4 ? std::stod(argv[4]) : 0.30;
        const std::size_t level =
            argc > 5 ? static_cast<std::size_t>(std::stoul(argv[5])) : 2;

        std::string surface_name = "cone";
        for (int i = 6; i + 1 < argc; ++i)
        {
            if (std::string(argv[i]) == "--surface")
                surface_name = argv[++i];
        }

        if (resolution < 2)
            throw std::invalid_argument("resolution must be >= 2");
        if (!(extent > 0.0))
            throw std::invalid_argument("extent must be positive");
        if (!(radius > 0.0))
            throw std::invalid_argument("radius must be positive");
        if (level == 0)
            throw std::invalid_argument("level must be >= 1");

        const auto surface = make_surface(surface_name, c);

        const iso::mc::Bounds bounds{
            {-extent, -extent, -extent},
            { extent,  extent,  extent}};

        iso::mc::Options global_options;
        global_options.nx = resolution;
        global_options.ny = resolution;
        global_options.nz = resolution;
        global_options.isovalue = 0.0;

        const auto global = iso::mc::extract(
            [&](double x, double y, double z)
            {
                return surface->eval(x, y, z);
            },
            bounds,
            global_options);

        const SingularityDetector detector(*surface, 1e-2, 1e-3);
        iso::mc::Point center{0.0, 0.0, 0.0};
        bool found = false;
        double best_gradient = 1e300;
        SingularityType best_type = SingularityType::Regular;

        for (const auto& p : global.vertices)
        {
            const pmp::Point point(
                static_cast<float>(p[0]),
                static_cast<float>(p[1]),
                static_cast<float>(p[2]));
            const auto classification = detector.classifyPoint(point);

            if (classification.type == SingularityType::Regular)
                continue;
            if (std::abs(classification.function_value) > 1e-5)
                continue;

            if (!found || classification.gradient_norm < best_gradient)
            {
                found = true;
                best_gradient = classification.gradient_norm;
                center = p;
                best_type = classification.type;
            }
        }

        if (!found)
        {
            std::cout
                << "No on-surface singularity was detected for surface "
                << surface_name << ".\n";
            return 0;
        }

        LocalMarchingCubes local_mc;
        LocalMarchingCubes::Options patch_options;
        patch_options.radius = radius;
        patch_options.level = level;
        patch_options.padding_cells = 1;
        patch_options.sampling_bounds = bounds;
        patch_options.base_resolution = resolution;
        patch_options.isovalue = 0.0;

        const auto patch = local_mc.extract(
            *surface,
            center,
            patch_options);

        iso::stitch::LocalVoxelStitcher stitcher;
        iso::stitch::LocalVoxelStitcher::Options stitch_options;
        stitch_options.radius = radius;
        stitch_options.use_sphere_seam = true;

        const auto stitched = stitcher.stitch(
            global,
            patch.mesh,
            center,
            stitch_options);

        std::cout
            << "========================================================\n"
            << "  LOCAL VOXEL PATCH -> GLOBAL STITCH\n"
            << "========================================================\n"
            << "Surface             : " << surface_name << '\n'
            << "c                   : " << c << '\n'
            << "Global resolution   : " << resolution << "^3\n"
            << "Local level         : " << level << '\n'
            << "Local patch         : "
            << patch.resolution_x << 'x'
            << patch.resolution_y << 'x'
            << patch.resolution_z << '\n'
            << "Radius              : " << radius << '\n'
            << "Singularity type    : "
            << (best_type == SingularityType::NonDegenerateSingular
                    ? "NonDegenerateSingular"
                    : "DegenerateSingular")
            << '\n'
            << "Patch vertices      : " << patch.mesh.vertices.size() << '\n'
            << "Patch triangles     : " << patch.mesh.triangles.size() << '\n'
            << "Stitched            : "
            << (stitched.stitched ? "YES" : "NO") << '\n'
            << "Seam loops          : " << stitched.seam_loops << '\n'
            << "Seam vertices       : " << stitched.seam_vertices << '\n'
            << "Seam triangles      : " << stitched.seam_triangles << '\n'
            << "Max seam distance   : "
            << stitched.max_seam_vertex_distance << '\n'
            << "Global triangles removed: "
            << stitched.global_triangles_removed << '\n'
            << "Local triangles inserted: "
            << stitched.local_triangles_inserted << '\n'
            << "Final vertices      : "
            << stitched.topology.vertices << '\n'
            << "Final triangles     : "
            << stitched.topology.triangles << '\n'
            << "Final edges         : "
            << stitched.topology.edges << '\n'
            << "Final boundary edges: "
            << stitched.topology.boundary_edges << '\n'
            << "Final nonmanifold edges: "
            << stitched.topology.nonmanifold_edges << '\n'
            << "Final components    : "
            << stitched.topology.connected_components << '\n';

        std::filesystem::create_directories("outputs");
        const std::string filename =
            "outputs/stitched_voxel_" + surface_name +
            "_c_" + std::to_string(c) +
            "_r_" + std::to_string(resolution) +
            "_l_" + std::to_string(level) + ".obj";

        const pmp::SurfaceMesh final_mesh =
            to_pmp_mesh(stitched.mesh);
        pmp::write(final_mesh, filename);

        std::cout << "Saved: " << filename << '\n';

        if (no_gui(argc, argv))
            return stitched.stitched ? 0 : 2;

        iso::MyViewer window(
            "Stitched Local Voxel MC",
            1024,
            768);

        window.load_mesh(filename.c_str());
        window.set_voxel_grid(
            pmp::Point(
                static_cast<float>(bounds.min[0]),
                static_cast<float>(bounds.min[1]),
                static_cast<float>(bounds.min[2])),
            pmp::Point(
                static_cast<float>(bounds.max[0]),
                static_cast<float>(bounds.max[1]),
                static_cast<float>(bounds.max[2])),
            static_cast<int>(resolution),
            static_cast<int>(resolution),
            static_cast<int>(resolution));

        return window.run();
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
}
