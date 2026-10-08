#include "ImplicitSurface.hpp"
#include "LocalMarchingCubes.hpp"
#include "LocalUnfolder.hpp"
#include "LocalVoxelStitcher.hpp"
#include "MarchingCubes.hpp"
#include "SingularityDetector.hpp"
#include "app/MyViewer.h"

#include <pmp/io/io.h>
#include <pmp/surface_mesh.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

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
        vertices.push_back(
            mesh.add_vertex(
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

bool has_flag(int argc, char** argv, const std::string& flag)
{
    for (int i = 1; i < argc; ++i)
    {
        if (std::string(argv[i]) == flag)
            return true;
    }

    return false;
}

std::string option_value(
    int argc,
    char** argv,
    const std::string& option,
    const std::string& fallback)
{
    for (int i = 1; i + 1 < argc; ++i)
    {
        if (std::string(argv[i]) == option)
            return argv[i + 1];
    }

    return fallback;
}

std::string singularity_type_name(SingularityType type)
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

struct SingularityCandidate
{
    bool found = false;
    iso::mc::Point center{0.0, 0.0, 0.0};
    double gradient_norm =
        std::numeric_limits<double>::infinity();
    double function_value = 0.0;
    SingularityType type = SingularityType::Regular;
};

SingularityCandidate find_surface_singularity(
    const iso::mc::Mesh& mesh,
    const SingularityDetector& detector)
{
    SingularityCandidate candidate;

    for (const auto& p : mesh.vertices)
    {
        const pmp::Point point(
            static_cast<float>(p[0]),
            static_cast<float>(p[1]),
            static_cast<float>(p[2]));

        const auto result = detector.classifyPoint(point);

        if (result.type == SingularityType::Regular)
            continue;

        if (std::abs(result.function_value) > 1e-5)
            continue;

        if (!candidate.found ||
            result.gradient_norm < candidate.gradient_norm)
        {
            candidate.found = true;
            candidate.center = p;
            candidate.gradient_norm = result.gradient_norm;
            candidate.function_value = result.function_value;
            candidate.type = result.type;
        }
    }

    return candidate;
}

void print_topology(
    const std::string& label,
    const iso::stitch::TopologyStatistics& topology)
{
    std::cout
        << label << '\n'
        << "  vertices           : " << topology.vertices << '\n'
        << "  triangles          : " << topology.triangles << '\n'
        << "  edges              : " << topology.edges << '\n'
        << "  boundary edges     : " << topology.boundary_edges << '\n'
        << "  boundary vertices  : " << topology.boundary_vertices << '\n'
        << "  non-manifold edges : " << topology.nonmanifold_edges << '\n'
        << "  components         : "
        << topology.connected_components << '\n';
}

std::filesystem::path output_path(
    const std::filesystem::path& dir,
    const std::string& prefix,
    const std::string& surface,
    double c,
    std::size_t resolution,
    std::size_t level)
{
    std::ostringstream c_string;
    c_string << std::fixed << std::setprecision(6) << c;

    return dir /
        (prefix + "_" +
         surface + "_c_" +
         c_string.str() +
         "_r_" + std::to_string(resolution) +
         "_l_" + std::to_string(level) +
         ".obj");
}

void save_mesh(
    const iso::mc::Mesh& mesh,
    const std::filesystem::path& filename)
{
    const auto pmp_mesh = to_pmp_mesh(mesh);
    pmp::write(pmp_mesh, filename.string());

    std::cout << "Saved: " << filename << '\n';
}

void write_report(
    const std::filesystem::path& filename,
    const std::string& surface_name,
    double c,
    std::size_t resolution,
    double extent,
    double radius,
    std::size_t level,
    const SingularityCandidate& singularity,
    const iso::stitch::TopologyStatistics& global_topology,
    const LocalMarchingCubes::Result& patch,
    const iso::stitch::TopologyStatistics& patch_topology,
    const iso::stitch::StitchResult* stitched,
    const std::chrono::duration<double, std::milli>& global_duration)
{
    std::ofstream report(filename);
    if (!report)
        return;

    report << "Implicit isosurface experiment\n";
    report << "==============================\n\n";

    report << "Surface: " << surface_name << '\n';
    report << "c: " << std::setprecision(17) << c << '\n';
    report << "Global resolution: " << resolution << "^3\n";
    report << "Extent: " << extent << '\n';
    report << "Local radius: " << radius << '\n';
    report << "Local level: " << level << "\n\n";

    report << "Global MC\n";
    report << "---------\n";
    report << "Extraction time ms: "
           << global_duration.count() << '\n';
    report << "Vertices: " << global_topology.vertices << '\n';
    report << "Triangles: " << global_topology.triangles << '\n';
    report << "Edges: " << global_topology.edges << '\n';
    report << "Boundary edges: "
           << global_topology.boundary_edges << '\n';
    report << "Non-manifold edges: "
           << global_topology.nonmanifold_edges << '\n';
    report << "Components: "
           << global_topology.connected_components << "\n\n";

    report << "Singularity\n";
    report << "-----------\n";
    report << "Found: " << (singularity.found ? "yes" : "no") << '\n';
    report << "Type: "
           << singularity_type_name(singularity.type) << '\n';
    report << "Position: "
           << singularity.center[0] << ' '
           << singularity.center[1] << ' '
           << singularity.center[2] << '\n';
    report << "Gradient norm: "
           << singularity.gradient_norm << '\n';
    report << "Function value: "
           << singularity.function_value << "\n\n";

    report << "Local voxel MC\n";
    report << "--------------\n";
    report << "Resolution: "
           << patch.resolution_x << 'x'
           << patch.resolution_y << 'x'
           << patch.resolution_z << '\n';
    report << "Vertices: "
           << patch.mesh.vertices.size() << '\n';
    report << "Triangles: "
           << patch.mesh.triangles.size() << '\n';
    report << "Boundary edges: "
           << patch.boundary_edges << '\n';
    report << "Boundary vertices: "
           << patch.boundary_vertices << '\n';
    report << "Non-manifold edges: "
           << patch.nonmanifold_edges << '\n';
    report << "Extraction time ms: "
           << patch.extraction_time_ms << '\n';
    report << "Topology edges: "
           << patch_topology.edges << '\n';
    report << "Topology components: "
           << patch_topology.connected_components << "\n\n";

    if (stitched)
    {
        report << "Stitching\n";
        report << "---------\n";
        report << "Success: "
               << (stitched->stitched ? "yes" : "no") << '\n';
        report << "Global triangles removed: "
               << stitched->global_triangles_removed << '\n';
        report << "Local triangles inserted: "
               << stitched->local_triangles_inserted << '\n';
        report << "Seam loops: "
               << stitched->seam_loops << '\n';
        report << "Seam vertices: "
               << stitched->seam_vertices << '\n';
        report << "Seam triangles: "
               << stitched->seam_triangles << '\n';
        report << "Max seam vertex distance: "
               << stitched->max_seam_vertex_distance << '\n';
        report << "Final vertices: "
               << stitched->topology.vertices << '\n';
        report << "Final triangles: "
               << stitched->topology.triangles << '\n';
        report << "Final edges: "
               << stitched->topology.edges << '\n';
        report << "Final boundary edges: "
               << stitched->topology.boundary_edges << '\n';
        report << "Final non-manifold edges: "
               << stitched->topology.nonmanifold_edges << '\n';
        report << "Final components: "
               << stitched->topology.connected_components << '\n';
    }
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        if (has_flag(argc, argv, "--help"))
        {
            std::cout
                << "Usage:\n"
                << "  Isosurfaces [c] [resolution] [extent] [radius] [level]\n"
                << "              [--surface cone|quartic|sphere]\n"
                << "              [--output DIR]\n"
                << "              [--no-gui]\n"
                << "              [--no-stitch]\n";
            return 0;
        }

        const double c =
            argc > 1 ? std::stod(argv[1]) : 0.0;

        const std::size_t resolution =
            argc > 2
                ? static_cast<std::size_t>(std::stoul(argv[2]))
                : 65;

        const double extent =
            argc > 3 ? std::stod(argv[3]) : 2.0;

        const double radius =
            argc > 4 ? std::stod(argv[4]) : 0.30;

        const std::size_t level =
            argc > 5
                ? static_cast<std::size_t>(std::stoul(argv[5]))
                : 2;

        const std::string surface_name =
            option_value(
                argc,
                argv,
                "--surface",
                "cone");

        const std::filesystem::path output_dir =
            option_value(
                argc,
                argv,
                "--output",
                "outputs");

        const bool no_gui =
            has_flag(argc, argv, "--no-gui");

        const bool no_stitch =
            has_flag(argc, argv, "--no-stitch");

        if (resolution < 2)
            throw std::invalid_argument(
                "resolution must be >= 2");

        if (!(extent > 0.0))
            throw std::invalid_argument(
                "extent must be positive");

        if (!(radius > 0.0))
            throw std::invalid_argument(
                "radius must be positive");

        if (level == 0)
            throw std::invalid_argument(
                "level must be >= 1");

        if (surface_name != "cone" &&
            surface_name != "quartic" &&
            surface_name != "sphere")
        {
            throw std::invalid_argument(
                "Unknown surface. Use cone, quartic or sphere.");
        }

        std::filesystem::create_directories(output_dir);

        const auto surface =
            make_surface(surface_name, c);

        const iso::mc::Bounds bounds{
            {-extent, -extent, -extent},
            { extent,  extent, extent}};

        iso::mc::Options global_options;
        global_options.nx = resolution;
        global_options.ny = resolution;
        global_options.nz = resolution;
        global_options.isovalue = 0.0;

        const auto global_start =
            std::chrono::steady_clock::now();

        const iso::mc::Mesh global =
            iso::mc::extract(
                [&](double x, double y, double z)
                {
                    return surface->eval(x, y, z);
                },
                bounds,
                global_options);

        const auto global_end =
            std::chrono::steady_clock::now();

        const auto global_duration =
            std::chrono::duration<double, std::milli>(
                global_end - global_start);

        const auto global_topology =
            iso::stitch::compute_topology(global);

        std::cout
            << "========================================================\n"
            << "  IMPLICIT ISOSURFACE EXPERIMENT\n"
            << "========================================================\n"
            << "Surface             : " << surface_name << '\n'
            << "c                   : " << c << '\n'
            << "Global resolution   : " << resolution << "^3\n"
            << "Extent              : +/-" << extent << '\n'
            << "Local radius        : " << radius << '\n'
            << "Local level         : " << level << '\n'
            << "Output directory    : " << output_dir << '\n'
            << '\n';

        std::cout
            << "GLOBAL MARCHING CUBES\n"
            << "  extraction time    : "
            << global_duration.count() << " ms\n";

        print_topology("  topology:", global_topology);
        std::cout << '\n';

        const SingularityDetector detector(
            *surface,
            1e-2,
            1e-3);

        const auto singularity =
            find_surface_singularity(
                global,
                detector);

        if (!singularity.found)
        {
            std::cout
                << "No on-surface singularity was detected.\n";

            const auto global_filename =
                output_dir /
                (
                    "global_" +
                    surface_name +
                    "_c_" +
                    std::to_string(c) +
                    "_r_" +
                    std::to_string(resolution) +
                    ".obj");

            save_mesh(global, global_filename);

            if (no_gui)
                return 0;

            iso::MyViewer window(
                "Implicit Isosurface",
                1024,
                768);

            window.load_mesh(
                global_filename.string().c_str());

            window.set_marching_cubes(
                [&](double x, double y, double z)
                {
                    return surface->eval(x, y, z);
                },
                bounds,
                resolution,
                resolution,
                resolution,
                0.0);

            return window.run();
        }

        std::cout
            << "ON-SURFACE SINGULARITY\n"
            << "  position            : ("
            << singularity.center[0] << ", "
            << singularity.center[1] << ", "
            << singularity.center[2] << ")\n"
            << "  type                : "
            << singularity_type_name(
                   singularity.type)
            << '\n'
            << "  gradient norm       : "
            << singularity.gradient_norm << '\n'
            << "  function value      : "
            << singularity.function_value << "\n\n";

        /*
         * Local conforming baseline.
         *
         * We retain the full LevelStatistics reporting here, but do not
         * access a nonexistent LevelStatistics::mesh member. The refined
         * mesh storage remains owned by the LocalUnfolder result.
         */
        LocalUnfolder::Options local_options;
        local_options.radius = radius;
        local_options.levels = level;
        local_options.region_mode =
            LocalUnfolder::RegionMode::TopologicalBfs;
        local_options.sampling_bounds = bounds;
        local_options.base_resolution = resolution;
        local_options.isovalue = 0.0;

        LocalUnfolder unfolder;

        const auto baseline_start =
            std::chrono::steady_clock::now();

        const auto baseline =
            unfolder.refine(
                global,
                *surface,
                singularity.center,
                local_options);

        const auto baseline_end =
            std::chrono::steady_clock::now();

        std::cout
            << "LOCAL CONFORMING REFINEMENT BASELINE\n"
            << "  levels returned     : "
            << baseline.levels.size() << "\n\n";

        for (std::size_t i = 0;
             i < baseline.levels.size();
             ++i)
        {
            const auto& result =
                baseline.levels[i];

            std::cout
                << "  Level " << (i + 1) << '\n'
                << "    selected faces      : "
                << result.selected_faces << '\n'
                << "    sphere faces        : "
                << result.selected_faces_sphere << '\n'
                << "    topological faces   : "
                << result.selected_faces_topological << '\n'
                << "    selection changed   : "
                << (result.selection_changed ? "YES" : "NO")
                << '\n'
                << "    output vertices     : "
                << result.output_vertices << '\n'
                << "    output triangles    : "
                << result.output_triangles << '\n'
                << "    interface edges     : "
                << result.interface_edges << '\n'
                << "    interface vertices  : "
                << result.interface_vertices << '\n'
                << "    projection failures : "
                << result.projection_failures << '\n'
                << "    scalar refinement   : "
                << result.scalar.refinement_factor << '\n'
                << "    scalar voxels       : "
                << result.scalar.selected_voxels << '\n'
                << "    scalar samples      : "
                << result.scalar.corner_value_samples << '\n'
                << "    scalar spacing      : "
                << result.scalar.spacing << '\n'
                << "    scalar min          : "
                << result.scalar.min_value << '\n'
                << "    scalar max          : "
                << result.scalar.max_value << '\n'
                << "    scalar mean         : "
                << result.scalar.mean_value << '\n';
        }

        const auto baseline_ms =
            std::chrono::duration<double, std::milli>(
                baseline_end - baseline_start)
                .count();

        std::cout
            << "  total baseline time  : "
            << baseline_ms << " ms\n\n";

        /*
         * Independently extracted local voxel MC patch.
         */
        LocalMarchingCubes local_mc;

        LocalMarchingCubes::Options patch_options;
        patch_options.radius = radius;
        patch_options.level = level;
        patch_options.padding_cells = 1;
        patch_options.sampling_bounds = bounds;
        patch_options.base_resolution = resolution;
        patch_options.isovalue = 0.0;

        const auto patch =
            local_mc.extract(
                *surface,
                singularity.center,
                patch_options);

        const auto patch_topology =
            iso::stitch::compute_topology(
                patch.mesh);

        std::cout
            << "LOCAL VOXEL MARCHING CUBES\n"
            << "  resolution          : "
            << patch.resolution_x << 'x'
            << patch.resolution_y << 'x'
            << patch.resolution_z << '\n'
            << "  vertices             : "
            << patch.mesh.vertices.size() << '\n'
            << "  triangles            : "
            << patch.mesh.triangles.size() << '\n'
            << "  boundary edges       : "
            << patch.boundary_edges << '\n'
            << "  boundary vertices    : "
            << patch.boundary_vertices << '\n'
            << "  non-manifold edges   : "
            << patch.nonmanifold_edges << '\n'
            << "  extraction time      : "
            << patch.extraction_time_ms << " ms\n";

        print_topology(
            "  topology:",
            patch_topology);

        std::cout << '\n';

        const auto local_filename =
            output_path(
                output_dir,
                "local_voxel",
                surface_name,
                c,
                resolution,
                level);

        save_mesh(
            patch.mesh,
            local_filename);

        iso::stitch::StitchResult stitched;
        bool have_stitched_result = false;

        if (!no_stitch)
        {
            iso::stitch::LocalVoxelStitcher stitcher;

            iso::stitch::LocalVoxelStitcher::Options stitch_options;
            stitch_options.radius = radius;
            stitch_options.use_sphere_seam = true;

            const auto stitch_start =
                std::chrono::steady_clock::now();

            stitched =
                stitcher.stitch(
                    global,
                    patch.mesh,
                    singularity.center,
                    stitch_options);

            const auto stitch_end =
                std::chrono::steady_clock::now();

            const double stitch_ms =
                std::chrono::duration<double, std::milli>(
                    stitch_end - stitch_start)
                    .count();

            have_stitched_result = true;

            std::cout
                << "\nSTITCHED GLOBAL + LOCAL VOXEL MC\n"
                << "  stitched            : "
                << (stitched.stitched ? "YES" : "NO")
                << '\n'
                << "  stitch time         : "
                << stitch_ms << " ms\n"
                << "  global triangles removed: "
                << stitched.global_triangles_removed << '\n'
                << "  local triangles inserted: "
                << stitched.local_triangles_inserted << '\n'
                << "  seam loops           : "
                << stitched.seam_loops << '\n'
                << "  seam vertices        : "
                << stitched.seam_vertices << '\n'
                << "  seam triangles       : "
                << stitched.seam_triangles << '\n'
                << "  max seam distance    : "
                << stitched.max_seam_vertex_distance << '\n';

            print_topology(
                "  final topology:",
                stitched.topology);

            const auto stitched_filename =
                output_path(
                    output_dir,
                    "stitched_voxel",
                    surface_name,
                    c,
                    resolution,
                    level);

            save_mesh(
                stitched.mesh,
                stitched_filename);

            const auto report_filename =
                output_dir /
                (
                    "report_" +
                    surface_name +
                    "_c_" +
                    std::to_string(c) +
                    "_r_" +
                    std::to_string(resolution) +
                    "_l_" +
                    std::to_string(level) +
                    ".txt");

            write_report(
                report_filename,
                surface_name,
                c,
                resolution,
                extent,
                radius,
                level,
                singularity,
                global_topology,
                patch,
                patch_topology,
                &stitched,
                global_duration);

            std::cout
                << "Report saved: "
                << report_filename << '\n';
        }

        if (no_gui)
        {
            if (have_stitched_result)
                return stitched.stitched ? 0 : 2;

            return 0;
        }

        std::filesystem::path viewer_filename;

        if (have_stitched_result &&
            stitched.stitched)
        {
            viewer_filename =
                output_path(
                    output_dir,
                    "stitched_voxel",
                    surface_name,
                    c,
                    resolution,
                    level);
        }
        else
        {
            viewer_filename =
                output_dir /
                (
                    "global_" +
                    surface_name +
                    "_c_" +
                    std::to_string(c) +
                    "_r_" +
                    std::to_string(resolution) +
                    ".obj");

            save_mesh(
                global,
                viewer_filename);
        }

        iso::MyViewer window(
            have_stitched_result &&
                    stitched.stitched
                ? "Stitched Local Voxel MC"
                : "Implicit Isosurface",
            1024,
            768);

        window.load_mesh(
            viewer_filename.string().c_str());

        window.set_marching_cubes(
            [&](double x, double y, double z)
            {
                return surface->eval(x, y, z);
            },
            bounds,
            resolution,
            resolution,
            resolution,
            0.0);

        return window.run();
    }
    catch (const std::exception& e)
    {
        std::cerr
            << "Error: "
            << e.what()
            << '\n';

        return 1;
    }
}
