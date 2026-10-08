#include "ImplicitSurface.hpp"
#include "LocalMarchingCubes.hpp"
#include "LocalUnfolder.hpp"
#include "LocalVoxelStitcher.hpp"
#include "MarchingCubes.hpp"
#include "SingularityDetector.hpp"

#include <algorithm>
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
#include <vector>

namespace
{

struct SweepConfig
{
    double c_min = -0.30;
    double c_max = 0.30;
    double c_step = 0.01;
    std::vector<std::size_t> resolutions{65, 129};
    std::vector<double> radii{0.20, 0.30};
    std::vector<std::size_t> levels{1, 2};
    double extent = 2.0;
    std::string surface = "cone";
    std::filesystem::path output = "outputs/parameter_sweep.csv";
};

struct SingularityCandidate
{
    bool found = false;
    iso::mc::Point center{0.0, 0.0, 0.0};
    double gradient_norm = std::numeric_limits<double>::infinity();
    SingularityType type = SingularityType::Regular;
};

std::string csv_number(double value)
{
    if (!std::isfinite(value))
        return "";

    std::ostringstream out;
    out << std::setprecision(17) << value;
    return out.str();
}

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
        "Unknown surface. Use cone, quartic, sphere, or all.");
}

std::vector<std::string> selected_surfaces(const std::string& name)
{
    if (name == "all")
        return {"cone", "quartic", "sphere"};
    return {name};
}

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
            candidate.type = result.type;
        }
    }

    return candidate;
}

iso::stitch::TopologyStatistics topology(
    const iso::mc::Mesh& mesh)
{
    return iso::stitch::compute_topology(mesh);
}

void write_header(std::ofstream& csv)
{
    csv <<
        "surface,c,resolution,extent,radius,levels,level,"
        "global_vertices,global_triangles,global_edges,global_boundary_edges,"
        "global_boundary_vertices,global_nonmanifold_edges,global_components,"
        "singularity_found,singularity_type,selection_changed,"
        "sphere_faces,topological_faces,selected_faces,baseline_vertices,"
        "baseline_triangles,baseline_interface_edges,baseline_interface_vertices,"
        "projection_failures,scalar_refinement_factor,scalar_selected_voxels,"
        "scalar_corner_samples,scalar_spacing,scalar_min,scalar_max,scalar_mean,"
        "local_patch_resolution_x,local_patch_resolution_y,local_patch_resolution_z,"
        "local_patch_vertices,local_patch_triangles,local_patch_boundary_edges,"
        "local_patch_boundary_vertices,local_patch_nonmanifold_edges,local_patch_ms,"
        "stitched,stitched_vertices,stitched_triangles,stitched_edges,"
        "stitched_boundary_edges,stitched_boundary_vertices,"
        "stitched_nonmanifold_edges,stitched_components,global_triangles_removed,"
        "local_triangles_inserted,seam_loops,seam_vertices,seam_triangles,"
        "max_seam_vertex_distance\n";
}

void write_empty_row(
    std::ofstream& csv,
    const std::string& surface_name,
    double c,
    std::size_t resolution,
    double extent,
    double radius,
    std::size_t levels,
    std::size_t level,
    const iso::stitch::TopologyStatistics& global_topology,
    const SingularityCandidate& singularity)
{
    csv << surface_name << ','
        << csv_number(c) << ','
        << resolution << ','
        << csv_number(extent) << ','
        << csv_number(radius) << ','
        << levels << ','
        << level << ','
        << global_topology.vertices << ','
        << global_topology.triangles << ','
        << global_topology.edges << ','
        << global_topology.boundary_edges << ','
        << global_topology.boundary_vertices << ','
        << global_topology.nonmanifold_edges << ','
        << global_topology.connected_components << ','
        << (singularity.found ? 1 : 0) << ','
        << (singularity.found ? "singular" : "none") << ","
        << ",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\n";
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        SweepConfig config;

        for (int i = 1; i < argc; ++i)
        {
            const std::string arg = argv[i];

            auto require_value = [&](const char* option) -> std::string
            {
                if (i + 1 >= argc)
                    throw std::invalid_argument(
                        std::string("Missing value for ") + option);
                return argv[++i];
            };

            if (arg == "--surface")
                config.surface = require_value("--surface");
            else if (arg == "--c-min")
                config.c_min = std::stod(require_value("--c-min"));
            else if (arg == "--c-max")
                config.c_max = std::stod(require_value("--c-max"));
            else if (arg == "--c-step")
                config.c_step = std::stod(require_value("--c-step"));
            else if (arg == "--extent")
                config.extent = std::stod(require_value("--extent"));
            else if (arg == "--resolutions")
            {
                config.resolutions.clear();
                std::stringstream ss(require_value("--resolutions"));
                std::string item;
                while (std::getline(ss, item, ','))
                    config.resolutions.push_back(
                        static_cast<std::size_t>(std::stoul(item)));
            }
            else if (arg == "--radii")
            {
                config.radii.clear();
                std::stringstream ss(require_value("--radii"));
                std::string item;
                while (std::getline(ss, item, ','))
                    config.radii.push_back(std::stod(item));
            }
            else if (arg == "--levels")
            {
                config.levels.clear();
                std::stringstream ss(require_value("--levels"));
                std::string item;
                while (std::getline(ss, item, ','))
                    config.levels.push_back(
                        static_cast<std::size_t>(std::stoul(item)));
            }
            else if (arg == "--output")
                config.output = require_value("--output");
            else if (arg == "--quick")
            {
                config.c_min = -0.10;
                config.c_max = 0.10;
                config.c_step = 0.02;
                config.resolutions = {33, 65};
                config.radii = {0.20};
                config.levels = {1};
            }
            else if (arg == "--help")
            {
                std::cout
                    << "IsosurfacesSweep options:\n"
                    << "  --surface cone|quartic|sphere|all\n"
                    << "  --c-min VALUE --c-max VALUE --c-step VALUE\n"
                    << "  --resolutions N,N,...\n"
                    << "  --radii R,R,...\n"
                    << "  --levels L,L,...\n"
                    << "  --extent VALUE\n"
                    << "  --output FILE\n"
                    << "  --quick\n";
                return 0;
            }
            else
            {
                throw std::invalid_argument(
                    "Unknown argument: " + arg);
            }
        }

        if (!(config.c_step > 0.0))
            throw std::invalid_argument("c-step must be positive");
        if (!(config.extent > 0.0))
            throw std::invalid_argument("extent must be positive");
        if (config.resolutions.empty() ||
            config.radii.empty() ||
            config.levels.empty())
        {
            throw std::invalid_argument(
                "resolutions, radii and levels must not be empty");
        }

        for (std::size_t resolution : config.resolutions)
        {
            if (resolution < 2)
                throw std::invalid_argument("All resolutions must be >= 2");
        }
        for (double radius : config.radii)
        {
            if (!(radius > 0.0))
                throw std::invalid_argument("All radii must be positive");
        }

        const auto output_parent = config.output.parent_path();
        if (!output_parent.empty())
            std::filesystem::create_directories(output_parent);

        std::ofstream csv(config.output);
        if (!csv)
            throw std::runtime_error(
                "Could not open output CSV: " + config.output.string());

        csv << std::fixed;
        write_header(csv);

        std::cout
            << "Parameter sweep\n"
            << "  surfaces     : " << config.surface << '\n'
            << "  c range      : [" << config.c_min << ", "
            << config.c_max << "] step " << config.c_step << '\n'
            << "  resolutions  :";
        for (std::size_t n : config.resolutions)
            std::cout << ' ' << n;
        std::cout << "\n  radii        :";
        for (double r : config.radii)
            std::cout << ' ' << r;
        std::cout << "\n  levels       :";
        for (std::size_t l : config.levels)
            std::cout << ' ' << l;
        std::cout << "\n  output       : " << config.output << "\n\n";

        std::size_t run_count = 0;

        for (const std::string& surface_name : selected_surfaces(config.surface))
        {
            for (double c = config.c_min;
                 c <= config.c_max + 0.5 * config.c_step;
                 c += config.c_step)
            {
                const double current_c =
                    std::round(c / config.c_step) * config.c_step;

                for (std::size_t resolution : config.resolutions)
                {
                    const auto surface =
                        make_surface(surface_name, current_c);

                    const SingularityDetector detector(
                        *surface,
                        1e-2,
                        1e-3);

                    const iso::mc::Bounds bounds{
                        {-config.extent, -config.extent, -config.extent},
                        { config.extent,  config.extent,  config.extent}};

                    iso::mc::Options options;
                    options.nx = resolution;
                    options.ny = resolution;
                    options.nz = resolution;
                    options.isovalue = 0.0;

                    const iso::mc::Mesh global = iso::mc::extract(
                        [&](double x, double y, double z)
                        {
                            return surface->eval(x, y, z);
                        },
                        bounds,
                        options);

                    const auto global_topology = topology(global);
                    const SingularityCandidate singularity =
                        find_surface_singularity(global, detector);

                    std::cout
                        << '[' << surface_name << "] c="
                        << current_c
                        << " resolution=" << resolution
                        << " global V/T="
                        << global.vertices.size() << '/'
                        << global.triangles.size()
                        << " components="
                        << global_topology.connected_components;

                    if (!singularity.found)
                    {
                        std::cout << " no on-surface singularity\n";

                        write_empty_row(
                            csv,
                            surface_name,
                            current_c,
                            resolution,
                            config.extent,
                            0.0,
                            0,
                            0,
                            global_topology,
                            singularity);
                        ++run_count;
                        continue;
                    }

                    std::cout
                        << " singularity at ("
                        << singularity.center[0] << ','
                        << singularity.center[1] << ','
                        << singularity.center[2] << ")\n";

                    for (double radius : config.radii)
                    {
                        for (std::size_t requested_levels : config.levels)
                        {
                            LocalUnfolder::Options local_options;
                            local_options.radius = radius;
                            local_options.levels = requested_levels;
                            local_options.region_mode =
                                LocalUnfolder::RegionMode::TopologicalBfs;
                            local_options.sampling_bounds = bounds;
                            local_options.base_resolution = resolution;
                            local_options.isovalue = 0.0;

                            LocalUnfolder unfolder;
                            const auto baseline = unfolder.refine(
                                global,
                                *surface,
                                singularity.center,
                                local_options);

                            LocalMarchingCubes local_mc;
                            iso::stitch::LocalVoxelStitcher stitcher;

                            for (std::size_t level = 1;
                                 level <= requested_levels;
                                 ++level)
                            {
                                const auto baseline_level =
                                    baseline.levels.size() >= level
                                        ? &baseline.levels[level - 1]
                                        : nullptr;

                                LocalMarchingCubes::Options patch_options;
                                patch_options.radius = radius;
                                patch_options.level = level;
                                patch_options.padding_cells = 1;
                                patch_options.sampling_bounds = bounds;
                                patch_options.base_resolution = resolution;
                                patch_options.isovalue = 0.0;

                                const auto patch = local_mc.extract(
                                    *surface,
                                    singularity.center,
                                    patch_options);

                                iso::stitch::LocalVoxelStitcher::Options stitch_options;
                                stitch_options.radius = radius;
                                stitch_options.use_sphere_seam = true;

                                const auto stitched = stitcher.stitch(
                                    global,
                                    patch.mesh,
                                    singularity.center,
                                    stitch_options);

                                const std::size_t sphere_faces =
                                    baseline_level
                                        ? baseline_level->selected_faces_sphere
                                        : 0;
                                const std::size_t topo_faces =
                                    baseline_level
                                        ? baseline_level->selected_faces_topological
                                        : 0;
                                const std::size_t selected_faces =
                                    baseline_level
                                        ? baseline_level->selected_faces
                                        : 0;

                                csv
                                    << surface_name << ','
                                    << csv_number(current_c) << ','
                                    << resolution << ','
                                    << csv_number(config.extent) << ','
                                    << csv_number(radius) << ','
                                    << requested_levels << ','
                                    << level << ','
                                    << global_topology.vertices << ','
                                    << global_topology.triangles << ','
                                    << global_topology.edges << ','
                                    << global_topology.boundary_edges << ','
                                    << global_topology.boundary_vertices << ','
                                    << global_topology.nonmanifold_edges << ','
                                    << global_topology.connected_components << ','
                                    << 1 << ','
                                    << (singularity.type == SingularityType::NonDegenerateSingular
                                            ? "NonDegenerateSingular"
                                            : "DegenerateSingular") << ','
                                    << (baseline_level && baseline_level->selection_changed ? 1 : 0) << ','
                                    << sphere_faces << ','
                                    << topo_faces << ','
                                    << selected_faces << ','
                                    << (baseline_level ? baseline_level->output_vertices : 0) << ','
                                    << (baseline_level ? baseline_level->output_triangles : 0) << ','
                                    << (baseline_level ? baseline_level->interface_edges : 0) << ','
                                    << (baseline_level ? baseline_level->interface_vertices : 0) << ','
                                    << (baseline_level ? baseline_level->projection_failures : 0) << ','
                                    << (baseline_level ? baseline_level->scalar.refinement_factor : 0) << ','
                                    << (baseline_level ? baseline_level->scalar.selected_voxels : 0) << ','
                                    << (baseline_level ? baseline_level->scalar.corner_value_samples : 0) << ','
                                    << (baseline_level ? csv_number(baseline_level->scalar.spacing) : "") << ','
                                    << (baseline_level ? csv_number(baseline_level->scalar.min_value) : "") << ','
                                    << (baseline_level ? csv_number(baseline_level->scalar.max_value) : "") << ','
                                    << (baseline_level ? csv_number(baseline_level->scalar.mean_value) : "") << ','
                                    << patch.resolution_x << ','
                                    << patch.resolution_y << ','
                                    << patch.resolution_z << ','
                                    << patch.mesh.vertices.size() << ','
                                    << patch.mesh.triangles.size() << ','
                                    << patch.boundary_edges << ','
                                    << patch.boundary_vertices << ','
                                    << patch.nonmanifold_edges << ','
                                    << csv_number(patch.extraction_time_ms) << ','
                                    << (stitched.stitched ? 1 : 0) << ','
                                    << stitched.mesh.vertices.size() << ','
                                    << stitched.mesh.triangles.size() << ','
                                    << stitched.topology.edges << ','
                                    << stitched.topology.boundary_edges << ','
                                    << stitched.topology.boundary_vertices << ','
                                    << stitched.topology.nonmanifold_edges << ','
                                    << stitched.topology.connected_components << ','
                                    << stitched.global_triangles_removed << ','
                                    << stitched.local_triangles_inserted << ','
                                    << stitched.seam_loops << ','
                                    << stitched.seam_vertices << ','
                                    << stitched.seam_triangles << ','
                                    << csv_number(stitched.max_seam_vertex_distance)
                                    << '\n';

                                std::cout
                                    << "  radius=" << radius
                                    << " level=" << level
                                    << " selection_changed="
                                    << (baseline_level && baseline_level->selection_changed ? "YES" : "NO")
                                    << " scalar_mean="
                                    << (baseline_level ? baseline_level->scalar.mean_value : 0.0)
                                    << " stitched="
                                    << (stitched.stitched ? "YES" : "NO")
                                    << "\n";

                                ++run_count;
                            }
                        }
                    }
                }
            }
        }

        std::cout
            << "\nSweep rows written: " << run_count << '\n'
            << "CSV: " << config.output << '\n';

        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
}
