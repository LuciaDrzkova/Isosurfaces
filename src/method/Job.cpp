// One complete run (see Job.h). This is the best place to see how the pieces
// fit together.

#include "method/Job.h"

#include "method/Experiments.h"

#include <pmp/io/io.h>

#include <chrono>
#include <sstream>
#include <stdexcept>

using namespace std;

namespace iso {

JobReport run_job(const JobSettings& settings)
{
    const ImplicitSurface* surface = find_surface(settings.function);
    if (!surface)
        throw invalid_argument("unknown function '" + settings.function + "'");

    const Experiment* experiment = nullptr;
    if (!settings.experiment.empty())
    {
        experiment = find_experiment(settings.experiment);
        if (!experiment)
            throw invalid_argument("unknown experiment '" +
                                   settings.experiment + "'");
    }

    pmp::SurfaceMesh mesh;
    pmp::read(mesh, settings.input);
    if (!mesh.is_triangle_mesh())
        throw invalid_argument(settings.input + " is not a triangle mesh");

    const float angle = settings.remesh.min_angle;
    JobReport report;
    report.input = compute_mesh_metrics(mesh, *surface, angle);

    project_onto_surface(mesh, *surface);
    report.projected = compute_mesh_metrics(mesh, *surface, angle);

    const auto start = chrono::steady_clock::now();
    report.collapses_per_iteration =
        remesh(mesh, *surface, settings.remesh).collapses_per_iteration;
    report.seconds =
        chrono::duration<double>(chrono::steady_clock::now() - start).count();

    if (experiment)
    {
        ostringstream log;
        experiment->run(mesh, *surface, log);
        report.experiment_log = log.str();
    }

    report.output = compute_mesh_metrics(mesh, *surface, angle);
    pmp::write(mesh, settings.output);
    return report;
}

} // namespace iso
