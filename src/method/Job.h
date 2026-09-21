// One complete run, in this order: read the mesh, project it onto the surface,
// remesh it, run the experiment (optional), write the result.
// Both the window and the command line call run_job().

#pragma once

#include "method/MeshMetrics.h"
#include "method/Remesher.h"

#include <string>
#include <vector>

namespace iso {

//! Everything needed for one run: files, surface and remeshing parameters.
struct JobSettings
{
    std::string input;
    std::string output = "output.obj";
    std::string function; //!< name of a built-in surface, see builtin_surfaces()
    std::string experiment; //!< optional, name from experiments(); "" = none
    RemeshOptions remesh;
};

struct JobReport
{
    MeshMetrics input;     //!< mesh as read from disk
    MeshMetrics projected; //!< after projecting every vertex onto the surface
    MeshMetrics output;    //!< after remeshing (this is what is saved)
    std::vector<int> collapses_per_iteration;
    double seconds = 0; //!< remeshing time only
    std::string experiment_log; //!< text written by the experiment, if any
};

//! Reads settings.input, projects it onto the surface, remeshes it and writes
//! settings.output.
//! \throws std::invalid_argument for bad settings or a non-triangle mesh,
//!         pmp::IOException if a file cannot be read or written.
JobReport run_job(const JobSettings& settings);

} // namespace iso
