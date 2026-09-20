// The list of research experiments. An experiment is a function that runs on the
// remeshed mesh, see Singularities.cpp for a starting point.

#pragma once

#include "method/ImplicitSurface.h"

#include <pmp/surface_mesh.h>

#include <functional>
#include <ostream>
#include <string>
#include <vector>

namespace iso {

//! A research experiment that runs on the remeshed mesh (see run_job()).
//! It may modify the mesh; whatever it writes to \p log is shown in the GUI
//! and printed by the command line tool.
//!
//! To add one: write a function with this signature in its own .cpp, add it to
//! src/Experiments.cpp and to isosurfaces_core in CMakeLists.txt.
using ExperimentFunction = std::function<void(
    pmp::SurfaceMesh& mesh, const ImplicitSurface& surface, std::ostream& log)>;

struct Experiment
{
    std::string name; //!< used on the command line (--experiment <name>)
    std::string description;
    ExperimentFunction run;
};

//! All registered experiments, in a stable order.
const std::vector<Experiment>& experiments();

//! Look up an experiment by name; nullptr if there is none.
const Experiment* find_experiment(const std::string& name);

} // namespace iso
