// Register new experiments in the list below.

#include "method/Experiments.h"

#include "method/Singularities.h"

namespace iso {

// Register new experiments here.
const std::vector<Experiment>& experiments()
{
    static const std::vector<Experiment> all = {
        {"singularities",
         "Starter for the singularities research (src/Singularities.cpp)",
         analyze_singularities},
    };
    return all;
}

const Experiment* find_experiment(const std::string& name)
{
    for (const auto& e : experiments())
        if (e.name == name)
            return &e;
    return nullptr;
}

} // namespace iso
