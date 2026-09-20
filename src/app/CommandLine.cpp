#include "app/CommandLine.h"

#include "method/Experiments.h"
#include "method/ImplicitSurface.h"

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace std;

namespace iso {

namespace {

// "  name          " padded so the descriptions line up
string padded(const string& name)
{
    return "  " + name + string(name.size() < 16 ? 16 - name.size() : 1, ' ');
}

void print_functions(ostream& out)
{
    for (const auto& s : builtin_surfaces())
        out << padded(s.name()) << s.description() << '\n';
}

void print_experiments(ostream& out)
{
    for (const auto& e : experiments())
        out << padded(e.name) << e.description << '\n';
}

double parse_number(const string& option, const string& text)
{
    char* end = nullptr;
    const double value = strtod(text.c_str(), &end);
    if (text.empty() || *end != '\0')
        throw invalid_argument("invalid number for " + option + ": '" + text +
                               "'");
    return value;
}

} // namespace

void print_usage(ostream& out)
{
    out << "Usage: Isosurfaces [input-mesh] [options]\n"
           "\n"
           "Projects a triangle mesh onto an implicit surface and coarsens it\n"
           "by collapsing small and skinny triangles.\n"
           "\n"
           "Without arguments a window opens where the mesh, surface and\n"
           "parameters can be chosen. With --no-gui the mesh and --function are\n"
           "required and the result is only written to a file.\n"
           "\n"
           "Options:\n"
           "  -f, --function <name>   implicit surface to project onto\n"
           "  -o, --output <file>     output mesh (default: output.obj)\n"
           "  -n, --iterations <N>    number of remeshing passes (default: 6)\n"
           "  -a, --angle <degrees>   collapse triangles with an angle below this\n"
           "                          (default: 25)\n"
           "  -d, --divisor <x>       collapse triangles smaller than\n"
           "                          median area / x (default: 1.5)\n"
           "  -e, --experiment <name> run an experiment on the result\n"
           "      --no-gui            do not open a window\n"
           "      --list-functions    list the built-in implicit surfaces\n"
           "      --list-experiments  list the available experiments\n"
           "  -h, --help              show this help\n";
}

bool parse_command_line(int argc, char** argv, CommandLine& cl)
{
    JobSettings& job = cl.job;

    for (int i = 1; i < argc; ++i)
    {
        const string arg = argv[i];

        // The text after an option, e.g. "sphere" in "--function sphere"
        auto value = [&]() -> string {
            if (i + 1 >= argc)
                throw invalid_argument("missing value for " + arg);
            return argv[++i];
        };

        if (arg == "-h" || arg == "--help")
        {
            print_usage(cout);
            return false;
        }
        else if (arg == "--list-functions")
        {
            print_functions(cout);
            return false;
        }
        else if (arg == "--list-experiments")
        {
            print_experiments(cout);
            return false;
        }
        else if (arg == "-f" || arg == "--function")
        {
            job.function = value();
            cl.function_given = true;
        }
        else if (arg == "-o" || arg == "--output")
            job.output = value();
        else if (arg == "-e" || arg == "--experiment")
            job.experiment = value();
        else if (arg == "-n" || arg == "--iterations")
            job.remesh.iterations = static_cast<int>(parse_number(arg, value()));
        else if (arg == "-a" || arg == "--angle")
            job.remesh.min_angle = static_cast<float>(parse_number(arg, value()));
        else if (arg == "-d" || arg == "--divisor")
            job.remesh.area_divisor =
                static_cast<float>(parse_number(arg, value()));
        else if (arg == "--no-gui")
            cl.window = false;
        else if (!arg.empty() && arg[0] == '-')
            throw invalid_argument("unknown option " + arg);
        else if (!cl.input_given)
        {
            job.input = arg;
            cl.input_given = true;
        }
        else
            throw invalid_argument("unexpected argument '" + arg + "'");
    }

    // Check the values now, so mistakes show up before anything starts
    if (!job.function.empty() && !find_surface(job.function))
    {
        cerr << "Unknown function '" << job.function << "'. Available:\n";
        print_functions(cerr);
        throw runtime_error("unknown function");
    }
    if (!job.experiment.empty() && !find_experiment(job.experiment))
    {
        cerr << "Unknown experiment '" << job.experiment << "'. Available:\n";
        print_experiments(cerr);
        throw runtime_error("unknown experiment");
    }
    if (job.remesh.iterations < 0)
        throw invalid_argument("--iterations must not be negative");
    if (job.remesh.area_divisor <= 0.0f)
        throw invalid_argument("--divisor must be positive");
    if (!cl.window && !(cl.input_given && cl.function_given))
        throw invalid_argument("--no-gui needs an input mesh and --function");
    return true;
}

} // namespace iso
