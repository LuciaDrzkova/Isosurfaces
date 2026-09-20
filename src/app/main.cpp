// Where the program starts. Two ways to run it:
//   - with a window (default): choose everything in a dialog
//   - with --no-gui: everything comes from the command line, result goes to a file

#include "app/CommandLine.h"
#include "app/RemeshViewer.h"
#include "app/Settings.h"
#include "method/Job.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace std;

namespace {

// --no-gui: run once and print the statistics
int run_without_window(const iso::JobSettings& job)
{
    const auto report = iso::run_job(job);

    iso::print_metrics(cout, "INPUT MESH", report.input);
    iso::print_metrics(cout, "PROJECTED MESH", report.projected);
    for (size_t i = 0; i < report.collapses_per_iteration.size(); ++i)
        cout << "Iteration " << i + 1 << ": "
             << report.collapses_per_iteration[i] << " collapses\n";
    cout << "Remeshing took " << report.seconds << " s\n";
    if (!report.experiment_log.empty())
        cout << "=== EXPERIMENT " << job.experiment << " ===\n"
             << report.experiment_log;
    iso::print_metrics(cout, "OUTPUT MESH", report.output);
    cout << "Wrote " << job.output << '\n';
    return EXIT_SUCCESS;
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        iso::CommandLine cl;

        // The window starts from the last used settings; anything typed on the
        // command line overrides them. Without a window the saved settings are
        // ignored, so scripts always do the same thing.
        const auto settings_file = iso::default_settings_path();
        const bool no_window = any_of(argv + 1, argv + argc, [](const char* a) {
            return string(a) == "--no-gui";
        });
        if (!no_window)
            cl.job = iso::load_settings(settings_file);

        if (!iso::parse_command_line(argc, argv, cl))
            return EXIT_SUCCESS;

        if (!cl.window)
            return run_without_window(cl.job);

        iso::RemeshViewer viewer(cl.job, settings_file);
        // Only a complete command line skips the dialog, saved values don't
        viewer.start(cl.input_given && cl.function_given);
        return viewer.run();
    }
    catch (const invalid_argument& e)
    {
        cerr << "Error: " << e.what() << "\n\n";
        iso::print_usage(cerr);
        return EXIT_FAILURE;
    }
    catch (const exception& e)
    {
        cerr << "Error: " << e.what() << '\n';
        return EXIT_FAILURE;
    }
}
