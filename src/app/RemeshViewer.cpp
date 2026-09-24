// The window (see RemeshViewer.h). Only user-interface code lives here; the
// actual work is done by run_job() in method/Job.cpp.

#include "app/RemeshViewer.h"

#include "method/Experiments.h"
#include "app/Settings.h"

#include <imgui.h>

#include <algorithm>
#include <cstdio>
#include <iomanip>
#include <limits>
#include <sstream>

using namespace std;
namespace fs = filesystem;

namespace iso {

namespace {

template <size_t N>
void set_text(char (&buffer)[N], const string& text)
{
    snprintf(buffer, N, "%s", text.c_str());
}

// Formats via printf so the table cells stay compact ("1.2e-07" etc.)
template <typename... Args>
string format(const char* fmt, Args... args)
{
    char buffer[64];
    snprintf(buffer, sizeof buffer, fmt, args...);
    return buffer;
}

string exact_metric(double value)
{
    ostringstream stream;
    stream << setprecision(numeric_limits<double>::max_digits10) << value;
    return stream.str();
}

} // namespace

// ---------------------------------------------------------------------------
// Viewer
// ---------------------------------------------------------------------------

RemeshViewer::RemeshViewer(const JobSettings& initial, fs::path settings_file)
    : MeshViewer("Isosurfaces", 1200, 800),
      remesh_(initial.remesh),
      settings_file_(std::move(settings_file))
{
    set_draw_mode("Smooth Shading");
    set_text(input_, initial.input);
    set_text(output_, initial.output);

    const auto& surfaces = builtin_surfaces();
    for (size_t i = 0; i < surfaces.size(); ++i)
        if (surfaces[i].name() == initial.function)
            function_index_ = static_cast<int>(i);

    const auto& all = experiments();
    for (size_t i = 0; i < all.size(); ++i)
        if (all[i].name == initial.experiment)
            experiment_index_ = static_cast<int>(i);
}

RemeshViewer::~RemeshViewer()
{
    save_current_settings();
}

void RemeshViewer::save_current_settings() const
{
    if (settings_file_.empty())
        return;

    // Store absolute paths so they still work when launched from another folder
    JobSettings settings = current_settings();
    error_code ec;
    for (string* path : {&settings.input, &settings.output})
    {
        const fs::path absolute = fs::absolute(*path, ec);
        if (!path->empty() && !ec)
            *path = absolute.string();
    }
    save_settings(settings_file_, settings);
}

void RemeshViewer::start(bool run_immediately)
{
    if (run_immediately && input_[0] && function_index_ >= 0)
        launch_job();
    else
        open_setup_ = true;
}

JobSettings RemeshViewer::current_settings() const
{
    JobSettings s;
    s.input = input_;
    s.output = output_;
    if (function_index_ >= 0)
        s.function = builtin_surfaces()[function_index_].name();
    if (experiment_index_ >= 0)
        s.experiment = experiments()[experiment_index_].name;
    s.remesh = remesh_;
    return s;
}

void RemeshViewer::launch_job()
{
    error_.clear();
    save_current_settings();
    job_output_ = output_;
    job_ = async(launch::async,
                 [settings = current_settings()] { return run_job(settings); });
}

void RemeshViewer::collect_finished_job()
{
    if (!busy() || job_.wait_for(chrono::seconds(0)) != future_status::ready)
        return;

    try
    {
        report_ = job_.get();
        print_metrics(cout, "INPUT MESH", report_->input);
        print_metrics(cout, "PROJECTED MESH", report_->projected);
        for (size_t i = 0; i < report_->collapses_per_iteration.size(); ++i)
            cout << "Iteration " << i + 1 << ": "
                 << report_->collapses_per_iteration[i] << " collapses\n";
        cout << "Remeshing took " << exact_metric(report_->seconds) << " s\n";
        if (!report_->experiment_log.empty())
            cout << "=== EXPERIMENT " << current_settings().experiment << " ===\n"
                 << report_->experiment_log;
        print_metrics(cout, "OUTPUT MESH", report_->output);
        cout << "Wrote " << job_output_ << '\n';
        load_mesh(job_output_.c_str()); // GL work must happen on this thread
    }
    catch (const exception& e)
    {
        report_.reset();
        error_ = e.what();
        open_setup_ = true; // let the user fix the settings
    }
}

void RemeshViewer::drop(int count, const char** paths)
{
    if (count > 0 && !busy())
    {
        set_text(input_, paths[0]);
        open_setup_ = true;
    }
}

void RemeshViewer::process_imgui()
{
    collect_finished_job();
    MeshViewer::process_imgui(); // "Mesh Info"

    if (ImGui::CollapsingHeader("Remeshing", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::BeginDisabled(busy());
        if (ImGui::Button("Setup..."))
            open_setup_ = true;
        ImGui::EndDisabled();

        if (busy())
            ImGui::TextUnformatted("Working...");
        else if (!error_.empty())
            ImGui::TextColored(ImVec4(1, 0.4f, 0.4f, 1), "%s", error_.c_str());
        draw_report();
    }

    draw_setup_dialog();
}

void RemeshViewer::draw_report()
{
    if (!report_)
        return;

    const ImGuiTableFlags flags =
        ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg;
    if (ImGui::BeginTable("stats", 4, flags))
    {
        ImGui::TableSetupColumn("");
        ImGui::TableSetupColumn("Input");
        ImGui::TableSetupColumn("Projected");
        ImGui::TableSetupColumn("Output");
        // pmp uses black text on a light panel, but the header background would
        // stay dark, so give it a light one
        ImGui::PushStyleColor(ImGuiCol_TableHeaderBg,
                              ImVec4(0.78f, 0.87f, 0.98f, 1.0f));
        ImGui::TableHeadersRow();
        ImGui::PopStyleColor();

        auto row = [&](const char* label, auto&& cell) {
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            ImGui::TextUnformatted(label);
            for (const MeshMetrics* m :
                 {&report_->input, &report_->projected, &report_->output})
            {
                ImGui::TableNextColumn();
                ImGui::TextUnformatted(cell(*m).c_str());
            }
        };

        row("Vertices", [](const MeshMetrics& m) { return format("%zu", m.n_vertices); });
        row("Triangles", [](const MeshMetrics& m) { return format("%zu", m.n_faces); });
        row("Avg |f|", [](const MeshMetrics& m) { return format("%.2e", m.avg_distance); });
        row("Max |f|", [](const MeshMetrics& m) { return format("%.2e", m.max_distance); });
        row("Skinny tris", [](const MeshMetrics& m) { return format("%zu", m.faces_with_small_angle); });
        row("Regular verts", [](const MeshMetrics& m) { return format("%.1f%%", m.regular_vertex_percentage); });
        ImGui::EndTable();
    }

    const int collapses = [&] {
        int total = 0;
        for (int c : report_->collapses_per_iteration)
            total += c;
        return total;
    }();
    ImGui::Text("%d collapses in %.3f s", collapses, report_->seconds);
    if (!report_->experiment_log.empty())
    {
        ImGui::SeparatorText("Experiment");
        ImGui::TextWrapped("%s", report_->experiment_log.c_str());
    }
    ImGui::TextDisabled("Saved to %s", job_output_.c_str());
}

void RemeshViewer::draw_setup_dialog()
{
    const float s = imgui_scaling();

    if (open_setup_)
    {
        ImGui::OpenPopup("Setup");
        open_setup_ = false;
    }
    ImGui::SetNextWindowPos(ImGui::GetMainViewport()->GetCenter(),
                            ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
    if (!ImGui::BeginPopupModal("Setup", nullptr,
                                ImGuiWindowFlags_AlwaysAutoResize))
        return;

    const auto& surfaces = builtin_surfaces();

    ImGui::BeginDisabled(busy());

    // --- Files -------------------------------------------------------------
    ImGui::SeparatorText("Input mesh");
    ImGui::SetNextItemWidth(400 * s);
    ImGui::InputText("##input", input_, sizeof input_);
    ImGui::SameLine();
    if (ImGui::Button("Browse..."))
        input_browser_.open(fs::path(input_).parent_path());
    ImGui::TextDisabled("Triangle mesh (.obj .off .stl .ply ...). "
                        "You can also drop a file onto the window.");

    // --- Surface -----------------------------------------------------------
    ImGui::SeparatorText("Implicit surface");
    ImGui::SetNextItemWidth(250 * s);
    const char* preview =
        function_index_ >= 0 ? surfaces[function_index_].name().c_str()
                             : "Select a surface...";
    if (ImGui::BeginCombo("##function", preview))
    {
        for (int i = 0; i < static_cast<int>(surfaces.size()); ++i)
            if (ImGui::Selectable(surfaces[i].name().c_str(), i == function_index_))
                function_index_ = i;
        ImGui::EndCombo();
    }
    if (function_index_ >= 0)
        ImGui::TextDisabled("f = 0:  %s", surfaces[function_index_].description().c_str());

    // --- Parameters --------------------------------------------------------
    ImGui::SeparatorText("Remeshing");
    ImGui::SetNextItemWidth(200 * s);
    ImGui::InputInt("Iterations", &remesh_.iterations);
    remesh_.iterations = clamp(remesh_.iterations, 0, 1000);
    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("Outer passes; the target triangle area is re-measured on each.");

    // A slider for quick changes and a number field to type an exact value
    auto slider_with_input = [&](const char* label, float* value, float min,
                                 float max, const char* tooltip,
                                 ImGuiSliderFlags flags = 0) {
        ImGui::PushID(label);
        ImGui::SetNextItemWidth(200 * s);
        ImGui::SliderFloat("##slider", value, min, max, "", flags); // value is shown in the field
        const bool over_slider = ImGui::IsItemHovered();
        ImGui::SameLine();
        ImGui::SetNextItemWidth(80 * s);
        ImGui::InputFloat("##input", value, 0.0f, 0.0f, "%.2f");
        const bool over_input = ImGui::IsItemHovered();
        ImGui::SameLine();
        ImGui::TextUnformatted(label);
        if (over_slider || over_input || ImGui::IsItemHovered())
            ImGui::SetTooltip("%s", tooltip);
        ImGui::PopID();
    };

    slider_with_input("Min angle (deg)", &remesh_.min_angle, 0.0f, 60.0f,
                      "Triangles with an angle below this are collapsed.");
    slider_with_input("Area divisor", &remesh_.area_divisor, 0.1f, 10.0f,
                      "Triangles smaller than (median area / divisor) are "
                      "collapsed.\nLower values coarsen the mesh more.",
                      ImGuiSliderFlags_Logarithmic);

    // Typed values may be outside the slider's range, but not nonsensical
    remesh_.min_angle = clamp(remesh_.min_angle, 0.0f, 180.0f);
    remesh_.area_divisor = clamp(remesh_.area_divisor, 0.01f, 1000.0f);

    // --- Experiment ----------------------------------------------------------
    ImGui::SeparatorText("Experiment (optional)");
    const auto& all_experiments = experiments();
    ImGui::SetNextItemWidth(250 * s);
    const char* experiment_preview =
        experiment_index_ >= 0 ? all_experiments[experiment_index_].name.c_str()
                               : "None";
    if (ImGui::BeginCombo("##experiment", experiment_preview))
    {
        if (ImGui::Selectable("None", experiment_index_ < 0))
            experiment_index_ = -1;
        for (int i = 0; i < static_cast<int>(all_experiments.size()); ++i)
            if (ImGui::Selectable(all_experiments[i].name.c_str(),
                                  i == experiment_index_))
                experiment_index_ = i;
        ImGui::EndCombo();
    }
    if (experiment_index_ >= 0)
        ImGui::TextDisabled("%s", all_experiments[experiment_index_].description.c_str());

    ImGui::SeparatorText("Output mesh");
    ImGui::SetNextItemWidth(400 * s);
    ImGui::InputText("##output", output_, sizeof output_);
    ImGui::SameLine();
    if (ImGui::Button("Browse...##output"))
        output_browser_.open(fs::path(output_).parent_path());

    ImGui::EndDisabled();

    // --- Actions -----------------------------------------------------------
    ImGui::Separator();
    if (busy())
        ImGui::TextUnformatted("Working...");
    else if (!error_.empty())
        ImGui::TextColored(ImVec4(1, 0.4f, 0.4f, 1), "%s", error_.c_str());

    const bool ready = input_[0] && output_[0] && function_index_ >= 0;
    ImGui::BeginDisabled(busy() || !ready);
    if (ImGui::Button("Run", ImVec2(100 * s, 0)))
    {
        launch_job();
        dialog_running_ = true;
    }
    ImGui::EndDisabled();
    ImGui::SameLine();
    ImGui::BeginDisabled(busy());
    if (ImGui::Button("Close", ImVec2(100 * s, 0)))
        ImGui::CloseCurrentPopup();
    ImGui::EndDisabled();

    // Close on success; on failure stay open and show the error
    if (dialog_running_ && !busy())
    {
        dialog_running_ = false;
        if (error_.empty())
            ImGui::CloseCurrentPopup();
    }

    // The file browsers open on top of this dialog
    string picked;
    if (input_browser_.draw(s, picked))
        set_text(input_, picked);
    if (output_browser_.draw(s, picked))
    {
        // Change the folder, keep the file name ("output.obj" if none)
        string name = fs::path(output_).filename().string();
        if (name.empty())
            name = "output.obj";
        set_text(output_, (fs::path(picked) / name).string());
    }

    ImGui::EndPopup();
}

} // namespace iso
