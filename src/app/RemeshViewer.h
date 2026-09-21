// The window: 3D view of the mesh, the setup dialog and the statistics panel.

#pragma once

#include "app/FileBrowser.h"
#include "method/Job.h"

#include <pmp/visualization/mesh_viewer.h>

#include <filesystem>
#include <future>
#include <optional>
#include <string>
#include <vector>

namespace iso {

//! pmp's mesh viewer plus a setup dialog: choose the input mesh, output file,
//! implicit surface and remeshing parameters, run, and see the result.
//! The job runs on a worker thread so the window stays responsive.
class RemeshViewer : public pmp::MeshViewer
{
public:
    //! \p settings_file is where the dialog's values are saved (on every run
    //! and when the window closes); leave empty to disable saving.
    RemeshViewer(const JobSettings& initial,
                 std::filesystem::path settings_file);
    ~RemeshViewer() override;

    //! Runs \p initial immediately if \p run_immediately (and it is complete),
    //! otherwise opens the setup dialog on the first frame.
    void start(bool run_immediately);

protected:
    void process_imgui() override;
    void drop(int count, const char** paths) override;

private:
    void draw_setup_dialog();
    void draw_report();

    void launch_job();
    void collect_finished_job();
    bool busy() const { return job_.valid(); }
    JobSettings current_settings() const;
    void save_current_settings() const;

    // Text buffers edited by ImGui
    char input_[1024] = "";
    char output_[1024] = "output.obj";
    int function_index_ = -1;
    int experiment_index_ = -1; //!< -1 = none
    RemeshOptions remesh_;

    bool open_setup_ = false;
    bool dialog_running_ = false; //!< the dialog's Run button started the job
    FileBrowser input_browser_{FileBrowser::Mode::OpenFile};
    FileBrowser output_browser_{FileBrowser::Mode::Folder};
    std::filesystem::path settings_file_;

    std::future<JobReport> job_;
    std::string job_output_; //!< output path of the job in flight
    std::optional<JobReport> report_;
    std::string error_;
};

} // namespace iso
