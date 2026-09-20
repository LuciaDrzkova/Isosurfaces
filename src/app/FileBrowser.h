// A small "choose a file" dialog drawn inside the window (no system dialog).

#pragma once

#include <filesystem>
#include <string>
#include <vector>

namespace iso {

class FileBrowser
{
public:
    //! Show the dialog (from the next draw() call), starting in \p start_dir.
    void open(const std::filesystem::path& start_dir);

    //! Call every frame. Returns true in the frame where the user picked a
    //! file, and stores its path in \p picked. \p scale is the window's UI scale.
    bool draw(float scale, std::string& picked);

private:
    struct Entry
    {
        std::string name;
        bool is_dir;
    };

    void go_to(const std::filesystem::path& target);

    std::filesystem::path dir_;
    std::vector<Entry> entries_; //!< sub-folders and mesh files of dir_
    std::string error_;
    bool open_requested_ = false;
};

} // namespace iso
