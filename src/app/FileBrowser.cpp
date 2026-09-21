#include "app/FileBrowser.h"

#include <imgui.h>

#include <algorithm>
#include <cctype>
#include <cstdlib>

using namespace std;
namespace fs = filesystem;

namespace iso {

namespace {

string lower(string s)
{
    transform(s.begin(), s.end(), s.begin(),
              [](unsigned char c) { return static_cast<char>(tolower(c)); });
    return s;
}

// Only files pmp can read are shown
bool is_mesh_file(const fs::path& p)
{
    static const char* extensions[] = {".obj", ".off", ".stl", ".ply",
                                       ".xyz", ".agi", ".pmp"};
    const string ext = lower(p.extension().string());
    return any_of(begin(extensions), end(extensions),
                  [&](const char* e) { return ext == e; });
}

fs::path home_directory()
{
    if (const char* home = getenv("HOME"))
        return home;
    if (const char* profile = getenv("USERPROFILE"))
        return profile;
    return fs::current_path();
}

} // namespace

FileBrowser::FileBrowser(Mode mode) : mode_(mode) {}

void FileBrowser::open(const fs::path& start_dir)
{
    error_code ec;
    fs::path p = start_dir.empty() ? fs::current_path()
                                   : fs::absolute(start_dir, ec);
    if (ec || !fs::is_directory(p, ec))
        p = fs::current_path();
    go_to(p);
    open_requested_ = true;
}

void FileBrowser::go_to(const fs::path& target)
{
    error_ = "";

    error_code ec;
    fs::directory_iterator it(
        target, fs::directory_options::skip_permission_denied, ec);
    if (ec)
    {
        error_ = ec.message(); // keep showing the previous folder
        return;
    }

    vector<Entry> found;
    for (const auto& e : it)
    {
        const string name = e.path().filename().string();
        if (name.empty() || name[0] == '.') // hidden
            continue;
        error_code type_ec;
        const bool is_dir = e.is_directory(type_ec);
        if (is_dir || is_mesh_file(e.path()))
            found.push_back({name, is_dir});
    }

    // Folders first, then alphabetical
    sort(found.begin(), found.end(), [](const Entry& a, const Entry& b) {
        if (a.is_dir != b.is_dir)
            return a.is_dir;
        return lower(a.name) < lower(b.name);
    });
    entries_ = std::move(found);
    dir_ = target;
}

bool FileBrowser::draw(float scale, string& picked)
{
    const bool folders = mode_ == Mode::Folder;
    const char* title = folders ? "Choose output folder" : "Choose mesh";

    if (open_requested_)
    {
        ImGui::OpenPopup(title);
        open_requested_ = false;
    }

    ImGui::SetNextWindowSize(ImVec2(520 * scale, 0), ImGuiCond_Appearing);
    if (!ImGui::BeginPopupModal(title, nullptr,
                                ImGuiWindowFlags_AlwaysAutoResize))
        return false;

    bool chosen = false;
    fs::path next_dir; // where to go after this frame's clicks

    ImGui::TextWrapped("%s", dir_.string().c_str());
    if (ImGui::Button("Up") && dir_.has_parent_path())
        next_dir = dir_.parent_path();
    ImGui::SameLine();
    if (ImGui::Button("Home"))
        next_dir = home_directory();
    if (!error_.empty())
        ImGui::TextColored(ImVec4(1, 0.4f, 0.4f, 1), "%s", error_.c_str());

    if (ImGui::BeginChild("files", ImVec2(0, 300 * scale), true))
    {
        for (const auto& e : entries_)
        {
            const string label = (e.is_dir ? "[dir]  " : "       ") + e.name;
            if (folders && !e.is_dir)
            {
                ImGui::TextDisabled("%s", label.c_str()); // shown for context
                continue;
            }
            if (!ImGui::Selectable(label.c_str(), false))
                continue;
            if (e.is_dir)
                next_dir = dir_ / e.name;
            else
            {
                picked = (dir_ / e.name).string();
                chosen = true;
            }
        }
        if (entries_.empty())
            ImGui::TextDisabled("(no folders or mesh files here)");
    }
    ImGui::EndChild();

    if (!next_dir.empty())
        go_to(next_dir);

    const bool cancelled = ImGui::Button("Cancel");
    if (folders)
    {
        // Bottom right, on the same row as Cancel
        const float width = ImGui::CalcTextSize("Use this folder").x +
                            2 * ImGui::GetStyle().FramePadding.x;
        ImGui::SameLine();
        ImGui::SetCursorPosX(ImGui::GetWindowContentRegionMax().x - width);
        if (ImGui::Button("Use this folder"))
        {
            picked = dir_.string();
            chosen = true;
        }
    }

    if (cancelled || chosen)
        ImGui::CloseCurrentPopup();
    ImGui::EndPopup();
    return chosen;
}

} // namespace iso
