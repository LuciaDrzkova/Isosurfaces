// Remembers the last used settings (see Settings.h).

#include "app/Settings.h"

#include "method/Experiments.h"
#include "method/ImplicitSurface.h"

#include <cstdlib>
#include <fstream>

using namespace std;
namespace fs = filesystem;

namespace iso {

namespace {

string trim(const string& s)
{
    const char* ws = " \t\r\n";
    const size_t begin = s.find_first_not_of(ws);
    if (begin == string::npos)
        return "";
    return s.substr(begin, s.find_last_not_of(ws) - begin + 1);
}

bool parse_double(const string& text, double& out)
{
    char* end = nullptr;
    const double value = strtod(text.c_str(), &end);
    if (text.empty() || *end != '\0')
        return false;
    out = value;
    return true;
}

} // namespace

fs::path default_settings_path()
{
    fs::path base;
#if defined(_WIN32)
    if (const char* appdata = getenv("APPDATA"))
        base = appdata;
#elif defined(__APPLE__)
    if (const char* home = getenv("HOME"))
        base = fs::path(home) / "Library" / "Application Support";
#else
    if (const char* xdg = getenv("XDG_CONFIG_HOME"); xdg && *xdg)
        base = xdg;
    else if (const char* home = getenv("HOME"))
        base = fs::path(home) / ".config";
#endif
    if (base.empty())
        return {};
    return base / "Isosurfaces" / "settings.ini";
}

JobSettings load_settings(const fs::path& file)
{
    JobSettings s;
    ifstream in(file);
    string line;
    while (getline(in, line))
    {
        const size_t eq = line.find('=');
        if (eq == string::npos)
            continue;
        const string key = trim(line.substr(0, eq));
        const string value = trim(line.substr(eq + 1));

        double number = 0;
        if (key == "input")
            s.input = value;
        else if (key == "output" && !value.empty())
            s.output = value;
        else if (key == "function" && find_surface(value))
            s.function = value;
        else if (key == "experiment" && find_experiment(value))
            s.experiment = value;
        else if (key == "iterations" && parse_double(value, number) &&
                 number >= 0 && number <= 1000)
            s.remesh.iterations = static_cast<int>(number);
        else if (key == "min_angle" && parse_double(value, number) &&
                 number >= 0 && number <= 180)
            s.remesh.min_angle = static_cast<float>(number);
        else if (key == "area_divisor" && parse_double(value, number) &&
                 number > 0)
            s.remesh.area_divisor = static_cast<float>(number);
    }
    return s;
}

bool save_settings(const fs::path& file, const JobSettings& s)
{
    if (file.empty())
        return false;

    error_code ec;
    fs::create_directories(file.parent_path(), ec);

    ofstream out(file);
    if (!out)
        return false;
    out << "# Isosurfaces: last used settings (edited by the app)\n"
        << "input=" << s.input << '\n'
        << "output=" << s.output << '\n'
        << "function=" << s.function << '\n'
        << "experiment=" << s.experiment << '\n'
        << "iterations=" << s.remesh.iterations << '\n'
        << "min_angle=" << s.remesh.min_angle << '\n'
        << "area_divisor=" << s.remesh.area_divisor << '\n';
    return static_cast<bool>(out);
}

} // namespace iso
