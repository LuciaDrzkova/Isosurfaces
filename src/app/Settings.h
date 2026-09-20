// Remembers the last used settings between program starts.

#pragma once

#include "method/Job.h"

#include <filesystem>

namespace iso {

//! Per-user settings file (created on first save), e.g.
//! ~/Library/Application Support/Isosurfaces/settings.ini on macOS,
//! %APPDATA%\Isosurfaces\settings.ini on Windows and
//! ~/.config/Isosurfaces/settings.ini elsewhere.
//! Empty if the user's home directory cannot be determined.
std::filesystem::path default_settings_path();

//! Reads the last-used settings. A missing or unreadable file, malformed lines
//! and unknown surface/experiment names all fall back to the defaults, so this
//! never throws.
JobSettings load_settings(const std::filesystem::path& file);

//! Writes \p settings to \p file, creating its directory. Returns false on
//! failure (which callers may ignore: it only costs the convenience).
bool save_settings(const std::filesystem::path& file,
                   const JobSettings& settings);

} // namespace iso
