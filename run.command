#!/bin/bash
# Rebuilds what changed, then starts the program.
# macOS: double-click this file.  Terminal: ./run.command
# Anything after the name is passed on, e.g. ./run.command assets/sphere.obj -f sphere

cd "$(dirname "$0")" || exit 1

# Where cmake usually lives: Homebrew (Apple Silicon / Intel), the CMake.app from cmake.org
export PATH="/opt/homebrew/bin:/usr/local/bin:/Applications/CMake.app/Contents/bin:$PATH"

if ! command -v cmake > /dev/null; then
    echo "CMake is not installed (or not on the PATH)."
    echo
    echo "Install it, then run this file again:"
    echo "  - with Homebrew:  brew install cmake      (Homebrew: https://brew.sh)"
    echo "  - or download it: https://cmake.org/download/"
    echo "    (open CMake.app once, then Tools > How to Install For Command Line Use)"
    echo
    read -r -p "Press Enter to close. "
    exit 1
fi

if cmake -S . -B build && cmake --build build --config Release --parallel; then
    ./build/Isosurfaces "$@"
else
    echo
    echo "Build failed - see the messages above."
    read -r -p "Press Enter to close. "
    exit 1
fi
