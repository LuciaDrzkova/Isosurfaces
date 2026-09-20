# Isosurfaces

A C++17 application for remeshing a mesh onto an implicit surface. It uses
[pmp-library](https://github.com/pmp-library/pmp-library) for mesh processing
and viewing, and [autodiff](https://github.com/autodiff/autodiff) for the
implicit-surface calculations. Both libraries are downloaded automatically
by CMake during the first configuration.

## Requirements

- CMake 3.16 or newer
- A C++17 compiler
- Git (required by CMake to download the dependencies)
- OpenGL support

On Windows, install **Build Tools for Visual Studio** with the **Desktop
development with C++** workload, then install CMake and Git.

On macOS, install Xcode Command Line Tools and CMake. On Linux, install a
C++17 compiler, CMake, Git, and the system OpenGL development packages.

## Input mesh

Place the input mesh at:

```text
assets/wavy_network.obj
```

The program currently reads this file automatically. CMake copies the contents
of `assets/` next to the executable after each build.

## Build

Run these commands from the project directory:

### Windows

```powershell
cmake -S . -B build
cmake --build build --config Release
```

### macOS or Linux

```bash
cmake -S . -B build
cmake --build build --config Release
```

The first build may take a few minutes while CMake downloads and builds the
dependencies.

## Run

### Windows

```powershell
.\build\Release\Isosurfaces.exe
```

### macOS or Linux

```bash
./build/Isosurfaces
```

The application opens a viewer, processes `wavy_network.obj`, and writes the
result to `output.obj` in the same directory as the executable.

## Configuration

The current implicit surface and processing parameters are defined near the
top of `src/main.cpp`. Rebuild the project after changing them.
