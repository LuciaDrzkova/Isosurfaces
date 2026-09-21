# Isosurfaces

Remeshing a triangle mesh onto an implicit surface `f(x, y, z) = 0`.

The program takes a triangle mesh that roughly follows an implicit surface,
projects every vertex exactly onto the surface, and then coarsens and
regularises the mesh by collapsing small and skinny triangles while keeping all
vertices on the surface. Derivatives of `f` are computed exactly with automatic
differentiation. It has an interactive window and a command line mode.

[pmp-library](https://github.com/pmp-library/pmp-library) (mesh processing and
viewer) and [autodiff](https://github.com/autodiff/autodiff) are downloaded and
built automatically; you do not need to install them.

## Installation

### macOS

1. **Open Terminal** (press `Cmd + Space`, type *Terminal*, press Enter).

2. **Install the Xcode Command Line Tools** (the C++ compiler and Git):

   ```bash
   xcode-select --install
   ```

   Click *Install* in the dialog that appears and wait until it finishes.

3. **Install CMake.** The easiest way is [Homebrew](https://brew.sh). Install it
   by pasting the command from that page, then run:

   ```bash
   brew install cmake
   ```

   Alternatively, download the macOS `.dmg` from
   [cmake.org/download](https://cmake.org/download/), drag CMake into
   *Applications*, open it once, and choose *Tools → How to Install For Command
   Line Use*.

4. **Get the project.** Either download it from GitHub (*Code → Download ZIP*) and
   unzip it, or clone it:

   ```bash
   git clone https://github.com/LuciaDrzkova/Isosurfaces
   ```

5. **Build it.** In Terminal, go to the project folder and run:

   ```bash
   cd path/to/Isosurfaces
   cmake -S . -B build
   cmake --build build --config Release --parallel
   ```

   The first build downloads and compiles the dependencies and takes a few
   minutes. Later builds take seconds.

6. **Run it:**

   ```bash
   ./build/Isosurfaces
   ```

### Windows

1. **Install Visual Studio Build Tools** (the C++ compiler). Download them from
   [visualstudio.microsoft.com/visual-cpp-build-tools](https://visualstudio.microsoft.com/visual-cpp-build-tools/).
   In the installer, tick **Desktop development with C++** and click *Install*.
   (If you already have Visual Studio with that workload, skip this step.)

2. **Install CMake** from [cmake.org/download](https://cmake.org/download/)
   (*Windows x64 Installer*). In the installer, choose **Add CMake to the system
   PATH**.

3. **Install Git** from [git-scm.com/downloads](https://git-scm.com/downloads)
   with the default options.

4. **Get the project.** Either download it from GitHub (*Code → Download ZIP*) and
   unzip it, or clone it:

   ```powershell
   git clone https://github.com/LuciaDrzkova/Isosurfaces
   ```

5. **Open PowerShell** (Start menu → *PowerShell*) **after** the installations, so
   it finds the new programs, and build the project:

   ```powershell
   cd path\to\Isosurfaces
   cmake -S . -B build
   cmake --build build --config Release --parallel
   ```

   The first build downloads and compiles the dependencies and takes a few
   minutes. Later builds take seconds.

6. **Run it:**

   ```powershell
   .\build\Release\Isosurfaces.exe
   ```

### Rebuild and run in one step

After the first build, double-click **`run.command`** (macOS) or **`run.bat`**
(Windows) in the project folder. It rebuilds whatever changed and starts the
program. From a terminal: `./run.command` or `.\run.bat`.

### Optional: Visual Studio Code

With [Visual Studio Code](https://code.visualstudio.com/) and its *C/C++* and
*CMake Tools* extensions, open the project folder and press
`Cmd/Ctrl + Shift + B` to rebuild and start the program in one step, or `F5` to
run it under the debugger.

### If something goes wrong

- **`cmake` or `git` is not found:** close and reopen the terminal after
  installing, so it picks up the new programs.
- **Windows: "no CMAKE_CXX_COMPILER found":** the *Desktop development with C++*
  workload from step 1 is missing. Re-run the Visual Studio installer and add it.
- **The configure step fails while downloading:** check your internet
  connection and run the same command again.

## Using the program

Start the program and a setup dialog opens.

1. **Input mesh:** click *Browse...* and choose a triangle mesh (`.obj`, `.off`,
   `.stl`, `.ply`), or drop a file onto the window. A sample is included:
   `assets/sphere.obj`.
2. **Implicit surface:** choose the surface your mesh approximates. For the
   sample, choose `sphere`.
3. **Parameters:**

   | Parameter | Meaning |
   |---|---|
   | Iterations | Number of remeshing passes. |
   | Min angle | Triangles with an angle below this (degrees) are collapsed. |
   | Area divisor | Triangles with an area below *median area / divisor* are collapsed. Lower values give a coarser mesh. |

4. **Output mesh:** where the result is saved (default `output.obj`). Click
   *Browse...* to choose a different folder.
5. Click **Run**.

The result appears in the window, and the side panel compares the input mesh with
the output mesh (the remeshed result, which is saved): vertex and triangle
counts, how far the vertices are from the surface (`|f|`), and how many skinny
triangles remain. The console additionally prints the mesh right after the
projection step.
*Setup...* in the side panel reopens the dialog, and the program remembers your
last choices. In the 3D view, drag to rotate and scroll to zoom; press `?` to
list all key bindings.

For the sample, an area divisor of about `1.0` shows the effect well (642 to
about 300 vertices).

### Command line

The same can be run without a window, for example for batch processing:

```bash
./build/Isosurfaces assets/sphere.obj --function sphere --divisor 1.0 --no-gui -o result.obj
```

```text
Isosurfaces [input-mesh] [options]

  -f, --function <name>    implicit surface to project onto
  -o, --output <file>      output mesh (default: output.obj)
  -n, --iterations <N>     number of remeshing passes (default: 6)
  -a, --angle <degrees>    minimum angle (default: 25)
  -d, --divisor <x>        area divisor (default: 1.5)
      --no-gui             do not open a window
      --list-functions     list the available implicit surfaces
  -h, --help               show all options
```

## Method

1. **Projection.** Every vertex `p` is moved onto the surface with Newton steps
   `p ← p − f(p) ∇f(p) / |∇f(p)|²` until `|f(p)|` is negligible.
2. **Collapsing.** For every triangle that is small (area below *median /
   divisor*) or skinny (an angle below *min angle*), the shortest edge is
   collapsed. The remaining vertex moves to the centroid of its neighbours and is
   projected back onto the surface. Vertices next to a collapse are left alone
   until the next sweep, and boundary vertices are never moved.
3. **Repetition.** Step 2 is repeated until no more triangles qualify, for the
   chosen number of iterations; the median area is re-measured on each.

Limitations:

- Projection is local, so the input mesh should already be a rough
  approximation of the surface.
- The input must be a triangle mesh.
- An area divisor below 1 makes every triangle a candidate, so a closed mesh can
  shrink to almost nothing. Values between `1.0` and `2.0` are a good start.
- At singular points of the surface (`∇f = 0`) the projection is undefined.
