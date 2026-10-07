# Isosurfaces

Research code for studying topology changes and local mesh refinement of an implicit surface

\[
f(x,y,z) = x^2 + y^2 - z^2 - c.
\]

The current executable uses an in-project **classic 256-case Marching Cubes** implementation. It samples the implicit function on a regular grid, extracts the isosurface, reports mesh topology, classifies critical mesh vertices from the analytic/numerical gradient and Hessian, groups coincident singular representatives geometrically, and can apply conforming local triangle refinement around a detected singularity.

For the parameter family above:

- `c < 0`: two-sheet hyperboloid
- `c = 0`: double cone with a non-degenerate Morse saddle at the origin
- `c > 0`: one-sheet hyperboloid

The finite sampling box may clip the surface, so boundary edges are expected in those experiments.

## Build

The project uses CMake and downloads the `pmp-library` dependency automatically.

```bash
cmake -S . -B build
cmake --build build --config Release --parallel
```

## Run

The command line uses positional arguments:

```text
Isosurfaces [c] [resolution] [extent] [local-radius] [local-levels] [--no-gui]
```

Defaults:

```text
c             = 0
resolution    = 129
extent        = 2
local-radius  = 0.3
local-levels  = 2
```

Examples:

```bash
./build/Isosurfaces 0 129 2 0.3 2 --no-gui
./build/Isosurfaces -0.25 129 2 --no-gui
./build/Isosurfaces 0.25 129 2 --no-gui
```

A resolution such as `129` is useful for the cone experiment because the grid then contains the origin exactly when sampling `[-2,2]^3`.

Outputs are written under `outputs/` as PLY and OBJ files.

## What is measured

The console report includes vertex/edge/triangle counts, boundary edges and vertices, non-manifold edges and vertices, connected components, coincident vertex groups, and singularity-classification counts.

At `c=0`, the geometric singularity grouping can contain multiple mesh representatives at the same physical point. They are intentionally **not merged in the mesh**: they may belong to different connected components of the extracted surface. The report groups them geometrically while preserving their topology.

## Local refinement

The local stage currently performs **conforming triangle subdivision plus projection back to the implicit surface**. Region selection supports both:

- Euclidean sphere selection by face centroid; and
- topological/geodesic selection by a weighted BFS/Dijkstra traversal over the triangle adjacency graph.

The default is the topological mode. Both selection counts are reported so the geometric and topological definitions can be compared directly.

For every refinement level the program reports selected faces, split/interface edges and vertices, surface boundary edges and vertices, newly created triangles/vertices, and projection failures.

This stage is intentionally separate from the global voxel Marching Cubes extraction. A voxel-local Marching Cubes patch with exact interface stitching is a future refinement of the experimental pipeline rather than something silently claimed by the current implementation.

## Tests

CTest is enabled and the project builds a current-pipeline test executable:

```bash
ctest --test-dir build --output-on-failure
```

The tests cover cone singularity classification, hyperboloid regularity and Hessian eigenvalues, numerical derivative fallbacks, Euler characteristic of a closed Marching Cubes sphere, and local-region/interface/boundary statistics.

## Source layout

```text
include/
  ImplicitSurface.hpp       implicit surface family and exact cone derivatives
  MarchingCubes.hpp         classic 256-case Marching Cubes
  SingularityDetector.hpp   regular / non-degenerate / degenerate classification
  GradientHessian.hpp       safe numerical derivative helper
  LocalUnfolder.hpp         local region selection and conforming refinement
src/
  main.cpp                  current research executable
  app/MyViewer.h            PMP viewer wrapper

tests/
  tests.cpp                 current-pipeline tests
```

Legacy MeshLab/remeshing entry points are not part of the current CMake target. The active experiment is driven by `src/main.cpp`.
