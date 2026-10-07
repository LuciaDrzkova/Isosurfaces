# Isosurfaces

Research code for studying topology changes and local refinement of an implicit surface.

The current executable uses an in-project classic 256-case Marching Cubes implementation to extract

\[
f(x,y,z)=x^2+y^2-z^2-c.
\]

The pipeline is:

1. global voxel sampling and classic 256-case Marching Cubes,
2. mesh topology diagnostics,
3. critical-point/singularity classification from the gradient and Hessian,
4. geometric grouping of coincident singular mesh representatives,
5. conforming local triangle subdivision with projection back to the implicit surface, and
6. an independent local voxel Marching Cubes extraction around the detected singularity.

The local voxel extractor is the next-stage voxel-based refinement path. It deliberately runs as a separate local patch at this stage. The existing `LocalUnfolder` output remains the conforming triangle-subdivision baseline used for visualization and transition-stitching measurements. Direct replacement of that baseline by a locally re-extracted voxel patch is the subsequent stitching step.

For the parameter family above:

- `c < 0`: two-sheet hyperboloid
- `c = 0`: double cone with a non-degenerate Morse saddle at the origin
- `c > 0`: one-sheet hyperboloid

Because the global implicit surface is clipped to a finite Marching Cubes box, its outer boundary edges are expected. These are different from the local refinement interface between selected and untouched surface regions.

## Build

```bash
cmake -S . -B build
cmake --build build --config Release --parallel
```

Tests are enabled through CTest:

```bash
ctest --test-dir build --output-on-failure
```

## Run

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

A resolution such as `129` is useful for the cone experiment because the grid then contains the origin exactly on `[-2,2]^3`.

## Local refinement measurements

The program distinguishes:

- **global surface boundary**: edges where the extracted surface is clipped by the finite sampling box;
- **refinement interface**: the interface between refined and untouched portions of the surface;
- **local voxel MC patch boundary**: the artificial boundary of the independently extracted local patch.

The local voxel Marching Cubes stage reports its aligned local box, refinement factor, coarse/refined spacing, local grid resolution, patch topology, and extraction time for every level.

The conforming local subdivision stage supports two region definitions: Euclidean sphere selection and a multi-source topological BFS seeded from all mesh representatives of the geometric singularity. Both counts are reported; the default mode is topological BFS. The stage reports refinement-interface counts, new vertices/triangles, projection failures, scalar-value diagnostics, and final mesh topology.

The two local methods are intentionally kept separate so that the thesis can compare voxel sampling against the conforming triangle baseline before implementing direct voxel-patch replacement and transition stitching.
