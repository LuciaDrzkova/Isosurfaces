# Isosurfaces

Research code for studying topology changes and local refinement of an implicit surface.

The project now has three related local-refinement paths:

1. global voxel sampling + classic 256-case Marching Cubes,
2. `LocalUnfolder` conforming triangle subdivision as the baseline,
3. independent local voxel Marching Cubes followed by geometric seam stitching back into the global mesh.

The original parameter family is

\[
f(x,y,z)=x^2+y^2-z^2-c.
\]

For this family:

- `c < 0`: two-sheet hyperboloid,
- `c = 0`: double cone with a non-degenerate Morse saddle at the origin,
- `c > 0`: one-sheet hyperboloid.

The finite Marching Cubes sampling box intentionally creates an outer surface boundary. The voxel-patch stitcher removes a local spherical region, inserts the independently extracted finer patch, and bridges the two local boundary loops/loops with a conforming triangle strip. The original `LocalUnfolder` result is preserved as the comparison baseline.

A second non-quadric family is available:

\[
f(x,y,z)=x^2+y^2-z^2+\alpha(x^4+y^4+z^4)-c,
\qquad \alpha=0.05.
\]

This exercises the generic numerical gradient/Hessian paths in `SingularityDetector` and `LocalUnfolder` while retaining a non-degenerate saddle at the origin for `c=0`.

## Build

```bash
cmake -S . -B build
cmake --build build --config Release --parallel
ctest --test-dir build --output-on-failure
```

## Existing viewer

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

## Stitched local voxel patch viewer

The new `IsosurfacesStitched` executable runs the independent local voxel MC extraction and stitches the resulting patch back into the global MC mesh.

```bash
./build/IsosurfacesStitched 0 65 2 0.3 2
```

Non-quadric family:

```bash
./build/IsosurfacesStitched 0 65 2 0.3 2 --surface quartic
```

Disable the GUI and use topology output only:

```bash
./build/IsosurfacesStitched 0 65 2 0.3 2 --no-gui
```

The resulting OBJ is written to `outputs/stitched_voxel_*.obj`.

## Parameter sweep

`IsosurfacesSweep` reuses one global MC extraction per `(surface,c,resolution)` and records one CSV row for every `(radius,requested-levels,level)` combination. Each row contains the unified global topology report, the sphere-vs-topological `selection_changed` flag, `LocalUnfolder::ScalarStatistics`, local voxel patch statistics, and the stitched topology/seam statistics.

Default research sweep:

```bash
./build/IsosurfacesSweep
```

The output is:

```text
outputs/parameter_sweep.csv
```

A smaller smoke sweep:

```bash
./build/IsosurfacesSweep --quick
```

Fine sweep around the topology transition:

```bash
./build/IsosurfacesSweep \
  --surface cone \
  --c-min -0.30 \
  --c-max 0.30 \
  --c-step 0.01 \
  --resolutions 65,129 \
  --radii 0.20,0.30 \
  --levels 1,2 \
  --output outputs/parameter_sweep_cone.csv
```

Run the same experiment for the non-quadric family:

```bash
./build/IsosurfacesSweep \
  --surface quartic \
  --c-min -0.30 \
  --c-max 0.30 \
  --c-step 0.01 \
  --resolutions 65,129 \
  --radii 0.20,0.30 \
  --levels 1,2 \
  --output outputs/parameter_sweep_quartic.csv
```

Run the regular closed-surface control:

```bash
./build/IsosurfacesSweep --surface sphere --c-min 0 --c-max 0 --c-step 1
```

For the cone experiment, a resolution such as `129` places the origin exactly on the sampling grid for `[-2,2]^3`.

## Local-refinement measurements

The research data distinguishes:

- global surface boundary: clipping against the finite MC sampling box,
- conforming refinement interface: the `LocalUnfolder` refined/untouched interface,
- local voxel patch boundary: the artificial boundary of an independently extracted local patch,
- stitched seam: the bridge between the global outside mesh and the local voxel replacement.

`selection_changed` is computed by comparing the Euclidean sphere face selection against the multi-source topological BFS selection used by `LocalUnfolder`.

`ScalarStatistics` report refinement factor, spacing, selected voxels, scalar corner samples, scalar range, and scalar mean for each refinement level.
