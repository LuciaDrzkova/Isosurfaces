ISOSURFACES RESEARCH UPDATE
===========================

What this update changes
------------------------
1. Replaces include/ImplicitSurface.hpp with one consistent catalogue of nine
   CLI surface names. It includes aliases needed by the branch's existing tests.
2. Adds include/SurfaceFactory.hpp so `main.cpp` and the batch sweep use exactly
   the same list/factory instead of maintaining divergent copies.
3. Generalizes include/MeshGeometryCheck.hpp. The original cone still has its
   exact analytic distance overload; arbitrary surfaces use damped Newton
   projection along a numerical gradient. Failed projections remain failures
   and are included in the output counts.
4. Adds include/DualContouring.hpp as a uniform-grid QEF/Dual Contouring baseline.
   This is a compact reference implementation, not an octree/manifold-guaranteed
   variant. Singular/ambiguous configurations may remain non-manifold.
5. Replaces src/research_sweep.cpp. It sweeps every surface, compares classic MC
   and Dual Contouring, and records LocalUnfolder as a secondary comparison plus
   the proposed LocalMarchingCubes + LocalVoxelStitcher method. It writes
   outputs/parameter_sweep.csv in long format, with global and near-singularity
   geometric errors, projection failures, topology, selector choices, and stitch status.
6. Replaces CMakeLists.txt to include a second comparison-test executable while
   preserving the project's current executable targets.
7. Adds tests/comparison_tests.cpp and docs/METHOD_AND_LITERATURE.md.
8. patches/update_main.py applies the small necessary integrations to the
   existing 1,300-line src/main.cpp without replacing or deleting its viewer,
   exporters, and existing diagnostics. It creates a .bak file before editing.
9. tools/summarize_sweep.py summarizes the generated CSV by surface and method.

Install
-------
Copy the files from this package into the matching locations in your repository.
Replace CMakeLists.txt, include/ImplicitSurface.hpp, include/MeshGeometryCheck.hpp,
include/DualContouring.hpp, include/SurfaceFactory.hpp, and src/research_sweep.cpp.
Add tests/comparison_tests.cpp and the two files in docs/ (`METHOD_AND_LITERATURE.md`
and `SURVEY_NEAR_SINGULARITIES.md`). Then from the repository root run:

    python patches/update_main.py

If your project root does not already contain the `patches` folder, copy that
folder too. The script intentionally fails instead of guessing if the current
main.cpp differs from the expected branch source. Review the `.bak` if needed.

Build and tests
---------------
Run:

    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build --config Release
    ctest --test-dir build --output-on-failure

The repository needs network access on first CMake configure if its pmp-library
FetchContent dependency is not already cached.

Run a quick batch check
-----------------------
    ./build/IsosurfacesSweep --quick --output outputs/quick_sweep.csv

Summarize the generated results (per-surface, stitch rate and local triangle-RMS
reduction versus MC for Dual Contouring, LocalUnfolder, voxel patch and stitching):

    python tools/summarize_sweep.py outputs/quick_sweep.csv --markdown

Run a larger `c=0` singularity family comparison:

    ./build/IsosurfacesSweep --surface all --c-min 0 --c-max 0 --c-step 1 \
      --resolutions 33,65,97 --radii 0.20,0.30 --levels 1,2 \
      --output outputs/singularities.csv

Run the parameter sweep for the full family and range:

    ./build/IsosurfacesSweep --surface all --c-min -0.30 --c-max 0.30 \
      --c-step 0.05 --resolutions 33,65,97 --radii 0.20,0.30 \
      --levels 1,2 --output outputs/parameter_sweep.csv

Notes on interpreting results
-----------------------------
* For singular surfaces, MC/DC records are repeated per local radius and level so
  their local-neighborhood error is comparable against that exact configuration.
  Their extracted mesh/topology is identical across those repeated records.
* A `crossing-planes` stitch failure is recorded as a result, not hidden. Its
  singular set is the whole Z-axis, so a point-centered circular seam is not a
  model of the entire singular set. This is the intended boundary-case test for
  the current point-centered proposed method; a line-aware patch is a separate
  algorithmic extension, not something the implementation should falsely mark
  as successful.
* A generic distance estimate is not an exact global closest-point distance.
  Read its projection-failure columns together with its RMS/p95 error.
* The batch sweep now detects origin-centered singularities before inspecting
  extracted mesh vertices, avoiding the common failure where an even resolution
  does not put the singular point on a sampled grid vertex.
* The local voxel patch method is the proposed contribution. `LocalUnfolder` is
  reported under `local_unfolder_secondary`; do not frame both as co-equal main
  contributions in the method chapter.
* Do not claim one class of surfaces is “handled faithfully” until the CSV is
  generated, the resolutions/radii are compared, and representative meshes are
  inspected. The docs provide a recommended results table and literature notes.

Validation caveat
-----------------
The complete repository could not be cloned from this execution environment,
so these files have not been built against the full pmp dependency tree here.
The new comparison tests and sweep source passed C++17 syntax checks using compatible
interface stubs. Standalone Dual Contouring + generic-distance smoke tests generated
non-empty meshes for all nine supplied surfaces at N=25; the sphere case had 1,160
vertices, 2,316 triangles, and a converged 0.2-unit projection distance. The CSV
writer and summarizer were also checked on synthetic output, including 40 fields per
row and local-RMS comparisons. The stubs intentionally do not implement the real MC
extractor, so they cannot validate the full runtime test suite. Run the full commands
above on the repository to confirm integration. Do not treat the stub checks as a
substitute for the full project build/tests.
