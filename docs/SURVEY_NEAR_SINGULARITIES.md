# Survey notes: isosurface extraction near singularities

This is a code-grounded survey draft for the project. It separates what the literature says about extraction methods from hypotheses that should only be claimed after running the project sweep.

## 1. Why singularities are a difficult test

For a regular implicit surface \(S=\{x\in\mathbb{R}^3:f(x)=0\}\), the gradient is nonzero at every point in the region of interest. The implicit-function theorem then locally represents the zero set as a smooth surface. At a singular point, \(f=0\) and \(\nabla f=0\), so this local graph guarantee fails. A grid cell may intersect several sheets, and the topology encoded by corner signs may not determine the intended local arrangement without additional tests. A singular set can also be non-isolated: for `crossing-planes`, \(f(x,y,z)=x^2-y^2\) has \(\nabla f=(2x,-2y,0)\), and both vanish along the entire Z-axis at level zero.

This project should therefore distinguish at least two independent questions:

1. **Geometric accuracy:** are extracted vertices and triangle interiors close to the target zero set?
2. **Topological fidelity:** are connectivity, boundary structure, and manifoldness consistent with the intended local surface arrangement?

Neither metric substitutes for the other. A mesh can lie close to the zero set while connecting the wrong sheets, and a mesh can have clean edge counts but still approximate the geometry poorly.

## 2. Classic Marching Cubes (implemented baseline)

Lorensen and Cline's classic Marching Cubes (MC) samples the scalar field on a regular grid, classifies each cube from its corner signs, and uses an edge-case table with linear interpolation to place vertices. It is simple, fast, and deterministic for a fixed grid. The implementation in this repository explicitly uses the classic 256-case table rather than MC33/Lewiner.

Near singularities, the important limitation is not that MC is universally incorrect; rather, classic case-table choices can become dependent on exact grid alignment and ambiguous configurations. If the singular point lies exactly on grid vertices or edges, multiple local sheets can meet at the same sampled position. The chosen triangulation can have different component counts or non-manifold vertices after coincident positions are welded. Rotating the surface relative to the grid is therefore a meaningful test, not a cosmetic one.

**Project role:** `classic_mc` is the baseline and the source mesh for the local patch integration. Record its original topology before comparing refinements.

## 3. Topology-resolving Marching Cubes variants (survey comparison)

Lewiner et al. describe an implementation of Marching Cubes 33 with ambiguity resolution and topological guarantees for the interpolated cell data. This differs from simply changing the sampling resolution: the extractor changes how ambiguous cases are resolved. It is a useful future baseline because it can help separate a limitation caused by the classic case table from one caused by the underlying scalar grid or the mathematical singularity itself.

This variant is **not implemented in this repository update**. Do not describe it in the thesis as an experimentally compared method until a compatible implementation is added and tested against the same field samples and bounds.

## 4. Marching Tetrahedra (literature alternative)

Marching Tetrahedra decomposes each cube into tetrahedra and extracts the zero set tetrahedron by tetrahedron. Its simple local cases avoid some of the ambiguities of a cube case table, but the chosen tetrahedralization can introduce directional bias, and triangle quality can be less favorable. Treece, Prager, and Gee's paper on regularised Marching Tetrahedra discusses triangle-shape issues and subsequent improvements.

This is another reasonable comparison candidate, but Dual Contouring is the implemented second extractor in this package because its QEF vertex placement offers a more different geometric mechanism than another edge-interpolation method. Marching Tetrahedra remains part of the literature survey, not a measured baseline here.

## 5. Dual Contouring (implemented comparison baseline)

Ju, Losasso, Schaefer, and Warren's Dual Contouring uses Hermite data (edge intersections and normals) and a quadratic error function (QEF) to position one dual vertex in each active cell. The approach can preserve sharp features better than methods whose vertices are constrained to grid edges when the relevant feature information is present.

A basic one-vertex-per-cell formulation is not topology-proof. If a cell contains multiple sheets or a topologically ambiguous configuration, a single dual vertex can merge them or lead to non-manifold connectivity. Later work specifically addresses topology ambiguity in Dual Contouring. These limitations make DC a useful comparison for this project's intersecting-sheet cases, but not a ground-truth oracle.

**Project role:** `dual_contouring` is a compact uniform-grid QEF reference implementation. It uses numerically estimated gradients at edge crossings and regularizes the QEF; it is not the full adaptive octree algorithm and should be labelled accordingly in the thesis.

## 6. LocalUnfolder (secondary project-specific comparison)

`LocalUnfolder` subdivides selected faces/edges of the existing triangle mesh and projects new points toward the implicit surface. Its region can be chosen geometrically or by a topological breadth-first selection. This lets the project compare a mesh-conforming refinement path with a voxel re-extraction path while keeping the outer mesh mostly unchanged.

It is a useful secondary comparison, but this thesis should not present it as a second primary contribution. It acts on the existing triangulation and its subdivision transitions; it does not independently reconstruct a local scalar grid in the way the proposed method does.

## 7. Proposed method: local voxel re-extraction plus stitching

The proposed contribution is `LocalMarchingCubes` followed by `LocalVoxelStitcher`: extract a denser local grid from the original implicit field around a detected singularity, then replace the corresponding region of the global MC mesh and match the cut loops. This gives an independently sampled local surface rather than just subdividing old triangles.

The assumptions must be stated honestly:

- The current local patch is centered at a point, and its interface is a sphere-cut seam.
- Stitching requires compatible boundary loops. A mismatch is recorded as a failure and the baseline global mesh is retained; the result must not be reported as a successful refinement in that case.
- `crossing-planes` has a line singularity, so a point-centered patch is not a model of the entire singular set. Its failure is a meaningful boundary case and indicates that a line-aware patch/seam strategy would be a separate extension.
- `tilted-cone` helps test grid orientation. A one-component mesh or a zero changed-selection count is a signal to inspect, not by itself proof that the continuous surface was misdefined.

## 8. Experimental design and claim discipline

Use at least two resolutions and multiple local radii. Include the surface families `cone`, `elliptic-cone`, `rotated-cone`, `x-cone`, `y-cone`, `tilted-cone`, `crossing-planes`, `quartic`, and the regular `sphere` control. Sweep both `c=0` and nonzero `c` values: the cone-family origin is singular only when the shifted level passes through the critical point at `c=0`.

For each method, report:

- global and local vertex-distance mean/RMS/p95/max;
- global and local triangle-sample distance mean/RMS/p95/max;
- projection attempts and failures for both sample sets;
- connected components, boundary edges, and non-manifold edges;
- stitch success, seam loops, and maximum matched seam-vertex distance, where applicable;
- selector face counts/change indicator for `LocalUnfolder`.

The generic distance in this project is a damped numerical projection to \(f=0\), not a guaranteed global closest-point distance. Its stopping criterion uses the estimated spatial Newton step rather than only \(|f|\), because a quadratic field can have a small residual at a point that is still measurably off-surface. Near a singular set the gradient can vanish, so projection failures must be shown and discussed. Keep the exact analytic double-cone metric as a separate validation for `cone, c=0`.

Do not state in advance that every isolated non-degenerate singularity works or that every degenerate one fails. Use the CSV to characterize the observed behavior and mention the range of resolutions, radii, levels, and parameter values tested. Visual inspection is still needed to diagnose geometry or connectivity failures that aggregate metrics cannot explain.

## References

1. W. E. Lorensen and H. E. Cline. “Marching Cubes: A High Resolution 3D Surface Construction Algorithm.” *Computer Graphics* 21(4), 163–169, 1987. DOI: https://doi.org/10.1145/37402.37422.
2. T. Ju, F. Losasso, S. Schaefer, and J. Warren. “Dual Contouring of Hermite Data.” *ACM Transactions on Graphics* 21(3), 339–346, 2002. DOI: https://doi.org/10.1145/566654.566586.
3. T. Lewiner, H. Lopes, A. W. Vieira, and G. Tavares. “Efficient Implementation of Marching Cubes' Cases with Topological Guarantees.” *Journal of Graphics Tools* 8(2), 1–15, 2003. DOI: https://doi.org/10.1080/10867651.2003.10487582.
4. G. M. Treece, R. W. Prager, and A. H. Gee. “Regularised marching tetrahedra: improved iso-surface extraction.” *Computers & Graphics* 23(4), 583–598, 1999. DOI: https://doi.org/10.1016/S0097-8493(99)00076-X.
5. Y. Zhang and J. Qian. “Dual Contouring for domains with topology ambiguity.” *Computer Methods in Applied Mechanics and Engineering* 217–220, 34–45, 2012. DOI: https://doi.org/10.1016/j.cma.2012.01.004.
6. T. Rashid, S. Sultana, and M. A. Audette. “2-manifold surface meshing using dual contouring with tetrahedral decomposition.” *Advances in Engineering Software* 102, 83–96, 2016. DOI: https://doi.org/10.1016/j.advengsoft.2016.09.008.
