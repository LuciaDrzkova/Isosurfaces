# Method framing and literature notes

## Recommended thesis framing

**Primary contribution:** local voxel re-extraction with `LocalMarchingCubes`, followed by `LocalVoxelStitcher` integration into the global mesh. Evaluate whether this produces a geometrically closer and topologically acceptable representation near each singularity class.

**Classical baseline 1:** classic Marching Cubes on the complete uniform grid (`iso::mc::extract`). This is the non-adaptive reference.

**Classical baseline 2:** the included uniform-grid Dual Contouring implementation (`iso::dc::extract`). It places one QEF-derived vertex in each active cell and connects cells around sign-changing grid edges. This is a compact comparison implementation, not the full octree simplification method from the original paper, and it does not guarantee manifold topology in ambiguous cells.

**Secondary project-specific comparison:** `LocalUnfolder`, the conforming subdivision/projection pipeline. It is useful as a comparison but should not be described as the main contribution if the thesis claims the voxel patch method.

## Surface family

| CLI name | Equation at `c=0` | Singular set at `c=0` | Role |
|---|---|---|---|
| `cone` | `x^2+y^2-z^2=0` | isolated origin | circular, grid-aligned baseline |
| `elliptic-cone` | `x^2/2.25+y^2/0.5625-z^2=0` | isolated origin | anisotropic cone |
| `rotated-cone` | rotated elliptic cone | isolated origin | in-plane rotation |
| `x-cone` | `-x^2+y^2/2.25+z^2/0.5625=0` | isolated origin | axis swapped to X |
| `y-cone` | `x^2/2.25-y^2+z^2/0.5625=0` | isolated origin | axis swapped to Y |
| `tilted-cone` | rotated elliptic cone in 3D | isolated origin | off-grid principal axes |
| `crossing-planes` | `(x-y)(x+y)=0` | whole Z-axis | degenerate singular set |
| `quartic` | `x^2+y^2-z^2+0.05(x^4+y^4+z^4)=0` | isolated origin | non-quadratic control |
| `sphere` | `x^2+y^2+z^2-1=0` | none on the level set | regular control |

For `sphere`, `c` is ignored because its radius is a separate fixed constructor parameter. The sweep therefore writes one sphere experiment at `c=0` rather than duplicate rows for every `c` value.

## What can and cannot be concluded from the tool

The sweep records topology (connected components, boundary edges and non-manifold edges), generic surface-distance estimates, projection-failure counts, selector diagnostics, local patch statistics, and stitch success. It does not automatically prove that a singularity is faithfully represented. Conclusions should be derived from the generated CSV and inspection of representative meshes, across at least two resolutions and several local radii.

The generic geometric metric projects sample points toward `f=0` using damped Newton steps along the numerical gradient. Its convergence test uses the estimated spatial Newton step rather than only the scalar residual: for quadratic fields, `f` can be tiny while a point remains measurably off-surface. It is more reusable than a cone-specific closed form, but it is **not guaranteed to find the global Euclidean closest point**. Where the gradient is nearly zero or the iteration cannot converge, a failure is counted rather than assigned zero error. Report the failure rate alongside the error values. In particular, points near the singular Z-axis of `crossing-planes` require careful interpretation.

Do not write the expected classification as a result before measuring it. Use a table with at least: method, surface, c, resolution, local radius, refinement level, stitch success, vertex error RMS/p95, triangle-sample error RMS/p95, connected components, boundary edges, non-manifold edges, and projection failures. A failure case is part of the result if it is documented reproducibly.

## References for the survey chapter

1. W. E. Lorensen and H. E. Cline. “Marching Cubes: A High Resolution 3D Surface Construction Algorithm.” *Computer Graphics* 21(4), 163–169, 1987. DOI: https://doi.org/10.1145/37402.37422. Marching Cubes uses a cell case table and linearly interpolated grid-edge intersections. [Bibliographic record](https://cir.nii.ac.jp/crid/1361981470015626752?lang=en)
2. T. Ju, F. Losasso, S. Schaefer, and J. Warren. “Dual Contouring of Hermite Data.” *ACM Transactions on Graphics* 21(3), 339–346, 2002. DOI: https://doi.org/10.1145/566654.566586. The original method uses Hermite intersection/normal data and quadratic error functions (QEFs). [Bibliographic record](https://profiles.wustl.edu/en/publications/dual-contouring-of-hermite-data/)
3. T. Rashid, S. Sultana, and M. A. Audette. “2-manifold surface meshing using dual contouring with tetrahedral decomposition.” *Advances in Engineering Software* 102, 83–96, 2016. DOI: https://doi.org/10.1016/j.advengsoft.2016.09.008. This is a useful cautionary reference because ordinary dual contouring is not generally guaranteed to produce watertight two-manifold output. [Publisher record](https://www.sciencedirect.com/science/article/pii/S0965997816304537)

These references support the algorithm descriptions; a broader literature survey should also cover topology-correct Marching Cubes variants (e.g. MC33/Lewiner), adaptive contouring and signed-distance projection, with a discussion of why a second baseline is informative rather than assuming Dual Contouring necessarily wins on singularities.
