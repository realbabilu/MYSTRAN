# ANS8BDG6 v4 — CQUAD8 quadratic shell

**Ary · Independent Researcher · Personal Research** · MYSTRAN v18.00.a · 19 October 2026

> Research documentation. Full derivation and primary references: [master formulation](CQUAD8_CTRIA6_QUADRATIC_SHELL_FORMULATIONS.md). Exact Fortran source selection, sign conventions and validation branch must be checked before asserting implementation identity.

## Element-specific formulation (preserved from earlier research manuscript)

### 7.4 ANS8BDG6 v4 — enhanced ANS quadrilateral

The normal-strain components employ two longitudinal sides and three points per side (six per component). Thus \(\bar\varepsilon_{rr}\) and \(\bar\varepsilon_{ss}\) use the \(P_iQ_j\) interpolation above, while in-plane shear uses four \((\pm1/\sqrt3,\pm1/\sqrt3)\) samples and bilinear reconstruction. Curvature and midside deviation determine whether the membrane projection is active; a flat straight-sided element retains direct membrane strains. Transverse shear uses the alternative four- and six-point fields of §8.

The drilling variable is weakly constrained by the physical in-plane spin:
\[
\Pi_d=\frac12\int_A k_tGh\left(\theta_n-\tfrac12(v_{,x}-u_{,y})\right)^2dA,\quad k_t=0.1\ \text{(selected reference parameter)}.
\]
**Adopted:** β/δ sampling philosophy and enhanced transverse shear collocation from Jung–Han (2013), plus the drilling-spin energy idea. **Modified:** modified-field shape interpolation, conditional membrane tying, blended shear and six-DOF small-strain implementation. **Not reproduced:** a proven exact implementation of the paper's full \(\beta\delta\gamma^*_6\) interpolation or its nonlinear solution history.

## Applicable verification and limitations

See [reference-anchored comparison](../../shell_benchmark/comparison_all.md), including 2-001 patch checks, 2-005 C/D convergence, 2-008 modal, 2-013 thermal stress, 2-014 thermal bending, 2-016 and 2-017 buckling. A solver pairing is not a reference accuracy proof. `2-009` foundation remains excluded.

### Mystran Benchmark: ANS8BDG6 versus Nastran CQUAD8

The available paired runs use Nastran **CQUAD8** as the like-for-like reference solver. Results are reported against the independent case reference, with Nastran and MYSTRAN errors kept as separate quantities. The complete case tables and qualification notes are in the [reference-anchored comparison archive](../../shell_benchmark/comparison_all.md); no additional solver run is implied here.

- **Modal:** frequency reference is `1 / period`; compare ratio `f_solver/f_ref` and absolute frequency errors separately.
- **Thermal:** local radial–tangential axes and resultants must match CSI definitions.
- **Buckling:** preserve signed eigenvalues, prestress load basis and mode pairing; do not assume mode indices are interchangeable.
- **Solver caveats:** report AUTOSPC/BAILOUT conditions where present.


## Current checkout recovery note — 19 October 2026

Commit1f97ea99 makes SIMOQ8 same-point tensor recovery explicit while preserving the OP2 center basis, and corrects basic/local basis handling in quadratic GPFORCE. Stiffness and pressure kernels are unchanged; applied loads and displacements remain unchanged. This does not remove transverse-shear, curved-shell or buckling-normalization limitations. See [reference comparison](../../shell_benchmark/comparison_all.md).
