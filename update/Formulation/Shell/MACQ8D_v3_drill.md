# MACQ8D v3_drill — CQUAD8 quadratic shell

**Ary · Independent Researcher · Personal Research** · MYSTRAN v18.00.a · 19 October 2026

> Research documentation. Full derivation and primary references: [master formulation](CQUAD8_CTRIA6_QUADRATIC_SHELL_FORMULATIONS.md). Exact Fortran source selection, sign conventions and validation branch must be checked before asserting implementation identity.

## Element-specific formulation (preserved from earlier research manuscript)

### 7.6 MacNeal Q8 v3 — eight-node modified interpolation with drilling

A geometry-dependent interpolation \(\widetilde N_I\) replaces the standard Q8 field for the displacement/strain approximation while preserving the eight external nodes. The surface geometry is still reconstructed from the physical Q8 nodal positions. From the modified interpolation form \(\widetilde B_m,\widetilde B_b,\widetilde B_s\) and independent drilling penalty \(\widetilde B_d\), assemble
\[
K^{\mathrm{MN8}}=\int_A h\widetilde B_m^TC_m\widetilde B_m+\frac{h^3}{12}\widetilde B_b^TC_m\widetilde B_b+h\widetilde B_s^TC_s\widetilde B_s+h c_d\widetilde B_d^T\widetilde B_d\,dA.
\]
**Adopted:** MacNeal–Harder (1992) motivation for improved eight-node field interpolation. **Modified:** separate drilling constraint and six-DOF shell extension. **Not reproduced:** claim of equivalence to all original nine-node formulations or general nonlinear constitutive theory.

## Applicable verification and limitations

See [reference-anchored comparison](../../shell_benchmark/comparison_all.md), including 2-001 patch checks, 2-005 C/D convergence, 2-008 modal, 2-013 thermal stress, 2-014 thermal bending, 2-016 and 2-017 buckling. A solver pairing is not a reference accuracy proof. `2-009` foundation remains excluded.

### Mystran Benchmark: MACQ8D versus Nastran CQUAD8

The available paired runs use Nastran **CQUAD8** as the like-for-like reference solver. Results are reported against the independent case reference, with Nastran and MYSTRAN errors kept as separate quantities. The complete case tables and qualification notes are in the [reference-anchored comparison archive](../../shell_benchmark/comparison_all.md); no additional solver run is implied here.

- **Modal:** frequency reference is `1 / period`; compare ratio `f_solver/f_ref` and absolute frequency errors separately.
- **Thermal:** local radial–tangential axes and resultants must match CSI definitions.
- **Buckling:** preserve signed eigenvalues, prestress load basis and mode pairing; do not assume mode indices are interchangeable.
- **Solver caveats:** report AUTOSPC/BAILOUT conditions where present.


## Current checkout recovery note — 19 October 2026

Commit1f97ea99 makes SIMOQ8 same-point tensor recovery explicit while preserving the OP2 center basis, and corrects basic/local basis handling in quadratic GPFORCE. Stiffness and pressure kernels are unchanged; applied loads and displacements remain unchanged. This does not remove transverse-shear, curved-shell or buckling-normalization limitations. See [reference comparison](../../shell_benchmark/comparison_all.md).
