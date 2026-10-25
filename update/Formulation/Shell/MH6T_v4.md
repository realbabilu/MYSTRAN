# MH6T v4 — CTRIA6 quadratic shell

**Ary · Independent Researcher · Personal Research** · MYSTRAN v18.00.a · 19 October 2026

> Research documentation. Full derivation and primary references: [master formulation](CQUAD8_CTRIA6_QUADRATIC_SHELL_FORMULATIONS.md). Exact Fortran source selection, sign conventions and validation branch must be checked before asserting implementation identity.

## Element-specific formulation (preserved from earlier research manuscript)

### 7.8 MH6T v4 — triangular modified-field MacNeal family

Using P2 geometry, form direct bending from interpolated rotations and reconstruct in-layer/shear strains by a weighted projection over separately specified triangular sample points:
\[
\bar B_\alpha(r,s)=\sum_{j=1}^{m}H_j(r,s)B^{n}_\alpha(r_j,s_j),\qquad K=\int_{\Omega_\triangle}(h\bar B_m^TC_m\bar B_m+\tfrac{h^3}{12}B_b^TC_mB_b+h\bar B_s^TC_s\bar B_s+h c_dB_d^TB_d)J_A\,dr\,ds.
\]
**Adopted:** modified/assumed-strain treatment of thin-shell locking. **Modified:** triangular tying and six-DOF drilling. **Not reproduced:** equivalence of this research interpolation to any single unmodified published MH6T element.

## Applicable verification and limitations

See [reference-anchored comparison](../../shell_benchmark/comparison_all.md), including 2-001 patch checks, 2-005 C/D convergence, 2-008 modal, 2-013 thermal stress, 2-014 thermal bending, 2-016 and 2-017 buckling. A solver pairing is not a reference accuracy proof. `2-009` foundation remains excluded.

### Mystran Benchmark: MH6T versus Nastran CQUAD8

The available paired runs use Nastran **CQUAD8** as the like-for-like reference solver. Results are reported against the independent case reference, with Nastran and MYSTRAN errors kept as separate quantities. The complete case tables and qualification notes are in the [reference-anchored comparison archive](../../shell_benchmark/comparison_all.md); no additional solver run is implied here.

- **Modal:** frequency reference is `1 / period`; compare ratio `f_solver/f_ref` and absolute frequency errors separately.
- **Thermal:** local radial–tangential axes and resultants must match CSI definitions.
- **Buckling:** preserve signed eigenvalues, prestress load basis and mode pairing; do not assume mode indices are interchangeable.
- **Solver caveats:** report AUTOSPC/BAILOUT conditions where present.


## Current checkout recovery note — 19 October 2026

Commit1f97ea99 makes SIMOQ8 same-point tensor recovery explicit while preserving the OP2 center basis, and corrects basic/local basis handling in quadratic GPFORCE. Stiffness and pressure kernels are unchanged; applied loads and displacements remain unchanged. This does not remove transverse-shear, curved-shell or buckling-normalization limitations. See [reference comparison](../../shell_benchmark/comparison_all.md).
