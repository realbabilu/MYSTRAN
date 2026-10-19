# SIMOT6 v2 — CTRIA6 quadratic shell

**Ary · Independent Researcher · Personal Research** · MYSTRAN v18.00.a · 19 October 2026

> Research documentation. Full derivation and primary references: [master formulation](CQUAD8_CTRIA6_QUADRATIC_SHELL_FORMULATIONS.md). Exact Fortran source selection, sign conventions and validation branch must be checked before asserting implementation identity.

## Element-specific formulation (preserved from earlier research manuscript)

### 7.2 SIMOT6 v2 — six-node triangular shell

Using the P2 triangle basis \(N_I(L_1,L_2,L_3)\) of §3.2 and \(\Omega_\triangle=\{r\ge0,s\ge0,r+s\le1\}\), construct directors independently for six nodes and the covariant surface gradient \(\mathbf g_r,\mathbf g_s\). The strain–rotation operators are evaluated in the moving physical tangent basis and the triangular stiffness is
\[
 K^{\mathrm{SIMOT6}}=\int_{\Omega_\triangle}\left[hB_m^TC_mB_m+\frac{h^3}{12}B_b^TC_mB_b+hB_s^TC_sB_s+h c_dB_d^TB_d\right]J_A\,dr\,ds.
\]
Where projected rows are selected, they must use a *triangular* polynomial space and triangle sample positions; Q8 tensor-product tying is not transferable. **Adopted:** stress-resultant director theory and objective frame requirements. **Modified:** P2 triangle interpolation, six-DOF representation and drilling constraint. **Not reproduced:** exact finite-rotation geometric Simo tangent. SIMOT6 is an independent formulation, not a corner-reduced SIMOQ8.

## Applicable verification and limitations

See [reference-anchored comparison](../../shell_benchmark/comparison_all.md), including 2-001 patch checks, 2-005 C/D convergence, 2-008 modal, 2-013 thermal stress, 2-014 thermal bending, 2-016 and 2-017 buckling. A solver pairing is not a reference accuracy proof. `2-009` foundation remains excluded.

- **Modal:** frequency reference is `1 / period`; compare ratio `f_solver/f_ref` and absolute frequency errors separately.
- **Thermal:** local radial–tangential axes and resultants must match CSI definitions.
- **Buckling:** preserve signed eigenvalues, prestress load basis and mode pairing; do not assume mode indices are interchangeable.
- **Solver caveats:** report AUTOSPC/BAILOUT conditions where present.


## Current checkout recovery note — 19 October 2026

Commit1f97ea99 makes SIMOQ8 same-point tensor recovery explicit while preserving the OP2 center basis, and corrects basic/local basis handling in quadratic GPFORCE. Stiffness and pressure kernels are unchanged; applied loads and displacements remain unchanged. This does not remove transverse-shear, curved-shell or buckling-normalization limitations. See [reference comparison](../../shell_benchmark/comparison_all.md).
