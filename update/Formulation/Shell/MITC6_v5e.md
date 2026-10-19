# MITC6 v5e — CTRIA6 quadratic shell

**Ary · Independent Researcher · Personal Research** · MYSTRAN v18.00.a · 19 October 2026

> Research documentation. Full derivation and primary references: [master formulation](CQUAD8_CTRIA6_QUADRATIC_SHELL_FORMULATIONS.md). Exact Fortran source selection, sign conventions and validation branch must be checked before asserting implementation identity.

## Element-specific formulation (preserved from earlier research manuscript)

### 7.7 MITC6 v5e — enhanced six-node triangular mixed element

On the reference triangle, form the covariant membrane matrix \(B^n_m(r,s)\), and covariant transverse shear matrix \(B^n_s(r,s)\). A projection of a sampled scalar strain row \(b\) to an affine field is determined from \(\mathbf p=(1,r,s)^T\):
\[
\mathbf A_{j,:}=[1,r_j,s_j],\qquad \mathbf C=(A^TWA)^{-1}A^TW\begin{bmatrix}b(r_1,s_1)\\\vdots\\b(r_m,s_m)\end{bmatrix},\qquad \bar b(r,s)=\mathbf p^T\mathbf C.
\]
The enhanced tying configuration changes selected internal sampling locations relative to the preceding research version; the operator must be checked on the triangle itself and cannot be inferred from a quadrilateral MITC rule. Compute physical membrane and transverse shear tensors using the surface metric; bending and drilling are separately integrated. **Adopted:** Kim–Bathe triangular MITC motivation and strain tying. **Modified:** affine reconstruction, particular tying positions and six-DOF drilling. **Not reproduced:** a claim that the research v5e sample selection was published unchanged in Kim–Bathe (2009).

## Applicable verification and limitations

See [reference-anchored comparison](../../shell_benchmark/comparison_all.md), including 2-001 patch checks, 2-005 C/D convergence, 2-008 modal, 2-013 thermal stress, 2-014 thermal bending, 2-016 and 2-017 buckling. A solver pairing is not a reference accuracy proof. `2-009` foundation remains excluded.

- **Modal:** frequency reference is `1 / period`; compare ratio `f_solver/f_ref` and absolute frequency errors separately.
- **Thermal:** local radial–tangential axes and resultants must match CSI definitions.
- **Buckling:** preserve signed eigenvalues, prestress load basis and mode pairing; do not assume mode indices are interchangeable.
- **Solver caveats:** report AUTOSPC/BAILOUT conditions where present.


## Current checkout recovery note — 19 October 2026

Commit1f97ea99 makes SIMOQ8 same-point tensor recovery explicit while preserving the OP2 center basis, and corrects basic/local basis handling in quadratic GPFORCE. Stiffness and pressure kernels are unchanged; applied loads and displacements remain unchanged. This does not remove transverse-shear, curved-shell or buckling-normalization limitations. See [reference comparison](../../shell_benchmark/comparison_all.md).
