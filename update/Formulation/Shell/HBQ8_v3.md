# HBQ8 v3 — CQUAD8 quadratic shell

**Ary · Independent Researcher · Personal Research** · MYSTRAN v18.00.a · 19 October 2026

> Research documentation. Full derivation and primary references: [master formulation](CQUAD8_CTRIA6_QUADRATIC_SHELL_FORMULATIONS.md). Exact Fortran source selection, sign conventions and validation branch must be checked before asserting implementation identity.

## Element-specific formulation (preserved from earlier research manuscript)

### 7.5 HBQ8 Kikuchi v3 — modified-field hybrid quadrilateral

The adopted field interpolation \(\widetilde N_I(r,s;\mathbf x)\) is geometry-sensitive (Kikuchi-type variant), giving physical surface gradients \(\nabla_\Gamma\widetilde N_I\); however surface area is integrated from the standard Q8 geometric map. Membrane, bending and drilling energies are evaluated from separate differential rows, with the selected shear assumed field from §8. Its total stiffness retains the four mechanical energies of §5 and the drilling energy of §7.4, with its own drilling coefficient. The geometry switch is
\[
\chi_c=\begin{cases}1,&\theta_{\max}>0.01\ \mathrm{rad}\quad\text{or}\quad\delta_{mid}>10^{-3},\\0,&\text{otherwise},\end{cases}
\quad\delta_{mid}=\max_e\frac{\left\|\mathbf x_{m_e}-\tfrac12(\mathbf x_{a_e}+\mathbf x_{b_e})\right\|}{\|\mathbf x_{b_e}-\mathbf x_{a_e}\|}.
\]
Here \(\theta_{\max}=\max_I\arccos(\mathbf d_I\cdot\bar{\mathbf d})\) with normalized average director \(\bar{\mathbf d}\). For \(\chi_c=0\), direct membrane field; for \(\chi_c=1\), covariantly tied membrane field. This is a **formulation criterion**, not a programming routine. **Adopted:** HBQ8 shell-energy decomposition and Kikuchi-type modified interpolation. **Modified:** switching and shear blending. **Not reproduced:** full geometrically nonlinear shell tangent.

## Applicable verification and limitations

See [reference-anchored comparison](../../shell_benchmark/comparison_all.md), including 2-001 patch checks, 2-005 C/D convergence, 2-008 modal, 2-013 thermal stress, 2-014 thermal bending, 2-016 and 2-017 buckling. A solver pairing is not a reference accuracy proof. `2-009` foundation remains excluded.

### Mystran Benchmark: HBQ8 versus Nastran CQUAD8

The available paired runs use Nastran **CQUAD8** as the like-for-like reference solver. Results are reported against the independent case reference, with Nastran and MYSTRAN errors kept as separate quantities. The complete case tables and qualification notes are in the [reference-anchored comparison archive](../../shell_benchmark/comparison_all.md); no additional solver run is implied here.

- **Modal:** frequency reference is `1 / period`; compare ratio `f_solver/f_ref` and absolute frequency errors separately.
- **Thermal:** local radial–tangential axes and resultants must match CSI definitions.
- **Buckling:** preserve signed eigenvalues, prestress load basis and mode pairing; do not assume mode indices are interchangeable.
- **Solver caveats:** report AUTOSPC/BAILOUT conditions where present.


## Current checkout recovery note — 19 October 2026

Commit1f97ea99 makes SIMOQ8 same-point tensor recovery explicit while preserving the OP2 center basis, and corrects basic/local basis handling in quadratic GPFORCE. Stiffness and pressure kernels are unchanged; applied loads and displacements remain unchanged. This does not remove transverse-shear, curved-shell or buckling-normalization limitations. See [reference comparison](../../shell_benchmark/comparison_all.md).
