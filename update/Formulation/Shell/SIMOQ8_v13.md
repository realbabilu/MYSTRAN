# SIMOQ8 v13 — CQUAD8 quadratic shell

**Ary · Independent Researcher · Personal Research** · MYSTRAN v18.00.a · 19 October 2026

> Research documentation. Full derivation and primary references: [master formulation](CQUAD8_CTRIA6_QUADRATIC_SHELL_FORMULATIONS.md). Exact Fortran source selection, sign conventions and validation branch must be checked before asserting implementation identity.

## Element-specific formulation (preserved from earlier research manuscript)

### 7.1 SIMOQ8 v13 — eight-node serendipity shell

The surface uses the Q8 serendipity functions of §3.1. Nodal directors \(\mathbf d_I\) define the rotation-induced through-thickness displacement \(\mathbf u(z)=\sum_I N_I[\mathbf u_I+z(\boldsymbol\theta_I\times\mathbf d_I)]\). The physical strain operators follow the covariant surface map \(\mathbf g_\alpha=\mathbf x_{,\alpha}\), with tensor transformation to an orthonormal tangent frame. The material stiffness is
\[
 K^{\mathrm{SIMOQ8}}=\int_A\!\left[h\bar B_m^TC_m\bar B_m+\frac{h^3}{12}B_b^TC_mB_b+h\bar B_s^TC_s\bar B_s+h\,c_dB_d^TB_d\right]dA.
\]
An optional assumed-strain branch replaces directly computed membrane or shear strain rows by projected rows \(\bar B_\alpha=\sum_j L_j^{(\alpha)} B_\alpha(\xi_j)\), whereas bending and drilling are treated separately. An optional incompatible/enhanced strain contributes an internal field \(\varepsilon^{E}=G\alpha\), with the condensed stiffness \(K_{qq}^{\star}=K_{qq}-K_{q\alpha}K_{\alpha\alpha}^{-1}K_{\alpha q}\) **only when that branch is active and nonsingular**. This is not a claim that the entire nonlinear Simo–Fox–Rifai algorithm is reproduced. **Adopted:** director kinematics, stress-resultant decomposition, local covariance. **Modified:** Q8 shape basis, conditional assumed/enhanced strains and drilling penalty. **Not reproduced:** full geometrically exact time integration, finite rotational tangent, large-strain plasticity and thickness stretching.

## Applicable verification and limitations

See [reference-anchored comparison](../../shell_benchmark/comparison_all.md), including 2-001 patch checks, 2-005 C/D convergence, 2-008 modal, 2-013 thermal stress, 2-014 thermal bending, 2-016 and 2-017 buckling. A solver pairing is not a reference accuracy proof. `2-009` foundation remains excluded.

- **Modal:** frequency reference is `1 / period`; compare ratio `f_solver/f_ref` and absolute frequency errors separately.
- **Thermal:** local radial–tangential axes and resultants must match CSI definitions.
- **Buckling:** preserve signed eigenvalues, prestress load basis and mode pairing; do not assume mode indices are interchangeable.
- **Solver caveats:** report AUTOSPC/BAILOUT conditions where present.


## Current checkout recovery note — 19 October 2026

Commit1f97ea99 makes SIMOQ8 same-point tensor recovery explicit while preserving the OP2 center basis, and corrects basic/local basis handling in quadratic GPFORCE. Stiffness and pressure kernels are unchanged; applied loads and displacements remain unchanged. This does not remove transverse-shear, curved-shell or buckling-normalization limitations. See [reference comparison](../../shell_benchmark/comparison_all.md).
