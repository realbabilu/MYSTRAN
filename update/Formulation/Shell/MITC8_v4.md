# MITC8 v4 — CQUAD8 quadratic shell

**Ary · Independent Researcher · Personal Research** · MYSTRAN v18.00.a · 19 October 2026

> Research documentation. Full derivation and primary references: [master formulation](CQUAD8_CTRIA6_QUADRATIC_SHELL_FORMULATIONS.md). Exact Fortran source selection, sign conventions and validation branch must be checked before asserting implementation identity.

## Element-specific formulation (preserved from earlier research manuscript)

### 7.3 MITC8 v4 — mixed covariant tensorial interpolation

The covariant in-layer strain is \(\varepsilon_{\alpha\beta}=\tfrac12(\mathbf g_\alpha\cdot\mathbf u_{,\beta}+\mathbf g_\beta\cdot\mathbf u_{,\alpha})\), with separate assumed membrane and transverse shear tensors. Curved surfaces use a metric-aware conversion into physical components rather than treating covariant strains as Cartesian strains. For membrane normal rows, define \(b=\sqrt{3/5}\), \(P_\pm(x)=(1\pm x)/2\), and \(Q_j(x)\), the quadratic cardinal functions at \((-b,0,b)\). Then
\[
\bar\varepsilon_{rr}(r,s)=\sum_{i\in\{-1,1\}}\sum_{j\in\{-b,0,b\}}P_i(r)Q_j(s)\varepsilon_{rr}(i,j),\quad
\bar\varepsilon_{ss}(r,s)=\sum_{i\in\{-1,1\}}\sum_{j\in\{-b,0,b\}}Q_j(r)P_i(s)\varepsilon_{ss}(j,i).
\]
For \(a=1/\sqrt3\), the mixed in-plane shear is \(\bar\varepsilon_{rs}=\sum_{i,j\in\{-a,a\}}L_i(r)L_j(s)\varepsilon_{rs}(i,j)\), with \(L_{\pm a}(x)=(a\pm x)/(2a)\). This membrane tying is activated on curved/warped geometry and otherwise the direct membrane interpolation is retained. The transverse shear operator is given by §8. Drilling stabilization and optional locally modified displacement fields are research adaptations. **Adopted:** MITC covariant component interpolation and physical strain transformation (Bathe–Dvorkin; Dvorkin). **Modified:** curvature-conditioned membrane tying, the slenderness-dependent shear blend and six-DOF drilling. **Not reproduced:** original five-DOF complete nonlinear tangent and unmodified historical MITC8 operators.

## Applicable verification and limitations

See [reference-anchored comparison](../../shell_benchmark/comparison_all.md), including 2-001 patch checks, 2-005 C/D convergence, 2-008 modal, 2-013 thermal stress, 2-014 thermal bending, 2-016 and 2-017 buckling. A solver pairing is not a reference accuracy proof. `2-009` foundation remains excluded.

- **Modal:** frequency reference is `1 / period`; compare ratio `f_solver/f_ref` and absolute frequency errors separately.
- **Thermal:** local radial–tangential axes and resultants must match CSI definitions.
- **Buckling:** preserve signed eigenvalues, prestress load basis and mode pairing; do not assume mode indices are interchangeable.
- **Solver caveats:** report AUTOSPC/BAILOUT conditions where present.


## Current checkout recovery note — 19 October 2026

Commit1f97ea99 makes SIMOQ8 same-point tensor recovery explicit while preserving the OP2 center basis, and corrects basic/local basis handling in quadratic GPFORCE. Stiffness and pressure kernels are unchanged; applied loads and displacements remain unchanged. This does not remove transverse-shear, curved-shell or buckling-normalization limitations. See [reference comparison](../../shell_benchmark/comparison_all.md).
