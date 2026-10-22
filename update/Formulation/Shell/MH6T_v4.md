# MH6T: Published Formulation versus Project Python v4 — MYSTRAN Linear Shell Research

**Author:** Ary  
**Date:** 23 October 2026  
**Target:** Published Li et al. (2015) versus project Python `MacNeal_MH6T_Tri_v4`
**Implementation baseline:** MH6T v3 is the project formulation previously ported to MYSTRAN; v4 is the latest supplied Python source for mathematical comparison. A v3 source/kernel was not supplied for direct auditing.  
**Python research reference:** `MacNeal_MH6T_Tri_v4` **(not automatically evidence of MYSTRAN v4 port)**  
**Source history:** `v1 → v2 → v3 (not supplied) → v4`  
**Evidence:** inspected uploaded Python v1, v2 and v4. Numerical claims must be kept separate from source-code evidence.

## Research scope: MYSTRAN linear SOL 101 first

**Priority:** implement and validate small-displacement **linear-elastic MYSTRAN SOL 101** with MH6T v4 as the baseline. The paper by Li et al. (2015) distinguishes assumed **linear membrane strain** ($\widetilde{\varepsilon}^L_m=P\alpha$, Eqs. 26, 34–38) from assumed **higher-order membrane strain** ($\widetilde{\varepsilon}^H_m=P\phi$, Eqs. 27, 44–46). The latter comes from a **geometrically nonlinear** strain contribution and is projected through a mixed variational condition. It is not needed to construct the small-displacement SOL 101 stiffness.

| Item in Li et al. (2015) | Current MYSTRAN SOL 101 decision |
|---|---|
| Nine line-integrated linear membrane constraints A–I (Eqs. 34–38) | **Use**, with documented v4 modifications |
| Six line-integrated transverse-shear constraints A–F (Eqs. 39–43) | **Use**, with documented v4 modifications |
| Linear bending and shear energy | **Use** in v4 elastic stiffness |
| Assumed higher-order membrane projection (Eqs. 27, 44–46) | **Exclude** from SOL 101; nonlinear research only |
| Co-rotational large rotations and elasto-plastic return mapping | **Out of scope** for SOL 101 |
| v5 projected higher-order membrane prototype | **Retain separately**, not the production linear kernel |

The v5 branch must **not** augment or replace `k_local()` in SOL 101. At the reference configuration its higher-order correction has zero linear tangent, so it is not expected to improve ordinary linear SAP2000/Nastran comparisons. **Linear eigen-buckling** based on initial-stress $K_G$ (SOL 105) is distinct from full nonlinear post-buckling and does not imply implementation of the paper's higher-order membrane field.

**Provenance:** Li et al. (2015), Section 3.2, Eqs. (26)–(28), (34)–(46), Appendix 1. DOI: https://doi.org/10.1007/s00466-015-1138-1.


## Paper-to-source verdict and scope boundary

The primary publication is **Li, Xiang, Izzuddin, Vu-Quoc, Zhuo & Zhang (2015), *A 6-node co-rotational triangular elasto-plastic shell element*, Computational Mechanics, DOI 10.1007/s00466-015-1138-1**. The published element has **five generalized nodal variables (30 element DOF)**; our research Python uses **six global DOF per node (36 element DOF)** and adds an artificial drilling constraint. They must not be described as identical elements. The paper cites MacNeal (1982) as the foundation for line-integrated assumed strains.

The exact adopted subset is Li et al.'s **linear membrane** affine field (paper Eqs. 26, 28, 34–38) and **transverse shear** affine field (Eqs. 31–32, 39–43), using the published sample coordinates in **Table 1 and Fig. 2**. These become 9×9 and 6×6 reconstruction systems in the Python source. By contrast, the **mixed higher-order membrane projection** (Eqs. 27, 44–46), co-rotational large rotations, material plasticity/return mapping, and nonlinear consistent tangent are **not used for MYSTRAN SOL 101**. Paper Eq. (46) is therefore not the source of the SOL 101 `Bm` operator.

**Differences worth highlighting against the paper:** Li Eq. (34) uses the chord vector divided by its full 3D squared norm for the member extension. Python v4 instead divides the full 3D chord numerator by the **squared *projected* centroid-plane chord length**, changing warped-element behavior. Li Eq. (39) uses a **subtriangle normal** `g_l` for the edge transverse-displacement term; Python v4 uses an **average endpoint director** `g_l=(n_i+n_j)/2`, with the project-specific rotation-vector formulation. These are **research modifications** with distinct kinematic invariance and accuracy questions. Li Eq. (47) supplies the standard plane stress/shear material operators; project code uses the isotropic specialized matrices and `k_shear=5/6`.

**Source audit hierarchy:** `v1` and `v2` files were inspected, `v3` is confirmed by project history as the previously ported MYSTRAN formulation but its exact source is absent here, and `v4` was inspected as the current Python research reference. Fortran/MYSTRAN parity and numerical SAP/Nastran benchmarks are delegated to Codex, not claimed here.

## Abstract

MH6T is a quadratic six-node Reissner–Mindlin-type shell research element using **line/edge-based assumed membrane and transverse shear strain reconstruction**, inspired by R. H. MacNeal's assumed-strain approach and the linear elastic operators presented by Li et al. (2015). A nine-coefficient membrane field is reconstructed from nine line constraints and a six-coefficient transverse shear field from six line constraints. Independent director-based bending and drilling operators are used. The latest supplied v4 extends the static kernel with consistent and HRZ-type diagonal mass, normal/directional pressure, gravity, thermal resultants and loading, and translational initial-stress geometric stiffness. These extensions are code-specific and must not be attributed to Li et al. as a verified reproduction. The claimed v3-to-v4 elastic-kernel identity is documented in the v4 header but cannot be independently proven without the absent v3 file.

## 1. Reference formulation: what was and was not adopted

| Feature | MacNeal (1982) / Li et al. (2015) | MH6T v4 implementation | Traceability |
|---|---|---|---|
| Assumed-strain principle | MacNeal (1982) | Yes | Conceptual foundation |
| Nine membrane edge samples A–I | Li et al., §3.2, Eqs. 26–38, Fig. 2a, Table 1, as identified in code | Yes, affine reconstruction | Adapted; **not** full paper algorithm |
| Six transverse shear samples A–F | Li et al., Eqs. 39–43, Fig. 2b | Yes, altered edge director operator | Modified |
| 6-node quadratic triangle | 6-node topology | Yes, P2 | Adopted |
| Co-rotational finite-rotation mechanics | Part of Li et al. | No | Outside linear implementation |
| Elasto-plasticity / return mapping | Part of Li et al. | No | Not implemented |
| Director-based bending | Not transcribed from Li's method | Yes | Derived from associated Simo-style shell approach, sign corrected |
| Drilling penalty | Code-specific | Yes, β = 0.02 by default | Added stabilization |
| 6 DOF per node | Compatible with solver interface | Yes, 36×36 elemental matrix | Implementation choice |
| Curved-shell handling | Li's flat sub-facet assumptions | Centroid frame + full 3D chord and averaged nodal director | Research modification; needs curved-mesh QA |
| Thermal, mass, pressure, geometric stiffness | Not source claims for the selected linear membrane/shear operators | Yes, v4 | Software extension, separately validated if numerical records exist |

**Important:** MH6T is **not Rezaiee T6/MITC6 point tying**. Its membrane and shear strains are reconstructed from edge/line constraints, not the Rezaiee in-plane/shear tying point interpolation. **The full published PDF `li2015.pdf` has now been supplied and inspected:** Sections 2–3, Fig. 2, Table 1, Eqs. (26)–(28), (34)–(47), and Appendix 1, Eqs. (55)–(62). The paper’s formulation and the Python approximation are not identical.

## 2. Element geometry, node order and interpolation

<img width="891" height="938" alt="nodes_t6" src="https://github.com/user-attachments/assets/ddab7954-d2d9-4031-8bd1-5ddedd1cdc80" />


Reference triangle: `r ≥ 0`, `s ≥ 0`, `r+s ≤ 1`. Nodes are 1=(0,0), 2=(1,0), 3=(0,1), 4=(1/2,0), 5=(1/2,1/2), 6=(0,1/2). In MYSTRAN/Nastran ordering, midnodes correspond to edges 1–2, 2–3 and 3–1.

Let `L1=1-r-s`, `L2=r`, and `L3=s`. The shape functions are

$$
\begin{aligned}
N_1&=L_1(2L_1-1),&N_2&=L_2(2L_2-1),&N_3&=L_3(2L_3-1),\\
N_4&=4L_1L_2,&N_5&=4L_2L_3,&N_6&=4L_3L_1.
\end{aligned}
$$

For nodal positions `x_I`, geometry and tangents are

$$\mathbf{x}(r,s)=\sum_{I=1}^{6}N_I\mathbf{x}_I,\qquad
\mathbf{g}_r=\sum_I N_{I,r}\mathbf{x}_I,\quad
\mathbf{g}_s=\sum_I N_{I,s}\mathbf{x}_I,\quad
J_A=\|\mathbf{g}_r\times\mathbf{g}_s\|.$$

Displacement order per node: `d_I = [u_x,u_y,u_z,θ_x,θ_y,θ_z]_I^T`; the assembled element vector has 36 entries. The code uses a centroid tangent plane `(E1,E2,E3)` for line reconstruction, and pointwise bases for bending/drilling calculations. A centroid-based orientation sign is selected once and applied to the nodal directors.

## 3. Nine-line assumed membrane strain reconstruction

<img width="891" height="938" alt="line_sampling_A_I" src="https://github.com/user-attachments/assets/5301bf61-97ca-4239-b185-4aa4fe862716" />


The 9 line groups are A=(1,4), B=(4,2), C=(2,5), D=(5,3), E=(3,6), F=(6,1), G=(6,4), H=(4,5), I=(5,6). Their associated sampling coordinates are A=(.25,0), B=(.75,0), C=(.75,.25), D=(.25,.75), E=(0,.75), F=(0,.25), G=(.25,.25), H=(.5,.25), I=(.25,.5).

Assume an affine membrane field in the fixed centroid-plane frame:

$$
\widehat{\boldsymbol\varepsilon}_m(r,s)=
\begin{bmatrix}1&r&s&0&0&0&0&0&0\\
0&0&0&1&r&s&0&0&0\\
0&0&0&0&0&0&1&r&s\end{bmatrix}
\underbrace{\begin{bmatrix}a_0&a_r&a_s&b_0&b_r&b_s&c_0&c_r&c_s\end{bmatrix}^{T}}_{\boldsymbol\alpha_m}.
$$

For line `(i,j)`, let `Δx_ij=x_j-x_i` (full 3D chord), `a_ij=(Δx·E1,Δx·E2)`, `ℓ_p²=a_ij·a_ij`, and directional cosines `(c_k,s_k)=a_ij/ℓ_p`. The **v4 modified** line-extension operator is

$$
\widehat e_k = \frac{\Delta\mathbf{x}_{ij}\cdot(\mathbf{u}_j-\mathbf{u}_i)}{\ell_p^2} = \mathbf{b}_k\mathbf{d}_e.
$$

It uses the *full 3D chord* in the numerator but the *projected squared length* in the denominator. This is not interchangeable with the earlier v1/v2 projected-chord numerator when the surface is curved. It is also not guaranteed to reproduce all curved-shell rigid modes without testing.

Directional interpolation at sample `(r_k,s_k)` gives

$$
\widehat e_k = [c_k^2,r_kc_k^2,s_kc_k^2,s_k^2,r_ks_k^2,s_ks_k^2,c_ks_k,r_kc_ks_k,s_kc_ks_k]\,\boldsymbol\alpha_m.
$$

Stack the nine constraints:

$$\boldsymbol\Gamma\boldsymbol\alpha_m=\mathbf{B}_{\rm line}\mathbf{d}_e,\quad
\boldsymbol\alpha_m=\boldsymbol\Gamma^{-1}\mathbf{B}_{\rm line}\mathbf{d}_e,\quad
\boxed{\mathbf B_m(r,s)=\mathbf P_m(r,s)\boldsymbol\Gamma^{-1}\mathbf B_{\rm line}}.$$

The source solves this linear system numerically (`np.linalg.solve`); it does not use explicit symbolic inversion. A singular or poorly conditioned `Γ` indicates a problematic mesh and should be reported rather than silently regularized.

## 4. Six-line assumed transverse shear reconstruction

The six sampling edges A–F follow the triangle boundary. Introduce

$$\widehat{\boldsymbol\gamma}_s(r,s)=\begin{bmatrix}1&r&s&0&0&0\\0&0&0&1&r&s\end{bmatrix}\boldsymbol\beta_s.$$

For edge `(i,j)`, let `L=||x_j-x_i||`, `\bar n_{ij}=(n_i+n_j)/2`, and `a_3=x_j-x_i`. The actual v4 code applies

$$
\widehat\gamma_k=\frac{\bar{\mathbf n}_{ij}\cdot(\mathbf u_j-\mathbf u_i)}{L}
+\frac{(\mathbf n_i\times\mathbf a_3)\cdot\boldsymbol\theta_i+(\mathbf n_j\times\mathbf a_3)\cdot\boldsymbol\theta_j}{2L}.
$$

Directional shear reconstruction uses the `6×6` matrix `Ω`, whose row at sample `(r_k,s_k)` is

$$\boldsymbol\Omega_k=[c_k,r_kc_k,s_kc_k,s_k,r_ks_k,s_ks_k],\quad c_k=\frac{\mathbf a_3\cdot\mathbf E_1}{L},\quad s_k=\frac{\mathbf a_3\cdot\mathbf E_2}{L}.$$

Thus

$$\boxed{\mathbf B_s(r,s)=\mathbf Q_s(r,s)\boldsymbol\Omega^{-1}\mathbf G_{\rm line}}.$$

**Revision distinction:** v1/v2 used subfacet normal `g_l` in the translation term, following their reading of Li et al.; v4 uses the *average nodal director* in both translational and rotational contributions to improve infinitesimal rigid-rotation compatibility on curved meshes. That is an MH6T research change, not a verbatim Li operator.

## 5. Bending, drilling and constitutive matrices

The bending operator `Bb` is constructed through director gradients. For node `I`, define `n_{,r}=Σ N_{I,r} n_I`, `n_{,s}=Σ N_{I,s} n_I`. The v4 rotational terms employ `g_r×n_I`, `g_s×n_I`, while translation terms employ the *opposite* director-gradient sign relative to v1/v2. For example, the covariant `rr` contributions are

$$\delta\rho_{rr}^{\rm rot}=N_{I,r}(\mathbf g_r\times\mathbf n_I)\cdot\delta\boldsymbol\theta_I,\qquad
\delta\rho_{rr}^{\rm trans}=-N_{I,r}\mathbf n_{,r}\cdot\delta\mathbf u_I.$$

The `ss` and `rs` terms are constructed analogously, then transformed into physical `xx,yy,xy` components by the pointwise contravariant/Cartesian mapping in the code. This paired sign change was introduced to cancel spurious bending under curved-shell infinitesimal rigid rotations.

The drilling operator is a penalty stabilization, not an independent constitutive shell strain:

$$\gamma_d=\mathbf B_d\mathbf d_e,\quad
\mathbf B_{d,I}^{(u)}=\tfrac12(N_{I,x}\mathbf e_2-N_{I,y}\mathbf e_1)^T,\quad
\mathbf B_{d,I}^{(\theta)}=-N_I\mathbf e_3^T.$$

The code projects geometric derivative coefficients to physical directions before assembling drilling; default drilling coefficient is `β_drill = 0.02`.

For isotropic elastic material, let `C = E/(1-ν²)[[1,ν,0],[ν,1,0],[0,0,(1-ν)/2]]`, `G=E/(2(1+ν))`, `C_s=(5/6)G I₂`, and `C_d=β_drill G`. A single average nodal thickness `\bar h` is used in the stiffness terms. The element stiffness is

$$
\boxed{K_e=\int_A\left[\bar h B_m^TCB_m+\frac{\bar h^3}{12}B_b^TCB_b+\bar h B_s^TC_sB_s+\bar h B_d^TC_dB_d\right]dA.}
$$

The four contributions are symmetric, provided the reconstruction systems are nonsingular. This is a 36×36 element matrix, assembled with the global six-DOF ordering (`T_matrix = I`).

## 6. Gaussian integration

| Integrated quantity | Active quadrature | Sum of reference-triangle weights |
|---|---|---:|
| Elastic membrane / bending / shear | 3 points | 1/2 |
| Drilling stiffness | 6 points | 1/2 |
| Mass (consistent or HRZ-lumped precursor) | 6 points | 1/2 |
| Pressure / traction / gravity, default `order=5` | 7 points | 1/2 |
| Thermal equivalent nodal forces | 3 points | 1/2 |
| Initial-stress stiffness, default `order=4` | 6 points | 1/2 |

<img width="1045" height="938" alt="gauss3" src="https://github.com/user-attachments/assets/8d623709-c0ee-4206-8840-867e47b0ce7a" />

<img width="891" height="938" alt="gauss6" src="https://github.com/user-attachments/assets/aab8557f-19a7-48d6-bcb3-4004239a0c14" />

<img width="912" height="938" alt="gauss7" src="https://github.com/user-attachments/assets/ea335300-a7d9-46b8-8e27-dbf38fcd8625" />

The seven-point option is available for certain v4 load/integration helpers;. The Gaussian integration points are different from the A–I line-sampling positions.

## 7. v4 multiphysics and stability operators

### 7.1 Consistent and HRZ diagonal mass

$$M_{IJ}^{tt}=\int_A\rho\bar h N_IN_J\,dA\,I_3,\qquad M_{IJ}^{rr}=\int_A\rho\frac{\bar h^3}{12}N_IN_J\,dA\,I_3.$$

The matrix retains all three rotational coordinates, including drilling rotation, as an implementation choice. The `consistent=False` branch diagonalizes each six-DOF component block by scaling its original diagonal to the block's algebraic total (`HRZ-type`). It should not be described as equal-mass nodal sharing.

### 7.2 Surface pressure, traction and gravity

$$\mathbf f_p=\int_A \mathbf N_u^T\mathbf t\,dA.$$

Default pressure acts along the stabilized `+e3` normal; a supplied `direction` applies a global directional traction per actual area. Gravity is implemented as `t=ρ \bar h g`. For pressure, scalar `p` and direction conventions must be recorded when comparing with Nastran `PLOAD4`.

### 7.3 Thermal

The implementation uses `ε_th = α ΔT [1,1,0]^T`, `κ_th = -α(dT/dz)[1,1,0]^T`, and equivalent force

$$\mathbf f_{th}=\int_A\left[B_m^T\bar h C\boldsymbol\varepsilon_{th}+B_b^T\frac{\bar h^3}{12}C\boldsymbol\kappa_{th}\right]dA.$$

Thermal force, membrane resultant and bending-moment recovery are implemented. The gradient direction and bending sign conventions must match the solver and the physical fiber side before quantitative comparison.

### 7.4 Initial-stress geometric stiffness

For membrane resultants `N0=[Nxx,Nyy,Nxy]` and tangent-plane shape gradients `∇N_I`,

$$K_{G,IJ}^{(u_a,u_a)}=\int_A (\nabla N_I)^T\begin{bmatrix}N_{xx}&N_{xy}\\N_{xy}&N_{yy}\end{bmatrix}(\nabla N_J)dA,\qquad a=x,y,z.$$

Only the **three translational blocks** are populated; no rotational initial-stress terms or full nonlinear consistent tangent are included. Eigen-buckling convention in source: `(K+λ K_G) φ = 0`, prestress tension positive. This is **linear bifurcation buckling**, not a postbuckling or large-rotation formulation.

## 8. Revision traceability: v1 → v2 → v4

| Topic | v1 | v2 | v4 (latest supplied) |
|---|---|---|---|
| Main assumed strain | 9 membrane + 6 shear line constraints | Same basis/operators | Same basis and sampling indices |
| Orientation | +Z tie-break at each point | Likewise | One centroid sign shared across directors and point bases |
| Membrane edge RHS | Projected tangent-plane chord | Projected tangent-plane chord | Full 3D chord numerator with projected squared-length denominator |
| Shear translation normal | Individual subfacet normal `g_l` | Individual subfacet normal `g_l` | Average endpoint director `(n_i+n_j)/2` |
| Shear rotation | Nodal directors | Nodal directors | Nodal directors, compatible with v4 translation operator |
| Bending rotation | `g×n` (corrected) | `g×n` | `g×n` |
| Bending translation | Positive director-gradient contribution | Positive | **Negative**, paired for curved rigid rotations |
| Drilling basis | Frozen centroid physical frame | Frozen centroid physical frame | Pointwise physical frame |
| Elastic Gauss | 3+6 | 3+6 | 3+6, unchanged |
| Mass | Consistent-only function present | Consistent-only function present | Consistent or HRZ diagonal |
| Pressure | Scalar normal force | Scalar normal force | Directional traction and configurable 3/6/7-point quadrature |
| Thermal | Not present | Not present | Equivalent load and resultant recovery |
| Geometric stiffness | Not present | Not present | Translational initial-stress `Kg` |

**MYSTRAN-port caveat:** MH6T **v3** was ported to MYSTRAN earlier in this project. The v4 file header claims that `k_local()` is unchanged from v3 and that v4 adds mass/thermal/load/buckling support. Because v3 Python and current MYSTRAN Fortran source were not supplied in this turn, identical stiffness/outputs across v3 and v4 cannot be independently guaranteed; Codex must verify directly. No reconstructed v3 source is substituted.

## 9. Numerical validation evidence and limits

This document audits **the supplied code**, not a fresh run of MSC Nastran or MYSTRAN. The v1/v2 source comments report historical 2-001 membrane/bending patch fixes and a twisted-beam warped-mesh discrepancy, but no full result provenance or corresponding OP2 data accompanies the three Python files. Therefore they are historical engineering notes, **not independently reproduced validation results**.

A previous project benchmark noted **ported MH6T v3** `T9` error approximately `+8.7818%`, while a separate 2-004 case showed material mesh sensitivity. These are contextual historical results, **not v4 accuracy claims**. The previously prepared generic CTRIA6 comparison table has blank entries for MH6T v4; those must not be replaced by Rezaiee or aggregate T6 results.

| Test group | Metric / check | v4 vs Nastran status |
|---|---|---|
| SAP 2-001 | Membrane and bending analytical patch / winding | Historical patch evidence only; v4 per-solver CSV pending |
| SAP 2-002–2-004 | Bending, curved response, twisted beam | Paired v4 CSV/OP2 needed |
| SAP 2-005–2-007 | Plate pressure, Scordelis–Lo, hemisphere | Per-model matching meshes and outputs needed |
| SAP 2-008 | Modes, mass consistency and MAC | Paired modal results needed |
| SAP 2-013 | Thermal resultants and temperature gradient | Paired stress/resultant outputs needed |
| SAP 2-016–2-017 | Signed buckling eigenvalues | Paired initial-stress and mode verification needed |

**Comparison protocol:** use physically identical coordinates, connectivity, shell thickness, director/winding sign, PSHELL `MID3`, constitutive/shear factor, SPCs, pressure orientation and mode/prestress conditions. Compare displacements in common global XYZ; for `GPSTRESS`, `GPSTRAIN`, membrane moments and transverse shear, convert element-local results to the same physical frame and fiber before computing percentages. Nastran stress/force data at absent sampling stations must remain missing, not interpolated and mislabeled as native OP2.

## 10. Conclusions and open technical questions

1. MH6T v4's mathematical identity is **line-based assumed membrane/shear strain**, with a 9×9 membrane and 6×6 shear reconstruction, rather than standard MITC6 point tying.
2. The major v4 change affecting elastic stiffness relative to v1/v2 is the treatment of curved geometry: full 3D chord contribution, average director shear translation, consistent normal-sign logic, corrected bending translation sign, and pointwise drilling basis.
3. Mass, pressure, thermal, and buckling operators are **implementations**, not automatically validated merely because methods exist in the Python class.
4. Particular follow-up checks should target element rigid-body rotation residuals under warped geometry, rank/condition of `Γ` and `Ω`, distorted patch consistency, curved-shell stress recovery, and mesh refinement in twisted-beam, Scordelis–Lo and signed eigen-buckling cases.

## References

1. **MacNeal, R. H. (1982).** Derivation of element stiffness matrices by assumed strain distributions. *Nuclear Engineering and Design*, 70(1), 3–12. https://doi.org/10.1016/0029-5493(82)90262-X.
2. **Li, Z., Xiang, Y., Izzuddin, B. A., Vu-Quoc, L., Zhuo, X., & Zhang, C. (2015).** A 6-node co-rotational triangular elasto-plastic shell element. *Computational Mechanics*, 55(5), 837–859. https://doi.org/10.1007/s00466-015-1138-1.
3. **Simo, J. C., & Fox, D. D. (1989).** On a stress resultant geometrically exact shell model. Part I: Formulation and optimal parametrization. *Computer Methods in Applied Mechanics and Engineering*, 72, 267–304. Used only for broad director-shell context, not as a verbatim source of MH6T's assumed-strain operator.
4. **Kim, D.-N., & Bathe, K.-J. (2009).** A triangular six-node shell element. *Computers & Structures*, 87, 1451–1460. Comparator for **MITC6**, not MH6T's line reconstruction. https://doi.org/10.1016/j.compstruc.2009.05.002.
5. **Computers and Structures, Inc.** *SAP2000 Verification Examples: Shells*, particularly Problems 2-001, 2-002, 2-004, 2-005, 2-006, 2-007, 2-008, 2-013, 2-016 and 2-017. Reference documents: `https://docs.csiamerica.com/manuals/sap2000/Verification/Analysis/Shells/Problem%202-001.pdf` (replace suffix by relevant ID; verify each published case and observable before claiming independent-reference correspondence). SAP2000 is developed by Computers and Structures, Inc. (CSI): https://www.csiamerica.com/products/sap2000.



