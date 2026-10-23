# Simo T6 Shell Element: CTRIA6 Shell 

**A technical research formulation and source-traceability report for MYSTRAN quadratic CTRIA6**  
**Author:** Ary  
**Date:** 23 October 2026  
**Formulation:** SIMOT6 v2; six-node quadratic one-director shell with conforming membrane and shear strains.  
**Numerical evidence:** existing formulation-specific MYSTRAN/Nastran CSV, compared separately against independent references; no new solver runs.  

## Abstract

This document describes the project's six-node, quadratic, linear-elastic triangular shell (SimoT6), its geometric interpolation and membrane/bending/transverse-shear/drilling operators, integration rules, auxiliary mass/pressure/thermal/geometric-stiffness routines, and the correspondence with the Simo–Fox–Rifai shell literature. It carefully separates **published shell-theory concepts** from **project-specific discretization and stabilization**, and presents existing formulation-specific benchmark results against independent references.

## 1. Formulation provenance and scope

The relevant published lineage is the geometrically exact one-director/stress-resultant shell formulation of Simo and Fox and its linear computational treatment by Simo, Fox and Rifai. **SIMOT6 v2 is a project-specific discretization, not a verbatim published 1993 six-node Simo element.** The SimoT6 name identifies the project formulation. The element employs standard isoparametric P2 displacement/geometry interpolation and directly evaluated conforming strains, rather than reproducing every mixed Hellinger–Reissner or nonlinear director-field treatment of the research series.

| Feature | Published Simo/Fox/Rifai lineage | SIMOT6 v2 formulation | Status |
|---|---|---|---|
| One-director shell kinematics | Fundamental continuum framework | Interpolated nodal directors; linearized rotation coupling | Adapted |
| Membrane strain | Geometric shell strain; mixed methods in published FE | Direct pointwise covariant `Bm` using P2 geometry | Adapted; not complete mixed scheme |
| Bending strain | Director/curvature mechanics | Director-gradient and rotational cross-product blocks | Adapted |
| Transverse shear | Shell shear-resultant kinematics; mixed interpolation also studied | Direct conforming `Bs`, **no MITC tying** | Adapted; locking risk |
| DOFs | Intrinsic tangent director rotations, 5 essential per node; six-parameter representations discussed | 6 global DOF/node = 36; drilling penalty | Project adaptation |
| Finite rotations / nonlinear response | Geometrically exact, nonlinear | No finite-rotation update or consistent nonlinear tangent | Not implemented |
| Drilling stabilization | Depends on discrete treatment | `beta_drill=0.02` penalty and 6-point integration | Project-specific |
| Membrane/bending/shear quadrature | Depends on published discrete element | 3 triangle points | Project-specific |
| Thermal, pressure, eigen, buckling | Auxiliary problems, not an equivalence criterion | Helper methods exist in source | Validation separate |

### 1.1 Nodal ordering and kinematics

<img width="1080" height="900" alt="simo_t6_nodes" src="https://github.com/user-attachments/assets/1142e1d2-7bc3-4188-9a7a-02d2bad62d51" />


Reference domain $\hat\Omega=\{(r,s):r\geq0,s\geq0,r+s\leq1\}$, with area $1/2$. Coordinates are $(0,0),(1,0),(0,1),(1/2,0),(1/2,1/2),(0,1/2)$ for nodes 1 to 6. Each node carries $\mathbf q_i=[u_x,u_y,u_z,\theta_x,\theta_y,\theta_z]^T$, so $\mathbf q_e\in\mathbb R^{36}$ with global nodal coordinates and an identity element-level coordinate mapping.

Define barycentric coordinates $L_1=1-r-s$, $L_2=r$, $L_3=s$. The six quadratic shape functions are

$$
N_1=L_1(2L_1-1),\quad N_2=L_2(2L_2-1),\quad N_3=L_3(2L_3-1),
$$
$$
N_4=4L_1L_2,\quad N_5=4L_2L_3,\quad N_6=4L_3L_1.
$$

The isoparametric mid-surface and displacement fields are $\mathbf x(r,s)=\sum_i N_i\mathbf x_i$ and $\mathbf u(r,s)=\sum_iN_i\mathbf u_i$. Surface tangents and the metric are

$$
\mathbf g_r=\mathbf x_{,r},\quad\mathbf g_s=\mathbf x_{,s},\quad a_{\alpha\beta}=\mathbf g_\alpha\cdot\mathbf g_\beta,\quad J_A=\|\mathbf g_r\times\mathbf g_s\|.
$$

At each node evaluate $\mathbf n_i\propto\mathbf g_r\times\mathbf g_s$. In v2 a **single centroid-derived sign** is applied to all nodal normals, to reduce winding-inconsistent normals. The local frame is $\mathbf e_1=\mathbf g_r/\|\mathbf g_r\|$, $\mathbf e_3=\mathrm{sign}_0(\mathbf g_r\times\mathbf g_s)/J_A$, $\mathbf e_2=\mathbf e_3\times\mathbf e_1$. The geometric conversion uses the metric inverse, not an explicit `2×2` Cartesian Jacobian.

## 2. Linear membrane formulation

The natural covariant membrane strains are defined via the linearized surface metric variation,

$$
\delta\varepsilon_{rr}=\mathbf g_r\cdot\mathbf u_{,r},\quad
\delta\varepsilon_{ss}=\mathbf g_s\cdot\mathbf u_{,s},\quad
\delta\varepsilon_{rs}=\frac12(\mathbf g_r\cdot\mathbf u_{,s}+\mathbf g_s\cdot\mathbf u_{,r}).
$$

For translational DOFs of node $I$, the corresponding natural rows are

$$
\mathbf B_{rr,I}=N_{I,r}\mathbf g_r^T,\quad
\mathbf B_{ss,I}=N_{I,s}\mathbf g_s^T,\quad
\mathbf B_{rs,I}=\frac12(N_{I,r}\mathbf g_s^T+N_{I,s}\mathbf g_r^T).
$$

Let $a^{\alpha\beta}$ denote the inverse surface metric, and set $c^\alpha_i=\mathbf g^\alpha\cdot\mathbf e_i$ for the physical tangent basis $i=1,2$. The tensor transformation converts the covariant components, with engineering shear third row:

$$
\varepsilon_{xx}=(c^r_1)^2\varepsilon_{rr}+(c^s_1)^2\varepsilon_{ss}+2c^r_1c^s_1\varepsilon_{rs},
$$
$$
\varepsilon_{yy}=(c^r_2)^2\varepsilon_{rr}+(c^s_2)^2\varepsilon_{ss}+2c^r_2c^s_2\varepsilon_{rs},
$$
$$
\gamma_{xy}=2\{c^r_1c^r_2\varepsilon_{rr}+c^s_1c^s_2\varepsilon_{ss}+(c^r_1c^s_2+c^s_1c^r_2)\varepsilon_{rs}\}.
$$

These relations assemble $\boldsymbol\varepsilon_m=B_m\mathbf q_e$, where $B_m\in\mathbb R^{3\times36}$; its rotational columns vanish. **Unlike MH6T or MITC6, no assumed membrane strain tying/projector is applied here.**

## 3. Linearized bending and director field

The director field uses vectors $\mathbf n_I$ and their natural gradients $\mathbf n_{,r}=\sum_IN_{I,r}\mathbf n_I$, $\mathbf n_{,s}=\sum_IN_{I,s}\mathbf n_I$. Each node contributes both a translational and a rotational curvature block. For example, the natural $rr$ curvature contributions in v2 are

$$
\delta\rho_{rr}^{(\theta_I)}=N_{I,r}\,[\mathbf g_r\times\mathbf n_I]\cdot\delta\boldsymbol\theta_I,
\qquad
\delta\rho_{rr}^{(u_I)}=-N_{I,r}\,\mathbf n_{,r}\cdot\delta\mathbf u_I.
$$

The $ss$ and symmetric $rs$ rows follow by exchanging $r,s$ and taking their half sums. Tensor conversion gives $B_b\in\mathbb R^{3\times36}$. This **sign pairing is specific to v2** and should not be attributed verbatim to an individual Simo equation without matching the paper's conventions.

## 4. Transverse shear and drilling

The shear director interpolant is **unnormalized**: $\mathbf t_0(r,s)=\sum_I N_I\mathbf n_I$. The covariant shear rows in v2 are

$$
\mathbf B_{rt,I}^{(u)}=N_{I,r}\mathbf t_0^T,\quad
\mathbf B_{st,I}^{(u)}=N_{I,s}\mathbf t_0^T,
$$
$$
\mathbf B_{rt,I}^{(\theta)}=N_I(\mathbf n_I\times\mathbf g_r)^T,\quad
\mathbf B_{st,I}^{(\theta)}=N_I(\mathbf n_I\times\mathbf g_s)^T.
$$

These rows are transformed to physical $xz$ and $yz$ shear components via the contravariant tangents to form $B_s\in\mathbb R^{2\times36}$. **Shear is conforming rather than tied:** inspect thin-shell locking and distorted/curved tests independently; matching one patch test does not establish uniform convergence.

For drilling stabilization, define the projected nodal derivatives $N_{I,x}$ and $N_{I,y}$ (obtained with the metric inverse and tangential projection). The formulation uses

$$
B_{d,I}^{(u)}=\frac12\big(N_{I,x}\mathbf e_2-N_{I,y}\mathbf e_1\big)^T,\qquad
B_{d,I}^{(\theta)}=-N_I\mathbf e_3^T,
$$

and $K_d=\int_{\Omega}\beta_d G h B_d^TB_d\,\mathrm dA$, with default $\beta_d=0.02$. The drilling treatment and its coefficient are **not a published Simo T6 specification**.

## 5. Elastic stiffness and integration

For an isotropic plane-stress shell,

$$
D=\frac{E}{1-\nu^2}\begin{bmatrix}
1 & \nu & 0 \\
\nu & 1 & 0 \\
0 & 0 & (1-\nu)/2
\end{bmatrix},\quad D_s=\kappa_s G I_2,\quad \kappa_s=5/6.
$$

The linear element stiffness is

$$
K_e=\int_{\Omega}\left[hB_m^TDB_m+\frac{h^3}{12}B_b^TDB_b+hB_s^TD_sB_s+\beta_dGhB_d^TB_d\right]\mathrm dA.
$$

The formulation uses **mean nodal thickness** $\bar h=(1/6)\sum_Ih_I$ rather than evaluating a variable interpolated $h(r,s)$. The integrations are grouped as follows.

| Contribution | Integration rule | Mathematical role |
|---|---|---|
| Membrane | 3 triangle points | Conforming covariant membrane strain, $B_m$ |
| Bending | 3 triangle points | Linearized director curvature, $B_b$ |
| Transverse shear | 3 triangle points | Conforming director shear, $B_s$; no MITC tying |
| Drilling stabilization | 6 triangle points | Penalty on relative normal rotation, $B_d$ |
| Consistent mass | 6 triangle points | Quadratic shape-function products and thickness moments |
| Pressure | Surface-load quadrature; order must match the load specification | Consistent traction integral over the actual surface |
| Thermal load | 3 triangle points | Membrane and bending thermal resultants |
| Initial-stress contribution | 3 triangle points in the restricted formulation | Translational geometric stiffness; scope qualified below |

<img width="1080" height="900" alt="simo_t6_gauss3" src="https://github.com/user-attachments/assets/ef9f4a2d-86af-4ac4-a529-db55a1d78bb0" />


<img width="1080" height="900" alt="simo_t6_gauss6" src="https://github.com/user-attachments/assets/f3ea42d9-b78a-4d21-8051-e36fc26159ac" />


The 3-point coordinates are $(1/6,1/6),(2/3,1/6),(1/6,2/3)$ with weights $1/6$ each. Six-point quadrature uses two symmetric triplets: $(a,a),(1-2a,a),(a,1-2a)$ with $a=0.445948490144588$ (weight $0.11169079483905$); and $(b,b),(1-2b,b),(b,1-2b)$ with $b=0.091576213509771$ (weight $0.054975871827661$).

## 6. Supplementary operators and limitations

**Mass:** v2 builds consistent translational and rotational mass blocks, $M_{ij}^{u}=\int\rho\bar hN_iN_j\,\mathrm dA$ and $M_{ij}^{\theta}=\int(\rho\bar h^3/12)N_iN_j\,\mathrm dA$. The direct shape-function lumping alternative integrates $N_i$, which for quadratic P2 corner functions can generate near-zero, cancellation-prone diagonal masses; **this is not a verified positive HRZ lumping scheme**. Validate independently before adopting it for eigen/SOL 103.

**Pressure:** use the consistent surface-traction integral. Pressure direction, winding, load sign and integration order must match the benchmark definition.

**Thermal:** the element builds an equivalent thermal force from isotropic $\varepsilon_{th}=\alpha\Delta T_0[1,1,0]^T$ and curvature input $\kappa_{th}=\alpha\Delta T_1[1,1,0]^T$; this is a project implementation whose gradient/sign interpretation must match the solver convention.

**Initial-stress geometric stiffness:** the restricted initial-stress approximation populates **only the global $U_z$ translational DOF rows/columns** with natural-gradient factors. A prescribed, physically justified membrane prestress is required; an arbitrary default compression is not a verified load basis. Consequently, it is **not a fully general, objective curved-shell $K_G$** and must not be treated as a verified Simo nonlinear tangent or production buckling formulation. Stress-based stability (SOL 105) is distinct from SOL 101 and geometrically nonlinear buckling.

## 7. Formulation evolution: v1p8 to v2

| Operator / control | v1p8 | v2 | Consequence / audit note |
|---|---|---|---|
| Director orientation | Nodal normals follow raw cross-product orientation | Single centroid sign applied at all six nodes | Improved director sign coherence for flipped winding; not a complete arbitrary-orientation proof |
| Bending rotation block | $\mathbf n_I\times\mathbf g_\alpha$ | $\mathbf g_\alpha\times\mathbf n_I$ | Rotation sign reversed |
| Bending translation block | Positive director-gradient term | Negative director-gradient term | Paired sign change, relevant to curved-shell rigid rotation |
| Shear director | $\mathbf t_0$ normalized pointwise | $\mathbf t_0=\sum N_I\mathbf n_I$ unnormalized | Restores consistent interpolation with the nodal director field |
| Drilling translation | Used metric-basis coefficients directly | Converted derivatives through local surface tangent projection | Corrects physical derivative interpretation |
| Drilling rotation | $-N_I\mathbf n_I$ | $-N_I\mathbf e_3(r,s)$ | Current local-point normal vs node normal: not algebraically identical on curved shell |
| Membrane | Pointwise covariant | Pointwise covariant | No MITC/assumed strain |
| Material matrix / 3+6 quadrature | As shown in earlier source | Same basic construction | Linear stiffness changes from bending/shear/drilling corrections |
| Mass / thermal / baseline geometric helper | Present | Present | Presence does not establish numerical validation |
| Pressure API | Local 6-point normal-force method | Delegates to shared pressure helper | Pressure conventions/availability depend on module |

The distinctions above describe v1p8 and v2; earlier revision history is outside this formulation summary.

## 8. Benchmark errors against independent references

The following figures use existing CSV data filtered by `formulation=SIMOT6`. They compare Nastran CTRIA6 and MYSTRAN SIMOT6 separately with the same independent reference. No solver was rerun. The formulation version recorded in a CSV is not by itself a proof of implementation identity.

For a nonzero target,

$$
e_{solver}(\%)=100\frac{|u_{solver}-u_{reference}|}{|u_{reference}|}.
$$

For a zero target, retain the absolute residual in physical units; never divide by zero. Bar heights show numerical error, not statistical confidence intervals.

### 8.1 CSI2-001: membrane and bending patch

Reference-bearing, finite paired rows from the production T6 master are used. Native element stress/resultant output and recovered grid values remain separate. Fiber stresses retain their signed analytic targets. Bars show the maximum pointwise error for each component/load case. Native transverse shear recovery remains subject to its recorded output-frame/recovery qualifications; no unconditional patch PASS is inferred.

![2-001 reference errors](<img width="2340" height="1620" alt="SIMOT6_2_001_reference_errors" src="https://github.com/user-attachments/assets/373a6970-1c55-4d17-aa0b-fac346b02a81" />
)

![2-001 absolute residuals at zero references](<img width="2340" height="630" alt="SIMOT6_2_001_zero_residuals" src="https://github.com/user-attachments/assets/23b8357d-9dd2-4f48-8efb-c7c911c92571" />
)


### 8.2 CSI2-005: plate convergence

Corrected C/D quarter-plate CSV results are used at nx=4,8,16,32, with Nastran and SIMOT6 on matching T6 meshes. Independent center-displacement magnitudes are12.97in (simple/uniform),16.95in (simple/point),2.60in (fixed/uniform),7.25in (fixed/point). Errors use displacement magnitude because these targets are magnitudes. Historical invalid-model master rows are excluded. Existing BAILOUT/AUTOSPC warnings and point-load T6 discrepancies remain visible.

![2-005 reference-error convergence](<img width="2160" height="1620" alt="SIMOT6_2_005_reference_convergence" src="https://github.com/user-attachments/assets/8c611db0-0666-4073-a5e1-6ade751e0416" />
)


### 8.3 CSI2-007: hemisphere convergence

The S2K-based thick hemisphere has an18-degree top opening. Independent loaded-point targets are $U_x=+0.0940$ft and $U_y=-0.0940$ft. Each curve uses the signed observation and target at the matching loaded location. The no-opening lower-bound0.0924ft value is not substituted. Meshes nx=4,8,16,32 belong to the same physical campaign.

![2-007 reference-error convergence](<img width="2160" height="900" alt="SIMOT6_2_007_reference_convergence" src="https://github.com/user-attachments/assets/19b15d15-efe9-47ca-9e44-3de2408cf7c6" />
)

## 9. Numerical provenance and limitations

Sources: `Shell_check_t6/production_20261018/exports/T6_comparison_master.csv`, `CSI_2005_quarter_20261019/release/paired.csv`, and `CSI_2007_s2k_comparison_20261020/convergence_Nastran_MYSTRAN_comparison.csv`. Source hashes are retained in `data_provenance.txt`. Missing results are not imputed. Grid-stress recovery is not substituted for native element output. Reference rounding, pressure conventions, support conditions and frame/fiber definitions remain part of interpreting errors. Agreement with Nastran alone is not a physical-accuracy criterion.

## 10. Conclusions

SIMOT6 v2 is a quadratic, linearized one-director shell with directly evaluated conforming membrane/shear strains, director-based curvature, and drilling stabilization. It is not an MITC tying element or a reproduction of the complete nonlinear/mixed Simo shell formulation. The benchmark figures separately quantify each solver's error to independent targets, with single-patch bars and genuine multimesh convergence curves. Stability, mass lumping, curved rigid-body invariance and recovery correctness require their own evidence.

## References

1. Simo, J. C., Fox, D. D. (1989). *On a stress resultant geometrically exact shell model. Part I: Formulation and optimal parametrization.* Computer Methods in Applied Mechanics and Engineering, 72, 267–304. DOI: https://doi.org/10.1016/0045-7825(89)90002-9.
2. Simo, J. C., Fox, D. D., Rifai, M. S. (1989). *On a stress resultant geometrically exact shell model. Part II: The linear theory; computational aspects.* Computer Methods in Applied Mechanics and Engineering, 73, 53–92. DOI: https://doi.org/10.1016/0045-7825(89)90098-4.
3. Simo, J. C., Fox, D. D., Rifai, M. S. (1990). *On a stress resultant geometrically exact shell model. Part III: Computational aspects of the nonlinear theory.* Computer Methods in Applied Mechanics and Engineering, 79, 21–70. DOI: https://doi.org/10.1016/0045-7825(90)90094-3.
4. MacNeal, R. H., Harder, R. L. (1985). *A proposed standard set of problems to test finite element accuracy.* Finite Elements in Analysis and Design, 1, 3–20. DOI: https://doi.org/10.1016/0168-874X(85)90003-4.
5. Computers and Structures, Inc. *SAP2000 / CSiBridge Software Verification—Shells, Problem 2-001 (prescribed displacement patch).* https://docs.csiamerica.com/manuals/csibridge/Verification/Analysis/Shells/Problem%202-001.pdf.
6. Computers and Structures, Inc. *Software Verification—Shells, Problem 2-005.* https://docs.csiamerica.com/manuals/csibridge/Verification/Analysis/Shells/Problem%202-005.pdf.
7. Computers and Structures, Inc. *Software Verification—Shells, Problem 2-007.* https://docs.csiamerica.com/manuals/csibridge/Verification/Analysis/Shells/Problem%202-007.pdf.

**Scientific attribution note:** The source’s year-based filename is a historical project designation. This documentation does **not** claim that a paper published in 1993 introduced the exact six-node discrete operators presented here.
