# Simo T6 Shell Element: CTRIA6 Shell 

**A technical research formulation and source-traceability report for MYSTRAN quadratic CTRIA6**  
**Author:** Ary  
**Date:** 23 October 2026  
**Version audited:** `Simo1993_Tri6_ShellElement_v2.py` (current); `Simo1993_Tri6_ShellElement_v1p8(1).py` (earlier).  
**Evidence classification:** source-code audit; historical project CSV excerpts; **not** a new MSC Nastran execution.  

## Abstract

This document describes the project's six-node, quadratic, linear-elastic triangular shell (SimoT6), its geometric interpolation and membrane/bending/transverse-shear/drilling operators, integration rules, auxiliary mass/pressure/thermal/geometric-stiffness routines, and the correspondence with the Simo–Fox–Rifai shell literature. It carefully separates **published shell-theory concepts** from **project-specific discretization and stabilization**, and separates existing Python numerical evidence from still-unconfirmed per-case Nastran comparisons. The immediate intended scope is **MYSTRAN SOL 101 linear statics**; modal, thermal and buckling integrations require separate verification.

## 1. Formulation provenance and scope

The relevant published lineage is the geometrically exact one-director/stress-resultant shell formulation of Simo and Fox and its linear computational treatment by Simo, Fox and Rifai. **This project Python v2 is NOT a direct implementation of a published 1993 six-node Simo element.** The class name is a project historical identifier. The source implements standard isoparametric P2 displacement/geometry interpolation and directly evaluated conforming strains, rather than reproducing every mixed Hellinger–Reissner or nonlinear director-field treatment of the research series.

| Feature | Published Simo/Fox/Rifai lineage | Project Python v2 (audited source) | Status |
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

![Nodal numbering for six-node triangle](figures/simo_t6_nodes.svg)

Reference domain $\hat\Omega=\{(r,s):r\geq0,s\geq0,r+s\leq1\}$, with area $1/2$. Coordinates are $(0,0),(1,0),(0,1),(1/2,0),(1/2,1/2),(0,1/2)$ for nodes 1 to 6. Each node carries $\mathbf q_i=[u_x,u_y,u_z,\theta_x,\theta_y,\theta_z]^T$, so $\mathbf q_e\in\mathbb R^{36}$ and `T_matrix()=I`.

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

At each node the code evaluates $\mathbf n_i\propto\mathbf g_r\times\mathbf g_s$. In v2 a **single centroid-derived sign** is applied to all nodal normals, to reduce winding-inconsistent normals. The local frame is $\mathbf e_1=\mathbf g_r/\|\mathbf g_r\|$, $\mathbf e_3=\operatorname{sign}_0(\mathbf g_r\times\mathbf g_s)/J_A$, $\mathbf e_2=\mathbf e_3\times\mathbf e_1$. In the current source, the geometric conversion uses the metric inverse, not an explicit `2×2` Cartesian Jacobian.

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

Let $a^{\alpha\beta}$ denote the inverse surface metric, and set $c^\alpha_i=\mathbf g^\alpha\cdot\mathbf e_i$ for the physical tangent basis $i=1,2$. The source's `_tensor_to_physical` converts the tensor components, with engineering shear third row:

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

The code creates director vectors $\mathbf n_I$ and their natural gradients $\mathbf n_{,r}=\sum_IN_{I,r}\mathbf n_I$, $\mathbf n_{,s}=\sum_IN_{I,s}\mathbf n_I$. Each node contributes both a translational and a rotational curvature block. For example, the natural $rr$ curvature contributions in source v2 are

$$
\delta\rho_{rr}^{(\theta_I)}=N_{I,r}\,[\mathbf g_r\times\mathbf n_I]\cdot\delta\boldsymbol\theta_I,
\qquad
\delta\rho_{rr}^{(u_I)}=-N_{I,r}\,\mathbf n_{,r}\cdot\delta\mathbf u_I.
$$

The $ss$ and symmetric $rs$ rows follow by exchanging $r,s$ and taking their half sums. Tensor conversion gives $B_b\in\mathbb R^{3\times36}$. This **sign pairing is specific to v2** and should not be attributed verbatim to an individual Simo equation without matching the paper's conventions.

## 4. Transverse shear and drilling

The code's shear director interpolant is **unnormalized**: $\mathbf t_0(r,s)=\sum_I N_I\mathbf n_I$. The covariant shear rows in v2 are

$$
\mathbf B_{rt,I}^{(u)}=N_{I,r}\mathbf t_0^T,\quad
\mathbf B_{st,I}^{(u)}=N_{I,s}\mathbf t_0^T,
$$
$$
\mathbf B_{rt,I}^{(\theta)}=N_I(\mathbf n_I\times\mathbf g_r)^T,\quad
\mathbf B_{st,I}^{(\theta)}=N_I(\mathbf n_I\times\mathbf g_s)^T.
$$

These rows are transformed to physical $xz$ and $yz$ shear components via the contravariant tangents to form $B_s\in\mathbb R^{2\times36}$. **Shear is conforming rather than tied:** inspect thin-shell locking and distorted/curved tests independently; matching one patch test does not establish uniform convergence.

For drilling stabilization, define the projected nodal derivatives $N_{I,x}$ and $N_{I,y}$ (code uses metric inverse and tangential projection). The source uses

$$
B_{d,I}^{(u)}=\frac12\big(N_{I,x}\mathbf e_2-N_{I,y}\mathbf e_1\big)^T,\qquad
B_{d,I}^{(\theta)}=-N_I\mathbf e_3^T,
$$

and $K_d=\int_{\Omega}\beta_d G h B_d^TB_d\,\mathrm dA$, with default $\beta_d=0.02$. The drilling treatment and its coefficient are **not identified as a published Simo T6 specification**.

## 5. Elastic stiffness and integration

For an isotropic plane-stress shell,

$$
D=\frac{E}{1-\nu^2}\begin{bmatrix}1&\nu&0\\\nu&1&0\\0&0&(1-\nu)/2\end{bmatrix},\quad D_s=\kappa_s G I_2,\quad \kappa_s=5/6.
$$

The linear element stiffness in `k_local()` is

$$
K_e=\int_{\Omega}\left[hB_m^TDB_m+\frac{h^3}{12}B_b^TDB_b+hB_s^TD_sB_s+\beta_dGhB_d^TB_d\right]\mathrm dA.
$$

The uploaded implementation uses **mean nodal thickness** $\bar h=(1/6)\sum_Ih_I$ rather than evaluating a variable interpolated $h(r,s)$. The integrations are grouped as follows.

| Stiffness block | Rule | Source method |
|---|---|---|
| Membrane | 3 points | `k_local`: `Bm` |
| Bending | 3 points | `k_local`: `Bb` |
| Transverse shear | 3 points | `k_local`: `Bs` |
| Drilling | 6 points | `k_local`: `Bdrill` |
| Consistent mass | 6 points | `m_local` |
| Pressure | Rule delegated to `shell_pressure.pressure_global` in v2 | `f_pressure` |
| Thermal vector | 3 points | `f_thermal` |
| Geometric stiffness helper | 3 points | `k_geometric` (restricted implementation) |

![Three-point Gauss quadrature](figures/simo_t6_gauss3.svg)

![Six-point Gauss quadrature](figures/simo_t6_gauss6.svg)

The 3-point coordinates are $(1/6,1/6),(2/3,1/6),(1/6,2/3)$ with weights $1/6$ each. Six-point quadrature uses two symmetric triplets: $(a,a),(1-2a,a),(a,1-2a)$ with $a=0.445948490144588$ (weight $0.11169079483905$); and $(b,b),(1-2b,b),(b,1-2b)$ with $b=0.091576213509771$ (weight $0.054975871827661$).

## 6. Supplementary operators and limitations

**Mass:** v2 builds consistent translational and rotational mass blocks, $M_{ij}^{u}=\int\rho\bar hN_iN_j\,\mathrm dA$ and $M_{ij}^{\theta}=\int(\rho\bar h^3/12)N_iN_j\,\mathrm dA$. Its `consistent=False` path integrates $N_i$ directly, which for quadratic P2 corner functions can generate near-zero, cancellation-prone diagonal masses; **this is not a verified positive HRZ lumping scheme**. Validate independently before adopting it for eigen/SOL 103.

**Pressure:** v2 calls `shell_pressure.pressure_global`; the dependency and default pressure orientation require verification on the same winding and load sign as Nastran.

**Thermal:** the element builds an equivalent thermal force from isotropic $\varepsilon_{th}=\alpha\Delta T_0[1,1,0]^T$ and curvature input $\kappa_{th}=\alpha\Delta T_1[1,1,0]^T$; this is a project implementation whose gradient/sign interpretation must match the solver convention.

**Initial-stress geometric stiffness:** the `k_geometric` helper populates **only the global $U_z$ translational DOF rows/columns** with natural-gradient factors. It also defaults to an arbitrary negative membrane stress when `stress_m=None`. Consequently, it is **not a fully general, objective curved-shell $K_G$** and must not be treated as a verified Simo nonlinear tangent or production buckling formulation. Stress-based stability (SOL 105) is distinct from SOL 101 and geometrically nonlinear buckling.

## 7. v1p8 → v2: code-proven revision log

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

The files supplied are `v1p8` and `v2`. Earlier v1–v1p7 historical evolution is not directly audited here.

## 8. Validation evidence — Python, Nastran, and MYSTRAN must remain distinct

### 8.1 Source attribution of CSV evidence

Two Library CSVs were located: `latest_displacement_summary.csv` (specific **Python SimoT6v2** records), and `solver_comparison_summary.csv` (pooled Nastran/MYSTRAN comparison with **no formulation key**). The latter is **not assignable to SimoT6**. The user-requested matched **SIMOT6 v2 / CTRIA6** per-case CSV for 2-001, 2-005 and 2-007 has **not been reliably identified** among the available files. The following table deliberately leaves unsupported cells blank.

### 8.2 SAP 2-001 — prescribed membrane and bending patch

The project history reports that SimoT6 v1p8 passed both membrane and bending patch with errors approximately 0.0000%. This is *historical Python patch-test evidence*, not a contemporaneous recovered **Nastran–SimoT6v2** matched-data table. Reproducible values per element/sample point, mesh ID, matching nodal output frame and GPSTRESS for this campaign remain to be attached.

| Case | Python SimoT6 | Nastran CTRIA6 | MYSTRAN SIMOT6 | Independent reference | Evidentiary status |
|---|---|---|---|---|---|
| 2-001 LC1 membrane patch | Historical v1p8 ~0.0000% error | — | — | prescribed field in CSI verification | Historical only; matched CSV missing |
| 2-001 LC2 bending patch | Historical v1p8 ~0.0000% error | — | — | constant plate bending field | Historical only; matched CSV missing |

### 8.3 SAP 2-005 — plate bending / pressure

The available **Python-only** CSV `latest_displacement_summary.csv` records `Simo1993_Tri6_ShellElement_v2`: plate $U_z=-13.230688877881176$, independent reference $U_z=-12.97$, relative magnitude error $+2.009937377649762\%$ (the supplied CSV uses this sign convention). The specific model's mesh, support variant and pressure convention must be verified before identifying it with a particular CSI 2-005 load case or a paired Nastran deck.

| Source | $U_z$ | Reference | Error (%) | Status |
|---|---:|---:|---:|---|
| Python SimoT6 v2, plate case | -13.23068888 | -12.97000000 | +2.00994 | Located in project CSV; exact CSI case mapping pending |
| MSC Nastran CTRIA6, matched case | — | — | — | Per-case result not found |
| MYSTRAN SIMOT6, matched case | — | — | — | Per-case result not found |

Do not substitute general solver-aggregate `2-005` RMS numbers here. For a valid 2-005 convergence comparison, use identical **quarter plate** models with $4\times4$, $8\times8$, $16\times16$, $32\times32$ meshes, matching concentrated/uniform load and clamped/simply-supported conditions; generate **line convergence plots** of $U_z$ and relative error by mesh refinement, separately distinguishing Nastran, MYSTRAN and Python.

### 8.4 SAP 2-007 — hemispherical shell benchmark

No **SIMOT6-specific** 2-007 Nastran/MYSTRAN matched values were verified in the accessible CSVs. The generic `solver_comparison_summary.csv` has a `2-007 displacement` aggregate, but it mixes formulation identity; assigning its 99.04% L2 difference to SimoT6 would be **incorrect**.

| Test | Python SimoT6 v2 | Nastran CTRIA6 | MYSTRAN SIMOT6 | Status |
|---|---|---|---|---|
| SAP 2-007 reference/check-point displacement | — | — | — | Need matched solver CSV / OP2 and exact node mapping |

For curved shells, distinguish common global XYZ displacements from frame-dependent membrane resultants, plate moments and fiber stress. Verify director/winding and nodal averaging separately; do not confuse contour-averaged GPSTRESS with per-element natural-coordinate recovery.

### 8.5 Other CSV observations (not substituted for requested cases)

The same Python-only CSV reports SimoT6 v2 displacement values for `linear_geometry` and `cylinder`: respectively $U_z=-0.3043374677479359$ and $-0.301721620617218$ with an independent reference $-0.3086$. These labels are not sufficient to relabel the cases SAP 2-007. They are only auxiliary curved-shell evidence.

## 9. Codex handoff: add verified 2-001 / 2-005 / 2-007 comparison

1. Resolve a precise mapping `formulation=SIMOT6` and `version=v2` (plus the actually ported Fortran selector/revision) and Nastran `CTRIA6` using matching mesh coordinates, connectivity, thickness, material, PLOAD4, constraints and analysis case.
2. Extract paired data from OP2, not F06 for large campaigns. Record solver/OP2 provenance, load case, target node, coordinate system and result type. For 2-001 include membrane stress resultant, bending moment and signed surface stresses; for 2-005 use center displacement and reaction/pressure resultant; for 2-007 use prescribed hemisphere check points and both Cartesian and element-frame recovery.
3. Recommended CSV schema: `benchmark,case,mesh,load_type,support,solver,element,version,subcase,quantity,component,probe_id,frame,reference,result,error_pct,source_op2,source_dat`.
4. Use one signed-error convention $e=100(u-u_{ref})/\max(|u_{ref}|,\varepsilon)$ and additionally absolute error and force equilibrium. Do not compare magnitudes alone where physical sign matters.
5. Complete the above tables with paired source rows, and produce line convergence figures for 2-005; no bar charts.
6. Assert regression of `k_local()`, symmetry, six rigid-body modes for the unconstrained geometric formulation, CW/CCW invariance and curved rigid rotation. Audit whether drilling penalty alters expected zero-energy modes and whether the 3-point shear term locks for thin curved shells.
7. Preserve independent CSV alongside Markdown in `update/Formulation/Shell/results/`; put SVG figures in `update/Formulation/Shell/figures/` and maintain relative paths.

## 10. Conclusions

The **Python v2 SimoT6** is a specific **linearized one-director P2 shell implementation** with conforming membrane/shear, director-based curvature and a drilling penalty. It is not a reproduction of the full Simo nonlinear or mixed shell finite element. Source-documented corrections since v1p8 are significant to winding, bending, shear and drilling; they deserve a distinct version identifier. A Python plate result near the CSI 2-005 reference is available, but the requested **solver-paired SIMOT6-vs-Nastran CSV** remains unconfirmed and is not manufactured here.

## References

1. Simo, J. C., Fox, D. D. (1989). *On a stress resultant geometrically exact shell model. Part I: Formulation and optimal parametrization.* Computer Methods in Applied Mechanics and Engineering, 72, 267–304. DOI: https://doi.org/10.1016/0045-7825(89)90002-9.
2. Simo, J. C., Fox, D. D., Rifai, M. S. (1989). *On a stress resultant geometrically exact shell model. Part II: The linear theory; computational aspects.* Computer Methods in Applied Mechanics and Engineering, 73, 53–92. DOI: https://doi.org/10.1016/0045-7825(89)90098-4.
3. Simo, J. C., Fox, D. D., Rifai, M. S. (1990). *On a stress resultant geometrically exact shell model. Part III: Computational aspects of the nonlinear theory.* Computer Methods in Applied Mechanics and Engineering, 79, 21–70. DOI: https://doi.org/10.1016/0045-7825(90)90094-3.
4. MacNeal, R. H., Harder, R. L. (1985). *A proposed standard set of problems to test finite element accuracy.* Finite Elements in Analysis and Design, 1, 3–20. DOI: https://doi.org/10.1016/0168-874X(85)90003-4.
5. Computers and Structures, Inc. *SAP2000 / CSiBridge Software Verification—Shells, Problem 2-001 (prescribed displacement patch).* https://docs.csiamerica.com/manuals/csibridge/Verification/Analysis/Shells/Problem%202-001.pdf.
6. Computers and Structures, Inc. *Software Verification—Shells, Problem 2-005.* https://docs.csiamerica.com/manuals/csibridge/Verification/Analysis/Shells/Problem%202-005.pdf.
7. Computers and Structures, Inc. *Software Verification—Shells, Problem 2-007.* https://docs.csiamerica.com/manuals/csibridge/Verification/Analysis/Shells/Problem%202-007.pdf.

**Scientific attribution note:** The source’s year-based filename is a historical project designation. This documentation does **not** claim that a paper published in 1993 introduced the exact six-node Python operators presented here.
