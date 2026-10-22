# Rezaiee–Pajand T6: Mathematical Formulation, Paper Traceability, and Revision Audit

**Author:** Ary  
**Project:** MYSTRAN — quadratic triangular shell (CTRIA6 / T6)  
**Date:** 23 October 2026  
**Implementation selected:** `Rezaiee2017_Tri6_v3.py` (**latest**). `Rezaiee2017_Tri6_Final(1).py` (**initial implementation**, despite the filename “Final”).  
**Scientific classification:** *Research adaptation of the linearized six-node assumed-strain shell; not a reproduction of the paper's full geometrically nonlinear formulation.*

## Abstract

This technical formulation describes a six-node quadratic triangular shell following the strain-assumption methodology of Rezaiee-Pajand, Arabi and Masoodi (published online 2017; journal 2018). It separates (i) concepts and equations sourced directly from the original paper, (ii) additional engineering provisions for a 36-DOF MYSTRAN-compatible element, and (iii) changes made between the original Python `Final` implementation and current `v3`. **“Final” is the historical first implementation; v3 is the authoritative current development.** The status of each theoretical feature and benchmark is recorded without attributing cross-formulation aggregates to Rezaiee.

## 1. Element topology and natural coordinates

<img width="1021" height="996" alt="nodal_numbering" src="https://github.com/user-attachments/assets/2c963946-fc9a-422a-b2e3-aec3399113d7" />

The standard triangular domain is `r ≥ 0`, `s ≥ 0`, `r+s ≤ 1`. Let `L1=1−r−s`, `L2=r`, `L3=s`. Local nodes 1,2,3 are corners, while nodes 4,5,6 are mid-edge (1–2, 2–3, 3–1). The **quadratic P2 interpolation**, adopted by the Python implementations and consistent with the paper's six-node interpolation, is

$$
\begin{aligned}
N_1&=L_1(2L_1-1), &N_2&=L_2(2L_2-1), &N_3&=L_3(2L_3-1),\\
N_4&=4L_1L_2,&N_5&=4L_2L_3,&N_6&=4L_3L_1.
\end{aligned}
$$

The shell midsurface mapping and covariant tangent vectors are

$$\mathbf{x}_0(r,s)=\sum_{I=1}^{6}N_I\mathbf{X}_I,\quad \mathbf{g}_r=\sum_I N_{I,r}\mathbf{X}_I,\quad\mathbf{g}_s=\sum_I N_{I,s}\mathbf{X}_I,\quad J_A=\|\mathbf{g}_r\times\mathbf{g}_s\|.$$

The six-node numbering and shape functions are retained from initial to v3.

## 2. Kinematics: published element versus MYSTRAN adaptation

**Paper, Eqs. (1)–(7), pp. 3–5:** The geometric shell is interpolated along a director using nodal thicknesses `a_I` and the through-thickness coordinate `t`:

$$\mathbf{x}(r,s,t)=\sum_{I=1}^6 N_I\mathbf{x}_I+\frac{t}{2}\sum_{I=1}^6a_I N_I\mathbf{V}_{nI}.$$

Each node has **five** independent variables, three translations and **two** director rotations (`α_I,β_I`), producing **30 DOF**. The paper includes the quadratic rotational increment

$$\Delta\mathbf{V}_{nI}=-\alpha_I\mathbf{V}_{2I}+\beta_I\mathbf{V}_{1I}-\frac12(\alpha_I^2+\beta_I^2)\mathbf{V}_{nI}.$$

**Python v3:** stores **six** global components `[u_x,u_y,u_z,θ_x,θ_y,θ_z]` per node, giving a **36×36** element matrix, and introduces a drilling stabilization term. This is a software-compatibility extension, **not** a drilling DOF proposed by Rezaiee-Pajand et al. Moreover v3 is a *linearized* strain/stiffness implementation: the paper's finite-rotation update, full total-Lagrangian tangent and nonlinear Newton–Raphson/GDCM procedure are not implemented. Do not equate v3's `k_local()` with the full nonlinear algorithm.

For the paper, Eqs. (8)–(11), the incremental Green–Lagrange strain has linear and nonlinear contributions:

$$\Delta E_{ij}=\frac12\left(\mathbf{u}_{,i}\cdot\mathbf{g}_j+\mathbf{g}_i\cdot\mathbf{u}_{,j}+\mathbf{u}_{,i}\cdot\mathbf{u}_{,j}\right)=e_{ij}+\eta_{ij}.$$

The v3 linear shell uses `B_m`, `B_b`, `B_s`, and `B_drill` instead of assembling this complete nonlinear expression.

## 3. Membrane assumed-strain interpolation

**Paper, Eqs. (12)–(15), (18), Figures 2–3:** Rather than directly use every displacement-derived membrane strain component at all points, the method samples selected **tying stations** and interpolates covariant strain components. The paper prescribes

$$\hat\varepsilon_{rr}=a_1+b_1r+c_1s,\qquad \hat\varepsilon_{ss}=a_2+b_2r+c_2s,$$
$$\hat\varepsilon_{qq}=a_3+b_3r+c_3(1-r-s),\qquad \hat\varepsilon_{rs}=\tfrac12(\hat\varepsilon_{rr}+\hat\varepsilon_{ss}-\hat\varepsilon_{qq}).$$

Here `q` is the third oblique direction; **`ε_qq=ε_rr+ε_ss−2ε_rs`** is the normalization used in v3. The initial `Final` code instead represented the `q` strain with an internally half-scaled definition; a change in that scale **alone** is not necessarily a constitutive bug, provided every conversion is consistent.

**Paper's special coordinates, Eq. (19):**

$$r_1=s_1=\tfrac12-\tfrac{1}{2\sqrt3},\quad r_2=s_2=\tfrac12+\tfrac{1}{2\sqrt3},\quad r_3=s_3=\tfrac13,\quad r_c=s_c=\tfrac{1}{\sqrt3}.$$

The **current v3 implementation** fits affine membrane operators through the following selected stations (as configured in source, not asserted mathematically identical to every published coefficient):

| Field | Three source v3 stations |
|---|---|
| `rr` | `(r1,0)`, `(r2,0)`, `(r1,1/√3)` |
| `ss` | `(0,r1)`, `(0,r2)`, `(1/√3,r1)` |
| `qq` | `(r2,r1)`, `(r1,r2)`, `(r1,r1)` |

These pointwise row fits remain to be proved *exactly equivalent* to the published closed-form coefficients of Eqs. (15) and centroid-reference transmission of Eqs. (12)–(13). The Python's output-point physical transformation is **not** the same as a demonstrated implementation of the paper's centroid-reference strain transmission.

## 4. Transverse shear assumed-strain interpolation

**Paper Eqs. (16)–(17)** specify coupled quadratic assumed shear fields:

$$\hat\varepsilon_{rt}=a_1+b_1r+c_1s+d_1rs+f_1s^2,$$
$$\hat\varepsilon_{st}=a_2+b_2r+c_2s+d_2rs+f_2r^2,\qquad d_1=-f_2,\quad d_2=-f_1.$$

The v3 code implements the equivalent **eight-coefficient polynomial structure** using coefficients `[a1,b1,c1,a2,b2,c2,d,e]`, with

$$\tilde\gamma_{rt}=a_1+b_1r+c_1s+d\,rs+e\,s^2,$$
$$\tilde\gamma_{st}=a_2+b_2r+c_2s-d\,r^2-e\,rs,$$

and solves an **8×8 tying-value system**, not the paper's literal coefficient-by-coefficient expressions. The v3 stations are `rt`: `(r1,0),(r2,0),(1/3,1/3)`; `st`: `(0,r1),(0,r2),(1/3,1/3)`; third shear combination `qt`: `(r2,r1),(r1,r2)`. The `qt` operator uses `(γ_st−γ_rt)/√2` in this reconstruction. **Caution:** sharing polynomial *form* does not prove identical quadrature, covariant transmission, or nonlinear response.

<img width="2064" height="1010" alt="tying_points" src="https://github.com/user-attachments/assets/e3f44659-ed85-4343-9f87-41964439781c" />

## 5. Membrane, bending, shear, drilling stiffness in Python v3

With isotropic elastic material, the code uses

$$\mathbf C_{mb}=\frac{E}{1-\nu^2}\begin{bmatrix}1&\nu&0\\\nu&1&0\\0&0&(1-\nu)/2\end{bmatrix},\qquad\mathbf C_s=\kappa G\mathbf I_2,$$

where `G=E/[2(1+ν)]`, `κ=5/6`, `β_drill=0.02`, and `h̄=(Σ h_I)/6`. The exact **implemented linearized stiffness** is

$$\boxed{\mathbf K_e=\int_A \left[\bar h\,\mathbf B_m^T\mathbf C_{mb}\mathbf B_m+\frac{\bar h^3}{12}\mathbf B_b^T\mathbf C_{mb}\mathbf B_b+\bar h\,\mathbf B_s^T\mathbf C_s\mathbf B_s+\bar h\,\beta_{drill}G\mathbf B_d^T\mathbf B_d\right]dA.}$$

Here **assumed strain** modifies the membrane/shear operators; `B_b` represents a separate director-based bending operator and `B_d` is a numerical drilling penalty. The paper supplies the assumed-strain method and nonlinear virtual-work framework; it **does not prescribe** the specific v3 drilling penalty `β=0.02`, six nodal DOFs, or the v3 reduced/selective quadrature arrangement.

## 6. Gauss-point locations: paper versus project code

<img width="1061" height="996" alt="gauss3" src="https://github.com/user-attachments/assets/0a6a7866-afbc-4a58-844a-139e4f97ade8" />

The v3 stiffness uses the 3-point triangle rule with stations `(1/6,1/6)`, `(2/3,1/6)`, `(1/6,2/3)` and weights `1/6` each (weights integrate the reference triangle of area `1/2`).



It uses the six-point degree-four triangle rule for drilling/mass, with `(a,a),(b,a),(a,b)` for `a=0.445948490144588`, `b=0.108103018168070` (each weight `0.11169079483905`), and `(c,c),(d,c),(c,d)` for `c=0.091576213509771`, `d=0.816847572980459` (each weight `0.054975871827661`). The *source paper*, Section 5, p. 9, explicitly uses **seven Gaussian points** in its nonlinear element integration. Thus **the displayed 3/6 Gauss diagrams depict project implementation only**, not the published 7-point pattern.

## 7. Version history — INITIAL `Final` → LATEST `v3`

| Feature | Initial `Rezaiee2017_Tri6_Final` | Latest `Rezaiee2017_Tri6_v3` | Source / rationale |
|---|---|---|---|
| Status | Initial prototype despite label “Final” | **Selected current implementation** | Explicit project version order |
| P2 geometry, nodes 1–6 | Present | Retained | Paper six-node element; standard quadratic interpolation |
| DOF | 6 per node | 6 per node | Both software adaptations (paper: 5/node) |
| Director orientation | Nodal normals computed independently | Centroid-based shared sign used consistently | Project correction for CW/CCW invariance; not a paper algorithm |
| Membrane tying layout | Earlier sample layout, including off-paper symmetry assumptions | Revised stations aligned to inspected Figure 2 | Paper-inspired reconstruction; exact coefficient equivalence pending |
| `q` membrane contraction | Half-scaled `qq` convention | `qq=rr+ss−2rs` | Consistent change of strain representation |
| Shear interpolation | Coupled eight-coefficient polynomial | Same polynomial family; revised station assignment | Eqs. (16)–(17) paper-inspired |
| Shear director | Pointwise interpolated director normalized in earlier version | **Unnormalized** interpolated director in translation operator | v3 internal kinematic consistency correction |
| Bending curvature | Earlier rotation cross-product and translation signs | Reversed rotational cross-product plus matched negative director-gradient translation block | Restores sampled infinitesimal curved rigid-motion behavior |
| Drilling | Earlier oblique-derivative treatment | Derivatives projected to physical directions with pointwise surface normal | Numerical stabilization/covariance fix |
| Quadrature | 3-point membrane/bending/shear; 6-point drilling | Same | Project choice; paper specifies 7 points |
| Nonlinear TL tangent / NR/GDCM | Not reproduced | Not reproduced | **Still excluded** despite paper emphasis |
| Pressure, mass | Source helper functionality varies between revisions | Mass 6-point; pressure delegates to `shell_pressure.pressure_global` | Engineering helper code, not validated as original paper equation |

### Precise progression from earliest source

The direct comparison is **initial Final → current v3**, while an intermediate v2 was used in debugging and testing. The previously documented v2 → v3 T9 improvement must **not be mislabelled Final → v3**, since no matched Final T9 baseline has been established for that entry. The initial `Final` has winding/patch issues documented in a prior separate test campaign; v3 correction is not proof of verbatim reproduction of the journal.

## 8. Numerical evidence and its limits

### 8.1 Project Python regression — latest v3 (not Nastran)

| Benchmark | Mesh or case | v3 response / reference | Error / qualification |
|---|---:|---:|---|
| T6 | N=16 | 1.009987 | +0.9987% |
| T7 | N=16 | 1.095289 | +9.5289% |
| T8 | N=16 | 1.002515 | +0.2515% |
| T9 | N=16 | 1.172421649 | **+17.2422%**, within project ±20% rule, fails ±10% |
| SAP 2-003 | LC1, Nx=64 | 0.999460 | −0.0540% (independent benchmark comparison) |
| SAP 2-004 | Fixed Uz, Nx=24 | 1.016236 | +1.6236% |
| SAP 2-004 | Fixed Uy, Nx=24 | 1.074576 | +7.4576% |

The underlying v3 verification document reports maximum rigid residual `1.38×10⁻¹⁵`, a six near-zero rigid-mode count in sampled models and a maximum sampled corner-permutation relative K discrepancy of `8.65×10⁻¹²`. These are *sampled code-verification tests*, not a universal mathematical proof.

### 8.2 Python intermediate v2 → v3 regression (not Final → v3)

| Test | v2 | v3 | Interpretation |
|---|---:|---:|---|
| T9, N=16 ratio | 0.465572 | 1.172422 | Significant response change; latest closer under ±20% criterion, but no demonstrated convergence |
| 2-004 Fixed Uz error | +124.8353% | +1.6236% | Major improvement in twisted-beam case |
| 2-004 Fixed Uy error | +96.4471% | +7.4576% | Major improvement |
| 2-002 case 4 error | −5.9795% | −10.2136% | **Regression/trade-off:** not every case improves |

### 8.3 Nastran CTRIA6 and MYSTRAN REZAIEE

The project **T6 comparison summary** dated 18 October 2026 records native OP2-based checks across four selectors (`SIMOT6`, `MITC6`, `MH6T`, `REZAIEE`), including flat patch tests and an extensive thick-shell campaign. It reports aggregate patch quality across those selectors and states that curved element force/stress frame equivalence is **pending**. The source available here does **not** give a filterable **Rezaiee-only** CSV with individual matched Nastran and MYSTRAN numerical outputs. Consequently **no Rezaiee-only solver percent differences are asserted** in this report. Full table construction must join equal case, physical model, node/element, station, load case, stress fiber, basis, units, and signed response; Nastran is not itself an independent CSI reference.

## 9. Published nonlinear capabilities deliberately NOT claimed for current v3

- Finite director rotations and `α²+β²` update as a nonlinear state variable (paper Eqs. 3–7)
- Green–Lagrange quadratic strain contribution and its consistent tangent (Eqs. 8–11, 20–21)
- Newton–Raphson or generalized displacement control method (paper Section 5)
- Post-buckling response curves, large-displacement equilibrium paths and the published nonlinear test reproductions
- Paper's exact centroid-reference covariant strain transmission (Eqs. 12–13)
- Seven-point integration of the original nonlinear paper

**Interpretation:** This project can legitimately cite the paper as inspiration for the assumed strain *shape* and tying layout while explicitly documenting independent v3 corrections. Compatibility with an MSC/Nastran CTRIA6 deck does not make the elements formulation-identical.

## 10. References and verification suite

1. Rezaiee-Pajand, M., Arabi, E., Masoodi, A. R. (2017 online; 2018 journal). **A triangular shell element for geometrically nonlinear analysis**. *Acta Mechanica*. DOI: https://doi.org/10.1007/s00707-017-1971-8. **Primary paper**, notably Eqs. (1)–(21), Figs. 1–3, Section 5.
2. Kim, D.-N. and Bathe, K.-J. (2009). **A triangular six-node shell element**. *Computers & Structures*, 87, 1451–1460. DOI: https://doi.org/10.1016/j.compstruc.2009.05.002.
3. Lee, P.-S. and Bathe, K.-J. (2004). **Development of MITC isotropic triangular shell finite elements**. *Computers & Structures*, 82, 945–962.
4. CSI, **SAP2000 Analysis Verification Manual**, Shell problems 2-001 through 2-017 (consult the matching published PDF and physical model for each case), e.g. https://docs.csiamerica.com/manuals/sap2000/Verification/Analysis/Shells/Problem%202-001.pdf and https://docs.csiamerica.com/manuals/sap2000/Verification/Analysis/Shells/Problem%202-004.pdf. **Separate verification cases, not citations from Rezaiee (2017).**
5. Project provenance: `Rezaiee2017_Tri6_Final(1).py` (initial), `Rezaiee2017_Tri6_v3.py` (latest), `Rezaiee_T6_v3_verification.md` (Python numerical diagnostics), `T6_comparison_summary.md` (multiformulation native OP2 audit).

---

*Editorial rule:* “From paper” means explicit in paper; “retained/adapted” means inspired by but modified; “added in MYSTRAN/Python” means absent from the paper; and “not implemented” must not be described as verified. All numeric results should identify solver and provenance.
