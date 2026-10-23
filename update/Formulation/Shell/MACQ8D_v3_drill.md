# MACQ8D: A Modified Eight-Node Flat-Shell Element with Independent Drilling Stabilization — formulation and Numerical Verification

**Technical formulation and benchmark report — MYSTRAN CQUAD8 (`MACQ8D`, v3 drilling)**  
**Author:** Ary  
**Revision:** 23 October 2026  
**Evidence:** Existing formulation notes; archived `Q8_comparison_master.csv` and audited corrected-quarter `paired.csv`. Numerical tables and plots here **reuse archived results, not fresh Nastran/MYSTRAN runs**. The project-specific element is not presented as a verbatim reproduction of the 1992 publication.

## Abstract

An eight-node flat-facet shell element with a MacNeal–Harder modified interpolation field and a weak, independent drilling-rotation stabilization is documented. The surface geometry is interpolated with standard Q8 serendipity functions while membrane, bending, and transverse shear strains employ modified field functions. Six nodal degrees of freedom yield 48 global degrees of freedom, and no internal degrees of freedom are condensed. The v3 drilling regularization penalizes only the independent normal rotation, in contrast to the alternative continuum-spin-coupled stabilization. Numerical evidence is limited here to separate solver errors against analytic2-001 targets and independent2-005 plate targets, using existing CSV results. Recorded status and solver warnings prevent an unconditional validation claim.

## 1. Origin and evolution of the formulation

| Research stage | Formulation change | Scope and qualification |
|---|---|---|
| Native MacNeal Q8 | Flat local orthonormal frame, modified displacement interpolation, standard geometry, full 3×3 integration | Drilling rotation untreated; eight elemental drilling zero modes in addition to the six rigid-body modes, as reported in source comments |
| v2 | Optional rotation-convention handling in historical comparison driver | **The uploaded base native class, not the v2 convention-switch class, is the parent of v3** |
| v2 drilling | Spin-coupled penalty: normal rotation minus in-plane continuum spin | Separate research alternative, **not the v3 drilling operator** |
| **v3 drilling / MACQ8D** | Independent penalty on local `θz` using modified field functions | Default `βd=9×10⁻⁶`; no coupling to in-plane translational gradients |

**Terminology.** The MacNeal–Harder field modification is associated with the 1992 *Eight Nodes or Nine?* paper; the **independent drilling penalty is an added project modification**, not an attributed result of that paper. Its description as “Wilson penalty” in the development file is not treated as proof of a particular published Wilson derivation.

## 2. Kinematics, reference domain and local frame

Let \((r,s)\in[-1,1]^2\) and let the corners be 1–4 in counterclockwise parametric order, with midside nodes 5–8 on edges 1–2, 2–3, 3–4 and 4–1. Each node carries

$$
\mathbf d_i=[u_i,v_i,w_i,\theta_{xi},\theta_{yi},\theta_{zi}]^T,\qquad \mathbf d_e\in\mathbb R^{48}.
$$

The element uses a **single local flat frame** constructed at its center. In the supplied base class, the local normal is derived from the standard Q8 surface tangents and then a global-Z sign preference is applied. The in-plane axes are obtained by orthonormalization. A block-diagonal transformation maps translations and rotations from global to local axes. Here, \(\mathbf I_8\) is the eight-node identity matrix and \(\otimes\) denotes the Kronecker product:

$$
\mathbf d_e^{\ell}=\mathbf T_e\mathbf d_e^g,\qquad
\mathbf T_e=\mathbf I_8\otimes\mathbf T_n,
$$

$$
\mathbf T_n=\begin{bmatrix}
\mathbf R&\mathbf 0_{3\times3}\\
\mathbf 0_{3\times3}&\mathbf R
\end{bmatrix},\qquad
\mathbf R=\begin{bmatrix}
\mathbf e_1^T\\
\mathbf e_2^T\\
\mathbf e_3^T
\end{bmatrix}.
$$

The z-offset of nodes from the flat reference plane is not included in the local in-plane metric, which is an important limitation for strongly warped geometry.


## 3. Standard geometry and modified displacement field

Physical coordinates are mapped by the standard eight-node interpolation:

$$
\mathbf x^h(r,s)=\sum_{i=1}^8 N_i^{(8)}(r,s)\mathbf x_i,
\qquad \mathbf J_0=\begin{bmatrix}x_{,r}&y_{,r}\\x_{,s}&y_{,s}\end{bmatrix}.
$$

Only the *kinematic field* is modified. Let \(N_9=(1-r^2)(1-s^2)\) and define the eight-vector \(\mathbf N^{(8)}(0,0)\). The node-9 location is

$$
\mathbf x_9=\sum_{i=1}^8 N_i^{(8)}(0,0)\mathbf x_i^{\ell}.
$$

Using re-centered local coordinates \(\mathbf x_i^\prime=\mathbf x_i^\ell-\mathbf x_9\), define the interpolation matrix \(\mathbf X\in\mathbb R^{8\times8}\) with node-wise rows

$$
\mathbf X_{i:}=[1,x_i^\prime,y_i^\prime,(x_i^\prime)^2,x_i^\prime y_i^\prime,(y_i^\prime)^2,r_i^2s_i,r_is_i^2].
$$

The eliminated central-node interpolation weights are

$$
\mathbf t^T=\mathbf e_1^T\mathbf X^{-1},\qquad
\widehat N_i(r,s)=N_i^{(8)}(r,s)+N_9(r,s)[t_i-N_i^{(8)}(0,0)].
$$

Thus, derivatives are

$$
\widehat N_{i,r}=N_{i,r}^{(8)}+N_{9,r}[t_i-N_i^{(8)}(0,0)],\quad
\widehat N_{i,s}=N_{i,s}^{(8)}+N_{9,s}[t_i-N_i^{(8)}(0,0)].
$$

The coefficients \(\mathbf t\) are geometry dependent but fixed per element, and \(\det\mathbf X\ne0\) is required. The physical derivatives use the **standard geometry** Jacobian, not a modified geometry Jacobian:

$$
\begin{bmatrix}\widehat N_{i,x}\\\widehat N_{i,y}\end{bmatrix}
=\mathbf J_0^{-1}\begin{bmatrix}\widehat N_{i,r}\\\widehat N_{i,s}\end{bmatrix}.
$$

This separation is essential for reproducing the project formulation.

## 4. Membrane, bending, transverse-shear and drilling operators

With membrane strain vector \(\boldsymbol\varepsilon_m=[\varepsilon_{xx},\varepsilon_{yy},\gamma_{xy}]^T\), the columns belonging to node \(i\) are

$$
\mathbf B_{m,i}=\begin{bmatrix}
\widehat N_{i,x}&0&0&0&0&0\\
0&\widehat N_{i,y}&0&0&0&0\\
\widehat N_{i,y}&\widehat N_{i,x}&0&0&0&0
\end{bmatrix}.
$$

Under the native source's curvature-sign convention and a winding-dependent factor \(\chi_e\in\{-1,+1\}\), the bending operator is

$$
\mathbf B_{b,i}=\chi_e\begin{bmatrix}
0&0&0&0&-\widehat N_{i,x}&0\\
0&0&0&\widehat N_{i,y}&0&0\\
0&0&0&\widehat N_{i,x}&-\widehat N_{i,y}&0
\end{bmatrix}.
$$

**Warning:** the source applies a winding-dependent sign to bending only. A proper orientation-invariance proof requires consistent physical rotation mapping and a corresponding tensor-recovery transformation. The global-Z preference is not frame-independent on closed shells.

The Mindlin transverse-shear operator is

$$
\mathbf B_{s,i}=\begin{bmatrix}
0&0&\widehat N_{i,x}&0&\widehat N_i&0\\
0&0&\widehat N_{i,y}&-\widehat N_i&0&0
\end{bmatrix},\quad
\boldsymbol\gamma_s=\mathbf B_s\mathbf d_e^\ell.
$$

The v3 drilling stabilization uses **only** the local drilling rotation:

$$
\boxed{\mathbf B_{d,i}^{\mathrm{v3}}=[0,0,0,0,0,\widehat N_i]},\qquad
\gamma_d=\sum_i \widehat N_i\theta_{zi}.
$$

Contrast the investigated, but **not selected**, spin-coupled v2-drill mode:

$$
\gamma_d^{\mathrm{spin}}=\theta_z-\frac12(v_{,x}-u_{,y}).
$$

Decoupled drilling removes coupling between drilling rotation and membrane-spin constraints but does **not** guarantee rigid-body rotation invariance under rotation about the element normal. Since a rigid-body normal rotation has nonzero \(\theta_z\), the v3 weak penalty introduces artificial energy into that mode. It should be regarded as a **weak stabilization**, not an exact continuum-consistent drilling field.

## 5. Constitutive law, energy and 3×3 quadrature

For a homogeneous isotropic material,

$$
\mathbf C_m=\frac{E}{1-\nu^2}\begin{bmatrix}1&\nu&0\\\nu&1&0\\0&0&(1-\nu)/2\end{bmatrix},\quad
\mathbf C_s=\kappa_s G\mathbf I_2,\quad G=\frac{E}{2(1+\nu)},\quad\kappa_s=5/6.
$$

The total local stiffness is

$$
\boxed{\mathbf K_e=\int_{A_e}\left[h\mathbf B_m^T\mathbf C_m\mathbf B_m+\frac{h^3}{12}\mathbf B_b^T\mathbf C_m\mathbf B_b+h\mathbf B_s^T\mathbf C_s\mathbf B_s+\beta_d Gh\mathbf B_d^T\mathbf B_d\right]dA}.
$$

The source takes \(h\) as the nodal-thickness mean. The default drilling coefficient is \(\beta_d=9\times10^{-6}\). Numerical integration employs 3×3 Gauss–Legendre points at \(\{-\sqrt{3/5},0,+\sqrt{3/5}\}^2\), with one-dimensional weights \(\{5/9,8/9,5/9\}\). The local integration Jacobian is \(|\det\mathbf J_0|\).

<img width="1045" height="1045" alt="image" src="https://github.com/user-attachments/assets/ba2e92fd-07ff-45bb-942c-9a89ad7fd1b6" />

*Figure 1. Q8 ordering and the nine full-integration points. Geometry and displacement fields are not to be confused.*

## 6. Stress and resultant recovery: primary acceptance gate

From the local displacement solution, recover strains and resultants at **specified points**:

$$
\boldsymbol\varepsilon_m=\mathbf B_m\mathbf d_e^\ell,\qquad
\boldsymbol\kappa=\mathbf B_b\mathbf d_e^\ell,\qquad
\boldsymbol\gamma_s=\mathbf B_s\mathbf d_e^\ell,
$$

$$
\mathbf N=h\mathbf C_m\boldsymbol\varepsilon_m,\quad
\mathbf M=\frac{h^3}{12}\mathbf C_m\boldsymbol\kappa,\quad
\mathbf Q=h\mathbf C_s\boldsymbol\gamma_s.
$$

At thickness coordinate \(z\in[-h/2,h/2]\), the local in-plane fiber stress is

$$
\boldsymbol\sigma_{\parallel}(z)=\mathbf C_m(\boldsymbol\varepsilon_m+z\boldsymbol\kappa)
=\frac{\mathbf N}{h}+\frac{12z}{h^3}\mathbf M.
$$

Stress must be compared in **matching component bases**, thickness sides, and unaveraged Gauss/centroid/corner positions. Global output is obtained through the proper tensor rotation; the local result cannot be directly compared with a differently oriented Nastran element axis. The v3 penalty has no physical stress resultant and should not be included in reported membrane stress. Drilling rank, membrane patch, bending patch and rigid-body rotation tests must be checked separately.

## 7. Numerical verification against independent references

**Scope:**2-001 patch and2-005 mesh convergence only. All numerical results come from existing CSV files; no solver was rerun. Nastran and MYSTRAN MACQ8D are evaluated separately against the same analytic/independent reference:

$$
e_{solver}(\%)=100\frac{|u_{solver}-u_{reference}|}{|u_{reference}|}.
$$

For a zero reference, use an absolute residual in physical units. Bar heights represent errors, not statistical uncertainty intervals. A missing numeric pair remains missing.

### 7.1 Single patch test2-001

The source is `Q8_comparison_master.csv`, filtered to `MACQ8D`, `v3_drill`,2-001, and populated analytic targets. Only finite populated Nastran/MYSTRAN pairs are plotted. The source status remains **PENDING**, so this is measured evidence rather than a certified validation PASS. Stress and resultants retain their own physical units and native recovery locations. Upper/lower fibers use their respective signed targets; each bar is the largest pointwise error for that component and load case.

<img width="2340" height="1620" alt="MACQ8D_2_001_reference_errors" src="https://github.com/user-attachments/assets/14632699-6327-44f0-a3fe-2d7fdbe5b0a5" />


The archived LC2 Nastran stress/moment rows show large errors against their recorded analytic targets (hundreds of percent). These values are retained without correction; source status PENDING and frame/fiber/reference consistency require review before any patch PASS. Small MACQ8D errors do not resolve that provenance issue.

Zero-target components are kept separate and never divided by zero:

<img width="2340" height="1890" alt="MACQ8D_2_001_zero_residuals" src="https://github.com/user-attachments/assets/28cb9e73-7cd1-4cef-a798-b7037279c732" />


Excluded **288** reference-bearing rows with absent/nonfinite solver pairs. Grid stress and fiber-derived grid resultants are recorded separately in [patch metrics](patch2001_reference_metrics.csv), not substituted for native element output.

| LC | Quantity | Component | Reference metric | Nastran | MACQ8D | Points per solver |
|---|---|---|---|---:|---:|---:|
| 1 | resultant | MXX | absolute residual (source_force) | 0 | 0 | 25 |
| 1 | resultant | MXY | absolute residual (source_force) | 0 | 0 | 25 |
| 1 | resultant | MYY | absolute residual (source_force) | 0 | 0 | 25 |
| 1 | resultant | NXX | max pointwise error (%) | 4.488917e-06 | 3.963558e-06 | 25 |
| 1 | resultant | NXY | max pointwise error (%) | 7.274771e-06 | 4.510622e-06 | 25 |
| 1 | resultant | NYY | max pointwise error (%) | 4.488917e-06 | 3.946794e-06 | 25 |
| 1 | resultant | QX | absolute residual (source_force/source_length) | 0 | 0 | 25 |
| 1 | resultant | QY | absolute residual (source_force/source_length) | 0 | 0 | 25 |
| 1 | stress | SXX | max pointwise error (%) | 3.431909e-06 | 3.051758e-06 | 50 |
| 1 | stress | SXY | max pointwise error (%) | 6.159678e-06 | 2.680598e-06 | 50 |
| 1 | stress | SXZ | absolute residual (source_force/source_length^2) | 0 | 0 | 50 |
| 1 | stress | SYY | max pointwise error (%) | 3.990957e-06 | 3.998855e-06 | 50 |
| 1 | stress | SYZ | absolute residual (source_force/source_length^2) | 0 | 0 | 50 |
| 1 | stress | SZZ | absolute residual (source_force/source_length^2) | 0 | 0 | 50 |
| 2 | resultant | MXX | max pointwise error (%) | 304.662 | 4.084285e-06 | 25 |
| 2 | resultant | MXY | max pointwise error (%) | 473.8092 | 5.885143e-06 | 25 |
| 2 | resultant | MYY | max pointwise error (%) | 175.9318 | 5.854351e-06 | 25 |
| 2 | resultant | NXX | absolute residual (source_force/source_length) | 0 | 0 | 25 |
| 2 | resultant | NXY | absolute residual (source_force/source_length) | 0 | 0 | 25 |
| 2 | resultant | NYY | absolute residual (source_force/source_length) | 0 | 0 | 25 |
| 2 | resultant | QX | absolute residual (source_force/source_length) | 0.0004218897 | 6.334189e-16 | 25 |
| 2 | resultant | QY | absolute residual (source_force/source_length) | 0.0001510312 | 1.277137e-15 | 25 |
| 2 | stress | SXX | max pointwise error (%) | 304.662 | 3.963558e-06 | 50 |
| 2 | stress | SXY | max pointwise error (%) | 473.8092 | 4.510622e-06 | 50 |
| 2 | stress | SXZ | absolute residual (source_force/source_length^2) | 0 | 0 | 50 |
| 2 | stress | SYY | max pointwise error (%) | 175.9318 | 3.946794e-06 | 50 |
| 2 | stress | SYZ | absolute residual (source_force/source_length^2) | 0 | 0 | 50 |
| 2 | stress | SZZ | absolute residual (source_force/source_length^2) | 0 | 0 | 50 |

### 7.2 Convergence2-005: corrected C/D plate references

The source is the audited quarter-model `CSI_2005_quarter_20261019/release/paired.csv`, filtered to Q8/MACQ8D. It supersedes historical invalid2-005 master-model comparisons. Four genuine mesh resolutions nx=4,8,16,32 are available. Independent center-displacement magnitudes: simply supported uniform12.97in, simply supported point16.95in, fixed uniform2.60in, fixed point7.25in. Each curve is a solver's absolute percentage error to its independent target, not its difference from the other solver. Native signed displacements are compared by magnitude because the listed targets are magnitudes.


<img width="2160" height="1620" alt="MACQ8D_2_005_reference_convergence" src="https://github.com/user-attachments/assets/1fe9dff8-b342-41e5-ac9e-c283dfb5f5cd" />


[All convergence values and warnings](convergence2005_reference_values.csv). Existing BAILOUT=-1/AUTOSPC qualifications remain; convergence does not erase mechanism or output-recovery qualifications. No unrelated benchmark results are inserted.


## 8. Verification hierarchy and trade-offs

- **Constant-strain and constant-curvature patches:** certify pointwise stress/moment, not merely displacement. Compare CW/CCW in the same global basis.
- **Drilling stabilization:** compare rank against the free 6-RBM expectation, and test a pure rigid-body normal rotation for penalty-induced residual. A very small penalty can control singular modes while still violating exact normal-spin invariance.
- **Distorted elements:** the field modification improves certain patches yet can be ill-conditioned when the \(8\times8\) polynomial interpolation matrix is nearly singular. Report its condition number.
- **Twisted/warped geometry:** single flat frame cannot reproduce pointwise curved director kinematics; this is a structural limitation rather than a recoverable sign bug.
- **Pressure:** use corrected C/D quarter-model convergence data in Section7; old invalid2-005 model comparisons are excluded.
- **Mesh convergence:** Section7 uses the same corrected2-005 physical cases at four mesh resolutions; each solver is compared against the independent target.

## 9. Reproduction and provenance

Numerical sources are the existing Q8 production master and corrected quarter-model comparison CSV. Their SHA256 hashes are listed in `data_provenance.txt`. The accompanying metrics preserve missing-data exclusions, physical units and warnings. No implementation source is included. Exact implementation parity requires a separate kernel audit. The supplied `macneal1992.pdf` is a bibliographic reference only and supplies no numerical results here.

## References

1. R. H. MacNeal and R. L. Harder (1992), “Eight Nodes or Nine?”, *International Journal for Numerical Methods in Engineering*, **33**, 1049–1058. [Publication identification](https://ntrs.nasa.gov/api/citations/19980019011/downloads/19980019011.pdf).
2. R. H. MacNeal and R. L. Harder (1985), “A proposed standard set of problems to test finite element accuracy,” *Finite Elements in Analysis and Design*, **1**, 3–20. DOI: [10.1016/0168-874X(85)90003-4](https://doi.org/10.1016/0168-874X(85)90003-4).
3. R. H. MacNeal (1992), “On the limits of finite element perfectability,” *International Journal for Numerical Methods in Engineering*, **35**, 1589–1601. DOI: [10.1002/nme.1620350804](https://doi.org/10.1002/nme.1620350804).
4. R. H. MacNeal, C. T. Wilson, R. L. Harder, and C. C. Hoff (1998), “The treatment of shell normals in finite element analysis,” *Finite Elements in Analysis and Design*, **30**, 235–242. DOI: [10.1016/S0168-874X(98)00035-3](https://doi.org/10.1016/S0168-874X(98)00035-3).
5. Computers and Structures, Inc., *SAP2000 Verification: Shell Analysis*, [Verification collection](https://docs.csiamerica.com/manuals/sap2000/Verification/Analysis/Shells/).
6. MYSTRAN development benchmarks, [`update/shell_benchmark`](https://github.com/realbabilu/MYSTRAN/tree/v18.00.a/update/shell_benchmark); supporting archived user dataset `Q8_CSV_delivery_20261018(1).zip`.

**Supplied bibliographic copy:** `D:/Downloads/macneal1992.pdf` (reference only; numerical results exclusively from CSV).
