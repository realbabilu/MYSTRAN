# CPENTA6 ABE: edge-tied assumed shear and condensed enhanced-strain wedge

**Technical formulation and numerical-validation note — MYSTRAN NEWSOLID (`PENTATYP,ABE`)**  
**Revision:** 19 October 2026 (GitHub math syntax corrected) | **Origin:** Ary, [project formulation note](https://github.com/realbabilu/MYSTRAN/blob/v18.00.a/update/Formulation/Solid/CPENTA6_ABE_Formulation_Ary.md)  
**Evidence status:** Figures and numerical tables reproduce reported Python-reference data from the linked original document. This revision does **not** claim new finite-element runs or a newly verified Fortran port.

## Abstract

ABE is an 18-DOF, six-node continuum prism using an edge-tied assumed transverse-shear operator and nine internal enhanced assumed strain (EAS) parameters grouped as **A (1), B (2), E (6)**. The internal parameters are statically condensed; **only three global translations per node** are retained. The method is a **project-specific combination**, not a separately published element. For SAP2000 verification 5-010, the reported `patch_consistency=0` buckling factor changes from **0.8912** at a 2×1×12 mesh to **0.8885** at 8×1×48, near Euler's **0.8883**. The same setting, however, fails the distorted constant-stress patch test 5-001; `patch_consistency=1` restores patch consistency but reintroduces thin-bending locking. The two settings must not be conflated.

## 1. Domain, interpolation, and degrees of freedom

The reference wedge is

$$
\widehat\Omega=\{(r,s,\zeta):r\ge0,\ s\ge0,\ r+s\le1,\ -1\le\zeta\le1\}.
$$

Set $L_1=1-r-s$, $L_2=r$, and $L_3=s$. The six interpolation functions are

$$
N_a(r,s,\zeta)=\frac{1-\zeta}{2}L_a(r,s),\quad
N_{a+3}(r,s,\zeta)=\frac{1+\zeta}{2}L_a(r,s),\quad a=1,2,3.
$$

The isoparametric geometry and displacement are

$$
\mathbf x(r,s,\zeta)=\sum_{a=1}^{6}N_a\mathbf x_a,
\qquad \mathbf u^h(r,s,\zeta)=\mathbf N\mathbf d,
\qquad \mathbf d\in\mathbb R^{18},
$$

$$
\mathbf N=\begin{bmatrix}N_1\mathbf I_3&N_2\mathbf I_3&\cdots&N_6\mathbf I_3\end{bmatrix}.
$$

The bottom nodes are (1,2,3) at $\zeta=-1$ and their corresponding upper nodes are (4,5,6) at $\zeta=+1$. The Jacobian convention in this report uses **rows** as physical coordinate derivatives:

$$
\mathbf J=\frac{\partial\mathbf x}{\partial(r,s,\zeta)}=
\begin{bmatrix}\mathbf x_{,r}^{T}\\\mathbf x_{,s}^{T}\\\mathbf x_{,\zeta}^{T}\end{bmatrix},
\qquad dV=\det\mathbf J\,dr\,ds\,d\zeta.
$$

Use an orientation with $\det\mathbf J>0$ at all quadrature points. The spatial gradient is $\nabla_{\mathbf x}N_a=\mathbf J^{-1}\nabla_{\xi}N_a$ with the corresponding row/column convention; equivalently, $\nabla_{\xi}N_a=\mathbf J\nabla_{\mathbf x}N_a$.

Let $\mathbf B$ denote the usual $6\times18$ engineering-strain operator in ordering

$$
\boldsymbol\varepsilon_V=[\varepsilon_{xx},\varepsilon_{yy},\varepsilon_{zz},\gamma_{xy},\gamma_{yz},\gamma_{zx}]^T,\qquad
\boldsymbol\varepsilon_V=\mathbf B\mathbf d.
$$

For isotropic linear elasticity, the constitutive relationship in engineering-shear Voigt notation is

$$
\boldsymbol{\sigma}_v = \mathbf{D}\boldsymbol{\varepsilon}_v
$$

where

$$
\mathbf{D} =
\begin{bmatrix}
\lambda+2\mu & \lambda & \lambda & 0 & 0 & 0 \\
\lambda & \lambda+2\mu & \lambda & 0 & 0 & 0 \\
\lambda & \lambda & \lambda+2\mu & 0 & 0 & 0 \\
0 & 0 & 0 & \mu & 0 & 0 \\
0 & 0 & 0 & 0 & \mu & 0 \\
0 & 0 & 0 & 0 & 0 & \mu
\end{bmatrix}
$$

The Lamé parameters are

$$
\lambda = \frac{E\nu}{(1+\nu)(1-2\nu)},
\qquad
\mu = \frac{E}{2(1+\nu)}.
$$
  

## 2. Frozen centroid covariant transformation

At $\boldsymbol\xi_0=(1/3,1/3,0)$, evaluate $\mathbf J_0$. A physical symmetric strain tensor $\boldsymbol\varepsilon$ maps to centroid-frame covariant components through

$$
\widehat{\boldsymbol\varepsilon}=\mathbf J_0\boldsymbol\varepsilon\mathbf J_0^{T},\qquad
\widehat{\boldsymbol\varepsilon}_V=
[\widehat\varepsilon_{rr},\widehat\varepsilon_{ss},\widehat\varepsilon_{\zeta\zeta},
2\widehat\varepsilon_{rs},2\widehat\varepsilon_{s\zeta},2\widehat\varepsilon_{r\zeta}]^T.
$$

Define $\mathbf T_0\in\mathbb R^{6\times6}$ by the tensor relation above, **including engineering shear factors**, then

$$
\widehat{\mathbf B}(\boldsymbol\xi)=\mathbf T_0\mathbf B(\boldsymbol\xi),\quad
\widehat{\mathbf D}_0=\mathbf T_0^{-T}\mathbf D\mathbf T_0^{-1},\quad
\boldsymbol\varepsilon_V=\mathbf T_0^{-1}\widehat{\boldsymbol\varepsilon}_V.
$$

Using the *constant* $\mathbf T_0$ in this transformation avoids inadvertently substituting a different local frame at each integration point. The mechanical energy is unchanged:

$$
\boldsymbol\varepsilon_V^T\mathbf D\boldsymbol\varepsilon_V
=\widehat{\boldsymbol\varepsilon}_V^T\widehat{\mathbf D}_0\widehat{\boldsymbol\varepsilon}_V.
$$

## 3. Edge-tied assumed transverse shear

For each quadrature plane $\zeta=\zeta_g$, evaluate the covariant compatible shear rows at triangle-edge midpoints:


$$
\begin{aligned}
\mathbf a(\zeta_g)&=\widehat{\mathbf B}_{r\zeta}(1/2,0,\zeta_g),\\
\mathbf b(\zeta_g)&=\widehat{\mathbf B}_{s\zeta}(0,1/2,\zeta_g),\\
\mathbf c(\zeta_g)&=\mathbf b-\mathbf a-
\left[\widehat{\mathbf B}_{s\zeta}-\widehat{\mathbf B}_{r\zeta}\right](1/2,1/2,\zeta_g).
\end{aligned}
$$


Here $\mathbf a$, $\mathbf b$, and $\mathbf c$ are **row vectors of length 18**; $\widehat{\mathbf B}_{r\zeta}$ and $\widehat{\mathbf B}_{s\zeta}$ denote rows 6 and 5 of $\widehat{\mathbf B}$, respectively. Interpolate them as

$$
\boxed{\overline{\mathbf B}_{r\zeta}(r,s,\zeta_g)=\mathbf a+\mathbf c\,s},
\qquad
\boxed{\overline{\mathbf B}_{s\zeta}(r,s,\zeta_g)=\mathbf b-\mathbf c\,r}.
$$

The other four rows remain compatible. Define $\overline{\mathbf B}$ as this six-row natural-frame matrix. The tied operator is only asserted at the quadrature planes used by this implementation; no higher-order $\zeta$ reconstruction is assumed.

The edge-tied construction follows the *isotropic MITC3-type* shear interpolation discussed by Lee and Bathe (2004), adapted here to the continuum wedge. It is **not** identical to a shell MITC3 element; solid normal strains and the nine local EAS parameters remain present.



For each quadrature plane $\zeta=\zeta_g$, evaluate the covariant compatible shear rows at triangle-edge midpoints.

```math
\mathbf{c}(\zeta_g)
=
\mathbf{b}(\zeta_g)
-
\mathbf{a}(\zeta_g)
-
\left[
\widehat{\mathbf{B}}_{s\zeta}
-
\widehat{\mathbf{B}}_{r\zeta}
\right]
\left(\frac{1}{2},\frac{1}{2},\zeta_g\right)
```

## 4. Nine local enhanced assumed strain parameters

With $\bar r=r-1/3$, $\bar s=s-1/3$, the natural strain is augmented as

$$
\boxed{\widehat{\boldsymbol\varepsilon}_V
=\overline{\mathbf B}\mathbf d+\mathbf G(\boldsymbol\xi)\boldsymbol\alpha},
\qquad
\mathbf G=\frac{\det\mathbf J_0}{\det\mathbf J(\boldsymbol\xi)}\mathbf M_{\xi},
\qquad \boldsymbol\alpha\in\mathbb R^9.
$$

Let $\mathbf e_i$ denote the $i$-th unit column in natural engineering-strain Voigt space. **Exact column definition and ordering:**

$$
\boxed{\mathbf M_{\xi}=\begin{bmatrix}
\zeta\mathbf e_3 &
\zeta\mathbf e_5 & \zeta\mathbf e_6 &
\bar r\mathbf e_1 & \bar s\mathbf e_1 &
\bar r\mathbf e_2 & \bar s\mathbf e_2 &
\bar r\mathbf e_4 & \bar s\mathbf e_4
\end{bmatrix}}.
$$

| Family | Parameters | Natural components | Interpretation |
|---|---:|---|---|
| **A** | 1 | $\zeta\mathbf e_3$ | Linear thickness-normal enhancement |
| **B** | 2 | $\zeta\mathbf e_5$, $\zeta\mathbf e_6$ | Transverse-shear enhancements |
| **E** | 6 | $(\bar r,\bar s)\otimes(\mathbf e_1,\mathbf e_2,\mathbf e_4)$ | In-plane extension/shear enhancements |

**Orthogonality.** Since every reference mode has zero volume integral,

$$
\int_{\Omega_e}\mathbf G\,dV
=\int_{\widehat\Omega}\frac{\det\mathbf J_0}{\det\mathbf J}\mathbf M_\xi
\det\mathbf J\,d\widehat V
=\det\mathbf J_0\int_{\widehat\Omega}\mathbf M_\xi\,d\widehat V=\mathbf0.
$$

This ensures the EAS modes alone are orthogonal to any *constant stress represented in the frozen natural frame*. It **does not** prove that the edge-tied displacement operator passes the patch test on distorted geometry (see Section 7).

### 4.1 Mixed stationarity and static condensation

For the linear strain-enhanced element, the local internal potential is

$$
\Pi_e(\mathbf d,\boldsymbol\alpha)
=\frac12\int_{\Omega_e}
(\overline{\mathbf B}\mathbf d+\mathbf G\boldsymbol\alpha)^T
\widehat{\mathbf D}_0
(\overline{\mathbf B}\mathbf d+\mathbf G\boldsymbol\alpha)\,dV
-\mathbf d^T\mathbf f_e.
$$

Its block Hessian is

$$
\begin{bmatrix}\mathbf K_{dd}&\mathbf K_{d\alpha}\\
\mathbf K_{\alpha d}&\mathbf K_{\alpha\alpha}\end{bmatrix}
=\int_{\Omega_e}
\begin{bmatrix}\overline{\mathbf B}^T\\\mathbf G^T\end{bmatrix}
\widehat{\mathbf D}_0
\begin{bmatrix}\overline{\mathbf B}&\mathbf G\end{bmatrix}\,dV.
$$

The enhanced stationarity equation reads

$$
\mathbf K_{\alpha d}\mathbf d+\mathbf K_{\alpha\alpha}\boldsymbol\alpha=\mathbf0,
\quad \boldsymbol\alpha=-\mathbf K_{\alpha\alpha}^{-1}\mathbf K_{\alpha d}\mathbf d.
$$

Hence the **18 × 18 condensed element stiffness** is

$$
\boxed{\mathbf K_e=\mathbf K_{dd}-
\mathbf K_{d\alpha}\mathbf K_{\alpha\alpha}^{-1}\mathbf K_{\alpha d}}.
$$

In code, **solve** $\mathbf K_{\alpha\alpha}\mathbf X=\mathbf K_{\alpha d}$ by a factorization rather than forming an inverse; verify its conditioning and reject/diagnose singular local enhanced systems. No new global DOFs are introduced.

## 5. Patch-consistency mean correction

Define the shear-only row difference $\Delta\mathbf B=\widehat{\mathbf B}-\overline{\mathbf B}$, where all non-transverse-shear rows of $\Delta\mathbf B$ are zero. Its volume mean is

$$
\langle\Delta\mathbf B\rangle
=\frac{1}{V_e}\int_{\Omega_e}(\widehat{\mathbf B}-\overline{\mathbf B})\,dV,
\qquad V_e=\int_{\Omega_e}dV.
$$

The corrected assumed operator for
$$
$f=\texttt{patch\_consistency}\in[0,1]$ 
$$
is

$$
\boxed{\overline{\mathbf B}_f=\overline{\mathbf B}
+f\langle\Delta\mathbf B\rangle}.
$$

The correction is applied **only to the two transverse-shear rows** and is constant within the element. At $f=1$,

$$
\int_{\Omega_e}(\overline{\mathbf B}_1-\widehat{\mathbf B})\,dV=\mathbf0.
$$

Combined with EAS zero-mean orthogonality, this provides the constant-stress consistency condition for the reported linear patch tests; the actual proof also assumes a valid mesh, exact Dirichlet setup, and correct element integration. At $f=0$, no correction is applied and **distorted patch consistency is not satisfied** in the reported tests. Neither $f=0$ nor $f=1$ is universally superior.

## 6. Quadrature, stress recovery, and other operators

Use the tensor product of a three-point triangular rule and two-point Gauss–Legendre integration through $\zeta$:

| Point | $r$ | $s$ | $\zeta$ | Reference weight |
|---|---:|---:|---:|---:|
| 1 | 1/6 | 1/6 | $-1/\sqrt3$ | 1/6 |
| 2 | 2/3 | 1/6 | $-1/\sqrt3$ | 1/6 |
| 3 | 1/6 | 2/3 | $-1/\sqrt3$ | 1/6 |
| 4 | 1/6 | 1/6 | $+1/\sqrt3$ | 1/6 |
| 5 | 2/3 | 1/6 | $+1/\sqrt3$ | 1/6 |
| 6 | 1/6 | 2/3 | $+1/\sqrt3$ | 1/6 |

Thus $\sum_g w_g=1$ (the reference prism has unit volume) and numerical integrals read

$$
\int_{\Omega_e}\mathbf H\,dV\approx
\sum_{g=1}^{6}w_g\,\det\mathbf J_g\,\mathbf H_g.
$$

After the displacement solution, recover each element's enhanced coordinates from the same operator used for $\mathbf K_e$:

$$
\boldsymbol\alpha_e=-\mathbf K_{\alpha\alpha}^{-1}\mathbf K_{\alpha d}\mathbf d_e,
\quad
\boldsymbol\varepsilon_V(\boldsymbol\xi_g)
=\mathbf T_0^{-1}[\overline{\mathbf B}_f(\boldsymbol\xi_g)\mathbf d_e
+\mathbf G(\boldsymbol\xi_g)\boldsymbol\alpha_e],
\quad
\boldsymbol\sigma_V=\mathbf D\boldsymbol\varepsilon_V.
$$

The original Python study uses standard CPENTA6 operators for the consistent/lumped mass and initial-stress geometric stiffness. This **is not** a consistent nonlinear EAS tangent; 5-010 validation must be interpreted as an eigen-buckling benchmark with the reported prestress assumptions. For an initial Cauchy prestress $\boldsymbol\sigma_0$, the usual small-strain geometric operator has component blocks

$$
(K_G)_{ai,bj}=\delta_{ij}\int_{\Omega_e}
(\nabla N_a)^T\boldsymbol\sigma_0(\nabla N_b)\,dV,
$$

subject to the precise signs and conventions of the MYSTRAN buckling eigenproblem.

## 7. Numerical verification and convergence (reported Python data)

The **mesh labels and ratios are reproduced from the original note**. Numerical data tables are retained below; convergence plots will be inserted after the illustrations are supplied. No new numerical solutions or Nastran runs are claimed.

### 7.1 SAP2000 5-010: first buckling factor

> **Figure reserved:** 5-010 buckling-factor convergence. Illustration will be inserted later.

**Planned Figure 1.** Connected-point convergence of the first buckling factor, using a logarithmic vertical axis.

| Mesh | Standard | EAS9 v2 | ABE $f=1$ | ABE $f=0$ |
|---|---:|---:|---:|---:|
| 2×1×12 | 20.980 | 7.680 | 3.888 | 0.8912 |
| 4×1×24 | 6.258 | 3.414 | 2.216 | 0.8890 |
| 8×1×48 | 2.388 | 1.714 | 1.322 | 0.8885 |

**Reference:** Euler buckling factor 0.8883; the original note also reports Nastran CPENTA6 factor 0.8910 for its audited mesh. These are **different kinds of reference**: Euler is a theoretical asymptote, while Nastran is a discrete FE result and should not be represented as a universal exact solution.

> **Figure reserved:** 5-010 relative buckling error. Illustration will be inserted later.

**Planned Figure 2.** Absolute relative error $100|\lambda_h/0.8883-1|$ (%) against the Euler factor on the three reported mesh levels. The logarithmic ordinate highlights the ABE $f=0$ result. Its relative errors are approximately **0.326%, 0.079%, 0.023%**.

### 7.2 SAP2000 5-001: distorted constant-stress patch

The reported model consists of seven irregular solids (14 wedges), $E=10^6$, $\nu=0.25$, with prescribed exact displacement

$$
\mathbf u(\mathbf x)=\frac1{2000}
\begin{bmatrix}2x+y+z\\x+2y+z\\x+y+2z\end{bmatrix}.
$$

This is an affine field, giving the uniform exact strain

$$
\boldsymbol\varepsilon_V^{\mathrm{exact}}=
\frac1{2000}[2,2,2,2,2,2]^T.
$$

| Formulation | Relative displacement error | Relative stress error |
|---|---:|---:|
| Standard CPENTA6 | $4.3\times10^{-16}$ | $1.8\times10^{-15}$ |
| EAS9 v2 | $4.3\times10^{-16}$ | $2.0\times10^{-15}$ |
| ABE $f=1$ | $2.0\times10^{-15}$ | $8.3\times10^{-14}$ |
| ABE $f=0.5$ | 0.059 | 0.19 |
| ABE $f=0.25$ | 0.079 | 0.27 |
| ABE $f=0$ | 0.091 | 0.33 |

> **Figure reserved:** 5-001 patch-consistency trade-off. Illustration will be inserted later.

**Planned Figure 3.** Observed patch-test accuracy versus mean-restoration parameter $f$. Values at $f=1$ are near machine precision; the graph uses a linear ordinate. **These data show that the high buckling accuracy at $f=0$ does not imply general patch consistency.**

### 7.3 Distorted patch refinement: nonzero consistency-error floor

> **Figure reserved:** distorted patch mesh-refinement. Illustration will be inserted later.

**Planned Figure 4.** Reported error under 10% relative interior-node distortion. The error changes with $n$, but does **not** tend to numerical zero over the tested mesh sequence. This is a consistency warning, not a definitive asymptotic-order estimate.

| Cube grid $n\times n\times n$ | $n=2$ | $n=3$ | $n=4$ | $n=6$ |
|---|---:|---:|---:|---:|
| ABE $f=0$ patch error | 0.0008 | 0.0027 | 0.0034 | 0.0033 |

At $n=6$, the reported errors under 10%, 20%, and 30% perturbations are 0.0033, 0.0072, and 0.012, respectively. Affine mapped/regular meshes passed near machine precision. The metric is reproduced as reported; its exact norm definition should be included in the benchmark driver before publication as an independent reproducible result.

### 7.4 Additional reported benchmarks

Values are ratios to the stated independent benchmark, and `IN / OUT` represent different response directions.

| Benchmark | EAS9 v2 | ABE $f=1$ | ABE $f=0$ | CPENTA15 |
|---|---|---|---|---|
| 5-003 curved beam A (6×1×1) | .821 / .785 | .850 / .758 | .851 / .757 | .859 / .926 |
| 5-003 curved beam B (90×1×1) | .940 / .895 | .978 / .874 | .978 / .874 | .998 / .954 |
| 5-004 twist 12×2×2 | .952 / .951 | .960 / .963 | .960 / .963 | .997 / 1.004 |
| 5-004 twist 24×4×2 | .960 / .963 | .962 / .969 | .962 / .969 | .998 / 1.003 |
| 5-008 cylinder $\nu=.3/.49/.499/.4999$ | .972/.976/.980/.981 | .971/.976/.980/.981 | .971/.976/.980/.981 | 1.013/1.003/.984/.875 |
| Cantilever 12×2×2 (bending/frequency/buckling) | .946/1.023/1.044 | .961/1.017/1.035 | .961/1.017/1.035 | .994/1.001/1.001 |

### 7.5 Rank and invariance gates

For an unconstrained linear elastic three-dimensional element with 18 DOFs, the expected unconstrained stiffness rank is **12** (six rigid-body modes). The source reports ABE rank 12 for regular and distorted single-element cases, a rotation residual near $4\times10^{-16}$, and cyclic node-renumbering differences $\le10^{-13}$. Its EAS9 v2 comparator differed by roughly 3% under two cyclic permutations. For reproducibility, always test rotated coordinates **and** cyclic renumbering with matching DOF permutations, and report the eigenvalue threshold used to decide numerical rank.

## 8. Orientation and known limitations

The original note's reported 2×1×12 buckling factors depend strongly on the prism-axis orientation:

| Prism axis | Standard | EAS9 v2 | ABE $f=1$ | ABE $f=0$ |
|---|---:|---:|---:|---:|
| $y$ (triangles in $xz$) | 20.98 | 7.680 | 3.888 | 0.8912 |
| $z$ (longitudinal axis) | 26.71 | 0.9456 | 0.9127 | 0.9127 |
| $x$ (triangles in bending plane) | 79.84 | 79.79 | 79.79 | 79.79 |

The $x$ orientation locks badly in all these wedge variants. ABE **cannot** reconstruct in-plane curvature that the one-layer linear triangle lacks. Strong warpage, highly distorted elements, nonlinear buckling tangents, complete stress-recovery validation, thermal benchmark 5-011, and broader reliability under arbitrary meshing remain open. The original note reports mass and $K_G$ inherited from standard CPENTA6 rather than fully re-derived from an enhanced nonlinear variational formulation.

**Implementation caution:** `patch_consistency=0` is a research setting with excellent *reported* thin-column buckling but a demonstrable patch-test failure. Treat its use on general distorted meshes as experimental. `patch_consistency=1` passes the reported patch but significantly overpredicts thin-column buckling. A practical MYSTRAN release should expose this trade-off explicitly and should not silently substitute either option for standard CPENTA6 across all engineering use cases.

## 9. Replication checklist

1. Rebuild $\mathbf T_0$ from the declared Jacobian row convention and test energy invariance.
2. Check quadrature weights and $\det\mathbf J_g>0$, local $\mathbf K_{\alpha\alpha}$ conditioning and rigid-body rank.
3. Verify shear row indices: $r\zeta\leftrightarrow6$, $s\zeta\leftrightarrow5$ in the declared order.
4. Use the **same** corrected $\overline{\mathbf B}_f$ for stiffness, condensation, and stress recovery.
5. Verify regular/affine and distorted patch tests separately, with $f=0$ and $f=1$.
6. Test 5-003 bending, 5-004 twisting, 5-008 near-incompressibility, and 5-010 buckling on a consistent mesh sequence.
7. Compare against independent SAP2000 results and Nastran results at the **same physical locations, load cases, coordinate system, and units**.
8. Archive mesh, displacements, stress extraction choices, eigenvalue/sign convention, solver version, and CSV results.

## References

1. P. S. Lee and K. J. Bathe (2004), “Development of MITC isotropic triangular shell finite elements,” *Computers & Structures* **82**, 945–962, DOI: [10.1016/j.compstruc.2004.02.004](https://doi.org/10.1016/j.compstruc.2004.02.004).
2. U. Andelfinger and E. Ramm (1993), “EAS-elements for two-dimensional, three-dimensional, plate and shell structures and their equivalence to HR-elements,” *International Journal for Numerical Methods in Engineering* **36**, 1311–1337.
3. J. C. Simo and M. S. Rifai (1990), “A class of mixed assumed strain methods and the method of incompatible modes,” *International Journal for Numerical Methods in Engineering* **29**, 1595–1638.
4. J. C. Simo and F. Armero (1992), “Geometrically non-linear enhanced strain mixed methods and the method of incompatible modes,” *International Journal for Numerical Methods in Engineering* **33**, 1413–1449.
5. R. Hauptmann and K. Schweizerhof (1998), “A systematic development of solid-shell element formulations for linear and non-linear analyses employing only displacement degrees of freedom,” *International Journal for Numerical Methods in Engineering* **42**, 49–69.
6. R. H. MacNeal and R. L. Harder (1985), “A proposed standard set of problems to test finite element accuracy,” *Finite Elements in Analysis and Design* **1**, 3–20.
7. Computers and Structures, Inc. (CSI), *SAP2000 Verification: Analysis / Solids*, [5-001](https://docs.csiamerica.com/manuals/sap2000/Verification/Analysis/Solids/Problem%205-001.pdf), [5-003](https://docs.csiamerica.com/manuals/sap2000/Verification/Analysis/Solids/Problem%205-003.pdf), [5-004](https://docs.csiamerica.com/manuals/sap2000/Verification/Analysis/Solids/Problem%205-004.pdf), [5-008](https://docs.csiamerica.com/manuals/sap2000/Verification/Analysis/Solids/Problem%205-008.pdf), [5-010](https://docs.csiamerica.com/manuals/sap2000/Verification/Analysis/Solids/Problem%205-010.pdf).
8. Original project report: [CPENTA6_ABE_Formulation_Ary.md (v18.00.a)](https://github.com/realbabilu/MYSTRAN/blob/v18.00.a/update/Formulation/Solid/CPENTA6_ABE_Formulation_Ary.md). Also consult its further references to Taylor–Beresford–Wilson (1976), Wilson et al. (1973), Dvorkin–Bathe (1984), Bathe–Dvorkin (1986), and Büchter–Ramm (1992).
