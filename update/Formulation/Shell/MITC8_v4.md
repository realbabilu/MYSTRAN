# An Enhanced Eight-Node Mixed-Interpolation Shell Element with Curvature-Adaptive Membrane Tying and Slenderness-Dependent Transverse-Shear Blending: Formulation and Verification of MITC8 v4

**Ary · Independent Researcher · Personal Research** · MYSTRAN v18.00.a · 19 October 2026
**Technical formulation and verification report**  
**Element family:** CQUAD8, eight-node degenerated curved shell  
**Revision:** v4 | 19 October 2026  

## Abstract

An enhanced eight-node shell formulation, here designated **MITC8 v4**, is presented for thin and moderately thick shell analysis. The element employs a three-dimensional isoparametric midsurface, an independent interpolated director field, and six degrees of freedom per node. The formulation combines a curvature-activated covariant membrane assumed-strain operator with a thickness-sensitive mixture of two transverse-shear tying operators. The latter addresses the artificial stiffening observed for a six-point-per-component shear interpolation in thin plates without completely sacrificing the coarse-mesh resistance supplied by that interpolation. Additional provisions include drilling stabilisation, consistent and positive diagonal mass, invariant pressure integration, membrane and flexural thermal loads, and translational initial-stress stiffness. Numerical evidence is organised into membrane/bending patch tests, curved beams, twisted members, plate-loading convergence, free thermal expansion, modal response, and linearised bifurcation of plates and columns. On a 16 × 16 adaptation of verification problem 2-016, the signed critical-stress errors decrease from +8.073% to −0.182% for shear and from +7.281% to −0.324% for bending when moving from v3 to v4. The precise family name does **not** imply identical reproduction of the classical Bucalem–Bathe MITC8 interpolation: the present shear-tying rules and their slenderness-dependent combination are a research modification. Applicability to arbitrary curved-shell buckling, nonlinear rotations, and anisotropic laminates remains open.

**Keywords:** eight-node shell; mixed interpolation; assumed natural strain; transverse shear locking; membrane tying; drilling degree of freedom; geometric stiffness; shell verification.

## 1. Introduction

Low-order and quadratic displacement-based Reissner–Mindlin shell elements may exhibit parasitic transverse-shear energy as the thickness-to-span ratio approaches zero. On distorted or doubly curved geometries, inconsistent representation of the covariant membrane strain can produce a second source of artificial stiffness. The present development retains the eight-node serendipity geometry and the 48 nodal degrees of freedom while altering the projection of membrane and transverse-shear strain fields. The intended improvements are: (i) membrane patch consistency on flat straight-sided elements; (ii) improved coarse-mesh curved response through selective membrane tying; (iii) reduced thin-plate shear locking; and (iv) a consistent set of operators for stiffness, internal-force recovery, pressure, thermal loading, mass and linear buckling.

The designation **MITC8 v4** is used to identify the tested variant, not to claim exact identity with every operation in the canonical MITC8 construction of Bucalem and Bathe. In particular, the six-point/two-direction `tensor6` operator, the four-point `bdg4` operator and their geometry-dependent mixture should be considered a modified covariant tying strategy.

## 2. Kinematic description and interpolation

### 2.1 Eight-node reference cell and degrees of freedom

Let $(r,s)\in[-1,1]^2$. The corner nodes occupy $(-1,-1),(1,-1),(1,1),(-1,1)$; midside nodes occupy $(0,-1),(1,0),(0,1),(-1,0)$. Their standard geometry functions are


```math
N_i=\tfrac14(1+r r_i)(1+s s_i)(r r_i+s s_i-1),\qquad i=1,\ldots,4,
```



```math
N_5=\tfrac12(1-r^2)(1-s),\quad N_6=\tfrac12(1+r)(1-s^2),
\quad N_7=\tfrac12(1-r^2)(1+s),\quad N_8=\tfrac12(1-r)(1-s^2).
```


Each node carries $\mathbf d_I=[u_x,u_y,u_z,\theta_x,\theta_y,\theta_z]_I^T$, and $\mathbf d_e\in\mathbb R^{48}$. Nodal translations and rotations are expressed in a common global Cartesian system; no element-wide frozen planar basis is imposed.

<img width="1448" height="1086" alt="image" src="https://github.com/user-attachments/assets/9b2eb1bb-0493-4e59-b1ca-7ac56da4c561" />

### 2.2 Midsurface, covariant basis and dual basis

The geometric midsurface is


```math
\mathbf x(r,s)=\sum_{I=1}^{8}N_I(r,s)\mathbf X_I,\qquad
\mathbf g_1=\mathbf x_{,r},\quad\mathbf g_2=\mathbf x_{,s},
```



```math
J_A=\|\mathbf g_1\times\mathbf g_2\|,\qquad
\mathbf n=\frac{\mathbf g_1\times\mathbf g_2}{J_A},\qquad dA=J_A\,dr\,ds.
```


With $g_{\alpha\beta}=\mathbf g_\alpha\cdot\mathbf g_\beta$, define $\mathbf g^\alpha=g^{\alpha\beta}\mathbf g_\beta$, where $[g^{\alpha\beta}]=[g_{\alpha\beta}]^{-1}$. At each integration point, $\mathbf e_1,\mathbf e_2,\mathbf n$ form an orthonormal tangent-normal triad. The two-by-two matrix $\mathbf C$ holds the covariant components of the orthonormal tangents; the strain transformation therefore changes over a curved element.

The reference surface functions $N_I$ determine geometry, while the active displacement/director interpolation uses the adopted quadratic assumed-field variant (Kikuchi by default, with a MacNeal–Harder alternative). No general identity between these field functions and standard serendipity functions is assumed.

### 2.3 Reissner–Mindlin displacement field

A convenient first-order approximation is


```math
\mathbf u(r,s,z)=\mathbf u_0(r,s)+z\,\boldsymbol\vartheta(r,s)\times\mathbf n_0(r,s),
\quad \mathbf u_0=\sum_I\widehat N_I\mathbf u_I,
```


where $\mathbf n_0$ is obtained from consistently oriented nodal directors and $\widehat N_I$ denote the active field-interpolation functions. This representation identifies the kinematic model; the discrete operators also include the derivatives of interpolated nodal directors on curved surfaces. The sign of the director cross product, bending curvature and physical bending-moment recovery must be kept consistent.

## 3. Generalised strain and elastic stiffness

Write the generalised mechanical strain as


```math
\boldsymbol\epsilon_m=\mathbf B_m\mathbf d_e,\qquad
\boldsymbol\kappa=\mathbf B_b\mathbf d_e,\qquad
\boldsymbol\gamma=\mathbf B_s\mathbf d_e,\qquad
\epsilon_d=\mathbf B_d\mathbf d_e.
```


For an isotropic material under plane stress,


```math
\mathbf C_m=\frac{E}{1-\nu^2}
\begin{bmatrix}1&\nu&0\\\nu&1&0\\0&0&(1-\nu)/2\end{bmatrix},
\quad
\mathbf C_s=\kappa_s G\mathbf I_2,\quad G=\frac{E}{2(1+\nu)}.
```


The shell tangent used in the reported linear calculations is


```math
\boxed{\mathbf K_e=\int_{A_e}\left[h\mathbf B_m^T\mathbf C_m\mathbf B_m+
\frac{h^3}{12}\mathbf B_b^T\mathbf C_m\mathbf B_b+
 h\mathbf B_s^T\mathbf C_s\mathbf B_s+
 k_d\mathbf B_d^T\mathbf B_d\right]dA.}
```


The drilling term is a numerical stabilisation; its exact dimensional coefficient is defined by the adopted penalty convention rather than by independent physical drilling rigidity. The default dimensionless drill-control parameter is increased from $10^{-4}$ in v3 to $10^{-2}$ in v4. The reported flat-element condition number consequently decreases from approximately $9.2\times10^6$ to $9.2\times10^4$, with negligible change in the cited patch and twisted-beam observables.

### 3.1 Covariant membrane strain

For an infinitesimal midsurface translation field, the covariant membrane components are


```math
\epsilon_{11}=\mathbf g_1\cdot\mathbf u_{0,r},\qquad
\epsilon_{22}=\mathbf g_2\cdot\mathbf u_{0,s},\qquad
\epsilon_{12}=\tfrac12(\mathbf g_1\cdot\mathbf u_{0,s}+\mathbf g_2\cdot\mathbf u_{0,r}).
```


The physical tangent strain is obtained by transforming the covariant tensor into the pointwise orthonormal tangent basis. Denote this linear transformation by $\mathcal T_m(\mathbf C)$. The direct membrane operator is $\mathbf B_m^{\rm dir}=\mathcal T_m\mathbf B_m^{\rm cov}$.

### 3.2 Curvature-activated membrane tying

In v4 the membrane projection is **selected by a geometric indicator** rather than applied unconditionally:


```math
\boxed{\mathbf B_m=\begin{cases}
\mathcal T_m\,\widetilde{\mathbf B}_m^{\rm cov},&\chi_e=1,\\
\mathbf B_m^{\rm dir},&\chi_e=0.
\end{cases}}
```


Here $\widetilde{\mathbf B}_m^{\rm cov}$ uses separate 6-point projections for $\epsilon_{11}$ and $\epsilon_{22}$, and a 2×2 projection for $\epsilon_{12}$. Its interpolation employs the two-point linear and three-point quadratic Lagrange polynomials defined below. The indicator is


```math
\chi_e=\mathbf 1\!\left[\max_I\cos^{-1}(\mathbf n_I\cdot\overline{\mathbf n})>\theta_0
\quad\text{or}\quad
\max_{\text{edges }(i,j,m)}\frac{\|\mathbf X_m-(\mathbf X_i+\mathbf X_j)/2\|}{\|\mathbf X_j-\mathbf X_i\|}>\delta_0\right],
```


where $\theta_0=0.01$ rad and $\delta_0=10^{-3}$. Consequently, a flat straight-sided patch employs direct membrane strain and avoids the patch-test error introduced by always-on tying. An optional in-plane strain smoothing/interpolation step is inactive by default because earlier sampling variants degraded the membrane patch test.

<img width="1448" height="1086" alt="image" src="https://github.com/user-attachments/assets/4bf22d0d-e0a6-4dc6-b4ec-a51529a03efb" />

## 4. Transverse-shear interpolation and the v4 modification

### 4.1 Covariant transverse shear

Let $\mathbf t=\sum_I\widehat N_I\mathbf n_I$ and let the rotation field be $\boldsymbol\theta=\sum_I\widehat N_I\boldsymbol\theta_I$. In a pointwise director construction, a representative covariant shear expression is


```math
\gamma_\alpha=\mathbf t\cdot\mathbf u_{0,\alpha}
+\mathbf g_\alpha\cdot(\boldsymbol\theta\times\mathbf t),\qquad \alpha=1,2.
```


The implemented discrete natural-shear rows use the nodal normals in the rotation contributions and the interpolated normal in the translational derivatives. Tangent-physical shear components are obtained using the pointwise dual-basis mapping, not by tying fixed Cartesian components.

### 4.2 Tying polynomials and sample coordinates

Put $a=1/\sqrt3$, $b=\sqrt{3/5}$. Define


```math
L_-(x)=\frac{a-x}{2a},\qquad L_+(x)=\frac{a+x}{2a},
```



```math
Q_-(x)=\frac{x(x-b)}{2b^2},\quad Q_0(x)=1-\frac{x^2}{b^2},\quad
Q_+(x)=\frac{x(x+b)}{2b^2},
```


and $H_-(x)=(1-x)/2$, $H_+(x)=(1+x)/2$. The six-point-per-component interpolation is


```math
\widetilde\gamma_1^{(6)}(r,s)=\sum_{i=\pm}\sum_{j\in\{-,0,+\}}L_i(r)Q_j(s)\gamma_1(r_i,s_j),
\quad (r_i,s_j)=(\pm a,\{-b,0,b\}),
```



```math
\widetilde\gamma_2^{(6)}(r,s)=\sum_{j\in\{-,0,+\}}\sum_{i=\pm}Q_j(r)L_i(s)\gamma_2(r_j,s_i),
\quad (r_j,s_i)=(\{-b,0,b\},\pm a).
```


The four-point-per-component alternative is


```math
\widetilde\gamma_1^{(4)}(r,s)=\sum_{i=\pm}\sum_{j=\pm}L_i(r)H_j(s)\gamma_1(\pm a,\pm1),
\quad
\widetilde\gamma_2^{(4)}(r,s)=\sum_{j=\pm}\sum_{i=\pm}H_j(r)L_i(s)\gamma_2(\pm1,\pm a).
```


These are **twelve** and **eight** scalar shear tying observations per element, respectively. Such sample counts are not a proof of locking; the locking inference below is supported by the documented numerical slenderness studies.

<img width="1448" height="1086" alt="image" src="https://github.com/user-attachments/assets/36f8fec0-2142-4823-b5b0-9b6947672c58" />

### 4.3 Slenderness-dependent mixed projection

The central v4 change is the continuous interpolation


```math
\boxed{\widetilde{\boldsymbol\gamma}^{\rm v4}=(1-w_e)
\widetilde{\boldsymbol\gamma}^{(6)}+w_e\widetilde{\boldsymbol\gamma}^{(4)},\qquad
w_e=\max\!\left(0,1-\frac{s_0}{s_e}\right),\quad
s_e=\frac{\sqrt{A_e}}{h},\quad s_0=15.}
```


The physical $\mathbf B_s$ is recovered through the pointwise transformation from covariant to orthonormal tangent components. For $s_e\le15$, only the six-point interpolation is used; in the thin limit, $w_e\to1$ and the four-point interpolation dominates. The interpolation is performed at the **strain-operator level**; it is not a mixture of separately assembled element stiffness matrices. Therefore, its shear energy generally includes cross terms in $w_e(1-w_e)$.

The threshold $s_0=15$ is a calibrated research choice from the available plate studies, **not** an analytically established universal stability bound. The mixing rule is a practical anti-locking modification, and must be reassessed for severely warped, distorted, anisotropic or highly curved thin shells.

## 5. Bending, drilling and invariance

The bending operator employs pointwise gradients of the field functions, the interpolated nodal normal derivatives, and rotation contributions constructed from the local tangent basis. Its consistent curvature measure includes terms induced by director variation, rather than assuming an element-wide constant normal. The drilling constraint couples the rotation about the local normal to an antisymmetric in-plane displacement gradient:


```math
\epsilon_d\simeq\frac12\left(\mathbf e_2\cdot\mathbf u_{0,1}
-\mathbf e_1\cdot\mathbf u_{0,2}\right)-\boldsymbol\theta\cdot\mathbf n.
```


All physical strains, stress resultants and bending moments are evaluated in the pointwise tangent basis. On reversing element winding, geometric and load normals must follow the new parameter orientation while identical physical boundary-value problems should recover equivalent global responses. The material and rotational conventions are to be checked at stress-resultant level, not only with displacement patches.

## 6. Consistent loading, thermal effects and mass

### 6.1 Normal pressure and arbitrary traction

The consistent force for a signed normal pressure is


```math
\boxed{\mathbf f_p=\int_{-1}^1\!\int_{-1}^1
\mathbf N_u^T(r,s)\,p(r,s)\,[\mathbf g_1\times\mathbf g_2]\,dr\,ds.}
```


When a fixed global unit direction $\mathbf d_p$ is specified, replace the vector cross product by $J_A\mathbf d_p$. Thus the prescribed directional traction is measured per actual midsurface area. Positive normal pressure follows the element node-ordering cross product, not a forced positive global-Z convention. Five-point Gauss integration per parametric direction is the recorded default for distributed pressure; field interpolation may be consistent with the active element displacement field or deliberately selected as standard geometry interpolation for controlled comparisons.

### 6.2 Membrane and thickness-gradient temperature

Let $\Delta T(z)=\Delta T_0+z\Delta T_1$, with $\Delta T_1$ in temperature per length. The equivalent thermal load is


```math
\boxed{\mathbf f_T=\int_{A_e}\left[
 h\mathbf B_m^T\mathbf C_m\alpha\Delta T_0\mathbf a_T+
\frac{h^3}{12}\mathbf B_b^T\mathbf C_m\alpha\Delta T_1\mathbf a_T
\right]dA,\quad \mathbf a_T=[1,1,0]^T.}
```


The temperature gradient and bending-resultant signs follow the adopted kinematic curvature convention; a physical top-versus-bottom temperature statement must be converted into that convention before comparing external solver moments.

### 6.3 Consistent and positive lumped mass

With $\rho$ the mass density, the translational consistent mass is


```math
\mathbf M_{uu}=\int_{A_e}\rho h\,\mathbf N_u^T\mathbf N_u\,dA.
```


The rotational mass uses the inertial density $\rho h^3/12$ and either the full spatial rotation projector $\mathbf I_3$ or the tangent projector $\mathbf I_3-\mathbf n\otimes\mathbf n$:


```math
\mathbf M_{\theta\theta}=\int_{A_e}\frac{\rho h^3}{12}\mathbf N_\theta^T\mathbf P_\theta\mathbf N_\theta\,dA.
```


Positive diagonal lumping rescales positive nodal quadratic weights so that the resultant nodal masses recover the integrated total. This avoids the negative corner masses possible with direct row-sum lumping of a serendipity interpolation. The three-component rotation choice represents a modelling convention rather than an assertion of physical transverse-normal rotary inertia in classical thin-shell theory.

## 7. Initial-stress stiffness and bifurcation

For in-plane prestress $\boldsymbol\sigma_m=\begin{bmatrix}\sigma_{11}&\sigma_{12}\\\sigma_{12}&\sigma_{22}\end{bmatrix}$ in the pointwise tangent frame and tangent derivatives $\nabla_s\widehat N_I$, the reported geometric stiffness is


```math
\boxed{[\mathbf K_G]_{Iu,Jv}=\int_{A_e}
 h\,(\nabla_s\widehat N_I)^T\boldsymbol\sigma_m
 (\nabla_s\widehat N_J)\,P_{uv}\,dA,}
```


where $\mathbf P=\mathbf I_3$ for all translational components, or $\mathbf P=\mathbf n\otimes\mathbf n$ for a bending-only normal-direction approximation. Rotation-related initial-stress terms and a complete nonlinear curved-shell tangent are **not** included. The linearised buckling equation is


```math
[\mathbf K+\lambda\mathbf K_G]\boldsymbol\phi=\mathbf0,
```


with the usual convention of compression-negative prestress in $\mathbf K_G$. The full three-translation option is essential for the in-plane/column buckling comparisons; restricting prestress stiffness to the surface normal would omit those mechanisms.

## 8. Formulation evolution and rationale

| Feature | Earlier v1p2 | v3 | v4, adopted research formulation |
|---|---|---|---|
| Geometry | Planar element frame | Pointwise curved 3-D midsurface | Preserved pointwise geometry |
| Nodal kinematics | Fixed-frame projections | Globally expressed 48 DOF; varying directors | Preserved |
| Interpolation | Standard eight-node field and legacy membrane smoothing | Kikuchi default; optional alternative field | Kikuchi default; ILS omitted |
| Membrane strain | Legacy smoothing that failed the documented patch | Direct covariant-to-physical strain | Direct on flat patches; 6+6+4 covariant tying on curved/warped elements |
| Shear | Four-point fixed-frame ANS-like interpolation | Six-point-per-component covariant tying | Slenderness-weighted six-/four-point covariant tying |
| Drilling parameter | $10^{-4}$ | $10^{-4}$ | $10^{-2}$ |
| Mass | Limited earlier implementation | Not provided | Consistent and positive diagonal |
| Pressure and thermal loads | Incomplete | Not provided | Consistent distributed loading and two thermal resultants |
| Initial-stress stiffness | Not provided | Not provided | Translational geometric stiffness |

The v3 improvement over v1p2 is chiefly geometric and kinematic; the key v4 changes are the selective membrane projection and mixed shear constraint, plus additional analysis operators. Consequently, the difference between v3 and v4 should not be attributed solely to the shear blend in every benchmark.

## 9. Verification evidence

### 9.1 Evidence conventions

The numerical tables distinguish **independent analytical/CSI-reference comparisons**, **same-environment formulation studies**, and **subsequent native MYSTRAN verification**. Values shown as signed percent error follow $100(q_h/q_{\rm ref}-1)$. A nearly zero error in a particular benchmark does not imply validation of every strain and resultant component.

### 9.2 Membrane/bending patches and fundamental shell cases

| Benchmark observable | v1p2 | v3 | v4 | Interpretation |
|---|---:|---:|---:|---|
| 2-001 membrane patch, displacement error (%) | 13.5 | approximately $2\times10^{-12}$ | approximately $4\times10^{-12}$ | Flat membrane consistency retained |
| 2-001 bending patch, displacement error (%) | 56.3 | approximately $3\times10^{-10}$ | approximately $4\times10^{-10}$ | Flat bending consistency retained |
| 2-003 curved-beam LC1, $n_x=2$, signed error (%) | −98.9 | −91.4 | −20.0 | Coarse-mesh curvature response improved |
| 2-003 curved-beam LC1, $n_x=8$, signed error (%) | −70.8 | −4.53 | −1.11 | Improvement from selective membrane tying |
| 2-003 curved-beam LC1, $n_x=32$, signed error (%) | −66.7 | −0.10 | −0.07 | Converged response |
| T9, $N=4/8/16$, response ratio | 0.37 / 0.87 / 1.07 | 0.82 / 0.88 / 1.00 | approximately 1.24 / 0.992 / 1.01 | Coarse-mesh overshoot remains |
| 2-004 twisted member, $n_x=24$, $U_z/U_y$ error (%) | −53 / −68 | +1.76 / +7.01 | +1.75 / +7.02 | No material regression from v3 |

<img width="1448" height="1086" alt="image" src="https://github.com/user-attachments/assets/ad2ad4d0-c747-47af-84c0-c0901fd3b640" />

<img width="1448" height="1086" alt="image" src="https://github.com/user-attachments/assets/7237ce6e-ba34-4293-9bad-d765813dda2e" />

### 9.3 Thin-plate locking and slenderness

A quarter simply supported square plate at $h=10^{-4}$ has the following reported central-displacement ratio to the independent thin-plate reference:

| Quarter-mesh | Earlier six-point-only v4 | Final v4, blended shear |
|---|---:|---:|
| 1 × 1 | 0.380 | 0.912 |
| 2 × 2 | 0.095 | 0.991 |
| 4 × 4 | 0.053 | 1.001 |
| 8 × 8 | 0.684 | 1.001 |

The non-monotonic earlier sequence confirms that mesh density alone cannot be assumed to cure the particular shear constraint. For a 4 × 4 quarter-plate sweep, the final mixed operator yields a ratio of approximately 1.001 even at thickness $10^{-4}$, whereas the six-point-only result falls to 0.053. On an elongated 2 × 10 plate the single-cell-per-quarter mixed result overshoots by approximately 35%, a clear counterexample to any blanket claim of superior coarse-mesh accuracy.

<img width="1672" height="941" alt="image" src="https://github.com/user-attachments/assets/3d6d1ffd-7145-48f1-b093-087305475bd7" />
**Figure 6 :** continuous convergence curves for central deflection ratio versus quarter-mesh refinement, two shear projections, both square and elongated panels.  
<img width="1448" height="1086" alt="image" src="https://github.com/user-attachments/assets/0b778f48-4b75-4215-804d-ea3f48211913" />
**Figure 7 :** deflection ratio versus thickness, fixed 4 × 4 plate mesh.

### 9.4 Adapted CSI problem 2-005: transverse pressure

The recorded study uses a quarter model of a 2 × 10 thin plate and a progressively refined Q8 mesh. Listed errors are relative to the recorded independent reference for the same support and loading case.

| Quarter-mesh density $N$ | 3 | 6 | 9 | 12 |
|---|---:|---:|---:|---:|
| Simple support, uniform pressure: error (%) | +0.90 | +0.10 | −0.03 | −0.07 |
| Clamped, uniform pressure: error (%) | +2.97 | +0.93 | +0.44 | +0.25 |
| Simple support, concentrated load: error (%) | −3.00 | −1.78 | −1.19 | −1.06 |
| Clamped, concentrated load: error (%) | −5.48 | −3.54 | −2.07 | −1.48 |

A later **native MYSTRAN** validation note reports 2-005 CQUAD8 at 16 × 16 with central $U_z=-12.97469330$ in, compared with $-12.97$ in, or +0.0362% signed error; the matched earlier native model gave $-13.51560497$ in (+4.2067%). This is a separate executable/deck campaign, not the same numerical series as the four-quarter-mesh study above. A Nastran solver value for that exact matched study is not supplied here.

<img width="1448" height="1086" alt="image" src="https://github.com/user-attachments/assets/ee2846ae-7583-4797-922d-2add745125d5" />

### 9.5 Dynamics and mass consistency

For a simply supported plate with $h/a=0.01$, the recorded first-frequency ratios $f_h/f_{\rm analytical}$ are:

| Mesh | Consistent mass | Positive diagonal mass |
|---|---:|---:|
| 4 × 4 | 1.0319 | 1.0295 |
| 8 × 8 | 1.00009 | 0.99996 |

Total translational mass ratios from planar, hyperbolic-paraboloid and cylindrical element checks were approximately 1.0000000, 0.99999997 and 0.99995. The residual discrepancy on the curved element was attributed in the source report to finite-order integration of a non-polynomial area Jacobian. A separate subsequent native implementation note reports six modal density/added-mass/diagonal-mass checks with maximum relative difference for the first five frequencies of $3.413\times10^{-7}$ against the independent element reference computation; individual cases and spectra are not reproduced in the available summary.

### 9.6 Thermal consistency

For a free plate under uniform temperature change, the displacement residual against the analytical expansion field is reported as $8.8\times10^{-17}$ (reference displacement scale $10^{-3}$); residual mechanical strain is approximately $1.5\times10^{-16}$. For a prescribed thickness gradient, the reported transverse-displacement and rotation residuals are $6.5\times10^{-18}$ and $1.2\times10^{-17}$. These are patch-level tests. The later MYSTRAN note additionally lists five native thermal checks (free and restrained expansion, two winding-related gradient patches and an annulus case) as completed, but does not provide a per-case numerical table. Curved thermal bending and thermal-prestress buckling remain outside this reported evidence.

### 9.7 Linear buckling: adapted CSI 2-016 and 2-017

For problem 2-016 at $h=0.01$, the reference **critical stresses** are $0.042501034$ for the shear condition and $0.116491057$ for bending. They are critical stresses, not raw eigenvalue multipliers; the recorded driver uses unit reference line load and $\sigma_{\rm cr}=\lambda/h$.

| Mesh | v4 shear critical stress | v4 bending critical stress | Shear error (%) | Bending error (%) |
|---|---:|---:|---:|---:|
| 8 × 8 | 0.04297455 | 0.11734247 | +1.11 | +0.73 |
| 16 × 16 | 0.04242379 | 0.11611342 | −0.182 | −0.324 |

At 16 × 16 the v3 shear and bending errors are +8.073% and +7.281%, respectively. The separate native MYSTRAN before/after campaign reports 2-016 16 × 16 final errors of −0.1817% and −0.3242%, consistent with the independent formulation reference. On 2-017 the first native column eigenfactors are 126.5028 (2 × 12) and 126.4854 (4 × 24), unchanged from its matched earlier run. This checks an in-plane geometric-stiffness mode, but does not validate rotational prestress coupling or curved-shell nonlinear stability.

<img width="1672" height="941" alt="image" src="https://github.com/user-attachments/assets/d6320221-a7d4-4d2b-ab6a-8ad318dcc9b8" />
**Figure 9 :** 2-016 critical shear and bending stresses versus mesh density, with independent reference lines and v3/v4 connected curves.  
<img width="1448" height="1086" alt="image" src="https://github.com/user-attachments/assets/7fa4537e-824a-41b6-9a5d-74393349b8f8" />
**Figure 10 :** 2-017 lowest column buckling factor versus refinement.

## 10. Discussion, limitations and reproducibility

The v4 modification reduces thin-plate shear locking while retaining the principal benefits of the curved-surface formulation introduced in v3. The change is neither purely geometric nor a universal four-point projection: it couples curvature-dependent membrane selection with continuous shear interpolation. The fact that v4 can overshoot on a single elongated quarter-panel and on the coarsest T9 mesh is important; broad claims of mesh-independent accuracy would be unsupported. A computed convergence trend should be reported for each load/support class, rather than substituting a single pass/fail metric.

The available evidence is strongest for flat-plate bending, consistent pressure forces, patch-type thermal response, selected modal checks, and thin-plate/column bifurcation. In contrast, the following remain **unverified or incomplete**: severely distorted high-curvature shell studies; rigorous thin-limit stability proof; laminate/orthotropic constitutive matrices; variable shell thickness; nonlinear finite rotations; full rotational/curvature initial-stress terms; thermal-prestress buckling; independently audited curved-shell fiber stresses and moments; and complete matched external-solver comparisons on every CSI case. The $s_0=15$ blending parameter is empirically selected and should be treated as a calibration parameter pending a broader convergence study.

The exact code revision, active field functions, nodal winding, density convention, pressure normal, reference prestress, integration order, mode correspondence and recovery basis should accompany each future comparison. **Reported native MYSTRAN evidence and standalone element studies must not be silently combined into a single matched Nastran comparison.**

## 11. Conclusions

1. A pointwise three-dimensional eight-node shell with curvature-selected covariant membrane projection retains flat membrane and bending patch consistency while improving the cited coarse curved-beam benchmarks.
2. The central v4 shear innovation is an operator-level, slenderness-dependent interpolation between six- and four-point-per-component covariant tying. Within the reported plate studies it materially reduces thin-plate locking, but a general mathematical locking-free guarantee has not been demonstrated.
3. Consistent pressure, temperature, mass and translational prestress operators enable static, modal and restricted eigen-buckling analyses within a common tangent-frame convention.
4. For the adapted 2-016 thin-plate benchmark, errors at 16 × 16 are −0.182% (shear) and −0.324% (bending), compared with +8.073% and +7.281% for the recorded v3 baseline.
5. Further independent tests are needed for curved-shell buckling, full stress-resultant frame consistency and highly distorted/slender shells before general engineering use can be claimed.

## References

1. Bathe, K.-J., and Dvorkin, E. N. (1985). A four-node plate bending element based on Mindlin/Reissner plate theory and a mixed interpolation. *International Journal for Numerical Methods in Engineering*, **21**, 367–383.
2. Dvorkin, E. N., and Bathe, K.-J. (1984). A continuum mechanics based four-node shell element for general nonlinear analysis. *Engineering Computations*, **1**, 77–88.
3. Bucalem, M. L., and Bathe, K.-J. (1993). Higher-order MITC general shell elements. *International Journal for Numerical Methods in Engineering*, **36**, 3729–3754. [Bibliographic details should be checked against the publisher's final record before journal submission.]
4. Hughes, T. J. R., and Brezzi, F. (1989). On drilling degrees of freedom. *Computer Methods in Applied Mechanics and Engineering*, **72**, 105–121.
5. MacNeal, R. H., and Harder, R. L. (1985). A proposed standard set of problems to test finite element accuracy. *Finite Elements in Analysis and Design*, **1**, 3–20.
6. Computers and Structures, Inc. (CSI). *SAP2000/CSiBridge Verification Manuals*, shell problems 2-001 through 2-017. Representative problem references: [2-001](https://docs.csiamerica.com/manuals/csibridge/Verification/Analysis/Shells/Problem%202-001.pdf), [2-003](https://docs.csiamerica.com/manuals/csibridge/Verification/Analysis/Shells/Problem%202-003.pdf), [2-004](https://docs.csiamerica.com/manuals/csibridge/Verification/Analysis/Shells/Problem%202-004.pdf), [2-005](https://docs.csiamerica.com/manuals/csibridge/Verification/Analysis/Shells/Problem%202-005.pdf), [2-016](https://docs.csiamerica.com/manuals/csibridge/Verification/Analysis/Shells/Problem%202-016.pdf), [2-017](https://docs.csiamerica.com/manuals/csibridge/Verification/Analysis/Shells/Problem%202-017.pdf).
7. Internal research records (2026): *MITC8 v1p2 versus v3 review*; *MITC8 v3 to v4 changes*; *MITC8v4–MITC8D formulations*; *SIMOQ8 V13 alignment and MITC8 V4 upgrade*. Internal validation reports; not independent peer-reviewed literature.

---

### Supplement: graphical material attached
| Figure | Intended content | Axes / annotation |
|---|---|---|
| 1 | CQUAD8 natural coordinate topology | Corner and midside numbering, positive normal |
| 2 | Membrane sampling pattern | 6+6+4 covariant observations |
| 3 | Transverse shear sample positions | Six- versus four-point-per-component patterns |
| 4 | 2-003 curved-beam convergence | Element count vs signed displacement error (%) |
| 5 | T9 curved-shell convergence | N vs response/reference |
| 6 | 2-005 load/support convergence | Quarter-mesh density vs signed displacement error (%) |
| 7 | Thickness-dependent locking | Thickness vs central-deflection ratio |
| 8 | 2-016 thin-plate critical stresses | Mesh refinement vs critical stress; two modes |
| 9 | 2-017 column stability | Panel refinement vs lowest eigenvalue |

**Data integrity statement:** The reported benchmark results were transcribed from the supplied formulation review and the available project reports. No fresh benchmark simulation was performed for this document. Unreported or non-matched Nastran values have not been manufactured.
