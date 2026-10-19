# ANS8–BDG6 v4: A Geometry-Adaptive Eight-Node Shell with Hybrid Assumed Strains and Thickness-Dependent Transverse-Shear Blending
**Ary · Independent Researcher · Personal Research** · MYSTRAN v18.00.a · 19 October 2026
**Technical formulation and numerical-verification manuscript | MYSTRAN NEWSHELL / `QUAD8TYP,ANS8BDG6`**  
**Revision:** 20 October 2026  
**Research status:** A project-specific hybrid inspired by Jung–Han (2013), Han et al. (2004), and quadratic-field corrections associated with Kikuchi/MacNeal–Harder. **Not** a verbatim implementation of the published Jung–Han 2013 element. Numerical values labeled *archived* are transcribed from the supplied ANS8BDG6 v4 validation record, not newly executed analyses. An independent source-to-source audit of all intermediate operators is still required for definitive publication.

## Abstract

This note defines the eight-node, six-degree-of-freedom-per-node ANS8–BDG6 v4 research shell and establishes the relationship between its assumed-natural-strain (ANS) membrane interpolation, selectively activated membrane tying, serendipity-field variant, and an aspect-ratio-dependent mixture of two transverse-shear interpolation patterns. The element combines the β and δ strain-sampling families associated with eight-node ANS shell literature with two different transverse-shear families, denoted here by *tensor6* and *BDG4*. The default variant applies a modified Kikuchi-type field, preserves a point-dependent output basis, and switches membrane tying on curved or warped geometries. Its transverse-shear blending function depends on the dimensionless characteristic length-to-thickness ratio, providing a continuous interpolation between the more constrained tensor6 and more flexible BDG4 shear operators. The presented benchmark record includes the CSI 2-005 thin-plate deflection, 2-016 out-of-plane plate buckling, 2-017 column buckling, modal parity, thermal loading, and direct comparison of production operators against the independent reference implementation. Crucially, these improvements carry nontrivial **patch-consistency, distortion, locking, recovery, and curved-shell limitations**; agreement between implementations is distinguished throughout from agreement with independent analytical or CSI references.

**Keywords:** quadratic shell; assumed natural strain; eight-node serendipity; membrane locking; transverse shear locking; Kikuchi-type interpolation; thickness-adaptive strain blending; geometric stiffness; MYSTRAN.

## 1. Kinematics and degrees of freedom

Let the parent surface occupy \$\\widehat\\Omega=[-1,1]^2\$, with natural coordinates \$\\boldsymbol\\xi=(r,s)\$. The eight interpolation nodes follow the conventional ordering: four corners \$(-1,-1),(1,-1),(1,1),(-1,1)\$ and midsides \$(0,-1),(1,0),(0,1),(-1,0)\$. The conventional Q8 corner shape functions are

\
```math
N_i(r,s)=\\frac14(1+r r_i)(1+s s_i)(r r_i+s s_i-1),\\qquad i=1,\\ldots,4,
\
```


with the remaining four

\
```math
\\begin{aligned}
N_5&=\\tfrac12(1-r^2)(1-s), & N_6&=\\tfrac12(1+r)(1-s^2),\\\\
N_7&=\\tfrac12(1-r^2)(1+s), & N_8&=\\tfrac12(1-r)(1-s^2).
\\end{aligned}
\
```


For each node, the global structural degrees of freedom are \$(u_x,u_y,u_z,\\theta_x,\\theta_y,\\theta_z)\$, giving \$\\mathbf d_e\\in\\mathbb R^{48}\$. Let \$\\mathbf Q_e=[\\mathbf e_1,\\mathbf e_2,\\mathbf e_3]\$ denote the element reference orthonormal frame and \$\\mathbf T_e=\\operatorname{diag}(\\mathbf Q_e^T,\\mathbf Q_e^T,\\ldots)\$ the six-DOF rotation operator (with separate three-dimensional blocks for translations and rotations). Then \$\\mathbf d_e^{\\ell}=\\mathbf T_e\\mathbf d_e^{g}\$. The geometric surface is mapped using \$\\mathbf x(r,s)=\\sum_I N_I(r,s)\\mathbf x_I\$, while the default modified interpolation may act on the kinematic field; **the actual modified shape functions must be distinguished from the standard geometry mapping**.

For local flat coordinates \$(x,y)\$, \$\\mathbf J=\\partial(x,y)/\\partial(r,s)\$. Its inverse transforms derivatives with a consistent row/column convention. The local constitutive matrices are

\
```math
\\mathbf C_m=\\frac{E}{1-\\nu^2}
\\begin{bmatrix}1&\\nu&0\\\\\\nu&1&0\\\\0&0&(1-\\nu)/2\\end{bmatrix},
\\quad \\mathbf C_s=\\kappa_s\\frac{E}{2(1+\\nu)}\\mathbf I_2,\\quad \\kappa_s=5/6.
\
```


## 2. Compatible and assumed strain operators

In first-order shear-deformation kinematics, let \$(u,v,w,\\theta_x,\\theta_y,\\theta_z)\$ be the local degrees of freedom, with the rotational sign convention

\
```math
\\boldsymbol\\gamma_s=
\\begin{bmatrix}w_{,x}+\\theta_y\\\\w_{,y}-\\theta_x\\end{bmatrix},
\\qquad
\\boldsymbol\\kappa=
\\begin{bmatrix}\\theta_{y,x}\\\\-\\theta_{x,y}\\\\\\theta_{y,y}-\\theta_{x,x}\\end{bmatrix}.
\
```


The compatible membrane strain is \$\\boldsymbol\\varepsilon_m=[u_{,x},v_{,y},u_{,y}+v_{,x}]^T\$. Corresponding discrete operators satisfy \$\\boldsymbol\\varepsilon_m=\\mathbf B_m\\mathbf d_e^\\ell\$, \$\\boldsymbol\\kappa=\\mathbf B_b\\mathbf d_e^\\ell\$, and \$\\boldsymbol\\gamma_s=\\mathbf B_s\\mathbf d_e^\\ell\$. With parametric tangent vectors \$\\mathbf g_1=\\mathbf x_{,r}\$, \$\\mathbf g_2=\\mathbf x_{,s}\$, define natural membrane tensor components

\
```math
\\widetilde\\varepsilon_{\\alpha\\beta}=\\tfrac12(\\mathbf g_\\alpha\\cdot\\mathbf u_{,\\beta}+\\mathbf g_\\beta\\cdot\\mathbf u_{,\\alpha}),
\\qquad \\alpha,\\beta\\in\\{1,2\\}.
\
```


Physical engineering-strain components are recovered by contravariant transformation; engineering shear factors of two must be included explicitly. If \$\\mathbf A=\\mathbf J^{-1}\$, \$\\varepsilon_{ij}=A_{\\alpha i}A_{\\beta j}\\widetilde\\varepsilon_{\\alpha\\beta}\$ under the declared mapping convention. **A physical tensor and its element-local coordinate components are not interchangeable.**

### 2.1 β membrane-normal interpolation and δ in-plane shear

Introduce \$a=1/\\sqrt3\$, \$b=\\sqrt{3/5}\$, and \$c=1\$. A β pattern has six sampling points for each natural membrane-normal component,

\
```math
\\mathcal P_{\\beta,11}=\\{(\\pm c,-b),(\\pm c,0),(\\pm c,b)\\},
\\quad
\\mathcal P_{\\beta,22}=\\{(-b,\\pm c),(0,\\pm c),(b,\\pm c)\\}.
\
```


The δ in-plane shear pattern samples at \$\\mathcal P_\\delta=\\{(\\pm a,\\pm a)\\}\$. With one-dimensional linear and quadratic Lagrange polynomials \$\\ell_i^{(1)}\$ and \$\\ell_j^{(2)}\$, the membrane operators are reconstructed in the general form

\
```math
\\overline{\\mathbf B}_{m,11}(r,s)=\\sum_{i=1}^2\\sum_{j=1}^3
\\ell_i^{(1)}(r;[-1,1])\\ell_j^{(2)}(s;[-b,0,b])
\\mathbf B_{m,11}(r_i,s_j),
\
```


\
```math
\\overline{\\mathbf B}_{m,22}(r,s)=\\sum_{i=1}^2\\sum_{j=1}^3
\\ell_i^{(1)}(s;[-1,1])\\ell_j^{(2)}(r;[-b,0,b])
\\mathbf B_{m,22}(r_j,s_i),
\
```


\
```math
\\overline{\\mathbf B}_{m,12}(r,s)=\\sum_{i=1}^{2}\\sum_{j=1}^{2}
\\ell_i^{(1)}(r;[-a,a])\\ell_j^{(1)}(s;[-a,a])
\\mathbf B_{m,12}(r_i,s_j).
\
```


These equations describe the **β–δ family**; the actual research element may construct covariant samples and map them to physical components before interpolation. These operations do not generally commute on distorted elements, so the order of the transformation is a validation-critical choice.

### 2.2 Selective membrane tying

The element uses a geometry-sensitive membrane selection. Let \$\\mathbf n_I\$ be the nodal directors, \$\\overline{\\mathbf n}\$ their normalized average, \$\\mathbf x_{m(e)}\$ the midpoint node of edge \$e=(i,j)\$, and \$L_e=\\|\\mathbf x_j-\\mathbf x_i\\|\$. Define

\
```math
\\theta_{\\max}=\\max_I\\arccos(\\mathbf n_I\\cdot\\overline{\\mathbf n}),
\\qquad
\\delta_{\\max}=\\max_e\\frac{\\|\\mathbf x_{m(e)}-\\tfrac12(\\mathbf x_i+\\mathbf x_j)\\|}{L_e}.
\
```


The adopted switching rule is

\
```math
\\boxed{\\mathbf B_m^{\\mathrm{active}}=
\\begin{cases}
\\overline{\\mathbf B}_m^{\\beta\\delta}, & \\theta_{\\max}>0.01\\;\\mathrm{rad}\\ \\text{or}\\ \\delta_{\\max}>10^{-3},\\\\
\\mathbf B_m^{\\mathrm{direct}}, & \\text{otherwise}.
\\end{cases}}
\
```


Thus flat, straight-sided patches retain compatible membrane strain, while curved/warped elements receive assumed strain. The switching rule is **heuristic and discontinuous** at the threshold; it is not a theorem of locking-free behavior. Variants forcing membrane tying on or off should remain exposed for verification.

## 3. Two alternative transverse-shear sampling spaces

The v4 research formulation retains three named transverse-shear choices: *tensor6*, *BDG4*, and *direct*. To prevent confusion with Jung and Han's published βδγ nomenclature, these are here explicitly treated as **project operator names**.

### 3.1 Four-sample BDG4 shear operator

For the classical γ-style four-sample pattern, the two components are tied independently:

\
```math
\\mathcal P_{\\gamma,13}=\\{(\\pm a,\\pm1)\\},
\\qquad
\\mathcal P_{\\gamma,23}=\\{(\\pm1,\\pm a)\\}.
\
```


For each component, let \$\\mathcal I_{4}\$ denote bilinear interpolation through its four points. Then

\
```math
\\mathbf B_s^{\\mathrm{BDG4}}=
\\begin{bmatrix}
\\mathcal I_4^{\\gamma,13}(\\mathbf B_{s,13}^{\\mathrm{direct}})\\\\
\\mathcal I_4^{\\gamma,23}(\\mathbf B_{s,23}^{\\mathrm{direct}})
\\end{bmatrix}.
\
```


### 3.2 Six-sample tensor6 shear operator

The alternate tensor6 operator uses six tying values **per transverse-shear component**, with mixed linear–quadratic tensor-product reconstruction:

\
```math
\\mathcal I_{6}[q](r,s)=\\sum_{i=1}^{2}\\sum_{j=1}^{3}
\\ell_i^{(1)}(r)\\ell_j^{(2)}(s)q(r_i,s_j)
\
```


for the first component, with the parametric directions interchanged for the second component. The exact collocation locations and transformation convention belong to the **v4 source dependency** and are not independently asserted here because the uploaded v4 wrapper delegates their construction to a shared neutral formulation. The constitutive role is nevertheless explicit: it produces \$\\mathbf B_s^{\\mathrm{tensor6}}\\in\\mathbb R^{2\\times48}\$. A later source-complete revision should reproduce its specific coordinate table instead of inferring it from the historical β–δ arrangement.

### 3.3 Thickness-adaptive shear blend (v4)

The characteristic aspect ratio is \$\\chi=\\sqrt{A_e}/h\$, where \$A_e\$ is the element's **standard-geometry** surface area. With transition parameter \$\\chi_0=15\$, define

\
```math
\\boxed{\\omega(\\chi)=\\max\\!\\left(0,1-\\frac{\\chi_0}{\\chi}\\right)},
\\qquad
\\boxed{\\mathbf B_s^{\\mathrm{auto}}=(1-\\omega)\\mathbf B_s^{\\mathrm{tensor6}}+\\omega\\mathbf B_s^{\\mathrm{BDG4}}}.
\
```


For \$\\chi\\le 15\$, the blend is entirely tensor6; as \$h\\to0\$, \$\\omega\\to1\$ and the shear operator approaches BDG4. Unlike switching the integrated shear stiffness directly, **the strain operators are blended before forming** \$\\mathbf B_s^T\\mathbf C_s\\mathbf B_s\$:

\
```math
\\mathbf K_s(\\omega)=\\int_A h\\,\\big[(1-\\omega)\\mathbf B_6+\\omega\\mathbf B_4\\big]^T\\mathbf C_s
\\big[(1-\\omega)\\mathbf B_6+\\omega\\mathbf B_4\\big]\\,dA.
\
```


This contains cross-coupling \$\\omega(1-\\omega)(\\mathbf B_6^T\\mathbf C_s\\mathbf B_4+\\mathbf B_4^T\\mathbf C_s\\mathbf B_6)\$, which is absent from an arithmetic combination of preintegrated stiffness matrices. The distinction is essential to reproduce v4 faithfully.

**Trade-off:** the archived v4 design notes report substantial thin-plate locking for pure tensor6 (approximately 5% of reference deflection at \$h=10^{-4}\$) and excess flexibility for pure BDG4 (approximately +27% deflection at 2×2 and +7% at 4×4). The auto blend improves the supplied plate test but is an **empirical thickness-sensitive choice**, not an unconditional proof of optimal convergence, patch consistency, or stability for arbitrary distortion.

![Shear blend coefficient](figures/fig03_shear_blend.png)

*Figure 3. Adopted shear blending weight versus normalized characteristic dimension \$\\sqrt{A_e}/h\$. This is a plot of the defining equation, not fitted numerical test data.*

## 4. Default quadratic field and bending/drilling strain

The default field variant is a **Kikuchi-type corrected Q8** rather than pure standard serendipity. The alternative MacNeal–Harder choice is a research configuration, not an independently validated native MYSTRAN selector. It is useful to express a modified field abstractly as

\
```math
\\widetilde N_I=N_I+\\Delta N_I,
\\qquad
\\sum_I\\Delta N_I=0,
\
```


subject to the constant- and linear-field completeness constraints appropriate to the correction. **This is a representation, not a claim that the correction \$\\Delta N_I\$ has been reconstructed from the supplied v4 wrapper**; publishing exact modified field coefficients requires the underlying v2/common source. The active membrane, bending, shear, and drilling matrices must be evaluated with compatible choices of interpolation, geometry, and local frame.

In the adopted local rotation convention, bending and drilling operators can be written

\
```math
\\mathbf B_b\\mathbf d_e^\\ell=
[\\theta_{y,x},-\\theta_{x,y},\\theta_{y,y}-\\theta_{x,x}]^T,
\
```


\
```math
\\mathbf B_d\\mathbf d_e^\\ell=
\\theta_z+\\tfrac12(v_{,x}-u_{,y}),
\
```


with the sign checked against the element's actual right-hand rotational convention. The default drilling parameter is \$k_t=0.1\$; the archived native formulation identifies \$\\beta_{\\mathrm{drill}}=k_t h\$. Let \$c_d\$ denote the implementation's fully dimensioned drilling energy coefficient, which must be audited so that no thickness factor is applied twice:

\
```math
\\mathbf K_d=\\int_A c_d\\,\\mathbf B_d^T\\mathbf B_d\\,dA.
\
```


## 5. Virtual work and discrete stiffness

For isotropic homogeneous thickness \$h\$, the discrete small-strain element energy is

\
```math
\\Pi_e=\\frac12\\int_A\\left[
 h\\,\\boldsymbol\\varepsilon_m^T\\mathbf C_m\\boldsymbol\\varepsilon_m
+\\frac{h^3}{12}\\boldsymbol\\kappa^T\\mathbf C_m\\boldsymbol\\kappa
+h\\,\\boldsymbol\\gamma_s^T\\mathbf C_s\\boldsymbol\\gamma_s
+c_d\\,\\gamma_d^2\\right]dA-\\mathbf d_e^T\\mathbf f_e.
\
```


The active 48×48 stiffness therefore has four contributions:

\
```math
\\boxed{\\mathbf K_e=\\int_{-1}^{1}\\!\\int_{-1}^{1}
\\Big(h\\mathbf B_m^T\\mathbf C_m\\mathbf B_m+
\\frac{h^3}{12}\\mathbf B_b^T\\mathbf C_m\\mathbf B_b+
 h\\mathbf B_s^T\\mathbf C_s\\mathbf B_s+
 c_d\\mathbf B_d^T\\mathbf B_d\\Big)|J|\\,dr\\,ds.}
\
```


The formulation uses full **3×3 Gauss–Legendre surface integration** for elastic stiffness, with abscissae \$0,\\pm\\sqrt{3/5}\$ and weights \$8/9,5/9,5/9\$. \$\\sum_{i,j}w_iw_j=4\$. Figure 1 identifies the nine integration points; Figure 2 identifies the classical β–δ–γ tying families, which are **not** to be confused with the tensor6 shared-operator locations.

![Q8 Gauss points](figures/fig01_q8_gauss.png)

*Figure 1. Node positions and nine Gauss–Legendre integration points in the Q8 parent plane.*

![ANS tying families](figures/fig02_ans_tying.png)

*Figure 2. The β normal-membrane, δ in-plane-shear, and γ four-point transverse-shear sampling sets, reflecting the family of interpolation patterns motivating the hybrid formulation.*

## 6. Mass, pressure, temperature, and initial-stress stiffness

### 6.1 Consistent and lumped mass

Let \$\\mathbf N_6(r,s)\$ map nodal six-DOF vectors to interpolated fields. A generic consistent inertial contribution is

\
```math
\\mathbf M_e=\\int_A\\mathbf N_6^T\\,\\operatorname{diag}
(\\rho h,\\rho h,\\rho h,\\rho h^3/12,\\rho h^3/12,I_{\\theta_3})\\,\\mathbf N_6\\,dA.
\
```


The native v4 record reports both consistent and positive HRZ-style lumped mass; default rotary inertia policy is **all**, whereas PSHELL nonstructural mass remains translational. Its precise drilling-rotary inertia scaling depends on that policy and should not be inferred from the standard Mindlin expression alone. Material density and nonstructural mass must be separately documented.

### 6.2 Consistent normal pressure and thermal resultants

For signed pressure \$p\$, using the connectivity normal \$\\mathbf n=\\mathbf g_1\\times\\mathbf g_2/\\|\\mathbf g_1\\times\\mathbf g_2\\|\$, the nodal load is

\
```math
\\mathbf f_e^p=\\int_A\\mathbf N_u^T p\\mathbf n\\,dA.
\
```


For free uniform temperature shift \$\\Delta T_0\$ and an independent through-thickness gradient parameter \$\\Delta T_1\$, introduce \$\\boldsymbol\\varepsilon_T=\\alpha\\Delta T_0[1,1,0]^T\$ and \$\\boldsymbol\\kappa_T=\\alpha\\Delta T_1[1,1,0]^T\$, so

\
```math
\\mathbf f_e^T=\\int_A\\left(h\\mathbf B_m^T\\mathbf C_m\\boldsymbol\\varepsilon_T+
\\frac{h^3}{12}\\mathbf B_b^T\\mathbf C_m\\boldsymbol\\kappa_T\\right)dA,
\
```


provided \$\\Delta T_1\$ uses the same sign convention as \$\\mathbf B_b\$. The archived implementation specifically documents an **opposite thermal gradient argument convention** between the native TEMPP1 mapping and its comparison source: this is not an elastic stiffness discrepancy.

### 6.3 Geometric stiffness for linear buckling

With physical membrane prestress \$\\boldsymbol\\sigma_0\$ and in-plane gradients \$\\nabla_sN_I\$, the translational initial-stress stiffness blocks have the familiar form

\
```math
(K_G)_{Ia,Jb}=\\delta_{ab}\\int_A h\\,(\\nabla_sN_I)^T\\boldsymbol\\sigma_0(\\nabla_sN_J)\\,dA,
\\qquad a,b=1,2,3.
\
```


The stress must be expressed in the same local tangent frame as the shape-gradient components. The recorded native port maintains 4×4 initial-stress integration, while the shared geometric API permits different orders. This **does not** constitute a fully consistent nonlinear shell tangent: rotational/director stress terms, curved-shell nonlinear terms, and arbitrary material layering have not been established. A typical linear buckling calculation solves \$(\\mathbf K+\\lambda\\mathbf K_G)\\boldsymbol\\phi=0\$, with sign and load normalization made explicit.

## 7. Pointwise recovery and physical stress verification

At a recovery point \$(r,s)\$, use **the same active membrane/shear choices as the stiffness formulation**, not a raw field accidentally substituted into a tied operator:

\
```math
\\mathbf N_m=h\\mathbf C_m\\mathbf B_m^{\\mathrm{active}}\\mathbf d_e,
\\qquad
\\mathbf M=\\frac{h^3}{12}\\mathbf C_m\\mathbf B_b\\mathbf d_e,
\\qquad
\\mathbf Q=h\\mathbf C_s\\mathbf B_s^{\\mathrm{active}}\\mathbf d_e.
\
```


The membrane stress and bending fiber stresses are

\
```math
\\boldsymbol\\sigma_m=\\mathbf C_m\\mathbf B_m^{\\mathrm{active}}\\mathbf d_e,
\\qquad
\\boldsymbol\\sigma_{\\pm}=\\boldsymbol\\sigma_m\\pm\\frac h2\\mathbf C_m\\mathbf B_b\\mathbf d_e.
\
```


For the physical global stress tensor,

\
```math
\\boxed{\\boldsymbol\\sigma^{g}(r,s)=\\mathbf Q(r,s)\\boldsymbol\\sigma^{\\ell}(r,s)\\mathbf Q(r,s)^T.}
\
```


If an element-center frame is used for output, its tensor must be transformed **from the pointwise basis exactly once**. Local element stresses, corner stresses, nodal averages, and global tensor components are separate observables; unaveraged stress at the actual benchmark location is the primary regression quantity. A zero theoretical stress requires an absolute-error metric, not a percentage.

## 8. Numerical evidence and reproducibility boundaries

The numerical tables below reproduce the **supplied ANS8BDG6 V4 validation record**. Original CSVs and detailed DAT/OP2/F06 reports are indexed by the upstream [shell benchmark library](https://github.com/realbabilu/MYSTRAN/tree/v18.00.a/update/shell_benchmark). The public benchmark index is reference-anchored and warns against invalid older comparisons. Where no matching ANS8BDG6 case CSV could be independently identified in this document-generation session, the result is explicitly treated as an **archived validation statement**, not a newly verified run.

### 8.1 CSI 2-005: uniform pressure on the corrected C plate

The corrected independent target is \$U_z=-12.97\$ in for the simply supported C uniform-load case, **not** the older invalid A/B benchmark. The native reported 16×16 Q8 model contains 833 nodes and 256 elements.

| Model | Center \$U_z\$ (in) | Signed error magnitude against 12.97 in |
|---|---:|---:|
| Original native baseline | −13.51560497 | 4.2067% |
| Native ANS8BDG6 v4 | −12.97469330 | 0.0362% |
| Independent v4-reference implementation | −12.97463019 | 0.0357% |

The comparison establishes a measurable accuracy improvement relative to the CSI target on this **specific** flat plate and support/load combination; it does not independently prove all curved-shell pressure loads are correct.

![Pressure comparison](figures/fig04_pressure_reference.png)

*Figure 4. Error relative to the corrected CSI 2-005 C uniform-pressure target. This is a before/after comparison at one mesh, **not** a mesh-convergence curve.*

### 8.2 CSI 2-016: thin-plate buckling, first reported factor

| Mode, thickness \$h\$ | Native baseline \$\\lambda_1\$ | Native v4 \$\\lambda_1\$ | Archived critical-stress error (baseline → v4) |
|---|---:|---:|---|
| Out-of-plane shear, 0.01 | 0.0004593225 | 0.0004242379 | +8.0733% → −0.1817% |
| Out-of-plane bending, 0.01 | 0.001249732 | 0.001161134 | +7.2814% → −0.3242% |
| Out-of-plane shear, 1.0 | 400.8203 | 400.8203 | Unchanged; not assessed against the thin plate target |
| Out-of-plane bending, 1.0 | 1086.480 | 1086.480 | Unchanged; not assessed against the thin plate target |

The published CSI thin-plate stress references are approximately \$\\tau_{cr}=0.0425\$ ksi and \$\\sigma_{cr}=0.1165\$ ksi, refined in the supplied driver to \$0.042501034\$ and \$0.116491057\$ ksi. **The values above are eigenvalue factors, not stress values;** conversion involves the prescribed prestress load basis and thickness. The report records the corrected reference interpretation, while the public library retains other pending Nastran and normalization flags.

![Buckling reference errors](figures/fig05_buckling_errors.png)

*Figure 5. Before/after signed critical-stress errors on the 0.01-thickness plate, from the archived ANS8 v4 validation record. Not a mesh-convergence plot.*

### 8.3 CSI 2-017 and additional solver parity

The first in-plane column factors were reported as unchanged by the v4 upgrade: **126.5028** on the 2×12 configuration and **126.4854** on the 4×24 configuration. The CSI independent force reference is \$P_{cr}=127.91\$ kip; *do not* compare it directly to an eigenvalue factor without verifying the applied load basis. Further native comparisons in the supplied record report:

| Verification | Archived result | Interpretation |
|---|---|---|
| 60 geometry/thickness/parameter settings × 18 points | Maximum operator absolute difference \$3.109\\times10^{-15}\$ | Native/reference **operator parity**, not absolute engineering accuracy |
| Consistent / HRZ mass parity | Maximum difference \$8.882\\times10^{-16}\$ | Matrix parity only |
| 2-005 pressure equilibrium | Maximum full-displacement per-component relative difference \$4.900\\times10^{-6}\$ | Native/reference match, thin-model sensitivity |
| Six SOL105 cases | First-three eigenfactor relative difference \$4.114\\times10^{-7}\$ | Native/reference solver parity |
| Six SOL103 checks | Five-frequency maximum relative difference \$3.413\\times10^{-7}\$ | Modal parity |
| Thermal free/restrained and gradient cases | Reported pass, including two windings | Archived physical regression; exact per-point CSV needed for further tables |

No new **Nastran-versus-ANS8 v4 at identical mesh and recovery location** figure is asserted here without a matching CSV and extraction protocol.

## 9. Patch-test and locking trade-offs: mandatory release criterion

A nominally accurate pressure or buckling result is **not sufficient** for a general-purpose Q8 element. Three independent tests are required:

1. **Constant membrane/bending patch** on an irregular five-element topology, with only perimeter values prescribed; assess the interior displacement and pointwise membrane stress/bending moment, in identical global axes. Reference CSI 2-001.
2. **Nonzero constant transverse shear** to ensure the shear operator produces physical \$\\mathbf Q\\ne0\$, and zero-shear quadratic bending to ensure fictitious transverse shear does not dominate.
3. **Thin-limit refinement and distortion sweeps** under force control: determine whether overconstraint induces shear/membrane locking, or underconstraint creates spurious flexibility or near-zero mechanisms. Compare eigen-rank and winding invariance in addition to displacement.

The v3 history reports that applying membrane tying uniformly can spoil the distorted flat patch, while disabling it for curved shells can lead to membrane locking. The v4 shear blend is an additional compromise between tensor6 locking and BDG4 excess compliance. **Passing the regular patch is not proof of passing a distorted patch**, and passing a displacement patch does not guarantee accurate recovered stress. The selection can be defensible for the tested regimes without being universally optimal.

The archived v4 memo reports several improvements but does **not** supply a complete independently audited distorted 2-001 stress/moment table for every v4 operator setting here; accordingly no universal 2-001 stress PASS is claimed. For comparison with the primary published paper, Jung and Han (2013) report normalized bending, shear, and tension results for their own βδγ choices, but their numbers must **not** be re-labeled as validation of the present Kikuchi–tensor6/BDG4 hybrid.

## 10. Limitations and reproducibility gates

- **Incomplete source attribution of some operators.** The supplied v4 class delegates strain/mass/pressure/buckling operators to `shell_q8_common` and inherits v2 through v3. The exact shared tensor6 coordinate table, modified-field coefficients, pointwise curved geometry conventions, and drilling factors should be checked against those implementation files before claiming source-exact completeness.
- **Curvature and stress.** Rotational/director geometric tangent terms are absent; pointwise stress axes and curved-shell pressure resultant invariance remain case-dependent. The archive acknowledges a physical first-moment issue for some modified-field curved pressure models.
- **Blending and objectivity.** The single scalar weight is geometry/thickness dependent and differentiable only away from its breakpoint; future nonlinear consistent tangents must account for this dependence if geometry updates change the blend.
- **Materials.** Orthotropic laminated composites, nonlinear constitutive response, and spatially variable thickness are not validated here; the original Jung–Han laminated-composite research must not be conflated with this isotropic native implementation.
- **Benchmark scope.** 2-009 is excluded by the public benchmark index; 2-015 lacks a confirmed independent reference; some 2-016/2-017 normalization and Nastran recovery states remain pending. Do not populate missing columns from other solvers.
- **Evidence levels.** Distinguish independent CSI reference agreement, Nastran agreement on exact matched mesh, MYSTRAN/reference operator parity, and isolated historical test values.

## 11. Figures supplied for insertion

1. `figures/fig01_q8_gauss.png` — Q8 node numbering and 3×3 Gauss points.
2. `figures/fig02_ans_tying.png` — β, δ, and γ classical tying patterns.
3. `figures/fig03_shear_blend.png` — dimensionless auto shear weight from the governing equation.
4. `figures/fig04_pressure_reference.png` — actual reported baseline/v4 pressure error versus CSI 2-005.
5. `figures/fig05_buckling_errors.png` — actual reported 2-016 thin plate critical-stress error.

## References

1. Jung, W.-Y., and Han, S.-C. (2013). **An 8-Node Shell Element for Nonlinear Analysis of Shells Using the Refined Combination of Membrane and Shear Interpolation Functions.** *Mathematical Problems in Engineering*, Article 276304. https://doi.org/10.1155/2013/276304 . **Primary theoretical basis for β/δ/γ sampling; not the exact v4 hybrid.**
2. Jung, W.-Y., and Han, S.-C. (2014). **Shear buckling responses of laminated composite shells using a modified 8-node ANS shell element.** *Composite Structures*, **109**, 119–129. https://doi.org/10.1016/j.compstruct.2013.10.055 . **Separate publication, not interchangeable with the 2013 Mathematical Problems in Engineering paper.**
3. Han, S.-C., Kim, K.-D., and Kanok-Nukulchai, W. (2004). Work on assumed-strain eight-node shell interpolation, *Structural Engineering and Mechanics*, **18**(6); confirm exact title and page range against the original article before formal submission.
4. MacNeal, R. H., and Harder, R. L. (1985). **A proposed standard set of problems to test finite element accuracy.** *Finite Elements in Analysis and Design*, **1**, 3–20.
5. Bathe, K.-J., and Dvorkin, E. N. (1986). **A formulation of general shell elements—the use of mixed interpolation of tensorial components.** *International Journal for Numerical Methods in Engineering*, **22**, 697–722. (Background on mixed/assumed interpolation.)
6. Hughes, T. J. R., and Brezzi, F. (1989). **On drilling degrees of freedom.** *Computer Methods in Applied Mechanics and Engineering*, **72**, 105–121.
7. Computers and Structures, Inc. **Shells Verification Problems**, particularly [2-001](https://docs.csiamerica.com/manuals/csibridge/Verification/Analysis/Shells/Problem%202-001.pdf), [2-005](https://docs.csiamerica.com/manuals/csibridge/Verification/Analysis/Shells/Problem%202-005.pdf), [2-016](https://docs.csiamerica.com/manuals/csibridge/Verification/Analysis/Shells/Problem%202-016.pdf), and [2-017](https://docs.csiamerica.com/manuals/csibridge/Verification/Analysis/Shells/Problem%202-017.pdf).
8. MYSTRAN project. [Quadratic shell benchmark library, v18.00.a](https://github.com/realbabilu/MYSTRAN/tree/v18.00.a/update/shell_benchmark); [reference-anchored comparison](https://github.com/realbabilu/MYSTRAN/blob/v18.00.a/update/shell_benchmark/comparison_all.md); [2-005 corrected benchmark](https://github.com/realbabilu/MYSTRAN/blob/v18.00.a/update/shell_benchmark/benchmark_2_005.md).
9. MYSTRAN project archive. **ANS8BDG6_V4_VALIDATION.md**, supplied validation record for the v4 Fortran port, with baseline and final executable/deck hashes and detailed case results.
