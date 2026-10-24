# κMITC3+: A Three-Node Triangular Shell with Allman Drilling Membrane, Mixed Shear, and Linear Composite Curvature

**Ary | 25 October 2026 | Mathematical Formulation and Numerical Assessment**

> **Scope.** A formulation-oriented research note based on Lee et al. (2014) and Nguyen (2024). Numerical comparisons describe shell-element performance and limitations; model adaptations are identified where relevant.

## Abstract

The κMITC3+ triangular shell combines an Allman-type membrane interpolation incorporating true drilling rotation, a membrane-spin rank constraint, MITC3+ assumed transverse shear and a first-order composite-curvature projection. The element has three corner nodes with six external degrees of freedom per node and two condensed internal rotational bubble variables. Its numerical properties are evaluated against a Lee-type MITC3+ element with consistent-spin stabilization using affine patch fields, torsion, plates, roof and hemisphere benchmarks, and in-plane column buckling. Enhanced in-plane accuracy and rapid buckling convergence are accompanied by pronounced underprediction on a coarse doubly curved hemisphere. The findings support a qualified assessment of benefits and limitations rather than a claim of universal superiority.

## 1. Introduction

Three-node shell elements enable efficient unstructured meshing but must address membrane overstiffness, transverse-shear locking, drilling rotational rank and sensitivity to facet curvature. The κMITC3+ formulation combines an Allman-like membrane field with a true nodal drilling rotation, mixed-interpolated transverse shear and a linear projected curvature field. Its response is compared with the Lee–Lee–Bathe MITC3+ formulation equipped with consistent-spin stabilization. The discussion concerns numerical behavior rather than implementation architecture.

## 2. Parent triangle, local basis and degrees of freedom

Let the reference triangle be

$$
T=\{(r,s): r\ge 0,\ s\ge 0,\ r+s\le 1\},\qquad
N_1=1-r-s,\quad N_2=r,\quad N_3=s.
$$

For physical coordinates $\mathbf x_i$, define a right-handed facet frame

$$
\mathbf e_x=\frac{\mathbf x_2-\mathbf x_1}{\lVert\mathbf x_2-\mathbf x_1\rVert},\qquad
\mathbf e_z=\frac{(\mathbf x_2-\mathbf x_1)\times(\mathbf x_3-\mathbf x_1)}{\lVert(\mathbf x_2-\mathbf x_1)\times(\mathbf x_3-\mathbf x_1)\rVert},\qquad
\mathbf e_y=\mathbf e_z\times\mathbf e_x.
$$

The external nodal vector is $\mathbf d_i=[u_i,v_i,w_i,\theta_{xi},\theta_{yi},\theta_{zi}]^T$, with $18$ element DOFs. Two relative internal rotation-bubble variables $\boldsymbol\beta\in\mathbb R^2$ are condensed locally. The polynomial bubble is $N_b=27rs(1-r-s)$; it does not create global nodes or DOFs.

![Three corner nodes and rotational bubble](figures/kmitc3_nodal_bubble.svg)

*Figure 1. Three-node connectivity and the internal rotation-bubble convention, expressed in a right-handed local frame.*

## 3. Constitutive operators

For isotropic plane stress,

$$
\mathbf C_0=\frac{E}{1-\nu^2}
\begin{bmatrix}1&\nu&0\\\nu&1&0\\0&0&(1-\nu)/2\end{bmatrix},\quad
\mathbf D_m=h\mathbf C_0,\quad
\mathbf D_b=\frac{h^3}{12}\mathbf C_0,\quad
\mathbf D_s=\kappa_s Gh\mathbf I_2,
$$

where $G=E/[2(1+\nu)]$ and the adopted formulation uses $\kappa_s=5/6$. Membrane, bending and shear contributions are assembled independently before internal-mode condensation.

### 3.1 Modified Allman membrane and actual drilling DOF

The three-node constant-strain-triangle (CST) membrane alone has constant derivatives and does not provide a genuine membrane role for $\theta_z$. The Nguyen-inspired branch adopts a modified Allman interpolation,

$$
\boxed{\mathbf B_m^{A}(r,s)=\mathbf B_{\mathrm{CST}}-\alpha_A\mathbf H_M(r,s)\mathbf G_M,\qquad \alpha_A=\tfrac23.}
$$

Here $\mathbf G_M$ maps cyclic differences of drilling rotations and $\mathbf H_M$ maps these rotations into membrane strains. The coefficient $\alpha_A$ controls the membrane interpolation independently of the consistent-spin stabilization parameter.

The complementary mean-spin constraint is a rank-controlling operator, not an added arbitrary numerical diagonal for $\theta_z$:

$$
\mathbf b_{\mathrm{Add},i}=
\begin{bmatrix}N_{i,y}/2&-N_{i,x}/2&0&0&0&1/3\end{bmatrix},\quad
\boxed{\mathbf K_{\mathrm{Add}}=GhA\,\mathbf b_{\mathrm{Add}}^T\mathbf b_{\mathrm{Add}}.}
$$

The block row above denotes the contribution from each node in the 18-column assembled operator. Its normalization, rigid-rotation annihilation and physical units are mandatory patch-test checks. Do **not** add a separate Hughes–Brezzi or Lee consistent-spin stabilization on top without an explicit independent derivation: the Allman and consistent-spin contributions are assessed as distinct variants.

### 3.2 Six interior tying points for transverse shear

The assumed shear operator $\widetilde{\mathbf B}_s(r,s)$ is formed from six interior locations A–F in the Lee MITC3+ family, retaining contributions from nodal rotations, transverse displacement and the rotational bubble. The numerical interior offset in the reconstruction is $d=10^{-4}$. The parent integration uses the symmetric three-point rule

$$
(r,s,w)=\left(\tfrac16,\tfrac16,\tfrac13\right),\ \left(\tfrac23,\tfrac16,\tfrac13\right),\ \left(\tfrac16,\tfrac23,\tfrac13\right),
\qquad \mathbf K_s=\sum_g A w_g\,\widetilde{\mathbf B}_s^T\mathbf D_s\widetilde{\mathbf B}_s.
$$

![Lee interior shear tying points](figures/kmitc3_shear_tying.svg)

*Figure 2. Lee-style interior shear samples; these are tying locations, not the integration rule.*

![Three symmetric Gaussian integration points](figures/kmitc3_gauss_3point.svg)

*Figure 3. Symmetric three-point triangular integration used by the current parent operator.*

The **constrained-displacement shear patch** ratio is $1.13408$, not unity. This limitation must be retained, and a separate equilibrated shear-traction/moment-gradient patch is still needed before evaluating whether it is a true consistency failure.

### 3.3 Hu–Washizu first-order composite-curvature projection

The compatible bending curvature $\boldsymbol\kappa=\mathbf B_b(r,s)\,[\mathbf d;\boldsymbol\beta]$ is projected onto the first-order polynomial basis $\mathbf f=[1,r,s]^T$. Define

$$
\mathbf F=\int_A\mathbf f\mathbf f^T\,dA,\qquad
\mathbf Q=\int_A\mathbf f\otimes\mathbf B_b\,dA,\qquad
\boxed{\widetilde{\mathbf B}_b(r,s)=\mathbf f(r,s)^T\mathbf F^{-1}\mathbf Q.}
$$

The projection integrals $\mathbf F$ and $\mathbf Q$ are evaluated using symmetric triangular quadrature. Divide the parent triangle around $G=(1/3,1/3)$ into $T_1=(G,1,2)$, $T_2=(G,2,3)$, $T_3=(G,3,1)$. If $G_a$ denotes the centroid of the subtriangle $T_a$, then

$$
\boxed{\mathbf K_b^{P1}=\sum_{a=1}^{3}\frac{A}{3}
\widetilde{\mathbf B}_b(G_a)^T\mathbf D_b\widetilde{\mathbf B}_b(G_a).}
$$

The seven-point quadrature sensitivity for constructing $\mathbf F,\mathbf Q$ is a **numerical integration variation**, not part of the baseline comparison. It changed 2-002 torsion from 1.12967 to about 1.04990 without correcting the coarse hemisphere n8 response (~0.2533); there is insufficient justification to tune the default based on one case.

## 4. Internal condensation, temperature and stress recovery

With $\mathbf q=[\mathbf d;\boldsymbol\beta]$, form the symmetric $20\times20$ local matrix

$$
\mathbf K_{20}=\mathbf K_m^{A}+\mathbf K_{\mathrm{Add}}+\mathbf K_s^{Lee}+\mathbf K_b^{P1}
=\begin{bmatrix}\mathbf K_{cc}&\mathbf K_{cb}\\\mathbf K_{bc}&\mathbf K_{bb}\end{bmatrix},
$$

$$
\boxed{\mathbf K_{18}=\mathbf K_{cc}-\mathbf K_{cb}\mathbf K_{bb}^{-1}\mathbf K_{bc}},\qquad
\boldsymbol\beta=-\mathbf K_{bb}^{-1}\mathbf K_{bc}\mathbf d.
$$

The thermal membrane vector must use the *same* Allman strain matrix:

$$
\mathbf f_T=\int_A\mathbf B_m^{A\,T}\mathbf D_m\,\boldsymbol\varepsilon_T\,dA,
\qquad \boldsymbol\varepsilon_T=\alpha_T\Delta T[1,1,0]^T.
$$

A manufactured free-expansion test reports a single-element free-expansion force residual of $4.38\times10^{-17}$, which is **not** a validated SAP 2-013 thermal-stress comparison. Local recovery reports membrane force $\mathbf N=\mathbf D_m(\mathbf B_m^A\mathbf d-\boldsymbol\varepsilon_T)$, bending moments $\mathbf M=\mathbf D_b\widetilde{\mathbf B}_b\mathbf q$, and shear resultants $\mathbf V=\mathbf D_s\widetilde{\mathbf B}_s\mathbf q$. Any nodal averaging across shell facets requires a deliberate output coordinate frame; do not average incompatible local tensors directly.

## 5. Numerical benchmark results

All values below are **calculated/reference ratios**, with target unity. Values originate from matched numerical research models; adaptations are distinguished from the original reference models.

| Benchmark and resolution | Lee + Spin | κMITC3+ Nguyen | Technical interpretation |
|---|---:|---:|---|
| 2-001 affine membrane/bending/twisting | 1.000 | 1.000 | Prescribed-field checks passed |
| 2-001 constrained-shear patch | 1.13408 | 1.13408 | Remains open |
| 2-002 LC2 in-plane, n24 | 0.50438 | **0.80004** | Allman improvement; remains too stiff |
| 2-002 LC4 torsion, n24 | **0.98446** | 1.12967 | Nguyen overshoots by about 13% |
| 2-002 LC5 in-plane, n24 | 0.50434 | **0.80028** | Allman improvement |
| 2-003 UY, n16 | 0.27621 | 0.68224 | Substantial curvature/membrane deficiency remains |
| 2-004 Fz, n24 | 0.99611 | 0.99221 | Comparable on this benchmark |
| 2-005 clamped point load, n24 | 0.94622 | **0.95778** | Modest improvement |
| 2-006 Scordelis roof, n12 | 0.95856 | **0.98895** | Nguyen closer |
| **2-007 supplied original CTRIA3 DAT, n8** | **0.99193** | **0.25567** | **Serious coarse curved-shell failure** |
| 2-009 Winkler k=800, n16 | 0.77304 | 0.93415 | Better, but n24 overshoots to 1.05373 |

The membrane/curvature benefit is genuine within these benchmark setups, but **does not establish global superiority**: the coarsest hemisphere model is considerably worse than Lee + Spin.

### 5.1 Explicit distinction for hemisphere 2-007

| Parent divisions | Triangles | Nguyen ratio | Model provenance |
|---:|---:|---:|---|
| 8 | 128 | **0.25567** | **Original reference input** |
| 16 | 512 | 0.81553 | Geometrically reconstructed refinement |
| 32 | 2,048 | 0.98735 | Geometrically reconstructed refinement |
| 64 | 8,192 | 0.99820 | Geometrically reconstructed refinement |

Although reconstructed fine meshes approach unity, their geometry, boundary discretization and loading have not been proven identical to original benchmark DATs. Thus **do not report the n64 result as a successful original-DAT convergence series**. The convergence graph explicitly distinguishes Allman-only and Allman plus composite curvature.

![Hemisphere refinement response](figures/kmitc3_2007_convergence.png)

*Figure 4. Allman and Allman + composite refinement; n8 uses the reference input, while refined meshes are reconstructed.*

### 5.2 In-plane column buckling: Example 2-017

A 12 × 100 in cantilever column (thickness 1 in, $E=3600$ ksi, $\nu=0.2$) was discretized using two triangular elements per quadrilateral parent cell. The active degrees of freedom were $U_x$, $U_z$ and $R_y$. The independent critical load is $127.91$ kip. The table reports the first buckling load calculated using a common uniform compressive membrane-force geometric stiffness. This is a controlled triangular-element comparison, **not** a replication of the original quadrilateral formulation or a claim of a separately recovered pre-buckling stress field.

| Parent mesh | Triangles | MITC3+ Lee + Spin (kip) | κMITC3+ Allman + Composite (kip) |
|---|---:|---:|---:|
| 2 × 12 | 48 | 267.043 | 142.320 |
| 4 × 24 | 192 | 162.468 | 130.618 |
| 6 × 36 | 432 | 142.751 | 128.522 |
| 8 × 48 | 768 | 135.813 | 127.791 |
| 12 × 72 | 1,728 | 130.841 | 127.269 |
| 16 × 96 | 3,072 | 129.097 | 127.086 |

![Convergence of the first in-plane buckling eigenvalue](figures/kmitc3_2017_mode1_convergence.png)

*Figure 5. Convergence of the first critical load; the horizontal reference denotes $127.91$ kip.*

At the finest mesh, the first six positive eigenvalues were:

| Mode | Lee + Spin (kip) | Allman + Composite (kip) |
|---:|---:|---:|
| 1 | 129.097 | 127.086 |
| 2 | 1,089.563 | 1,073.240 |
| 3 | 2,691.635 | 2,653.901 |
| 4 | 4,525.035 | 4,466.186 |
| 5 | 6,288.716 | 6,211.788 |
| 6 | 7,836.735 | 7,742.962 |

Because only in-plane degrees of freedom were active, this example isolates the membrane/interpolation advantage; it does not independently demonstrate the benefit of composite curvature in out-of-plane bending. Its first-mode agreement also does not resolve the coarse curved-shell limitation reported in Section 5.1.

## 6. Discussion: drilling, membrane behavior and curved-shell limitations

### 6.1 Why the Allman drilling mechanism is not Hughes–Brezzi penalty stabilization

The added nodal rotation $\theta_z$ does not enter an ordinary CST membrane. In a conventional rank-fix design, its energetic participation may be supplied only through a small penalty or an auxiliary spin stiffness. In κMITC3+, in contrast, the off-diagonal edge-difference operator $\mathbf G_M$ makes $\theta_z$ affect $\boldsymbol\varepsilon_m$ directly. The mean continuum-spin constraint $\mathbf K_{\mathrm{Add}}$ complements this membrane interpolation to avoid the otherwise deficient drilling rank. It should therefore not be interpreted as adding an unrelated duplicate penalty. The Allman formulation does not include an additional Lee-type consistent-spin penalty.

### 6.2 Why an enriched membrane can improve in-plane bending but harm a coarse hemisphere

The approximately 0.50-to-0.80 improvement in 2-002 LC2/LC5 is consistent with additional membrane kinematics reducing CST overstiffness. The deterioration in the coarse 2-007 hemisphere is equally important: curved response depends on the variation of the local shell basis across facets, the geometry discretization, drilling–membrane coupling and possible membrane/bending contamination. The available tests do not uniquely isolate which mechanism dominates. Since the n8 original DAT and n16/n32/n64 reconstructed meshes are not an identical-origin mesh family, the observed refinement cannot be used to declare that a curved-shell limitation is resolved.

### 6.3 Interpretation of the remaining transverse-shear patch discrepancy

A computed/reference value of 1.13408 for the constrained-displacement shear test is a measurable excess stiffness under that test construction, but it is not sufficient to prove failure of a properly equilibrated constant-shear traction patch. A more discriminating verification would enforce compatible displacement and rotation fields, evaluate internal force consistency, and introduce a linear moment gradient where the governing equilibrium equations demand it. This is a scientific limitation of the present evidence and is recorded without tuning the shear operator.

## 7. Conclusions

1. The modified Allman field introduces a membrane-energy contribution from nodal drilling rotations, supplemented by a mean-spin rank constraint.
2. A linear projected curvature field and interior MITC shear tying provide the bending and transverse-shear ingredients of the κMITC3+ family.
3. The numerical evidence shows improvements in membrane-dominated tests and in-plane buckling convergence, alongside torsional overshoot and severe loss of accuracy in a coarse doubly curved hemisphere.
4. The outcomes are benchmark-dependent; no claim of unconditional superiority over MITC3+ Lee + Spin is justified.

## 8. References

1. Lee, Y., Lee, P.-S., and Bathe, K.-J. (2014). *The MITC3+ shell element and its performance.* **Computers & Structures**, 138, 12–23. https://doi.org/10.1016/j.compstruc.2014.02.005.
2. Nguyen, S. H. (2024). *Linear composite curvature MITC3+ flat shell elements.* **Computers & Structures**, 295, 107272. https://doi.org/10.1016/j.compstruc.2024.107272.
