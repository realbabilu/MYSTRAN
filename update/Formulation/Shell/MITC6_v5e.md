# MITC6 Six-Node Triangular Shell — Research Formulation, Evolution, and MYSTRAN Validation Boundaries

**Author:** Ary  
**Date:** 23 October 2026  
**Status:** Source-audited technical research manuscript; not a claim of exact reproduction of the published MITC6 shell  
**Selected reproducible Python source:** `MITC6_Tri_v5e 
**MYSTRAN selector:** `PARAM,TRIA6TYP,MITC6`   


## Abstract

A 36-degree-of-freedom, six-node triangular shell is constructed with quadratic midsurface interpolation and mixed-interpolated covariant membrane and transverse-shear strains. The research lineage begins with an MITC-inspired edge scheme (v0), changes to an interior-sample affine membrane recovery and an eight-condition coupled shear projection (v1), corrects director/bending/drilling objectivity and internal kinematic consistency (v4), adds mass, thermal actions, surface loading, and linearized initial-stress stiffness (v5), and relocates two transverse-shear samples from an interior line to the true edge `r=0` (v5e). The element retains six conventional degrees of freedom per node, with no internal globally assembled unknowns. This report states the operator algebra, quadrature, force/mass definitions, implementation limitations, and source provenance separately from archived MYSTRAN validation claims.

## 1. Degrees of freedom, nodal ordering, and geometry

Each node has the *global* vector

$$
\mathbf d_I=[u_{xI},u_{yI},u_{zI},\theta_{xI},\theta_{yI},\theta_{zI}]^T, \qquad \mathbf d_e\in\mathbb R^{36}.
$$

The nodal order is **[corner 1, corner 2, corner 3, mid12, mid23, mid31]**, with natural triangle

$$
0\le r,\quad 0\le s,\quad r+s\le1; \qquad L_1=1-r-s,\; L_2=r,\; L_3=s.
$$

![MITC6 nodal numbering](figures/mitc6_nodes.svg)

The six quadratic shape functions are

$$
\begin{aligned}
N_1&=L_1(2L_1-1),&N_2&=L_2(2L_2-1),&N_3&=L_3(2L_3-1),\\
N_4&=4L_1L_2,&N_5&=4L_2L_3,&N_6&=4L_3L_1.
\end{aligned}
$$

The isoparametric midsurface, tangent vectors, surface Jacobian and first fundamental form are

$$
\mathbf x(r,s)=\sum_I N_I\mathbf x_I,\quad
\mathbf g_r=\sum_I N_{I,r}\mathbf x_I,\quad
\mathbf g_s=\sum_I N_{I,s}\mathbf x_I,
$$

$$
J_A=\|\mathbf g_r\times\mathbf g_s\|,\qquad
\mathbf A=\begin{bmatrix}\mathbf g_r\cdot\mathbf g_r&\mathbf g_r\cdot\mathbf g_s\\\mathbf g_s\cdot\mathbf g_r&\mathbf g_s\cdot\mathbf g_s\end{bmatrix}.
$$

The source computes contravariant coefficients from $\mathbf A^{-1}$ and evaluates a local orthonormal frame $\{\mathbf e_1,\mathbf e_2,\mathbf e_3\}$ at each sample. Director normals are evaluated at the six nodes. For its orientation policy, v4/v5/v5e chooses a **single element-wide sign** from the centroid normal, with a +global-Z tie break. This is an implementation convention, not a general guarantee of objective normal orientation on a closed curved surface. When `nodal_normals` are supplied, they must be coherent.

## 2. Natural covariant membrane and shear operators

Set $\mathbf B_{rr}^0,\mathbf B_{ss}^0,\mathbf B_{rs}^0\in\mathbb R^{1\times36}$ to the covariant membrane rows computed before tying. For node $I$, with translation entries occupying columns $6I:6I+3$ in zero-based indexing,

$$
\begin{aligned}
\mathbf B_{rr,I}^{0,u}&=N_{I,r}\mathbf g_r^T,\\
\mathbf B_{ss,I}^{0,u}&=N_{I,s}\mathbf g_s^T,\\
\mathbf B_{rs,I}^{0,u}&=\tfrac12(N_{I,r}\mathbf g_s^T+N_{I,s}\mathbf g_r^T).
\end{aligned}
$$

These rows have zero rotational coefficients. They encode covariant tensor components, with $\varepsilon_{rs}$ as tensor shear strain (not engineering shear).

The source defines an interpolated (not renormalized) director $\mathbf t=\sum_I N_I\mathbf n_I$. Its natural shear rows $\mathbf B_{rt}^0,\mathbf B_{st}^0$ have the node blocks

$$
\begin{aligned}
\mathbf B_{rt,I}^{0,u}&=N_{I,r}\mathbf t^T,&
\mathbf B_{rt,I}^{0,\theta}&=N_I(\mathbf n_I\times\mathbf g_r)^T,\\
\mathbf B_{st,I}^{0,u}&=N_{I,s}\mathbf t^T,&
\mathbf B_{st,I}^{0,\theta}&=N_I(\mathbf n_I\times\mathbf g_s)^T.
\end{aligned}
$$

The choice of an *unnormalized* $\mathbf t$ in **both** shear contributions is a v4 consistency correction, retained in v5 and v5e. Normalize only the local geometric normal as required for the local basis; do not normalize $\mathbf t$ in one term alone.

## 3. Membrane MITC-style projection: current uploaded v5e

Define the two sample parameters and centroid

$$
R_1=\tfrac12-\frac{1}{2\sqrt3}=0.2113248654,\quad
R_2=\tfrac12+\frac{1}{2\sqrt3}=0.7886751346,\quad
C=(\tfrac13,\tfrac13).
$$

For the **uploaded** v5e, membrane samples remain inherited from v1/v4. Three affine row fields are fitted independently:

$$
\widetilde{\mathbf B}_{rr}=\mathbf a_1+r\mathbf b_1+s\mathbf c_1,\qquad
\widetilde{\mathbf B}_{ss}=\mathbf a_2+r\mathbf b_2+s\mathbf c_2,
$$

$$
\widetilde{\mathbf B}_{qq}=\mathbf a_3+r\mathbf b_3+(1-r-s)\mathbf c_3,
\qquad \mathbf B_{qq}^0=\tfrac12(\mathbf B_{rr}^0+\mathbf B_{ss}^0)-\mathbf B_{rs}^0.
$$

For $rr$, fit values at $(R_1,R_1)$, $(R_2,R_1)$ and $C$. For $ss$, fit $(R_1,R_1)$, $(R_1,R_2)$ and $C$. For $qq$, fit $(R_2,R_1)$, $(R_1,R_2)$ and $C$. These are **interior-sample affine closures**, not edge sampling.

For each of the three groups, the coefficient vector is computed by exactly solving a $3\times3$ interpolation system:

$$
\begin{bmatrix}\mathbf a_k\\\mathbf b_k\\\mathbf c_k\end{bmatrix}
=\mathbf H_k^{-1}\begin{bmatrix}\mathbf b_k^0(P_1)\\\mathbf b_k^0(P_2)\\\mathbf b_k^0(C)\end{bmatrix},
$$

where rows of $\mathbf H_k$ are $(1,r,s)$ for $rr,ss$ and $(1,r,1-r-s)$ for $qq$. In each equation the right side consists of row vectors across all 36 DOFs. The remaining covariant in-plane shear row is

$$
\widetilde{\mathbf B}_{rs}=\tfrac12(\widetilde{\mathbf B}_{rr}+\widetilde{\mathbf B}_{ss})-\widetilde{\mathbf B}_{qq}.
$$

Let $a^{\alpha\beta}$ denote the inverse metric and define $c_i^{\alpha}=a^{\alpha\beta}(\mathbf g_\beta\cdot\mathbf e_i)$, with $i=1,2$ and $\alpha,\beta\in\{r,s\}$. The tensor-to-physical map is

$$
\begin{aligned}
\mathbf B_{m,xx}&=(c_1^r)^2\widetilde{\mathbf B}_{rr}+(c_1^s)^2\widetilde{\mathbf B}_{ss}+2c_1^rc_1^s\widetilde{\mathbf B}_{rs},\\
\mathbf B_{m,yy}&=(c_2^r)^2\widetilde{\mathbf B}_{rr}+(c_2^s)^2\widetilde{\mathbf B}_{ss}+2c_2^rc_2^s\widetilde{\mathbf B}_{rs},\\
\mathbf B_{m,xy}&=2[c_1^rc_2^r\widetilde{\mathbf B}_{rr}+c_1^sc_2^s\widetilde{\mathbf B}_{ss}+(c_1^rc_2^s+c_1^sc_2^r)\widetilde{\mathbf B}_{rs}].
\end{aligned}
$$

The final row stores **engineering** membrane shear strain $\gamma_{xy}=2\varepsilon_{xy}$, consistent with the isotropic plane-stress matrix below.

## 4. Coupled transverse-shear tying in v5e

![MITC6 v5e shear tying locations](figures/mitc6_tying.svg)

The two-component assumed covariant shear field is

$$
\begin{aligned}
\widetilde{\mathbf B}_{rt}(r,s)&=\mathbf a_1+\mathbf b_1r+\mathbf c_1s+\mathbf d\,rs+\mathbf e\,s^2,\\
\widetilde{\mathbf B}_{st}(r,s)&=\mathbf a_2+\mathbf b_2r+\mathbf c_2s-\mathbf d\,r^2-\mathbf e\,rs.
\end{aligned}
$$

There are eight unknown coefficient row vectors $[\mathbf a_1,\mathbf b_1,\mathbf c_1,\mathbf a_2,\mathbf b_2,\mathbf c_2,\mathbf d,\mathbf e]$. Eight sampling conditions define the interpolation:

| Observed covariant component | Sample 1 | Sample 2 |
|---|---|---|
| $rt$ | $(R_1,0)$ | $(R_2,0)$ |
| **$st$ — v5e corrected** | **$(0,R_1)$** | **$(0,R_2)$** |
| Common $rt$ and $st$ | $(1/3,1/3)$ | both components at the same point |
| $qt=(st-rt)/\sqrt2$ | $(R_2,R_1)$ | $(R_1,R_2)$ |

Write the observed eight-vector of 36-DOF rows as $\mathbf Y$ and the sampling/interpolation matrix as $\mathbf M\in\mathbb R^{8\times8}$. Then the coefficients follow the source algorithm:

$$
\mathbf C=\mathbf M^{-1}\mathbf Y,\qquad
\widetilde{\mathbf B}_{rt}=\mathbf p_{rt}(r,s)\mathbf C,\quad
\widetilde{\mathbf B}_{st}=\mathbf p_{st}(r,s)\mathbf C,
$$

$$
\mathbf p_{rt}=[1,r,s,0,0,0,rs,s^2],\qquad
\mathbf p_{st}=[0,0,0,1,r,s,-r^2,-rs].
$$

For $qt$ rows of $\mathbf M$, use $(\mathbf p_{st}-\mathbf p_{rt})/\sqrt2$, and sample $(\mathbf B^0_{st}-\mathbf B^0_{rt})/\sqrt2$. Convert to physical transverse shear using the same contravariant coefficients:

$$
\mathbf B_{s,x}=c_1^r\widetilde{\mathbf B}_{rt}+c_1^s\widetilde{\mathbf B}_{st},\qquad
\mathbf B_{s,y}=c_2^r\widetilde{\mathbf B}_{rt}+c_2^s\widetilde{\mathbf B}_{st}.
$$

**Exactly what changed v5 → v5e:** the only mathematical change in the two supplied files is replacing $st$ observations $(R_1,R_1)$ and $(R_1,R_2)$ with $(0,R_1)$ and $(0,R_2)$. The coefficient matrix $\mathbf M$ must be regenerated at these new points: moving the observations without updating its polynomial rows is incorrect. All other class logic is identical except the class name.

## 5. Bending and drilling stabilization

The v4 correction retained in v5/v5e defines bending through covariant curvature contributions, then transforms them to physical $xx,yy,xy$ components. At node $I$, let $\mathbf n_I$ be the nodal director and $\mathbf t_{,r}=\sum_J N_{J,r}\mathbf n_J$ (analogously in $s$). The covariant bending row contributions include

$$
\begin{aligned}
\mathbf B_{\rho,rr,I}^{\theta}&=N_{I,r}(\mathbf g_r\times\mathbf n_I)^T,&
\mathbf B_{\rho,rr,I}^{u}&=-N_{I,r}\mathbf t_{,r}^T,\\
\mathbf B_{\rho,ss,I}^{\theta}&=N_{I,s}(\mathbf g_s\times\mathbf n_I)^T,&
\mathbf B_{\rho,ss,I}^{u}&=-N_{I,s}\mathbf t_{,s}^T,\\
\mathbf B_{\rho,rs,I}^{\theta}&=\tfrac12[N_{I,r}(\mathbf g_s\times\mathbf n_I)+N_{I,s}(\mathbf g_r\times\mathbf n_I)]^T,\\
\mathbf B_{\rho,rs,I}^{u}&=-\tfrac12[N_{I,r}\mathbf t_{,s}+N_{I,s}\mathbf t_{,r}]^T.
\end{aligned}
$$

The crucial v4 repair is a **common curvature sign convention** in translation and rotational blocks. The element continues to use a spin-related drilling stabilization rather than a pure rotational penalty:

$$
\mathbf B_{d,I}^{u}=\tfrac12(N_{I,x}\mathbf e_2-N_{I,y}\mathbf e_1)^T,
\qquad \mathbf B_{d,I}^{\theta}=-N_I\mathbf e_3^T.
$$

Here $N_{I,x},N_{I,y}$ are physical tangent derivatives in the pointwise frame as computed by the source. The default numerical coefficient is $\beta_d=0.02$. The added energy is $\beta_dGh\int\mathbf B_d^T\mathbf B_d\,dA$. This is a continuum-spin-coupled penalty; it must not be misrepresented as a complete derivation of the original Hughes–Brezzi mixed variational formulation.

## 6. Stiffness and quadrature

With constant isotropic $E,\nu,h$, the in-plane constitutive matrix and the shear modulus are

$$
\mathbf C_m=\frac{E}{1-\nu^2}\begin{bmatrix}1&\nu&0\\\nu&1&0\\0&0&(1-\nu)/2\end{bmatrix},\quad
G=\frac{E}{2(1+\nu)},\quad k_s=\frac56.
$$

The source stiffness is

$$
\boxed{\mathbf K_e=\int_A\left[h\mathbf B_m^T\mathbf C_m\mathbf B_m+
\frac{h^3}{12}\mathbf B_b^T\mathbf C_m\mathbf B_b+
 k_sGh\mathbf B_s^T\mathbf B_s+
\beta_dGh\mathbf B_d^T\mathbf B_d\right]dA.}
$$

| Operator contribution | Integration in `k_local` |
|---|---|
| Membrane and bending | 3-point triangle |
| Transverse shear | 3-point triangle |
| Drilling | 6-point triangle |
| Mass, pressure in supplied v5e | 6-point triangle |
| Geometric stiffness | 6-point triangle |
| Thermal equivalent loads | 3-point triangle |

The three-point triangle positions are $(1/6,1/6)$, $(2/3,1/6)$, $(1/6,2/3)$ with weights $1/6$ each (summing to reference triangle area $1/2$). The six-point rule has three permutations of $(a,a)$, $(b,a)$ with $a=0.445948490144588$, $b=0.108103018168070$, and three permutations of $(c,c)$, $(d,c)$ with $c=0.091576213509771$, $d=0.816847572980459$, and associated weights $0.111690794839050$ and $0.054975871827661$. Their sum is also $1/2$.

**Implementation warning:** the original source uses scalar-average thickness `np.mean(self.h)` even when a six-value thickness array is provided; this is **not** a fully spatially variable-thickness shell formulation.

## 7. Consistent mass, HRZ mass, and NSM qualification

For material density $\rho$, the six-point consistent mass is

$$
(\mathbf M_e)_{I a,J b}=\delta_{ab}\int_A N_IN_J
\begin{cases}\rho h,&a=1,2,3,\\\rho h^3/12,&a=4,5,6\end{cases}dA.
$$

The v5 and v5e functions permit HRZ-type diagonalization per DOF family: diagonal entries of each consistent sub-block are scaled to match the corresponding block total. The supplied Python implementation can divide by zero when material density is zero, whereas the reported MYSTRAN adapter handles NSM-only models separately. **PSHELL NSM is not directly implemented in this standalone Python class**; native NSM rules must be checked in the MYSTRAN port.

## 8. Surface traction, normal pressure, and thermal actions

For normal pressure $p$ and fixed-direction traction $\mathbf q$ (units of force per physical surface area),

$$
\mathbf f_{I}^{p}=\int_A N_Ip\mathbf n\,dA,\qquad
\mathbf f_I^{q}=\int_A N_I\mathbf q\,dA.
$$

The source applies the v4/v5 oriented normal convention for $p$; its gravity function supplies $\mathbf q=\rho h\mathbf g$. All these loads create **translational nodal forces only**. Importantly, the *uploaded* v5e uses **six-point triangle pressure integration**, **not** the eight-per-direction Duffy integration reported in the separate enhanced MYSTRAN pressure-helper checkpoint. Treat that distinction as an unverified port difference.

For thermal expansion coefficient $\alpha$, uniform temperature change $\Delta T$, and through-thickness gradient $\Delta T_{,z}$ (positive along the source's chosen normal),

$$
\boldsymbol\varepsilon_{th}=\alpha\Delta T[1,1,0]^T,\qquad
\boldsymbol\kappa_{th}=-\alpha\Delta T_{,z}[1,1,0]^T,
$$

$$
\mathbf f_e^{th}=\int_A\left[h\mathbf B_m^T\mathbf C_m\boldsymbol\varepsilon_{th}+
\frac{h^3}{12}\mathbf B_b^T\mathbf C_m\boldsymbol\kappa_{th}\right]dA.
$$

The source membrane/moment recovery subtracts the thermal free strain and uses the corresponding bending sign. TEMPP1 physical top/bottom temperature interpretation and winding must be checked against the native solver's convention.

## 9. Initial-stress geometric stiffness and eigen buckling

The source recovers membrane stress resultants from the **active tied** membrane operator:

$$
\mathbf N_0=h\mathbf C_m[\mathbf B_m\mathbf d_e-\boldsymbol\varepsilon_{th}],\qquad
\mathbf S_0=\begin{bmatrix}N_{xx}&N_{xy}\\N_{xy}&N_{yy}\end{bmatrix}.
$$

Let $\mathbf G$ be the $2\times6$ matrix of physical tangent gradients of the **standard quadratic** shape functions. At each of the six integration points,

$$
\mathbf k_g=\mathbf G^T\mathbf S_0\mathbf G,\qquad
(K_G)_{I a,J b}=\delta_{ab}\int_A(k_g)_{IJ}\,dA \quad(a,b=1,2,3).
$$

The rotational geometric blocks are zero in this simplified linearized formulation. The corresponding eigenproblem is $(\mathbf K+\lambda\mathbf K_G)\boldsymbol\phi=0$ with tensile $\mathbf N_0$ positive. This is a **linear initial-stress buckling model**, not a complete finite-rotation or nonlinear shell tangent.

## 10. Version history: source-confirmed changes

| Revision | Main formulation state | Scope / caveat |
|---|---|---|
| `MITC6_Tri_v0.py` | Early MITC-inspired edge-tangential membrane and edge/coupled shear sampling | Source identifies a best-effort literature reconstruction; not a certified Kim–Bathe reproduction |
| `MITC6_Tri_v1.py` | Interior affine membrane sampling; eight-condition coupled shear interpolation; original `st` interior line | Kinematic inconsistencies persisted |
| `MITC6_Tri_v4.py` | One director orientation sign, same unnormalized interpolated director in both shear terms, bending translation-sign correction, point-normal-consistent drilling | Elastic-production baseline reported for native MYSTRAN before subsequent upgrade |
| `MITC6_Tri_v5.py` | Identical elastic stiffness operators to v4; consistent/HRZ mass, thermal, normal/directional pressure, gravity and geometric stiffness | Same interior `st` shear samples as v4 |
| **`MITC6_Tri_v5e.py`** | **Moves `st` shear observations from $(R_1,R_1),(R_1,R_2)$ to $(0,R_1),(0,R_2)$** | **All other functions match v5 aside from the class name** |

No changes in `MITC6_Tri_v5e.py` to the membrane interpolation are present in the uploaded source. Thus the phrase *“enhanced V5e”* is ambiguous unless accompanied by a source checksum and a definition of the active membrane operator.

## 11. Numerical validation: what the available evidence can and cannot establish

### 11.1 Archived native checkpoint, **not yet tied to this uploaded v5e source**

The earlier port summary supplied by the researcher records the following before/after **MYSTRAN V4 → enhanced checkpoint**. Because its deployed enhanced membrane differs from the uploaded source, the results must not be relabelled as head-to-head proof for this manuscript's Python v5e.

| CSI verification | Native V4 | Reported enhanced port | Independent reference / interpretation |
|---|---:|---:|---|
| 2-005 simply supported plate center $U_z$ | −13.06639671 | −12.97978592 | −12.9700 (reported comparison target) |
| 2-016 thin shear first eigenfactor, $h=0.01$ | 0.0004285837 | 0.0004234285 | Critical shear stress factor/$h$; reference 0.042501034 |
| 2-016 thin bending first eigenfactor, $h=0.01$ | 0.001173875 | 0.001158739 | Critical stress reference 0.116491057 |
| 2-017 column, 2×12 | 126.5009 | 126.5009 | Unchanged at F06 precision |
| 2-017 column, 4×24 | 126.4852 | 126.4852 | Unchanged at F06 precision |
| 2-006 linear roof $U_z$ | −0.30377501 | −0.30441040 | −0.3086 |
| 2-006 cylindrical roof $U_z$ | −0.28382701 | −0.30164146 | −0.3086 |

The checkpoint additionally reports 165-point operator comparisons ($6.141\times10^{-15}$ max difference), six SOL105 parity checks (relative difference $4.222\times10^{-7}$), six SOL103 checks (five frequencies within $3.670\times10^{-7}$), thermal and 36 pressure tests. These are **reported archival checks**, not newly executed in producing this manuscript, and they refer to the reported enhanced deployment rather than a separately established binary/source match to the uploaded `MITC6_Tri_v5e.py`.

### 11.2 Required independent comparisons and test definitions

| Benchmark | Primary observable | Release/interpretation gate |
|---|---|---|
| CSI 2-001 | Interior prescribed-field displacement, membrane stresses, moments and shear resultants | Compare same global tensor frame, local fiber and recovery point; inspect all components including zero-shear consistency |
| CSI 2-002 | Cantilever tip responses, six loading cases | No blanket Rz SPC if assessing drilling behavior |
| CSI 2-003 | Curved-beam $U_y,U_z$ convergence | Same midside-on-curve geometry and identical nodal loads |
| CSI 2-004 | Twisted-beam tip response | Distinguish chord-midpoint from true twisted midside surface |
| CSI 2-005 | Uniform/point-loaded quarter plate | Exact center, symmetry supports, and consistent element load; use connected line-convergence plots |
| CSI 2-006 | Scordelis–Lo roof vertical displacement | Verify normal vs global vertical gravity and projected-vs-surface load area |
| CSI 2-008 | Modal frequencies and shapes | Compare mass/NSM, consistent vs HRZ, mode MAC |
| CSI 2-013 | Free/restrained thermal and temperature gradients | Match thermal sign, winding, top/bottom fiber conventions |
| CSI 2-016 | Shear/bending plate buckling | Normalize first eigenfactor by actual prestress and thickness before comparing independent critical stress |
| CSI 2-017 | Column buckling first mode | Check actual force basis and matching eigenmode, not factor alone |

No numerical Nastran-v5e comparisons are supplied in this package. The [MYSTRAN comparison index](https://github.com/realbabilu/MYSTRAN/blob/v18.00.a/update/shell_benchmark/comparison_all.md) is a repository pointer, **not** evidence of a matching direct row unless its deck, formulation, observation and coordinate basis are traced.

## 12. MYSTRAN implementation notes and known limits

The reported selector `PARAM,TRIA6TYP,MITC6` identifies the active Fortran element family, but **not by itself** whether the compiled solver uses uploaded Python v5e interior membrane or the separately reported edge-membrane enhanced checkpoint. Verify the compiled `CTRIA6_MITC6.f90` operator and executable hash before describing production as mathematically identical to this source. The port summary reports original and upgraded binaries with SHA256 `80cccf87736c26a192beb675ca2f4e2cd8450b6c40c65e22af77e53bc7f0e714` and `13ce3108c61f9ac75c1c92ca576caebf58b3524bbabbd6d658d8faa287d28cae`, respectively; hashes establish archive identity, not the mathematical relationship to an unattested source.

The supplied code assumes homogeneous isotropic small-strain shell material, scalar-average thickness, standard six global nodal DOFs, and linearized membrane-prestress geometric stiffness. It does not prove complete behavior for anisotropic materials, layered sections, arbitrary nonlinear constitutive response or finite rotations. The `T_matrix` is the identity because its DOFs are global, while operator recovery uses pointwise physical tangent directions. Stress tensors must be rotated into basic axes **once** before solver-to-solver comparisons.

## 13. Comparison with published MITC6 formulations (source-to-literature audit)

**Scope and evidence.** The comparison below is made by examining the uploaded `MITC6_Tri_v5e.py` against the three full papers supplied by the researcher: Lee and Bathe (2004), Kim and Bathe (2009), and Lee and Bathe (2010). The purpose is to identify **mathematical equivalence and differences**, not to imply that the project's code is the authors' original implementation. Numerical comparisons reported in the articles are results of the published elements, **not measurements of the project v5e**.

### 13.1 Literature lineage and what is actually implemented

| Topic | Lee–Bathe 2004 | Kim–Bathe 2009 | Lee–Bathe 2010 | Uploaded Python `v5e` | Consequence |
|---|---|---|---|---|---|
| Main subject | MITC3, MITC6-a and MITC6-b isotropic triangular shells | Improved MITC6 shell addressing curved-shell instability | Quadratic MITC *plate* versus *shell* in plate bending | Linear research T6 shell with six global DOFs per node | The papers represent different baselines and must not be treated as interchangeable |
| Membrane interpolation | Three affine normal-strain fields for `rr`, `ss`, and `qq`, tied using edge-related samples and interior closure (2004, §3.3 and Fig. 6) | Retains polynomial forms but changes covariant basis used during interpolation (2009, §3) | Discusses the MITC6 shell membrane projection | **Interior-sample** affine fitting in `_compute_Bm` with centroid closure; not the original 2004 edge layout | Not an exact reproduction of published membrane tying |
| Shear interpolation | MITC6-a coupled quadratic/cross terms; three edge families and centroid conditions (2004, §4, Fig. 7a, Eqs. 36–39) | Keeps polynomial and tying layout, but evaluates components in an element-constant contravariant basis (2009, Eqs. 12–13) | Same MITC6 shell shear family studied numerically | Coupled eight-parameter fit of `rt`, `st` and `qt`; samples `st` at `(0,R1),(0,R2)` | Similar *polynomial family* and edge sampling to MITC6-a, not proof of operator equivalence |
| Covariant/contravariant basis | Uses the original geometry-dependent covariant tensor components | **Central improvement:** strain interpolation is referred to a basis fixed at the element centroid, while source strains retain variable geometry | Uses improved quadratic MITC shell in comparison | `_natural_membrane_rows`/`_natural_shear_rows` use pointwise `g_r,g_s`; `_covariant_maps(r,s)` maps at each evaluation point | **Centroid-constant contravariant interpolation is absent** from this Python implementation |
| Drilling | Published six-node shell is a continuum-mechanics-based shell with two in-surface director rotations | Improved element expressly avoids an artificial stabilization constant | Focuses on plate/shell quadratic bending | Uses six global DOFs per node including drilling and `beta_drill=0.02` penalty | Numerical stiffness is additionally affected by a *project-specific* drilling term |
| Transverse-shear correction | Depends on published shell convention | Depends on published shell convention | Explicitly states `kappa=1` for compared MITC plate and shell solutions (2010, §2) | `k_shear=5/6` | Literature plate-bending curves cannot be compared as exactly identical setups without changing this parameter |
| Modal/pressure/thermal/buckling | Not a source for the project's specific operators | Not a source for the project's specific operators | Plate-bending numerical study, not a code verification of such features | Added in v5 and inherited in v5e: six-point mass, HRZ, surface traction, thermal force, simplified initial-stress `Kg` | These are documented as implementation extensions, not claims directly from the papers |

### 13.2 Comparison of the shear approximation spaces

The uploaded Python formulation interpolates the two covariant shear rows through eight unknown coefficient rows:

$$
\begin{aligned}
\widetilde B_{rt}(r,s)&=a_1+b_1r+c_1s+d\,rs+e\,s^2,\\
\widetilde B_{st}(r,s)&=a_2+b_2r+c_2s-d\,r^2-e\,rs,\\
\widetilde B_{qt}(r,s)&=\frac{\widetilde B_{st}-\widetilde B_{rt}}{\sqrt2}.
\end{aligned}
$$

These are **operator rows**, rather than already-evaluated scalar strains: the coefficients are `1×36` vectors, solved from an `8×8` interpolation system. This polynomial structure agrees with the **family of coupled MITC6-a shear fields** discussed by Lee and Bathe (2004, Eq. 39), and reviewed by Kim and Bathe (2009, Eq. 11). However, equality of the polynomial space does not ensure equality of the operator: the component definitions, interpolation basis, exact tying points, and physical-coordinate transformation must also agree.

For the supplied `v5e`, the eight interpolation conditions are:

| Condition | Natural coordinate | Source sampled |
|---|---|---|
| `rt` 1, 2 | `(R1,0)`, `(R2,0)` | `_natural_shear_rows(...)[0]` |
| `st` 1, 2 | `(0,R1)`, `(0,R2)` | `_natural_shear_rows(...)[1]` |
| `rt`, `st` common | `(1/3,1/3)` (twice) | Both natural rows |
| `qt` 1, 2 | `(R2,R1)`, `(R1,R2)` | `(B_st-B_rt)/sqrt(2)` |

with $R_{1,2}=\tfrac12\mp(2\sqrt3)^{-1}$. The transition **v5 → v5e** changes only the two `st` observation coordinates. It does **not** implement Kim–Bathe's centroid-frozen contravariant components, nor does it replace the existing v4/v5 membrane operator.

### 13.3 Critical distinction: the 2009 basis correction

Kim and Bathe (2009, **page 1456, Eqs. 12–13**) retained the shapes of the earlier MITC6 strain polynomials and their tying pattern but replaced the *basis in which the components are interpolated*. Schematically, the updated approach defines centroid-frozen contravariant basis vectors

$$
\widehat{\mathbf g}^{\,i}(r,s,t)
=\mathbf g^{i}(1/3,1/3,t),
\qquad
\widehat e_{ij}(r,s,t)
=\varepsilon_{kl}(r,s,t)\bigl(\mathbf g^k\!\cdot\widehat{\mathbf g}^{\,i}\bigr)
\bigl(\mathbf g^l\!\cdot\widehat{\mathbf g}^{\,j}\bigr),
$$

followed by interpolation of these **hatted components** at the tying locations. The above tensor-change-of-basis identity is explanatory notation; use the paper's Eqs. (12)–(13) for the authors' precise through-thickness definitions.

In contrast, uploaded v5e forms the raw rows at the **sampling point** using `g_r(r,s), g_s(r,s)` and subsequently maps the tied rows with a pointwise orthonormal frame. A pointwise frame transformation after interpolation is *not* equivalent to tying components expressed in a fixed centroid contravariant basis. Thus it is technically incorrect to describe Python v5e as the **improved Kim–Bathe 2009 formulation**. The difference may be insignificant on straight, flat elements with constant bases, but can matter on irregular curved geometries; Kim and Bathe's discriminating hyperboloid examples motivated the change.

### 13.4 Membrane tying: Python is not the published edge layout

Both the 2004 MITC6-a scheme and the 2009 paper's retained membrane polynomial use

$$
\widetilde B_{rr}=a_1+b_1r+c_1s,\quad
\widetilde B_{ss}=a_2+b_2r+c_2s,\quad
\widetilde B_{qq}=a_3+b_3r+c_3(1-r-s),
$$

with $\widetilde B_{rs}=\tfrac12(\widetilde B_{rr}+\widetilde B_{ss})-\widetilde B_{qq}$. **Identical polynomial forms do not imply identical sampling.** Lee–Bathe (2004, **Fig. 6 and Eqs. 30–34**) develops edge-oriented strain samples. Uploaded v5e instead constructs affine fits from interior `(R1,R1)`, `(R2,R1)`, `(R1,R2)` and centroid conditions. This is inherited from v1/v4/v5 and is deliberately retained in the present baseline. No edge-based membrane variant is substituted in this manuscript.

### 13.5 Published tests versus source-specific validation

Kim–Bathe (2009, **Fig. 4/Table 1**) reports isotropy, absence of spurious zero-energy modes, and membrane/bending patch passes for the **published original and improved** elements. These outcomes cannot be transferred directly to uploaded v5e. The 2009 paper also emphasizes a hyperboloid that exhibited oscillatory displacement under the older formulation; this is an appropriate additional curved-shell discrimination test for the project. Lee–Bathe (2010, **Figs. 4–11**) shows that transverse-shear **stress recovery** can be much more difficult than convergence of deflections, and highlights sampling-location sensitivity. Consequently, MYSTRAN validation must check raw pointwise `Q_x,Q_y`, not displacement alone.

Recommended same-model tests for an eventual explicit implementation comparison: (i) a flat five-patch membrane/bending and nonzero-shear test with six rigid-body vectors; (ii) isosceles-triangle node-permutation/isotropy tests; (iii) distorted and curved MITC6 shells with different thicknesses; (iv) the 2009 clamped/free hyperboloid cases; and (v) identical plate meshes with `k_shear=5/6` versus `1`, with direct edge-midpoint shear recovery. Those are **proposed tests**, not newly executed results.

### 13.6 Engineering conclusion and publication wording

The project's **uploaded Python v5e** is best described as a *MITC6-inspired, quadratic six-node shell using a source-specific affine membrane tying projection and an eight-condition coupled assumed transverse-shear projection, augmented with project drilling, mass, loading, thermal, and initial-stress operators*. It is **not asserted** to be the complete improved Kim–Bathe 2009 shell. The primary source references support the MITC methodology and polynomial design, but do not substitute for operator-by-operator mathematical and numerical parity checks.

---

## 14. Conclusions

The delivered `MITC6_Tri_v5e.py` differs from `MITC6_Tri_v5.py` solely through the physically meaningful relocation of two **transverse shear** tying samples to the edge $r=0$. It retains the same interior/centroid-affine **membrane** recovery as v4/v5. Its principal mathematical identity is a source-specific MITC-inspired shell with two-component coupled shear projection, director-consistent bending/drilling, and supplementary linear mass/pressure/thermal/buckling operators. The distinct enhanced-MYSTRAN checkpoint may show improved benchmark behavior, but those numerical results cannot yet be attributed to the uploaded source without an exact operator match. Accordingly, no universal accuracy superiority over v5, SimoT6, or the independent REZAIEE element is asserted here.

## References

1. Kim, D.-N., and Bathe, K.-J. (2009). **A triangular six-node shell element.** *Computers & Structures* **87**, 1451–1460. DOI: [10.1016/j.compstruc.2009.05.002](https://doi.org/10.1016/j.compstruc.2009.05.002). Primary MITC6 motivation and original tying concept. Bibliographic page details should be cross-checked with the publisher's version before formal publication.
2. Lee, P.-S., and Bathe, K.-J. (2004). **Development of MITC isotropic triangular shell finite elements.** *Computers & Structures* **82**, 945–962. DOI: [10.1016/j.compstruc.2004.02.004](https://doi.org/10.1016/j.compstruc.2004.02.004). Source for MITC6-a, MITC6-b and original tying patterns.
3. Lee, P.-S., and Bathe, K.-J. (2010). **The quadratic MITC plate and MITC shell elements in plate bending.** *Advances in Engineering Software* **41**, 712–728. DOI: [10.1016/j.advengsoft.2009.12.011](https://doi.org/10.1016/j.advengsoft.2009.12.011). Plate/shell convergence and shear recovery comparison.
4. Bathe, K.-J., and Dvorkin, E. N. (1986). **A formulation of general shell elements—the use of mixed interpolation of tensorial components.** *International Journal for Numerical Methods in Engineering* **22**, 697–722. MITC methodology.
5. Hughes, T. J. R., and Brezzi, F. (1989). **On drilling degrees of freedom.** *Computer Methods in Applied Mechanics and Engineering* **72**, 105–121. Background on drilling, not direct proof that the coded penalty is the same mixed method.
6. Computers and Structures, Inc. **SAP2000/CSiBridge Verification Manual**, shell problems [2-001](https://docs.csiamerica.com/manuals/csibridge/Verification/Analysis/Shells/Problem%202-001.pdf), [2-005](https://docs.csiamerica.com/manuals/csibridge/Verification/Analysis/Shells/Problem%202-005.pdf), [2-006](https://docs.csiamerica.com/manuals/csibridge/Verification/Analysis/Shells/Problem%202-006.pdf), [2-016](https://docs.csiamerica.com/manuals/csibridge/Verification/Analysis/Shells/Problem%202-016.pdf), and [2-017](https://docs.csiamerica.com/manuals/csibridge/Verification/Analysis/Shells/Problem%202-017.pdf).
7. Primary project implementation archive: `MITC6_Tri_v0.py`, `MITC6_Tri_v1.py`, `MITC6_Tri_v4.py`, `MITC6_Tri_v5.py`, `MITC6_Tri_v5e.py` (source supplied to this study, 2026).
8. Internal MYSTRAN V5e upgrade checkpoint supplied by the researcher, October 2026; not independently verified against this uploaded source.

**Research attribution note:** Rezaiee-Pajand's separate triangular element is maintained as an independent implementation in this project and is **not represented as the origin or a synonym of the MITC6 v5e formulation**.
