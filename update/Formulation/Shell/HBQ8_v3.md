# HBQ8–Kikuchi v3: Modified Serendipity Surface Shell with Geometry-Adaptive Membrane Tying and Slenderness-Dependent Transverse-Shear Blending

**Ary · Independent Researcher · Personal Research** · MYSTRAN v18.00.a · 19 October 2026

**Mathematical formulation and numerical verification note — MYSTRAN quadratic shells**  
**Revision:** 20 October 2026 — explicit modified-field coefficients addendum  
**Evidence classification:** Mathematical operators audited from supplied HBQ8 base, Kikuchi/MacNeal v1, complete HBQ8 v2 and HBQ8 v3 source, shared Q8 operator module, and modified shape-function source; benchmark numbers reproduced from the user's reported runs unless explicitly noted. No new Python, MYSTRAN, or Nastran benchmark is claimed in this document.

## Abstract

A 48-degree-of-freedom, eight-node quadrilateral shell is described in which a Kikuchi-modified serendipity interpolation field is combined with a surface-covariant membrane operator, geometry-adaptive membrane assumed-strain interpolation, and slenderness-dependent transverse-shear blending. The formulation follows a four-stage research evolution: the original displacement-based HBQ8 flat shell, a Kikuchi/MacNeal modified-field substitution, a curved-surface extension with tensor-product assumed shear, and a revised v3 shell employing a controlled blend of six-point (`tensor6`) and four-point (`bdg4`) shear tying. The central numerical trade-off is explicit: six-point shear tying produces locking in very thin plates, whereas four-point tying reduces overconstraint but can lead to excessive flexibility at coarse resolution. For the reported simply supported quarter-plate pressure case with thickness $h=10^{-4}$, the $2\times2$ mesh center-deflection ratio changes from 0.095 in v2 to 0.990 in v3; the corresponding $4\times4$ ratios are 0.053 and 1.001, respectively. The modified element preserves the reported distorted-mesh displacement patch results, but full stress-recovery, curved geometric stiffness, and native MYSTRAN parity remain independent validation requirements.

**Keywords:** modified serendipity; Kikuchi interpolation; HBQ8; assumed natural strain; transverse shear locking; patch test; shell buckling; MYSTRAN.

## 1. Element genealogy, attribution and scope

The **implemented research lineage** is

$$
\text{HBQ8 base}\longrightarrow\text{KikuchiMacNeal--HBQ8 v1}
\longrightarrow\text{HBQ8--Kikuchi v2}
\longrightarrow\text{HBQ8--Kikuchi v3}.
$$

The first stage supplies a conventional six-degree-of-freedom-per-node membrane–Mindlin-bending–shear shell with a drilling penalty and $3\times3$ quadrature. In v1, the interpolating displacement field is replaced by **either** Kikuchi-modified Q8 **or** MacNeal–Harder, selected as alternative fields rather than sequential corrections. v2 adds a three-dimensional surface-tangent metric, covariant membrane/shear rows and a tensor-product tying operator. v3 retains this structure and introduces two independent *automatic* decisions: selective membrane tying on curved/warped geometries and a slenderness-dependent blend of transverse-shear interpolation patterns.

**Attribution boundary.** Darılmaz and Kumbasar (2006) proposed an eight-node *assumed-stress hybrid* shell based on the Hellinger–Reissner principle. The supplied `HBQ8_ShellElement` source is a **displacement-based** membrane/bending/shear assembly with a drilling penalty. No independent stress parameters or Hellinger–Reissner condensation appear in that source. Consequently the present v3 is a **project-developed hybrid inspired by several published ideas**, not a faithful implementation of the Darılmaz–Kumbasar variational element. Likewise, the thickness-dependent blend is not asserted to be a published Jung–Han or Kikuchi operator.

## 2. Q8 interpolation, geometry and 48 degrees of freedom

Let $(r,s)\in[-1,1]^2$, and let the four corner nodes be $(r_i,s_i)=(\pm1,\pm1)$ counter-clockwise. Standard Q8 geometric corner functions are

$$
N_i(r,s)=\tfrac14(1+rr_i)(1+ss_i)(rr_i+ss_i-1),\quad i=1,\ldots,4.
$$

With conventional numbering for midsides $(5,6,7,8)$,

$$
\begin{aligned}
N_5&=\tfrac12(1-r^2)(1-s), & N_6&=\tfrac12(1+r)(1-s^2),\\
N_7&=\tfrac12(1-r^2)(1+s), & N_8&=\tfrac12(1-r)(1-s^2).
\end{aligned}
$$

Define standard geometry and its covariant tangent vectors using **standard** $N_i$:

$$
\mathbf x(r,s)=\sum_{i=1}^8N_i(r,s)\mathbf x_i,\qquad
\mathbf g_\alpha^{(0)}=\sum_i N_{i,\alpha}\mathbf x_i,\qquad
dA_0=\|\mathbf g_1^{(0)}\times\mathbf g_2^{(0)}\|\,dr\,ds.
$$

Each node carries $\mathbf d_i=[u_i,v_i,w_i,\theta_{xi},\theta_{yi},\theta_{zi}]^T$ in the global Cartesian frame, hence $\mathbf d_e\in\mathbb R^{48}$. The modified displacement-interpolation functions $\widetilde N_i(r,s)$ are now defined explicitly from the supplied original Kikuchi/MacNeal field-function source. The field is modified, while **the geometry and surface-area Jacobian retain the standard Q8 interpolation**:

$$
\mathbf u^h(r,s)=\sum_{i=1}^8\widetilde N_i(r,s)\mathbf u_i,
\qquad
\boldsymbol\theta^h(r,s)=\sum_{i=1}^8\widetilde N_i(r,s)\boldsymbol\theta_i.
$$

### 2.1 Kikuchi–Okabe–Fujio correction: explicit coefficients

Write the complete nine-node biquadratic Lagrange interpolation on $[-1,1]^2$ as $L_i(r,s)$, $i=1,\ldots,9$, where nodes $1$–$8$ lie on the boundary and node $9$ lies at $(0,0)$. The center bubble is

$$
L_9(r,s)=(1-r^2)(1-s^2),\quad
L_{9,r}=-2r(1-s^2),\quad L_{9,s}=-2s(1-r^2).
$$

The eight retained Lagrange functions are tensor products of the one-dimensional basis

$$
\ell_{-1}(x)=\tfrac12x(x-1),\qquad
\ell_0(x)=1-x^2,\qquad
\ell_{+1}(x)=\tfrac12x(x+1),
\quad L_i=\ell_{r_i}(r)\ell_{s_i}(s).
$$

Rather than deleting the ninth shape function, Kikuchi-type interpolation redistributes its contribution to the eight boundary fields:

$$
\boxed{\widetilde N_i^{\mathrm K}(r,s)=L_i(r,s)+c_i L_9(r,s)},\qquad
\widetilde N_{i,\alpha}^{\mathrm K}=L_{i,\alpha}+c_i L_{9,\alpha}.
$$

For the physical **corner** coordinates $(x_i,y_i)$ expressed in one planar element basis, with cyclic indices $i=1,\ldots,4$, define

$$
D_i=(x_{i+1}-x_i)(y_{i-1}-y_i)
-(x_{i-1}-x_i)(y_{i+1}-y_i),
$$

where indices wrap cyclically, and let $k=i+2$, $m=i-1$ (modulo four). The supplied closed-form construction uses

$$
\boxed{c_i=-\frac14+\frac{D_i-D_k}{8(D_i+D_k)}},\qquad
\boxed{c_{i+4}=\frac12+\frac{D_m-D_i}{4(D_i+D_k)}}.
$$

The denominator $D_i+D_k$ must be nonzero. These expressions are the implemented *Kikuchi field coefficients*; they are not used as extra geometrical nodes or additional independent global DOFs. With these coefficients, interpolation identities must be checked for each supported quadrilateral, particularly severe distortion and winding reversal.

### 2.2 Alternative MacNeal–Harder field (not an additional sequential correction)

The source also supplies a selectable eight-node interpolation constructed using an eliminated ninth-node constraint. First determine the physical center point from the **standard Q8 geometry interpolation**:

$$
\mathbf x_9=\sum_{i=1}^8N_i^{\mathrm Q8}(0,0)\,\mathbf x_i,\qquad
\bar{\mathbf x}_i=\mathbf x_i-\mathbf x_9.
$$

Here $\mathbf x_i=(x_i,y_i)$ are **planar projected field coordinates**. Assemble the $8\times8$ interpolation matrix

$$
\mathbf X_{i,:}=\left[1,\bar x_i,\bar y_i,\bar x_i^2,
\bar x_i\bar y_i,\bar y_i^2,r_i^2s_i,r_is_i^2\right],
\qquad \mathbf T_{:}=(\mathbf X^{-1})_{1,:}.
$$

Then, with $a_i=N_i^{\mathrm Q8}(0,0)$,

$$
\boxed{\widetilde N_i^{\mathrm{MH}}=N_i^{\mathrm Q8}
+L_9(T_i-a_i)},\qquad
\widetilde N_{i,\alpha}^{\mathrm{MH}}=N_{i,\alpha}^{\mathrm Q8}
+L_{9,\alpha}(T_i-a_i).
$$

This alternative requires a nonsingular $\mathbf X$ and recentering at the eliminated node. In the selected **HBQ8–Kikuchi v3** branch the default field is Kikuchi; MacNeal–Harder remains an explicitly different field variant.

### 2.3 Separation of modified fields and standard geometry

For clarity, the implementation uses the following distinct interpolants:

$$
\mathbf x^{\mathrm{geom}}(r,s)=\sum_i N_i^{\mathrm Q8}\mathbf x_i,
\qquad
\mathbf u^{\mathrm{field}}(r,s)=\sum_i\widetilde N_i\mathbf u_i.
$$

Therefore

$$
\mathbf g^{\mathrm{geom}}_\alpha=\sum_i N_{i,\alpha}^{\mathrm Q8}\mathbf x_i,
\qquad dA=\left\|\mathbf g^{\mathrm{geom}}_1\times
\mathbf g^{\mathrm{geom}}_2\right\|\,dr\,ds.
$$

The field-tangent metric introduced subsequently in Section 3 is a separate construction using $\widetilde N_{i,\alpha}$. Mixing these two measures without acknowledging the distinction changes the discrete method.

The standard geometry integration surface and modified-field strain metric are intentionally different in v2/v3. Modified shape functions follow Kikuchi, Okabe and Fujio (1999), whose design motivation is to improve polynomial reproduction in distorted bilinear quadrilaterals; the precise shell extension used here is project specific.


<img width="1470" height="1260" alt="fig01_gauss_Q8" src="https://github.com/user-attachments/assets/116946a0-4a5f-425d-9c92-beeff96b36a6" />

**Figure 1.** Standard Q8 connectivity and $3\times3$ Gauss–Legendre integration stations $r,s\in\{-\sqrt{3/5},0,+\sqrt{3/5}\}$, weights $(5/9,8/9,5/9)$ per direction. These are integration sites, **not** ANS tying sites.

## 3. Modified-field surface metric and local basis

Following v2, the *field* tangents and their metric are

$$
\mathbf g_\alpha=\sum_i\widetilde N_{i,\alpha}\mathbf x_i,\qquad
\mathbf a=[\mathbf g_1;\mathbf g_2],\qquad
\mathbf G=\mathbf a\mathbf a^T,\qquad
\mathbf g^\alpha=\sum_\beta G^{\alpha\beta}\mathbf g_\beta.
$$

The pointwise orthonormal basis is constructed as

$$
\mathbf e_1=\frac{\mathbf g_1}{\|\mathbf g_1\|},\qquad
\mathbf n=\varsigma\frac{\mathbf g_1\times\mathbf g_2}{\|\mathbf g_1\times\mathbf g_2\|},\qquad
\mathbf e_2=\mathbf n\times\mathbf e_1.
$$

Here $\varsigma\in\{-1,+1\}$ is selected at the element center and retained throughout the element. Let

$$
C_{\alpha a}=\mathbf g^\alpha\cdot\mathbf e_a,\qquad
\widetilde N_{i,a}=\sum_\alpha C_{\alpha a}\widetilde N_{i,\alpha},\quad a=1,2.
$$

The practical significance is **frame coherence:** both strain reconstruction and output transformation must use consistent pointwise frames. The element stiffness is accumulated with **standard geometric surface measure $dA_0$**, not a silently substituted $\|\mathbf g_1\times\mathbf g_2\|$.

## 4. Membrane operator and selective tying

Define the covariant membrane strain operator acting on $\mathbf d_e$ by

$$
\widehat\varepsilon_{11}=\mathbf g_1\cdot\mathbf u^h_{,1},\qquad
\widehat\varepsilon_{22}=\mathbf g_2\cdot\mathbf u^h_{,2},\qquad
\widehat\varepsilon_{12}=\tfrac12(\mathbf g_1\cdot\mathbf u^h_{,2}+\mathbf g_2\cdot\mathbf u^h_{,1}).
$$

Thus the corresponding **raw** matrix rows, acting only on translational DOFs, read

$$
\begin{aligned}
\widehat{\mathbf B}_{m,11}^{(i)}&=\widetilde N_{i,1}\mathbf g_1^T,\\
\widehat{\mathbf B}_{m,22}^{(i)}&=\widetilde N_{i,2}\mathbf g_2^T,\\
\widehat{\mathbf B}_{m,12}^{(i)}&=\tfrac12(\widetilde N_{i,1}\mathbf g_2^T+\widetilde N_{i,2}\mathbf g_1^T).
\end{aligned}
$$

The raw/tied covariant tensor is mapped to the pointwise orthonormal frame by the $2\times2$ transformation $\mathbf C$:

$$
\begin{bmatrix}\varepsilon_{11} & \varepsilon_{12}\\\varepsilon_{12}&\varepsilon_{22}\end{bmatrix}
=\mathbf C^T
\begin{bmatrix}\widehat\varepsilon_{11}&\widehat\varepsilon_{12}\\\widehat\varepsilon_{12}&\widehat\varepsilon_{22}\end{bmatrix}\mathbf C,
\qquad
\boldsymbol\varepsilon_m=[\varepsilon_{11},\varepsilon_{22},2\varepsilon_{12}]^T.
$$

For the **tied** field, introduce $a=1/\sqrt3$, $b=\sqrt{3/5}$, the linear cardinal functions $L_{\pm}(x)$ on $(\!-a,+a)$ and quadratic cardinal functions $Q_{-,0,+}(x)$ on $(\!-b,0,+b)$:

$$
L_-(x)=\frac{a-x}{2a},\quad L_+(x)=\frac{a+x}{2a};\qquad
Q_-(x)=\frac{x(x-b)}{2b^2},\quad Q_0(x)=1-\frac{x^2}{b^2},\quad Q_+(x)=\frac{x(x+b)}{2b^2}.
$$

The tensor-product membrane assumed fields are

$$
\begin{aligned}
\overline{\mathbf B}_{11}(r,s)
&=\sum_{i\in\{-,+\}}\sum_{j\in\{-,0,+\}}
L_i(r)Q_j(s)\widehat{\mathbf B}_{11}(r_i,s_j),\\
\overline{\mathbf B}_{22}(r,s)
&=\sum_{i\in\{-,+\}}\sum_{j\in\{-,0,+\}}
Q_j(r)L_i(s)\widehat{\mathbf B}_{22}(r_j,s_i),\\
\overline{\mathbf B}_{12}(r,s)
&=\sum_{i,j\in\{-,+\}}L_i(r)L_j(s)\widehat{\mathbf B}_{12}(r_i,s_j).
\end{aligned}
$$

The **v3 geometry-dependent switch** applies tied membrane operators only when

$$
\boxed{
\max_i\cos^{-1}(\mathbf n_i\cdot\bar{\mathbf n})>0.01\ \mathrm{rad}
\quad\text{or}\quad
\max_{\mathrm{edges}}\frac{\|\mathbf x_{\mathrm{mid}}-(\mathbf x_A+\mathbf x_B)/2\|}{\|\mathbf x_B-\mathbf x_A\|}>10^{-3}.}
$$

On flat, straight-sided distorted patches, v3 retains the **untied raw membrane operator** by default, protecting membrane patch consistency. On curved/warped elements, the tied operator targets membrane locking. The thresholds are empirical switching parameters, not material constants or mathematical convergence guarantees.

## 5. Transverse shear: v2 locking mechanism and v3 blend

With interpolated director $\mathbf t=\sum_i\widetilde N_i\mathbf n_i$ and field tangents $\mathbf g_\alpha$, the natural shear operator of v2 is, for each node $i$,

$$
\boxed{
\widehat{\mathbf B}_{s,\alpha}^{(i)}=
\begin{bmatrix}
\widetilde N_{i,\alpha}\mathbf t^T &
\widetilde N_i(\mathbf n_i\times\mathbf g_\alpha)^T
\end{bmatrix},\quad \alpha=1,2.}
$$

These two covariant rows are transformed to the physical tangent components through $\mathbf B_s=\mathbf C^T\widehat{\mathbf B}_s$.

### 5.1 Original six-point tensor-product tying

For each transverse shear component the v2 operator samples at **six** stations, corresponding to **12 componentwise constraints**:

$$
\begin{aligned}
\widehat{\mathbf B}_{s,1}^{(6)}(r,s)
&=\sum_{i=\pm}\sum_{j=-,0,+}L_i(r)Q_j(s)
\widehat{\mathbf B}_{s,1}(r_i,s_j),\\
\widehat{\mathbf B}_{s,2}^{(6)}(r,s)
&=\sum_{i=\pm}\sum_{j=-,0,+}Q_j(r)L_i(s)
\widehat{\mathbf B}_{s,2}(r_j,s_i).
\end{aligned}
$$

The six-point shear interpolation is the operative distinction: the reported v2 thin-plate pressure results exhibit **locking despite the Kikuchi field**. Therefore this numerical defect belongs to the **shear tying pattern**, not merely to the polynomial field definition.

### 5.2 Four-point BDG tying

Replace the quadratic transverse interpolation $Q_j$ by endpoint interpolation $H_{\pm}(x)=(1\pm x)/2$:

$$
\begin{aligned}
\widehat{\mathbf B}_{s,1}^{(4)}(r,s)
&=\sum_{i=\pm}\sum_{j=\pm}L_i(r)H_j(s)
\widehat{\mathbf B}_{s,1}(r_i,\eta_j),\\
\widehat{\mathbf B}_{s,2}^{(4)}(r,s)
&=\sum_{i=\pm}\sum_{j=\pm}H_j(r)L_i(s)
\widehat{\mathbf B}_{s,2}(\eta_j,s_i),
\quad\eta_j\in\{-1,1\}.
\end{aligned}
$$

The reduced number of shear sampling constraints avoids the severe v2 thin-shell overconstraint but may produce **coarse-mesh excess flexibility**. This trade-off motivated the blended v3 construction.

### 5.3 Slenderness-dependent v3 operator

For standard geometric area $A_e=\int_{\Omega_e}dA_0$, define

$$
\chi_e=\frac{\sqrt{A_e}}{h},\qquad
\omega_e=\max\left(0,1-\frac{15}{\chi_e}\right),
$$

and the **source-exact strain interpolation blend**

$$
\boxed{\mathbf B_s^{\mathrm{v3}}(r,s)
=\mathbf C(r,s)^T\left[(1-\omega_e)\widehat{\mathbf B}_s^{(6)}(r,s)
+\omega_e\widehat{\mathbf B}_s^{(4)}(r,s)\right].}
$$

Crucially, this is **blending the shear operators before squaring them in the energy**, not a linear interpolation of stiffness matrices:

$$
\mathbf K_s^{\mathrm{v3}}=h\int_{\Omega_e}
(\mathbf B_s^{\mathrm{v3}})^T\mathbf C_s\mathbf B_s^{\mathrm{v3}}\,dA_0
\ne(1-\omega_e)\mathbf K_s^{(6)}+\omega_e\mathbf K_s^{(4)}
\quad\text{in general}.
$$

This distinction is important because the mixed cross-term contributes to the stiffness and to the thin-plate behavior. Direct (`tying=False`), forced six-point, and forced four-point branches remain distinct research settings.

<img width="2730" height="840" alt="fig02_shear_tying" src="https://github.com/user-attachments/assets/dadaf7a0-d701-4be5-8724-f4e2a414bfe1" />

**Figure 2.** Shear sampling locations for the $r$-component and the combined six-point pattern. The $s$-component uses transposed coordinates, not coincident duplication.

<img width="1680" height="987" alt="fig03_blend_weight" src="https://github.com/user-attachments/assets/b862d05e-4672-4d35-a218-99a2baf048ea" />

**Figure 3.** Continuous blending parameter $\omega_e$ as a function of $\chi_e=\sqrt{A_e}/h$; $\omega_e=0$ for $\chi_e\le15$, while $\omega_e\to1$ for increasingly thin shells.

## 6. Bending, drilling and assembled stiffness

Let $\widetilde N_{i,a}$ denote gradients in the orthonormal pointwise basis and $n_{a}=\sum_i\widetilde N_{i,a}\mathbf n_i$. The **bending operator** is the v2/v3 surface-director linearization:

$$
\begin{aligned}
\mathbf B_{b,11}^{(i),u}&=-\widetilde N_{i,1}\mathbf n_{,1}^{T}, &
\mathbf B_{b,11}^{(i),\theta}&=\widetilde N_{i,1}(\mathbf e_1\times\mathbf n_i)^T,\\
\mathbf B_{b,22}^{(i),u}&=-\widetilde N_{i,2}\mathbf n_{,2}^{T}, &
\mathbf B_{b,22}^{(i),\theta}&=\widetilde N_{i,2}(\mathbf e_2\times\mathbf n_i)^T,\\
\mathbf B_{b,12}^{(i),u}&=-\widetilde N_{i,1}\mathbf n_{,2}^{T}-\widetilde N_{i,2}\mathbf n_{,1}^{T}, &
\mathbf B_{b,12}^{(i),\theta}&=\widetilde N_{i,1}(\mathbf e_2\times\mathbf n_i)^T+\widetilde N_{i,2}(\mathbf e_1\times\mathbf n_i)^T.
\end{aligned}
$$

With $\mathbf d_e$ in global Cartesian DOFs, the drilling penalty operator is

$$
\boxed{\mathbf B_d^{(i)}=
\begin{bmatrix}
\tfrac12\left(\widetilde N_{i,1}\mathbf e_2-\widetilde N_{i,2}\mathbf e_1\right)^T &
-\widetilde N_i\mathbf n^T
\end{bmatrix}.}
$$

Thus drilling is not introduced as a spurious membrane-strain term; it enters only through a separate stabilization energy. With $G=E/[2(1+\nu)]$, the isotropic plane-stress and transverse-shear constitutive matrices are

$$
\mathbf C_m=\frac{E}{1-\nu^2}
\begin{bmatrix}1&\nu&0\\\nu&1&0\\0&0&(1-\nu)/2\end{bmatrix},
\qquad \mathbf C_s=k_sG\mathbf I_2.
$$

The v3 elemental stiffness is

$$
\boxed{
\mathbf K_e=\int_{\Omega_e}\left[
 h\mathbf B_m^T\mathbf C_m\mathbf B_m
 +\frac{h^3}{12}\mathbf B_b^T\mathbf C_m\mathbf B_b
 +h\mathbf B_s^T\mathbf C_s\mathbf B_s
 +h\beta_dG\mathbf B_d^T\mathbf B_d
\right]dA_0.}
$$

Integration uses the $3\times3$ Gauss rule in Figure 1. The source default drilling penalty is $\beta_d=10^{-3}$ at the native base level. No EAS parameters or independent stress parameters are statically condensed in the inspected source. A free isolated 48-DOF element should be audited for six rigid-body modes and no spurious zero modes rather than assumed stable solely from stiffness symmetry.

## 7. Mass, pressure, thermal loads and geometric stiffness

### 7.1 Consistent and HRZ-lumped mass

The consistent translational and rotational inertial coefficients are $\rho h$ and $\rho h^3/12$ respectively. Using the modified field for mass interpolation,

$$
\mathbf M_e=\int_{\Omega_e}\mathbf N_u^T(\rho h\mathbf I_3)\mathbf N_u\,dA_0
+\int_{\Omega_e}\mathbf N_{\theta}^T
\frac{\rho h^3}{12}\mathbf P_n\mathbf N_{\theta}\,dA_0.
$$

The projector $\mathbf P_n$ is either $\mathbf I_3$ (all rotational inertias) or $\mathbf I_3-\mathbf n\otimes\mathbf n$ (tangential rotational inertia). The HRZ lumped mass distributes positive diagonal values using $\int\widetilde N_i^2dA_0$ normalized to exactly conserve the integrated component mass. This is a positive-diagonal redistribution, **not** row-sum lumping of sign-changing Q8 mass entries.

### 7.2 Surface pressure

For a scalar signed connectivity-normal pressure $p(r,s)$,

$$
\mathbf f_{p,i}=\int_{-1}^1\int_{-1}^1\widetilde N_i(r,s)
 p(r,s)\left[\mathbf g_1^{(0)}\times\mathbf g_2^{(0)}\right]dr\,ds.
$$

For fixed unit global direction $\hat{\mathbf q}$, replace the normal-area vector by $\hat{\mathbf q}\|\mathbf g_1^{(0)}\times\mathbf g_2^{(0)}\|$; an adjustable quadrature (default $5\times5$) is retained. The force is a **global translational nodal vector**, with no pressure-induced rotational force in this linear model. The alternative standard-field branch is an audit tool, not the consistent modified-field default.

### 7.3 Uniform and linearly varying temperature

For isotropic linear thermal expansion,

$$
\boldsymbol\varepsilon_T=\alpha\Delta T_0[1,1,0]^T,
\qquad
\boldsymbol\kappa_T=\alpha\Delta T_1[1,1,0]^T,
$$

and the consistent nodal thermal action is

$$
\mathbf f_T=\int_{\Omega_e}\left[
 h\mathbf B_m^T\mathbf C_m\boldsymbol\varepsilon_T
+\frac{h^3}{12}\mathbf B_b^T\mathbf C_m\boldsymbol\kappa_T
\right]dA_0.
$$

Here $\Delta T_1$ means a **through-thickness temperature gradient** (temperature per length), not top-minus-bottom temperature without thickness normalization.

### 7.4 Initial-stress geometric stiffness

For in-plane membrane stress
$\boldsymbol\sigma_m=[\sigma_{11},\sigma_{22},\sigma_{12}]^T$ in the pointwise tangent basis, set

$$
\mathbf S_m=h\begin{bmatrix}\sigma_{11}&\sigma_{12}\\\sigma_{12}&\sigma_{22}\end{bmatrix},
\qquad
H_{ij}=\int_{\Omega_e}(\nabla_t\widetilde N_i)^T\mathbf S_m(\nabla_t\widetilde N_j)dA_0.
$$

The initial-stress translational blocks are

$$
\boxed{(\mathbf K_G)_{iu,jv}=H_{ij}\delta_{uv},\quad u,v\in\{1,2,3\},}
$$

or a normal-only projection for plate-bending studies. The source does **not** include the complete rotational, director-curvature or nonlinear mixed formulation tangent. Hence its eigen-buckling validation is restricted to the reported initial-stress assumptions.

## 8. Stress and resultant recovery: required independent quality gates

At a given physical station $(r,s)$, the *same active* operators used for stiffness should be used for elemental recovery:

$$
\boldsymbol\varepsilon_m=\mathbf B_m\mathbf d_e,
\quad \boldsymbol\kappa=\mathbf B_b\mathbf d_e,
\quad\boldsymbol\gamma_s=\mathbf B_s\mathbf d_e.
$$

Then

$$
\mathbf N=h\mathbf C_m\boldsymbol\varepsilon_m,\qquad
\mathbf M=\frac{h^3}{12}\mathbf C_m\boldsymbol\kappa,\qquad
\mathbf Q=h\mathbf C_s\boldsymbol\gamma_s.
$$

The fiber stress at $z\in[-h/2,+h/2]$ is

$$
\boxed{\boldsymbol\sigma_{\parallel}(r,s,z)
=\mathbf C_m\left[\boldsymbol\varepsilon_m(r,s)+z\boldsymbol\kappa(r,s)\right],}
$$

with sign convention following the chosen outward director. For comparison with Nastran/MYSTRAN, transform the physical $2\times2$ tensor with the actual pointwise tangent basis; never reinterpret already-rotated values in an unrelated element-center basis. **Unaveraged** Gauss/center/corner/fiber values and **averaged** nodal values must be reported separately.

A successful displacement patch test does **not** establish stress recovery accuracy. Stress acceptance requires componentwise $\sigma_{11},\sigma_{22},\tau_{12}$ and $M_{11},M_{22},M_{12}$ checks, as well as appropriate transverse-shear resultant checks. The supplied benchmark summary does not furnish complete Nastran elementwise stress parity for HBQ8 v3; therefore no such parity is claimed.

## 9. Numerical verification and convergence (reported data)

**Evidence policy.** All numbers in this section reproduce the user's benchmark report in this conversation. They are not presented as newly executed numerical runs. Mesh counts and response definitions are retained where reported; missing reference values and missing solvers are explicitly identified.

### 9.1 Thin simply supported quarter-plate under uniform pressure

The reported response is the **central deflection divided by the analytical/reference value**, with $n$ denoting the quarter-plate mesh subdivision. The following data exhibit the locking-versus-flexibility trade-off.

| Thickness $h$ | $n$ | SIMO v13 | HBQ8 v2 | **HBQ8 v3** | MITC8 v4 |
|---|---:|---:|---:|---:|---|
| $10^{-2}$ | 1 / 2 / 4 | .820 / .978 / .996 | .409 / .815 / .997 | **.896 / .987 / 1.000** | Same as v3 (reported) |
| $10^{-3}$ | 1 / 2 / 4 | .820 / .977 / .995 | .380 / .130 / .765 | **.910 / .990 / 1.001** | Same as v3 (reported) |
| $10^{-4}$ | 1 / 2 / 4 / 8 | .820 / .977 / .995 / .999 | .379 / **.095 / .053** / .683 | **.911 / .990 / 1.001 / 1.001** | Same as v3 (reported) |


<img width="1890" height="1176" alt="fig04_thinplate_convergence" src="https://github.com/user-attachments/assets/e56e9e37-2c7a-47a4-9316-dddfbacf19b8" />

**Figure 4.** **Connected-point line convergence** of the reported quarter-plate response. Each marker corresponds to a supplied numerical value; line segments are guides only, not interpolated FE solutions. Nonmonotonic deterioration of v2 with plate thinning is a strong locking indicator. The plot does not establish stress convergence, which remains separately required.

### 9.2 Distorted membrane/bending patch 2-001

The user reports for v3: **membrane displacement error $4\times10^{-12}$**, **bending displacement error $4\times10^{-10}$**, with a reversed-winding error of **$7\times10^{-10}$** (numerical error convention as in the source benchmark report). These values were reportedly unchanged from v2 when flat-element membrane tying remained disabled. The result documents **displacement patch consistency**, not a direct elementwise membrane/bending stress or top/bottom fiber stress acceptance gate.

The deliberate trade-off is that curved-element membrane tying improves resistance to membrane locking but may alter the patch behavior if imposed on flat distorted elements. The default switch preserves direct membrane interpolation on flat straight-sided elements. It is not a universal proof of nonlinear patch consistency on arbitrary warped surfaces.

### 9.3 Curved beam 2-003 and pinched cylinder T9

| Case | Mesh | HBQ8 v2 | HBQ8 v3 | Quantity |
|---|---|---:|---:|---|
| 2-003 LC1 | $n_x=2$ | −91.1% | **−4.77%** | Signed response error |
| 2-003 LC1 | $n_x=8$ | −4.55% | **−0.13%** | Signed response error |
| 2-003 LC1 | $n_x=32$ | −0.14% | **−0.12%** | Signed response error |
| T9 | $N=4$ | .809 | **.966** | Displacement/reference |
| T9 | $N=8$ | .875 | **1.02** | Displacement/reference |
| T9 | $N=16$ | 1.00 | **1.02** | Displacement/reference |

<img width="1470" height="965" alt="fig05_2003_convergence" src="https://github.com/user-attachments/assets/b926cb0b-63dd-4c15-96bd-889897aa3444" />


**Figure 5.** Connected-point refinement of 2-003 LC1 for v3. The residual approximately 2% overshoot at T9 $N=16$ is **not erased**: the curved geometry membrane tying improves coarse meshes but can overshoot the reference at finer levels.

### 9.4 Twisted beam and plate bending

| Benchmark | Reported v3 result | Interpretation |
|---|---|---|
| 2-004, $n_x=24$ | +1.78% / +7.02% | Two reported displacement components; not a stress acceptance |
| T8, $N=2$ | .731 (v2 .664) | Improved coarse-mesh displacement ratio |
| T8, $N=4$ | .897 (v2 .864) | Improved, still short of 1.0 |

### 9.5 Buckling 2-016 and 2-017

| Case | Mesh | Reported HBQ8 v3 error | Caution |
|---|---|---:|---|
| 2-016 shear, $h=.01$ | $8\times8$ | +1.11% | Critical-stress error, not raw eigenvalue |
| 2-016 shear, $h=.01$ | $16\times16$ | −0.18% | Critical-stress error |
| 2-016 bending, $h=.01$ | $8\times8$ | +0.73% | Critical-stress error |
| 2-016 bending, $h=.01$ | $16\times16$ | −0.32% | Critical-stress error |
| 2-017 column | Mesh not specified in supplied summary | −1.10% | Exact mesh/reference required for strict parity |

Reported HBQ8 v3 buckling values agree with ANS8 v4 and MITC8 v4 for those settings; this is a **user-reported comparison**, not a claim of cross-solver stress-tangent equivalence or a new Nastran run. The translation-only $\mathbf K_G$ limitations in §7.4 remain in effect.

### 9.6 Operator and force consistency checks

The reported prior development checks include $\mathbf K_{\mathrm{v3}}=\mathbf K_{\mathrm{v2}}$ at old settings (tensor6, membrane tying off) for flat/cylindrical geometry, pressure force differences $\sim5\times10^{-17}$, and normalized mass ratio 1.000000. These are highly useful regression baselines but **were originally run with a replacement base**, and should be repeated with the complete original base and dependencies before making strong implementation-parity claims. The benchmark figures above have the same provenance qualifier.

## 10. Trade-offs, failure mechanisms and validation priorities

| Design choice | Gain | Risk / limitation | Required quality gate |
|---|---|---|---|
| Modified Kikuchi Q8 field | Explicit center-bubble redistribution via $c_i$ | Degenerate denominator or curvature/projection effects still require checks | Nodal interpolation, partition of unity, linear completeness, stress patch |
| `tensor6` shear | Richer interpolation constraints | Severe thin-plate shear locking | $h\to0$ pressure bending convergence |
| `bdg4` shear | Lower constraint count | Excess flexibility on coarse meshes | Low-$n$ pressure and shear force |
| v3 shear blend | Balances locking and coarse flexibility | Empirical threshold 15; mesh/geometry sensitivity | Multi-thickness, multi-mesh convergence, pressure/stress |
| Selective membrane tying | Large gains on curved shell T9 / curved beam 2-003 | Mild T9 overshoot; switching discontinuity | Flat distorted patch + warped/cylindrical regression |
| Modified-field consistent loads | Matching interpolation for pressure/thermal | Curved-surface resultant/first moment not universally proved | Force/torque identities; winding-normal audit |
| Translational $\mathbf K_G$ | Matches reported linear eigen-buckling cases | Missing rotational and curved nonlinear terms | 2-016/2-017 plus curved prebuckling stress |
| Pointwise stress recovery | Explicit local tensor basis | Output averaging and double-rotation hazards | 2-001 stress, LE5, Nastran/MYSTRAN OP2 parity |

The **decision rule** is not to select a formulation solely because it matches a single center displacement. A production candidate must pass (i) patch displacement **and stress**, (ii) shear nonzero regression, (iii) buckling/thermal/mass/pressure consistency, and (iv) convergence on curved or distorted meshes without unacceptable regression.

## 11. Reproduction and publication checklist

1. Reproduce the Kikuchi and MacNeal–Harder coefficient construction above; verify nodal interpolation, partition of unity, derivative sums, and patch consistency on flat and distorted meshes, including reversed winding.
2. Re-run single-element rigid-body rank, rigid rotations and reversed winding with the complete original class chain.
3. Test 2-001 membrane and bending **stress** at Gauss/center/corners/top-bottom fiber, not only displacement.
4. Reproduce plate-pressure $h=10^{-2},10^{-3},10^{-4}$ for $n=1,2,4,8$ with both `tensor6` and `bdg4` forced branches and `auto`.
5. Perform connected-point convergence curves for 2-003/2-004, T8/T9, 2-016/2-017.
6. Compare Nastran and MYSTRAN results on **identical Q8 geometry, BCs, load definitions and output coordinates**; distinguish Nastran element-centered stress, corner stress and nodal averaged output.
7. Audit nonzero transverse shear and thermal curvature sign; audit pressure integrated force and torque on curved geometry.
8. Archive all mesh, tolerance, eigenvalue, stress extraction and solver data. Report failures, including T9 overshoot and the thin-plate v2 nonmonotonic response.

## References

1. F. Kikuchi, M. Okabe, and H. Fujio (1999), “Modification of the 8-node serendipity element,” *Computer Methods in Applied Mechanics and Engineering* **179**(1–2), 91–109. DOI: [10.1016/S0045-7825(99)00031-6](https://doi.org/10.1016/S0045-7825(99)00031-6).
2. K. Darılmaz and N. Kumbasar (2006), “An 8-node assumed stress hybrid element for analysis of shells,” *Computers & Structures* **84**(29–30), 1990–2000. DOI: [10.1016/j.compstruc.2006.08.003](https://doi.org/10.1016/j.compstruc.2006.08.003). **Historical inspiration only; source implementation differs variationally.**
3. R. H. MacNeal and R. L. Harder (1985), “A proposed standard set of problems to test finite element accuracy,” *Finite Elements in Analysis and Design* **1**(1), 3–20. DOI: [10.1016/0168-874X(85)90003-4](https://doi.org/10.1016/0168-874X(85)90003-4).
4. P.-S. Lee and K. J. Bathe (2004), “Development of MITC isotropic triangular shell finite elements,” *Computers & Structures* **82**, 945–962. DOI: [10.1016/j.compstruc.2004.02.004](https://doi.org/10.1016/j.compstruc.2004.02.004). Methodological ANS/MITC context, not an identical Q8 construction.
5. J. C. Simo and M. S. Rifai (1990), “A class of mixed assumed strain methods and the method of incompatible modes,” *International Journal for Numerical Methods in Engineering* **29**, 1595–1638. Variational context; v3 does not employ the Simo–Rifai EAS static condensation.
6. Computers and Structures, Inc., CSI verification suite: [Shells 2-001](https://docs.csiamerica.com/manuals/csibridge/Verification/Analysis/Shells/Problem%202-001.pdf), [2-003](https://docs.csiamerica.com/manuals/csibridge/Verification/Analysis/Shells/Problem%202-003.pdf), [2-004](https://docs.csiamerica.com/manuals/csibridge/Verification/Analysis/Shells/Problem%202-004.pdf), [2-016](https://docs.csiamerica.com/manuals/csibridge/Verification/Analysis/Shells/Problem%202-016.pdf), [2-017](https://docs.csiamerica.com/manuals/csibridge/Verification/Analysis/Shells/Problem%202-017.pdf).
7. [MYSTRAN v18.00.a shell benchmark repository](https://github.com/realbabilu/MYSTRAN/tree/v18.00.a/update/shell_benchmark). Consult exact CSV case identity, solver, and mesh before treating a comparison as native solver evidence.
8. Project source lineage: `HBQ8_ShellElement`, `KikuchiMacNeal_HBQ8_v1`, `HBQ8_Kikuchi_v2`, `HBQ8_Kikuchi_v3`, `shell_q8_common`, and `macneal_kikuchi_shapefun` (supplied source files, October 2026).
ls are unchanged; applied loads and displacements remain unchanged. This does not remove transverse-shear, curved-shell or buckling-normalization limitations. See [reference comparison](../../shell_benchmark/comparison_all.md).
