# κMITC3+ for the MYSTRAN CTRIAR Element: Allman Drilling Membrane, Mixed Transverse Shear, and Linear Composite Curvature

**Ary | 24 October 2026 | Research formulation article |

> **Attribution and scope.** This technical manuscript presents the theoretical basis and project-specific formulation proposed for the `CTRIAR` shell element in MYSTRAN, based on the Nguyen (2024) κMITC3+ concept. The numerical evidence arises from the project's independent 

## Abstract

The project κMITC3+ is a three-node, six-DOF-per-node triangular shell intended for the MYSTRAN `CTRIAR` family. Its formulation combines a modified Allman in-plane interpolation that activates the physical drilling rotation, the associated membrane-rank correction, a six-interior-point transverse-shear interpolation inherited from the Lee–Lee–Bathe MITC3+ construction, and a first-order Hu–Washizu-compatible curvature projection evaluated over three centroid-based subtriangles. The global element has 18 degrees of freedom and two statically condensed internal rotational bubble variables. The present project implementation passes affine membrane, bending and twisting checks and has six rigid-body modes. Numerical improvements in membrane bending and several shell benchmarks are accompanied by a significant coarse curved-shell failure: for the supplied original SAP 2-007 hemisphere DAT, the n=8 response/reference ratio is **0.25567**, compared with **0.99193** for the Lee+Spin baseline. Consequently κMITC3+ is selected as the CTRIAR research replacement for MITC3+HB, **but is not yet validated for native production use**.

## 1. Introduction and positioning in MYSTRAN

Three-node shell elements are computationally attractive because they conform to unstructured triangular meshes and avoid additional midside nodes. Their principal difficulties are overstiff membrane behavior under in-plane bending, transverse-shear locking, insufficient drilling-rotation rank and geometric sensitivity on coarse curved surfaces. A conventional CST membrane enriches neither the in-plane displacement gradient nor the physical role of the sixth nodal rotational degree of freedom. A penalty applied solely to restore rank is mathematically distinct from a drilling rotation that participates in the membrane deformation field.

For the MYSTRAN research program, the κMITC3+ family is proposed to address these mechanisms through three mutually consistent blocks: modified Allman membrane, mixed transverse shear and projected composite bending curvature. The *CTRIAR* designation reflects the membrane drilling capability; Lee + Spin MITC3+ remains a distinct *CTRIA3* formulation. The purpose of the comparison below is to explain the formulation choice and its measurable benefits and limitations, rather than to present κMITC3+ as categorically superior.

| MYSTRAN card | Parameter | Selected formulation | Role |
|---|---|---|---|
| `CTRIA3` | `PARAM,TRIA3TYP,MITC3+` | Lee (2014) + project consistent-spin | Separate CTRIA3 implementation |
| **`CTRIAR`** | **`PARAM,TRIARTYP,KMITC3+`** | **Nguyen-inspired κMITC3+** | **Replaces active CTRIAR MITC3+HB once checked** |

The former **MITC3+HB** was principally an artificial drilling-rank-repair research variant. κMITC3+ introduces a membrane operator in which nodal drilling enters the strain field, plus a separate rank constraint derived from membrane continuum spin. This provides a theoretical justification for the CTRIAR identity. The command `PARAM,TRIARTYP,KMITC3+` is the **chosen MYSTRAN interface**; operational availability must be confirmed against the eventual native executable.

## 2. Parent triangle, local basis and degrees of freedom

Let the reference triangle be

\[
T=\{(r,s): r\ge 0,\ s\ge 0,\ r+s\le 1\},\qquad
N_1=1-r-s,\quad N_2=r,\quad N_3=s.
\]

For physical coordinates \(\mathbf x_i\), define a right-handed facet frame

\[
\mathbf e_x=\frac{\mathbf x_2-\mathbf x_1}{\lVert\mathbf x_2-\mathbf x_1\rVert},\qquad
\mathbf e_z=\frac{(\mathbf x_2-\mathbf x_1)\times(\mathbf x_3-\mathbf x_1)}{\lVert(\mathbf x_2-\mathbf x_1)\times(\mathbf x_3-\mathbf x_1)\rVert},\qquad
\mathbf e_y=\mathbf e_z\times\mathbf e_x.
\]

The external nodal vector is \(\mathbf d_i=[u_i,v_i,w_i,\theta_{xi},\theta_{yi},\theta_{zi}]^T\), with \(18\) element DOFs. Two relative internal rotation-bubble variables \(\boldsymbol\beta\in\mathbb R^2\) are condensed locally. The polynomial bubble is \(N_b=27rs(1-r-s)\); it does not create global nodes or DOFs.

![Three corner nodes and rotational bubble](../figures/mitc3plus_nodal_bubble.svg)

*Figure 1. Three-node connectivity and the internal rotation-bubble convention, adapted to the project's right-handed local frame.*

## 3. Constitutive operators

For isotropic plane stress,

\[
\mathbf C_0=\frac{E}{1-\nu^2}
\begin{bmatrix}1&\nu&0\\\nu&1&0\\0&0&(1-\nu)/2\end{bmatrix},\quad
\mathbf D_m=h\mathbf C_0,\quad
\mathbf D_b=\frac{h^3}{12}\mathbf C_0,\quad
\mathbf D_s=\kappa_s Gh\mathbf I_2,
\]

where \(G=E/[2(1+\nu)]\) and the *current project* uses \(\kappa_s=5/6\). Membrane, bending and shear contributions are assembled independently before internal-mode condensation.

### 3.1 Modified Allman membrane and actual drilling DOF

The three-node constant-strain-triangle (CST) membrane alone has constant derivatives and does not provide a genuine membrane role for \(\theta_z\). The Nguyen-inspired project branch adopts a modified Allman interpolation,

\[
\boxed{\mathbf B_m^{A}(r,s)=\mathbf B_{\mathrm{CST}}-\alpha_A\mathbf H_M(r,s)\mathbf G_M,\qquad \alpha_A=\tfrac23.}
\]

Here \(\mathbf G_M\) maps cyclic edge differences of drilling rotations and \(\mathbf H_M\) maps the resulting Allman edge modes into in-plane strains. Their exact components and permutation signs must be verified directly against `MITC3Plus_Allman_Nguyen2024_Experimental.py` and the referenced paper before Fortran porting. In particular, \(\alpha_A\) is a membrane coefficient, **not** Lee's `spin_alpha`.

The independent mean-spin constraint in the project is a rank-controlling operator, not an added arbitrary numerical diagonal for \(\theta_z\):

\[
\mathbf b_{\mathrm{Add},i}=
\begin{bmatrix}N_{i,y}/2&-N_{i,x}/2&0&0&0&1/3\end{bmatrix},\quad
\boxed{\mathbf K_{\mathrm{Add}}=GhA\,\mathbf b_{\mathrm{Add}}^T\mathbf b_{\mathrm{Add}}.}
\]

The block row above denotes the contribution from each node in the 18-column assembled operator. Its normalization, rigid-rotation annihilation and physical units are mandatory patch-test checks. Do **not** add a separate Hughes–Brezzi or Lee consistent-spin stabilization on top without an explicit independent derivation: the present Nguyen project class calls `keep_spin=False`.

### 3.2 Six interior tying points for transverse shear

The assumed shear operator \(\widetilde{\mathbf B}_s(r,s)\) is formed from six interior locations A–F in the Lee MITC3+ family, retaining contributions from nodal rotations, transverse displacement and the rotational bubble. The numerical interior offset in the reconstruction is \(d=10^{-4}\). The parent integration uses the symmetric three-point rule

\[
(r,s,w)=\left(\tfrac16,\tfrac16,\tfrac13\right),\ \left(\tfrac23,\tfrac16,\tfrac13\right),\ \left(\tfrac16,\tfrac23,\tfrac13\right),
\qquad \mathbf K_s=\sum_g A w_g\,\widetilde{\mathbf B}_s^T\mathbf D_s\widetilde{\mathbf B}_s.
\]

![Lee interior shear tying points](../figures/mitc3plus_lee2014_tying_points.svg)

*Figure 2. Lee-style interior shear samples; these are tying locations, not the integration rule.*

![Three symmetric Gaussian integration points](../figures/mitc3plus_gauss_3point.svg)

*Figure 3. Symmetric three-point triangular integration used by the current parent operator.*

The existing **constrained-displacement shear patch** ratio is \(1.13408\), not unity. This limitation must be retained, and a separate equilibrated shear-traction/moment-gradient patch is still needed before evaluating whether it is a true consistency failure.

### 3.3 Hu–Washizu first-order composite-curvature projection

The compatible bending curvature \(\boldsymbol\kappa=\mathbf B_b(r,s)\,[\mathbf d;\boldsymbol\beta]\) is projected onto the first-order polynomial basis \(\mathbf f=[1,r,s]^T\). Define

\[
\mathbf F=\int_A\mathbf f\mathbf f^T\,dA,\qquad
\mathbf Q=\int_A\mathbf f\otimes\mathbf B_b\,dA,\qquad
\boxed{\widetilde{\mathbf B}_b(r,s)=\mathbf f(r,s)^T\mathbf F^{-1}\mathbf Q.}
\]

In code, \(\mathbf F\) and \(\mathbf Q\) are assembled using three symmetric points, and a linear solve is used rather than explicitly forming \(\mathbf F^{-1}\). Divide the parent triangle around \(G=(1/3,1/3)\) into \(T_1=(G,1,2)\), \(T_2=(G,2,3)\), \(T_3=(G,3,1)\). If \(G_a\) denotes the centroid of the subtriangle \(T_a\), then

\[
\boxed{\mathbf K_b^{P1}=\sum_{a=1}^{3}\frac{A}{3}
\widetilde{\mathbf B}_b(G_a)^T\mathbf D_b\widetilde{\mathbf B}_b(G_a).}
\]

The optional seven-point quadrature for constructing \(\mathbf F,\mathbf Q\) is a **research sensitivity switch**, not the selected default. It changed 2-002 torsion from 1.12967 to about 1.04990 without correcting the coarse hemisphere n8 response (~0.2533); there is insufficient justification to tune the default based on one case.

## 4. Internal condensation, temperature and stress recovery

With \(\mathbf q=[\mathbf d;\boldsymbol\beta]\), form the symmetric \(20\times20\) local matrix

\[
\mathbf K_{20}=\mathbf K_m^{A}+\mathbf K_{\mathrm{Add}}+\mathbf K_s^{Lee}+\mathbf K_b^{P1}
=\begin{bmatrix}\mathbf K_{cc}&\mathbf K_{cb}\\\mathbf K_{bc}&\mathbf K_{bb}\end{bmatrix},
\]

\[
\boxed{\mathbf K_{18}=\mathbf K_{cc}-\mathbf K_{cb}\mathbf K_{bb}^{-1}\mathbf K_{bc}},\qquad
\boldsymbol\beta=-\mathbf K_{bb}^{-1}\mathbf K_{bc}\mathbf d.
\]

The thermal membrane vector must use the *same* Allman strain matrix:

\[
\mathbf f_T=\int_A\mathbf B_m^{A\,T}\mathbf D_m\,\boldsymbol\varepsilon_T\,dA,
\qquad \boldsymbol\varepsilon_T=\alpha_T\Delta T[1,1,0]^T.
\]

The current code reports a single-element free-expansion force residual of \(4.38\times10^{-17}\), which is **not** a validated SAP 2-013 thermal-stress comparison. Local recovery reports membrane force \(\mathbf N=\mathbf D_m(\mathbf B_m^A\mathbf d-\boldsymbol\varepsilon_T)\), bending moments \(\mathbf M=\mathbf D_b\widetilde{\mathbf B}_b\mathbf q\), and shear resultants \(\mathbf V=\mathbf D_s\widetilde{\mathbf B}_s\mathbf q\). Any nodal averaging across shell facets requires a deliberate output coordinate frame; do not average incompatible local tensors directly.

## 5. The validation against Lee + Spin

All values below are **calculated/reference ratios**, with target unity. They are supplied research runs, not native MYSTRAN parity tests.

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
| 8 | 128 | **0.25567** | **Original supplied DAT** |
| 16 | 512 | 0.81553 | Reconstructed refinement |
| 32 | 2,048 | 0.98735 | Reconstructed refinement |
| 64 | 8,192 | 0.99820 | Reconstructed refinement |

Although reconstructed fine meshes approach unity, their geometry, boundary discretization and loading have not been proven identical to original benchmark DATs. Thus **do not report the n64 result as a successful original-DAT convergence series**. The separate earlier Allman-only plot below illustrates the same qualitative coarse-mesh issue but **must not be relabeled Nguyen P1**.

![Historical Allman-only hemisphere refinement](../figures/mitc3plus_allman_2007_convergence.png)

*Figure 4. Historical Allman-only hemisphere study, included for context; numerical Nguyen P1 refinement is given accurately in the table above.*

### 5.2 Further checks and open coverage

- Adapted 2-008 triangular modal plate: five fine-mesh period ratios approximately 1.00025, 1.00937, 1.00339, 1.00634, and 1.00888. **Not exact CSI CTRIA3 native parity.**
- 2-010 pressurized cylinder: original SOL101 DAT/OP2 supplied, but independent Nguyen native validation is pending.
- 2-017 plate buckling: original SOL105 DAT/OP2 supplied, but Nguyen eigen/geometric-stiffness parity has not been demonstrated.
- 2-013 thermal stresses and additional 2-011/2-012/2-014/2-015/2-016 models remain unvalidated as complete official benchmarks.

## 6. Discussion: drilling, membrane behavior and curved-shell limitations

### 6.1 Why the Allman drilling mechanism is not Hughes–Brezzi penalty stabilization

The added nodal rotation \(\theta_z\) does not enter an ordinary CST membrane. In a conventional rank-fix design, its energetic participation may be supplied only through a small penalty or an auxiliary spin stiffness. In κMITC3+, in contrast, the off-diagonal edge-difference operator \(\mathbf G_M\) makes \(\theta_z\) affect \(\boldsymbol\varepsilon_m\) directly. The mean continuum-spin constraint \(\mathbf K_{\mathrm{Add}}\) complements this membrane interpolation to avoid the otherwise deficient drilling rank. It should therefore not be interpreted as adding an unrelated duplicate penalty. The present research reconstruction contains **no separate Lee `spin_alpha` contribution** in its Nguyen stiffness.

### 6.2 Why an enriched membrane can improve in-plane bending but harm a coarse hemisphere

The approximately 0.50-to-0.80 improvement in 2-002 LC2/LC5 is consistent with additional membrane kinematics reducing CST overstiffness. The deterioration in the coarse 2-007 hemisphere is equally important: curved response depends on the variation of the local shell basis across facets, the geometry discretization, drilling–membrane coupling and possible membrane/bending contamination. The available tests do not uniquely isolate which mechanism dominates. Since the n8 original DAT and n16/n32/n64 reconstructed meshes are not an identical-origin mesh family, the observed refinement cannot be used to declare that a native curved-shell limitation is resolved.

### 6.3 Interpretation of the remaining transverse-shear patch discrepancy

A computed/reference value of 1.13408 for the constrained-displacement shear test is a measurable excess stiffness under that test construction, but it is not sufficient to prove failure of a properly equilibrated constant-shear traction patch. A more discriminating verification would enforce compatible displacement and rotation fields, evaluate internal force consistency, and introduce a linear moment gradient where the governing equilibrium equations demand it. This is a scientific limitation of the present evidence and is recorded without tuning the shear operator.

### 6.4 Implementation-facing research considerations

The intended MYSTRAN deployment is a three-node `CTRIAR` with right-handed six-DOF external kinematics and a condensed two-variable bubble. Native implementation must respect local-to-global transformations, element winding, thermal membrane interpolation, pressure loads, mass, stress/resultant recovery, and geometric stiffness under buckling. In particular, production claims would require comparison against independent reference results for original 2-010 pressure and 2-017 buckling inputs. 

## 7. Conclusions

1. κMITC3+ offers a mechanistically justified drilling-active CTRIAR membrane formulation by combining modified Allman interpolation and a continuum-spin rank constraint, rather than retaining an arbitrary rank-repair penalty alone.
2. First-order projected composite curvature and interior transverse-shear tying form a mathematically explicit energy formulation with 18 external DOFs and two condensed internal bubble coordinates.
3. The result evidence demonstrates gains in 2-002 membrane-dominated responses, 2-005 plate bending and 2-006 roof deformation, but also a significant 2-002 torsion shift and severe **original n8 2-007** coarse-curvature under-response.


## 8. References and research provenance

1. Lee, Y., Lee, P.-S., and Bathe, K.-J. (2014). *The MITC3+ shell element and its performance.* **Computers & Structures**, 138, 12–23. DOI: [10.1016/j.compstruc.2014.02.005](https://doi.org/10.1016/j.compstruc.2014.02.005).
2. Nguyen, S. H. (2024). *Linear composite curvature MITC3+ flat shell elements.* **Computers & Structures**, 295, 107272. DOI: [10.1016/j.compstruc.2024.107272](https://doi.org/10.1016/j.compstruc.2024.107272).
3. Project research source: `MITC3Plus_Nguyen2024_CompositeCurvature.py`; precursor `MITC3Plus_Allman_Nguyen2024_Experimental.py`; comparison `MITC3Plus_Lee2014_ConsistentSpin.py`. The detailed numerical audit is documented separately in `MITC3Plus_Nguyen2024_CompositeCurvature_Validation.md`.

**Conclusion:** κMITC3+ is the chosen future CTRIAR research formulation. Its Allman drilling participation and P1 curvature projection are substantive upgrades over the former HB rank fix, but the documented coarse-mesh curvature deficiency and absent native parity prevent a production-ready claim as of 24 October 2026.
