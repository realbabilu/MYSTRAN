# CHEXA8 formulation versions: EAS9 and Tied-EAS9

## 1. Versions

| # | Version | Strain field | Difference from the previous version |
|---|---|---|---|
| 1 | Standard | $\varepsilon = B\mathbf u$ | Full 2×2×2 Gauss, no enhancement |
| 2 | EAS9, scalar scaling | $\varepsilon = B\mathbf u + \frac{\det J_0}{\det J}\hat M\boldsymbol\alpha$ | Nine enhanced modes, scaled by $\det J_0/\det J$ only, defined in global Cartesian axes |
| 3 | EAS9 (covariant) | $\varepsilon = B\mathbf u + \frac{\det J_0}{\det J}T_0\hat M\boldsymbol\alpha$ | Modes mapped to physical axes with the covariant transformation $T_0$ |
| 4 | EAS12 / EAS18 / EAS21 | as version 3 | $\hat M$ extended by 3 / 9 / 12 modes |
| 5 | Tied-EAS9 | $\varepsilon = \tilde B\mathbf u + M\boldsymbol\alpha$ | $B$ replaced by $\tilde B$: transverse shear $xz$, $yz$ sampled at $\xi=0$, $\eta=0$ |
| 6 | Tied-EAS9 + correction | $\varepsilon = \tilde B'\mathbf u + M\boldsymbol\alpha$ | Weighted mean of $\tilde B - B$ removed (affine and patch consistency) |
| 7 | Tied-EAS9 + correction + auto axis | as version 6 | $\zeta$ relabelled to the shortest element direction |

All versions: 8 nodes, 24 translational DOF, 2×2×2 Gauss, enhanced parameters condensed at element level.

## 2. Common definitions

Natural coordinates $\boldsymbol\xi=(\xi,\eta,\zeta)\in[-1,1]^3$, trilinear shape functions

$$N_k=\tfrac18(1+\xi_k\xi)(1+\eta_k\eta)(1+\zeta_k\zeta),\quad k=1,\dots,8.$$

Jacobian $J_{ai}=\partial x_i/\partial\xi_a=\sum_k \partial N_k/\partial\xi_a\,x_{k,i}$, physical gradient $\nabla N=J^{-1}\,\partial N/\partial\boldsymbol\xi$. $B$ is the 6×24 strain-displacement matrix in Voigt order $[xx,yy,zz,xy,yz,zx]$ (engineering shear). $D$ is isotropic elasticity with $\lambda=E\nu/((1+\nu)(1-2\nu))$ and $\mu=E/(2(1+\nu))$. Gauss points $\pm1/\sqrt3$, weights $w_g=\det J_g$, $J_0=J(0,0,0)$.

## 3. Formulations

### 3.1 Standard

$$K=\sum_g w_g\,B_g^TDB_g.$$

### 3.2 EAS9, scalar scaling

Nine zero-mean enhanced modes in natural strain order $[\xi\xi,\eta\eta,\zeta\zeta,\xi\eta,\eta\zeta,\zeta\xi]$:

$$\hat M=\begin{bmatrix}\xi&0&0&0&0&0&0&0&0\\0&\eta&0&0&0&0&0&0&0\\0&0&\zeta&0&0&0&0&0&0\\0&0&0&\xi&\eta&0&0&0&0\\0&0&0&0&0&\eta&\zeta&0&0\\0&0&0&0&0&0&0&\zeta&\xi\end{bmatrix},\qquad M(\boldsymbol\xi)=\frac{\det J_0}{\det J(\boldsymbol\xi)}\hat M(\boldsymbol\xi).$$

The matrices

$$K_{uu}=\sum_g w_gB_g^TDB_g,\quad K_{u\alpha}=\sum_g w_gB_g^TDM_g,\quad K_{\alpha\alpha}=\sum_g w_gM_g^TDM_g$$

are condensed:

$$K=K_{uu}-K_{u\alpha}K_{\alpha\alpha}^{-1}K_{u\alpha}^T.$$

Because $\hat M$ is used directly as a Cartesian strain field, $K$ depends on the global axis orientation.

### 3.3 EAS9, covariant

Each column of $\hat M$ is a covariant strain tensor $\hat E$ (shear components halved). It is mapped to physical axes with the centre Jacobian:

$$E=J_0^{-1}\hat E\,J_0^{-T},\qquad M(\boldsymbol\xi)=\frac{\det J_0}{\det J(\boldsymbol\xi)}\,T_0\hat M(\boldsymbol\xi),$$

where $T_0$ is the 6×6 matrix of this map in Voigt form $[E_{11},E_{22},E_{33},2E_{12},2E_{23},2E_{31}]$. Since $w_gM_g=\det J_0\,T_0\hat M_g$ and the monomials $\xi,\eta,\zeta$ sum to zero over the eight Gauss points, $\sum_g w_gM_g=0$. Enhanced strains are therefore orthogonal to constant stress, $K_{\alpha u}\mathbf u_{\rm aff}=0$ and $\boldsymbol\alpha=0$ for affine states (patch consistency). $K$ is invariant to rotation of the element.

Recovery: $\boldsymbol\alpha=K_{\alpha\alpha}^{-1}(F_\alpha-K_{\alpha u}\mathbf u)$, $\boldsymbol\varepsilon=B\mathbf u+M\boldsymbol\alpha$, $\boldsymbol\sigma=D(\boldsymbol\varepsilon-\boldsymbol\varepsilon_{th})$. Thermal load: $F_{th}=F_u-K_{u\alpha}K_{\alpha\alpha}^{-1}F_\alpha$. Consistent mass uses the mechanical interpolation only.

### 3.4 Extended mode sets

$$\text{EAS12}=\text{EAS9}+\{\eta\zeta\text{ in }\varepsilon_{\xi\xi},\ \zeta\xi\text{ in }\varepsilon_{\eta\eta},\ \xi\eta\text{ in }\varepsilon_{\zeta\zeta}\}$$

EAS18: $\xi,\eta,\zeta$ in each of the six strains. EAS21: EAS18 plus the three quadratic modes of EAS12. Mapping and condensation as in 3.3.

### 3.5 Tied-EAS9

Local frame from the centre Jacobian:

$$\mathbf n=\frac{J_0[2,:]}{\|J_0[2,:]\|},\quad \mathbf t=J_0[0,:]-(J_0[0,:]\cdot\mathbf n)\mathbf n,\ \ \hat{\mathbf t}=\frac{\mathbf t}{\|\mathbf t\|},\quad \mathbf b=\mathbf n\times\hat{\mathbf t},\quad R=[\hat{\mathbf t}\ \ \mathbf b\ \ \mathbf n].$$

$P_{xz}$ and $P_{yz}$ (6×6) take the $xz$ or $yz$ shear component of $R^T\varepsilon R$ and rotate it back. The tied operator is

$$\tilde B(\xi,\eta,\zeta)=B(\xi,\eta,\zeta)+P_{xz}\big[B(0,\eta,\zeta)-B(\xi,\eta,\zeta)\big]+P_{yz}\big[B(\xi,0,\zeta)-B(\xi,\eta,\zeta)\big].$$

$K$ is built as in 3.3 with $B$ replaced by $\tilde B$; $M$ is unchanged. The operator assumes $\zeta$ is the thickness direction.

### 3.6 Correction

Let $\Delta B=\tilde B-B$. Affine and patch consistency require

- (C1) $\tilde B\mathbf u_{\rm aff}=\varepsilon_0$: constant strain is reproduced;
- (C2) $\sum_g w_g\tilde B_g=\sum_g w_gB_g$: equilibrium for constant stress.

Tied-EAS9 satisfies C1 pointwise (single frame $R$) but not C2 on distorted elements. The correction

$$\tilde B'=\tilde B-\frac{\sum_g w_g\Delta B_g}{\sum_g w_g}$$

enforces C2 exactly and keeps C1, since $\Delta B_g\mathbf u_{\rm aff}=0$ at every point. On parallelepipeds the mean of $\Delta B$ is zero and the correction vanishes.

### 3.7 Auto axis

The mean edge length along $\xi$, $\eta$, $\zeta$ is computed; the shortest direction becomes $\zeta$ by a cyclic relabelling of the natural axes (handedness preserved). The EAS9 mode set is invariant under cyclic permutation. $K$ is computed on the permuted node order and mapped back.

## 4. Solution steps

**Element**

1. Read node coordinates $X$ (8×3), $E$, $\nu$; form $D$.
2. Compute $J_0$, $\det J_0>0$. Versions 3–7: form $T_0$. Version 7: find the shortest direction and permute nodes.
3. Versions 5–7: form $R$, $P_{xz}$, $P_{yz}$.
4. At each of the 8 Gauss points: compute $J$, $\det J>0$, $\nabla N$, $B$. Versions 5–7: also $B(0,\eta,\zeta)$, $B(\xi,0,\zeta)$, then $\tilde B$ and $\Delta B$. Versions 2–7: $M$.
5. Version 6–7: subtract the weighted mean of $\Delta B$ from $\tilde B$.
6. Accumulate $K_{uu}$, $K_{u\alpha}$, $K_{\alpha\alpha}$ with weight $\det J$ (using $\tilde B'$ in versions 5–7).
7. Condense: $K=K_{uu}-K_{u\alpha}K_{\alpha\alpha}^{-1}K_{u\alpha}^T$, symmetrise; version 7: map back to the original node order.

**Global**

8. Assemble element matrices into a sparse global $K$ (3 DOF per node).
9. Partition free and constrained DOF; build the load vector.
10. Solve $K_{ff}\mathbf u_f=F_f$ with sparse LU.
11. Average the displacement over the probe nodes and compare with the reference.

## 5. Benchmarks and results

### 5.1 Benchmarks

| Benchmark | SAP2000 problem | Model | Reference | SAP2000 published | Meshes |
|---|---|---|---|---|---|
| Pinched dome | 5-007 | $E=68{,}250{,}000$ kip/ft², $\nu=0.3$, $R=10$ ft, $t=0.04$ ft, 18° hole, quarter model, one element through thickness | $U_x(A)=0.094$ ft | 0.0701 ft (HEXA8 16×16×1), 0.093 ft (HEXA8 48×48×1) | N×N×1, N = 4, 8, 16, 24, 48 |
| Twisted beam | 5-004 | $E=29\times10^6$, $\nu=0.22$, $L=12$, 90° twist, unit tip force | IN 0.005429, OUT 0.001749 | IN 0.005322, OUT 0.001719 | 24×4×2 |
| Pinched cylinder | – | $R=300$, $L/2=300$, $t=3$, $E=3\times10^6$, $\nu=0.3$, $P=1$, octant, end diaphragm | $1.8248\times10^{-5}$ | – | N×N×1, N = 8, 16, 24 |
| Affine patch | – | 3×3×3 bricks, interior nodes perturbed 15% and 30% | exact affine field | – | 3×3×3 |

### 5.2 Pinched dome 5-007

![Pinched dome convergence](figures_en/fig1_dome.png)

Ratio $U_x(A)/0.094$ (error against reference).

| Mesh | Standard | EAS9, scalar scaling | EAS9 (covariant) | Tied-EAS9 | Tied-EAS9 + correction | Nastran CHEXA8 | MYSTRAN EAS9 v2 | SAP2000 |
|---|---|---|---|---|---|---|---|---|
| 4×4×1 | 0.0007 (99.93%) | 0.0024 (99.76%) | 0.0105 (98.95%) | 0.0415 (95.85%) | 0.0402 (95.98%) | 0.0408 (95.92%) | 0.0105 (98.95%) | – |
| 8×8×1 | 0.0027 (99.73%) | 0.0111 (98.89%) | 0.1625 (83.75%) | 0.7435 (25.65%) | 0.7361 (26.39%) | 0.7422 (25.78%) | 0.1625 (83.75%) | – |
| 16×16×1 | 0.0103 (98.97%) | 0.0416 (95.84%) | 0.7460 (25.40%) | 0.9836 (1.64%) | 0.9833 (1.67%) | 0.9851 (1.49%) | 0.7460 (25.40%) | 0.7457 (25.43%) |
| 24×24×1 | 0.0225 (97.75%) | 0.0867 (91.33%) | 0.9267 (7.33%) | 0.9924 (0.76%) | 0.9924 (0.76%) | 0.9932 (0.68%) | 0.9267 (7.33%) | – |
| 48×48×1 | – | 0.2785 (72.15%) | 0.9892 (1.08%) | 0.9947 (0.53%) | 0.9947 (0.53%) | – | – | 0.9894 (1.06%) |

### 5.3 Twisted beam 5-004

![Twisted beam](figures_en/fig2_beam.png)

| Source | IN (ft) | error IN | OUT (ft) | error OUT |
|---|---|---|---|---|
| SAP2000 | 0.0053220 | 1.97% | 0.0017190 | 1.72% |
| Nastran CHEXA8 | 0.0054018 | 0.50% | 0.0017470 | 0.11% |
| MYSTRAN EAS9 v2 | 0.0053959 | 0.61% | 0.0017462 | 0.16% |
| EAS9, scalar scaling | 0.0053847 | 0.82% | 0.0017342 | 0.85% |
| EAS9 (covariant) | 0.0054057 | 0.43% | 0.0017477 | 0.08% |
| Tied-EAS9 | 0.0054058 | 0.43% | 0.0017477 | 0.07% |
| Tied-EAS9 + correction | 0.0054060 | 0.42% | 0.0017476 | 0.08% |
| Tied-EAS9 + correction + auto axis | 0.0054060 | 0.42% | 0.0017476 | 0.08% |
| Reference | 0.005429 | – | 0.001749 | – |

### 5.4 Pinched cylinder

![Pinched cylinder convergence](figures_en/fig3_cylinder.png)

Ratio to reference.

| Version | 8×8×1 | 16×16×1 | 24×24×1 |
|---|---|---|---|
| Standard | 0.0693 | 0.1485 | 0.2341 |
| EAS9, scalar scaling | 0.1442 | 0.4274 | 0.6312 |
| EAS9 (covariant) | 0.4055 | 0.8360 | 0.9388 |
| Tied-EAS9 | 0.4939 | 0.9066 | 0.9685 |
| Tied-EAS9 + correction | 0.4939 | 0.9066 | 0.9685 |

### 5.5 Affine patch test and rotation invariance

![Patch test](figures_en/fig4_patch.png)

Relative error of free-node displacement. Thin shell: 3×3×0.1 box.

| Version | Cube, 15% | Cube, 30% | Thin shell, 15% | Thin shell, 30% | Rotation error of $K$ |
|---|---|---|---|---|---|
| EAS9, scalar scaling | 3.0e-16 | 4.6e-16 | 6.1e-16 | 7.0e-16 | 8.1e-02 |
| EAS9 (covariant) | 2.8e-16 | 3.6e-16 | 5.6e-16 | 5.5e-16 | 4.2e-16 |
| Tied-EAS9 | 4.8e-05 | 1.9e-04 | 2.6e-04 | 2.8e-04 | 4.2e-16 |
| Tied-EAS9 + correction | 2.9e-16 | 5.2e-16 | 5.4e-16 | 4.1e-16 | 6.3e-16 |
| Tied-EAS9 + correction + auto axis | 2.6e-16 | 3.6e-16 | 5.4e-16 | 4.1e-16 | 6.3e-16 |

### 5.6 Extended EAS mode sets

Ratio to reference.

| Version | Dome 8 | Dome 16 | Dome 24 | Cylinder 8 | Cylinder 16 | Cylinder 24 |
|---|---|---|---|---|---|---|
| EAS9 (covariant) | 0.1625 | 0.7460 | 0.9267 | 0.4055 | 0.8360 | 0.9388 |
| EAS12 | 0.1627 | 0.7467 | 0.9272 | 0.4065 | 0.8386 | 0.9419 |
| EAS18 | 0.2171 | 3.2157 | 15.7140 | 1.0741 | 5.7539 | 15.1887 |
| EAS21 | 0.2214 | 3.2377 | 15.8823 | 1.1582 | 5.8948 | 15.9126 |

### 5.7 Dependence on thickness axis

Dome 24×24×1, ratio to reference, with the thickness direction mapped to $\zeta$, $\eta$ or $\xi$ by relabelling the element nodes.

| Version | thickness = ζ | thickness = η | thickness = ξ |
|---|---|---|---|
| EAS9 (covariant) | 0.9267 | 0.9267 | 0.9267 |
| Tied-EAS9 + correction | 0.9924 | 0.9267 | 0.9269 |
| Tied-EAS9 + correction + auto axis | 0.9924 | 0.9924 | 0.9924 |

## 6. References

### 6.1 Journals

1. J. C. Simo and M. S. Rifai (1990), "A class of mixed assumed strain methods and the method of incompatible modes," *Int. J. Numer. Methods Eng.* 29, 1595–1638. DOI 10.1002/nme.1620290802.
2. U. Andelfinger and E. Ramm (1993), "EAS-elements for two-dimensional, three-dimensional, plate and shell structures and their equivalence to HR-elements," *Int. J. Numer. Methods Eng.* 36, 1311–1337.
3. R. L. Taylor, P. J. Beresford and E. L. Wilson (1976), "A non-conforming element for stress analysis," *Int. J. Numer. Methods Eng.* 10, 1211–1219.
4. J. C. Simo, F. Armero and R. L. Taylor (1993), "Improved versions of assumed enhanced strain tri-linear elements for 3D finite deformation problems," *Comput. Methods Appl. Mech. Eng.* DOI 10.1016/0045-7825(93)90215-J.
5. K. J. Bathe and E. N. Dvorkin (1985), "A four-node plate bending element based on Mindlin/Reissner plate theory and a mixed interpolation," *Int. J. Numer. Methods Eng.* 21, 367–383.
6. R. Hauptmann and K. Schweizerhof (1998), "A systematic development of 'solid-shell' element formulations for linear and non-linear analyses employing only displacement degrees of freedom," *Int. J. Numer. Methods Eng.* 42, 49–69.
7. S. Klinkel, F. Gruttmann and W. Wagner (1999), "A continuum based three-dimensional shell element for laminated structures," *Comput. Struct.* 71, 43–62.
8. L. Vu-Quoc and X. G. Tan (2003), "Optimal solid shells for non-linear analyses of multilayer composites. I. Statics," *Comput. Methods Appl. Mech. Eng.* 192, 975–1016.
9. S. Schwarze et al. (2009), "A reduced integration solid-shell finite element based on the EAS and the ANS concept — geometrically linear problems," *Int. J. Numer. Methods Eng.* DOI 10.1002/nme.2653.
10. T. Belytschko and L. P. Bindeman (1993), "Assumed strain stabilization of the eight node hexahedral element," *Comput. Methods Appl. Mech. Eng.* DOI 10.1016/0045-7825(93)90124-G.
11. M. A. Puso (2000), "A highly efficient enhanced assumed strain physically stabilized hexahedral element," *Int. J. Numer. Methods Eng.* 49(8), 1029–1064. DOI 10.1002/1097-0207(20001120)49:8<1029::AID-NME990>3.0.CO;2-3.
12. T. J. R. Hughes (1980), "Generalization of selective integration procedures to anisotropic and nonlinear media," *Int. J. Numer. Methods Eng.* 15, 1413–1418.
13. R. H. MacNeal and R. L. Harder (1985), "A proposed standard set of problems to test finite element accuracy," *Finite Elem. Anal. Des.* 1, 3–20.

### 6.2 SAP2000 tests

Computers and Structures, Inc., *SAP2000 Verification: Solids*, Problems 5-001, 5-004 and 5-007: https://docs.csiamerica.com/manuals/sap2000/Verification/Analysis/Solids/
