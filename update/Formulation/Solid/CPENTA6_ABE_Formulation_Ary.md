# CPENTA6 ABE: Edge-Tied Assumed Transverse Shear with Enhanced Strain Families for the Six-Node Continuum Wedge

**Ary** · Personal · 18 October 2026  \|  Technical note, MYSTRAN NEWSOLID (`PENTATYP,ABE`). Numerical results are from the Python reference implementation; the Fortran port is documented in `NEWSOLID_CONTINUUM_FORMULATIONS.md`.

## Abstract

ABE is the six-node prism (Nastran CPENTA, 18 translational DOF) obtained by combining (i) a MITC3-type assumed transverse-shear field tied at the three edge midpoints of the triangle, evaluated in a covariant frame frozen at the element centroid, with (ii) nine locally condensed enhanced-strain columns grouped in three families A, B and E. It was developed to remove the parasitic transverse shear that makes ordinary and EAS9 wedges lock in thin bending: on the SAP2000 5-010 column the first buckling factor falls from 20.98 (standard) and 7.68 (EAS9 v2) to 0.8912, against 0.8910 for Nastran CPENTA6 and 0.8883 from Euler theory. ABE is not a published element; its ingredients are taken from the assumed-natural-strain, enhanced-assumed-strain and solid-shell literature cited in Section 3. A single parameter, `patch_consistency`, selects between the best thin-bending behaviour (0) and an exact constant-stress patch test (1); the trade-off is structural and is quantified.

## 1. Introduction

A linear wedge combines a constant-strain triangle in its plane with a linear variation along the prism axis. When one element spans a thin direction (a plate-like layer), pure bending produces a transverse shear that the interpolation cannot cancel. On a triangle, the slope of the normal displacement is constant along an edge while the rotation (the thickness difference of the in-plane displacement) is linear, so the two agree only at the *edge midpoint*, not at the centroid. The discrepancy, of order κh/6, appears as shear strain and dominates the strain energy. Evaluating the accurate (CPENTA15) buckling mode with each wedge operator showed transverse shear energy of 25.8 (standard) and 15.8 (EAS9 v2) against a physical bending energy of 0.89.

EAS9 for the wedge defines modes such as <span class="m">r̄</span> on the (r,ζ) shear and <span class="m">s̄</span> on the (s,ζ) shear. That set is not closed under renumbering of the triangle nodes: the same wedge numbered (2,3,1) instead of (1,2,3) changes its stiffness by about 3%. ABE replaces it by families that are closed, and removes the parasitic shear with an assumed field that is intrinsic to the triangle.

## 2. Element formulation

### 2.1 Nodal numbering and natural coordinates

<div>
<img width="518" height="456" alt="image" src="https://github.com/user-attachments/assets/7fe68c51-d147-430b-94d1-6b5e8bce9061" />

[Figure 1. CPENTA6 grid order (Nastran): G1-G2-G3 on the bottom triangle (ζ = −1), G4-G5-G6 on the top triangle (ζ = +1), G4 above G1. Natural axes r (G1→G2), s (G1→G3), ζ (bottom to top) are right-handed; the volume is positive when G1→G2→G3 is counter-clockwise seen from the G4-G6 face.]

***Figure 1.** CPENTA6 grid order (Nastran): G1-G2-G3 on the bottom triangle (ζ = −1), G4-G5-G6 on the top triangle (ζ = +1), G4 above G1. Natural axes r (G1→G2), s (G1→G3), ζ (bottom to top) are right-handed; the volume is positive when G1→G2→G3 is counter-clockwise seen from the G4-G6 face.*


    L1 = 1 - r - s,   L2 = r,   L3 = s,      r,s ≥ 0,  r + s ≤ 1,   ζ in [-1, 1]
    N_i     = (1/2) L_i (1 - ζ)        i = 1,2,3     (bottom)
    N_(i+3) = (1/2) L_i (1 + ζ)        i = 1,2,3     (top)
    u^h = Σ N_a d_a ,   d_a = (u,v,w)_a ,   18 translational DOF per element

Enhanced parameters are element-local and eliminated by condensation, so the global system keeps three translations per grid. Engineering-shear Voigt order is used for physical quantities, <span class="m">ε = [εxx, εyy, εzz, γxy, γyz, γzx]</span>, with <span class="m">σ = Dε</span> for isotropic D (λ, μ).

### 2.2 Frozen covariant frame at the centroid

Let J<sub>0</sub> be the Jacobian (rows ∂x/∂ξ<sub>i</sub>) at the centroid (r,s,ζ) = (1/3,1/3,0), and let ε̂ = J<sub>0</sub> ε J<sub>0</sub><sup>T</sup> be the strain tensor in the covariant frame. The natural Voigt vector has the order used throughout this note:

    ε̂ = [ ε̂rr , ε̂ss , ε̂ζζ , 2ε̂rs , 2ε̂sζ , 2ε̂rζ ]^T = T0 ε        (columns e1 ... e6)
    B̂(r,s,ζ) = T0 B(r,s,ζ)         (6 x 18)          D̂0 = T0^-T D T0^-1
    ε = T0^-1 ε̂                                    (return to the physical frame)

T<sub>0</sub> is constant, so energy is evaluated as ∫ (·)<sup>T</sup> D̂<sub>0</sub> (·) dV with no point-dependent transformation. Because T<sub>0</sub> depends only on the element geometry, the element is frame invariant (K rotated by Q equals Q K Q<sup>T</sup> to 4·10<sup>-16</sup>).

### 2.3 Edge-tied assumed transverse shear

On every ζ Gauss plane ζ<sub>g</sub> the two covariant transverse-shear rows are replaced by a field fixed by their tangential values at the three edge midpoints (points t1, t2, t3 in Figure 2):

    a = B̂rζ(1/2, 0 , ζg)
    b = B̂sζ(0 , 1/2, ζg)
    c = b - a - [ B̂sζ - B̂rζ ](1/2, 1/2, ζg)

    B̄rζ(r,s,ζg) = a + c s
    B̄sζ(r,s,ζg) = b - c r              the other four rows keep their compatible values B̂

This is the isotropic MITC3 transverse-shear interpolation of Lee and Bathe (2004), Eqs. (22)–(25), applied to the wedge at each thickness level; the field is constant along each edge and its space is intrinsic to the triangle, hence independent of node numbering (their “spatial isotropy”). A common mistake is to tie at the centroid: that preserves the volume mean of the shear but also the parasitic offset (5-010 gives 4.04 to 5.11 instead of 0.93).

### 2.4 Enhanced strain families

The enhancement follows the EAS form of Simo and Rifai and the natural-frame, determinant-scaled form of Andelfinger and Ramm (their Eq. 17 and Eqs. 32–33):

    ε̂ = B̄ d + (det J0 / det J) M(r̄,s̄,ζ) α          r̄ = r - 1/3,  s̄ = s - 1/3

    family A (1 mode)   : ζ e3                              thickness / Poisson strain, linear in ζ
    family B (2 modes)  : ζ e5 , ζ e6                        transverse shear linear in ζ
    family E (6 modes)  : r̄ e1 , s̄ e1 , r̄ e2 , s̄ e2 , r̄ e4 , s̄ e4    in-plane strain linear in the triangle

    K_dd = Σ w B̄^T D̂0 B̄      K_dα = Σ w B̄^T D̂0 M      K_αα = Σ w M^T D̂0 M
    K    = K_dd - K_dα K_αα^-1 K_αd              α = - K_αα^-1 K_αd d       (w = Gauss weight x det J)

All nine columns have zero mean over the reference volume, and the factor det J<sub>0</sub>/det J makes ∫ M dV vanish, so the enhancement alone cannot disturb the constant-stress condition (orthogonality, Andelfinger–Ramm Eq. 4). The families are closed under renumbering: each is a tensor type (scalar, vector or symmetric tensor) times a function type, so a node permutation maps a family to itself. Family A is the single, in-plane-constant member of the thickness enhancement of the solid-shell element ANS3DEAS (Hauptmann and Schweizerhof 1998). Two further natural families, <span class="m">ζζ ~ (r̄,s̄)</span> and <span class="m">in-plane ~ ζ</span>, were tested and rejected: either one reduces the unit-element rank (to 9 or 10 of 12), producing a mechanism (5-010 factor 0).

### 2.5 Constant-stress consistency parameter

The tied shear field has a different volume mean from the compatible one, which is what removes the parasitic offset and also what the constant-stress patch test checks. The parameter `f = patch_consistency` in [0,1] restores a fraction of that mean:

    ΔB = Σ w ( B̂ - B̄ )[sζ, rζ] / Σ w          (constant over the element)
    B̄[sζ, rζ]  ←  B̄[sζ, rζ] + f ΔB

f = 1 makes ∫(B̄ − B̂) dV = 0, so nodal forces under constant stress equal those of the compatible element and the patch test is exact on any mesh. The price is the thin-bending behaviour (Section 6.1). f = 0, the operating point implemented in the Fortran port, gives the buckling accuracy but passes the patch test only on affine meshes.

### 2.6 Stress recovery, mass and initial stress

Stress is recovered pointwise from the total strain, <span class="m">ε = T0^-1 ( B̄ d + (det J0/det J) M α )</span>, <span class="m">σ = D ε</span> (Andelfinger–Ramm Remark 2a). Mass and the initial-stress (geometric) stiffness are inherited from the standard CPENTA6; the latter is not a tangent derived for the assumed-strain field, so buckling results are those of the physical-prestress operator (cf. Simo and Armero 1992 for the finite-strain extension, not implemented).

## 3. Source attribution by reference

**ABE is a project label, not a literature name.** A, B and E are the names of the enhanced-strain families defined in Section 2.4 (A: ζζ, B: shear, E: in-plane). No journal article proposes this exact combination; each ingredient comes from the sources below, and the combination was selected by the tests in Section 6.

| Reference                                                  | Contribution                                                                                                                                                                                                                          | Used in ABE                                                                      | Not used                                                                                       |
|------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------|
| Lee & Bathe (2004)                                         | Isotropic MITC triangle: constant-along-edge shear (MITC3, Eqs. 22–25), tying at edge centres, “spatial isotropy” (independence of node numbering), isotropic element test                                                            | Shear rows B̄rζ, B̄sζ (Section 2.3), per ζ plane; renumbering test                 | MITC6, in-plane/membrane tying (Eqs. 26–34), shell kinematics                                  |
| Andelfinger & Ramm (1993)                                  | EAS in the natural frame with M = (det J0/det J) T0<sup>-T</sup> Mξ (Eq. 17; 6×6 form Eqs. 32–33); 3D.EAS-9 (modes 25–33); orthogonality condition (Eq. 4); stress-recovery options; remark that ANS shear can still lock when warped | Frozen-frame T0, det-ratio scaling, condensation (Eqs. 8–12), pointwise recovery | HR equivalence, EAS-15/21/30, hexahedral mode set, shear-stiffness scaling for warped elements |
| Simo & Rifai (1990)                                        | EAS framework; incompatible modes as special case; patch-test (orthogonality) interpretation                                                                                                                                          | Enhanced-strain variational basis, zero-mean modes                               | Incompatible-mode forms                                                                        |
| Hauptmann & Schweizerhof (1998)                            | Solid-shell concept: ANS transverse shear plus a linear normal strain through the thickness to avoid thickness locking (§2.1.2, §3.2, ANS3DEAS §4.3.2)                                                                                | Rationale and form of family A                                                   | Plane-stress variants (ANS6z, ANS3DLr), quadratic w′ (ANS3Dq), non-linear algorithm            |
| Simo & Armero (1992)                                       | Geometrically non-linear EAS and incompatible modes                                                                                                                                                                                   | Context only                                                                     | All (finite strain, tangent Kg)                                                                |
| Taylor, Beresford & Wilson (1976); Wilson et al. (1973)    | Incompatible modes, constant-stress (Taylor) correction                                                                                                                                                                               | Context for the orthogonality parameter f                                        | Displacement modes                                                                             |
| MacNeal & Harder (1985); SAP2000 verification 5-001, 5-010 | Patch-test geometry and column-buckling benchmark                                                                                                                                                                                     | Gates and benchmarks                                                             | –                                                                                              |

## 4. Adopted design and measured rationale

Each choice below was made on the SAP2000 5-010 column (first buckling factor, Nastran CPENTA6 0.8910) and checked against rank and invariance gates.

| Variant (shear / enhancement)                                           | 5-010 factor  | Unit-element rank (12 expected)                                              |
|-------------------------------------------------------------------------|---------------|------------------------------------------------------------------------------|
| compatible shear, no enhancement (standard)                             | 20.98         | 12                                                                           |
| centroid-tied shear (one-point in-plane)                                | 4.04          | 11 (hourglass)                                                               |
| EAS-only, best admissible set from a greedy search with rank constraint | 4.07          | 12                                                                           |
| edge-tied shear, no enhancement                                         | 0.9287        | 12                                                                           |
| edge-tied shear + A                                                     | 0.8912        | 12                                                                           |
| edge-tied shear + A + B                                                 | 0.8912        | 12                                                                           |
| **edge-tied shear + A + B + E (ABE)**                                   | **0.8912**    | 12 (min nonzero eigenvalue 1–3·10<sup>-2</sup> of max, up to distortion 0.3) |
| ABE + family F or G                                                     | 0 (mechanism) | 9–10                                                                         |

Family A lowers 5-010 from 0.9287 to 0.8912 (thickness locking, Hauptmann–Schweizerhof §2.1.2). Family B is needed when the prism axis runs along a beam (twist and curved-beam meshes recover from 0.21/0.33 to 0.94/0.94 once B is added). Family E adds no 5-010 change but improves the curved beam (0.866 → 0.978 IN on mesh B). The evidence that edge midpoints, not the centroid, are the right tying points is the Rayleigh quotient of the accurate CPENTA15 mode: 0.89 with the edge-tied field against 10.8 with centroid-tied shear.

## 5. Coordinate systems, Gauss points and integration

<div>
<img width="795" height="530" alt="image" src="https://github.com/user-attachments/assets/d77b6ff9-4404-428d-97d8-e67c13f38257" />
[Figure 2. Triangle plane (r,s). Three Gauss points, the three edge-midpoint tying points used by the assumed shear (re-evaluated on each ζ plane), and the centroid at which the frame J 0 , T 0 and det J 0 are frozen.]

***Figure 2.** Triangle plane (r,s). Three Gauss points, the three edge-midpoint tying points used by the assumed shear (re-evaluated on each ζ plane), and the centroid at which the frame J<sub>0</sub>, T<sub>0</sub> and det J<sub>0</sub> are frozen.*


<div>
<img width="417" height="381" alt="image" src="https://github.com/user-attachments/assets/cb8b0b9e-c1d6-4bbd-8d02-a53a469007e7" />
[Figure 3. The six volume Gauss points: the 3-point triangle rule at the two Gauss planes ζ = ±1/√3. Tying points use the same two planes.]

***Figure 3.** The six volume Gauss points: the 3-point triangle rule at the two Gauss planes ζ = ±1/√3. Tying points use the same two planes.*


| Point | r   | s   | ζ             | weight (before det J) |
|-------|-----|-----|---------------|-----------------------|
| 1     | 1/6 | 1/6 | −0.5773502692 | 1/6                   |
| 2     | 2/3 | 1/6 | −0.5773502692 | 1/6                   |
| 3     | 1/6 | 2/3 | −0.5773502692 | 1/6                   |
| 4     | 1/6 | 1/6 | +0.5773502692 | 1/6                   |
| 5     | 2/3 | 1/6 | +0.5773502692 | 1/6                   |
| 6     | 1/6 | 2/3 | +0.5773502692 | 1/6                   |

The triangle weights sum to 1/2 (the reference area) and the ζ weights to 2, so the rule integrates the reference volume (1) exactly. Tying values are not extra integration points: B̂ is evaluated at t1, t2, t3 for each ζ<sub>g</sub> and interpolated to the Gauss points with the fields of Section 2.3.

## 6. Internal benchmarks: ABE against EAS9 v2

All values below come from the Python reference implementation and the benchmark scripts of this project (ratio to the independent value unless stated). EAS9 v2 is the frame-invariant EAS9 wedge with the same quadrature. “ABE f = 0” is the pure edge-tied field; “ABE f = 1” restores the volume mean (Section 2.5).

### 6.1 SAP2000 5-010 column buckling (12×1×100 in, 48 wedges)

<div>
<img width="726" height="430" alt="image" src="https://github.com/user-attachments/assets/beb4d7ce-7fcc-43cd-bb21-2b7ed252e7da" />

[Figure 4. First buckling factor of the 5-010 column; Nastran CPENTA6 value from the audit archive.]

***Figure 4.** First buckling factor of the 5-010 column; Nastran CPENTA6 value from the audit archive.*


| Mesh (h/t)    | standard | EAS9 v2 | ABE f = 1 | ABE f = 0 |
|---------------|----------|---------|-----------|-----------|
| 2×1×12 (5.0)  | 20.98    | 7.680   | 3.888     | 0.8912    |
| 4×1×24 (2.5)  | 6.258    | 3.414   | 2.216     | 0.8890    |
| 8×1×48 (1.25) | 2.388    | 1.714   | 1.322     | 0.8885    |

With f = 0 the factor converges to the Euler value (0.8883) from above; the other variants converge much more slowly. Orientation of the prism axis for the 2×1×12 mesh:

| Prism axis                          | standard | EAS9 v2 | ABE f = 1 | ABE f = 0 |
|-------------------------------------|----------|---------|-----------|-----------|
| y (triangles in xz, the 5-010 mesh) | 20.98    | 7.680   | 3.888     | 0.8912    |
| z (axis along the column)           | 26.71    | 0.9456  | 0.9127    | 0.9127    |
| x (triangles in the bending plane)  | 79.84    | 79.79   | 79.79     | 79.79     |

Orientation x locks for every wedge tested, ABE included: one constant-strain triangle across the thickness cannot represent in-plane bending, and neither assumed shear nor enhanced strain adds displacement capability.

### 6.2 SAP2000 5-001 constant-stress patch test (real model)

Seven irregular solids, E = 10<sup>6</sup>, ν = 0.25, joint 2 at x = 0.826 (14 wedges, prescribed boundary field u = (2x+y+z, x+2y+z, x+y+2z)/2000). Error is relative to the field magnitude.

| Element          | displacement | stress  | Result |
|------------------|--------------|---------|--------|
| CPENTA6 standard | 4.3e-16      | 1.8e-15 | pass   |
| EAS9 v2          | 4.3e-16      | 2.0e-15 | pass   |
| ABE, f = 1       | 2.0e-15      | 8.3e-14 | pass   |
| ABE, f = 0.5     | 5.9e-2       | 1.9e-1  | fail   |
| ABE, f = 0.25    | 7.9e-2       | 2.7e-1  | fail   |
| ABE, f = 0       | 9.1e-2       | 3.3e-1  | fail   |

For f = 0 the error is a floor rather than a coarse-mesh effect. On cube meshes with interior nodes displaced by a fixed fraction of the cell size, it is 3.3·10<sup>-3</sup>, 7.2·10<sup>-3</sup> and 1.2·10<sup>-2</sup> at 10%, 20% and 30% displacement on the 6×6×6 mesh, and it does not decrease with refinement (n = 2, 3, 4, 6: 8·10<sup>-4</sup>, 2.7·10<sup>-3</sup>, 3.4·10<sup>-3</sup>, 3.3·10<sup>-3</sup> at 10%). Regular and affine-mapped meshes pass at 4·10<sup>-16</sup>.

### 6.3 Other benchmarks (ratio to reference; IN / OUT where two numbers)

| Benchmark                                         | EAS9 v2             | ABE f = 1           | ABE f = 0           | CPENTA15              |
|---------------------------------------------------|---------------------|---------------------|---------------------|-----------------------|
| 5-003 curved beam A (6×1×1)                       | .821 / .785         | .850 / .758         | .851 / .757         | .859 / .926           |
| 5-003 curved beam B (90×1×1)                      | .940 / .895         | .978 / .874         | .978 / .874         | .998 / .954           |
| 5-004 twist 12×2×2                                | .952 / .951         | .960 / .963         | .960 / .963         | .997 / 1.004          |
| 5-004 twist 24×4×2                                | .960 / .963         | .962 / .969         | .962 / .969         | .998 / 1.003          |
| 5-008 cylinder ν = 0.3 / 0.49 / 0.499 / 0.4999    | .972/.976/.980/.981 | .971/.976/.980/.981 | .971/.976/.980/.981 | 1.013/1.003/.984/.875 |
| Cantilever 12×2×2: bending / frequency / buckling | .946/1.023/1.044    | .961/1.017/1.035    | .961/1.017/1.035    | .994/1.001/1.001      |

Apart from the thin-bending buckling case, f = 0 and f = 1 coincide to the third digit; the thin-wedge buckling and the patch test are where they differ.

### 6.4 Element-level gates

| Gate                                              | standard | EAS9 v2         | ABE     |
|---------------------------------------------------|----------|-----------------|---------|
| frame rotation, K vs Q K Q<sup>T</sup>            | 4e-16    | 4e-16           | 4e-16   |
| renumbering (2,3,1) and (3,1,2)                   | 2e-16    | 3.2e-2 / 3.0e-2 | ≤ 1e-13 |
| unit-element rank (regular, distortion 0.15, 0.3) | 12       | 12              | 12      |

## 7. Limitations and open items

- **Patch test and thin bending cannot both be exact.** Exact orthogonality to constant stress forces the volume mean of the tied shear back to the compatible one, and that mean is the thin-bending locking (5-010: 3.89 with f = 1). Intermediate f values are not useful (the patch error falls slowly while the buckling factor rises fast). Tying the covariant components in the local basis instead of the frozen one made the floor 2.5 times worse, and a distortion-only mean correction did not help.
- **Orientation x** (triangles in the bending plane, one element thick) locks like every wedge tested; it needs an enlarged displacement space.
- **One value of f per model:** prestress and buckling must use the same operator. A buckling load with f = 1 on a plate-like wedge mesh is about 3.9 times too high (unconservative).
- The initial-stress stiffness is the standard physical one and the 5-010 prestress was a uniform stress; recovery of the true SOL 101 prestress, SAP2000 5-011 (thermal) and the 5-007 dome were not run in the Python reference. The MYSTRAN dome results are in the NEWSOLID document.
- Warped or strongly distorted elements: Andelfinger and Ramm report that ANS-formulated transverse shear can still lock when warped; this was not studied for the wedge.

## References

1.  P. S. Lee and K. J. Bathe (2004), “Development of MITC isotropic triangular shell finite elements,” *Computers & Structures*, 82 (11–12), 945–962. doi:10.1016/j.compstruc.2004.02.004
2.  U. Andelfinger and E. Ramm (1993), “EAS-elements for two-dimensional, three-dimensional, plate and shell structures and their equivalence to HR-elements,” *International Journal for Numerical Methods in Engineering*, 36 (8), 1311–1337.
3.  J. C. Simo and M. S. Rifai (1990), “A class of mixed assumed strain methods and the method of incompatible modes,” *International Journal for Numerical Methods in Engineering*, 29, 1595–1638.
4.  J. C. Simo and F. Armero (1992), “Geometrically non-linear enhanced strain mixed methods and the method of incompatible modes,” *International Journal for Numerical Methods in Engineering*, 33, 1413–1449.
5.  R. Hauptmann and K. Schweizerhof (1998), “A systematic development of ‘solid-shell’ element formulations for linear and non-linear analyses employing only displacement degrees of freedom,” *International Journal for Numerical Methods in Engineering*, 42, 49–69.
6.  R. L. Taylor, P. J. Beresford and E. L. Wilson (1976), “A non-conforming element for stress analysis,” *International Journal for Numerical Methods in Engineering*, 10, 1211–1220 (as cited in ref. 2).
7.  E. L. Wilson, R. L. Taylor, W. P. Doherty and J. Ghaboussi (1973), “Incompatible displacement models,” in S. T. Fenves et al. (eds.), *Numerical and Computational Methods in Structural Mechanics*, Academic Press, 43–57 (as cited in ref. 2).
8.  E. N. Dvorkin and K. J. Bathe (1984), “A continuum mechanics based four-node shell element for general nonlinear analysis,” *Engineering Computations*, 1, 77–88; K. J. Bathe and E. N. Dvorkin (1986), “A formulation of general shell elements—the use of mixed interpolation of tensorial components,” *International Journal for Numerical Methods in Engineering*, 22, 697–722 (origin of ANS/MITC, as cited in ref. 1).
9.  N. Büchter and E. Ramm (1992), “3D-extension of nonlinear shell equations based on the enhanced assumed strain concept,” in C. Hirsch, J. Périaux and E. Oñate (eds.), *Computational Methods in Applied Sciences*, Elsevier, Amsterdam (as cited in ref. 5; origin of the thickness enhancement).
10. R. H. MacNeal and R. L. Harder (1985), “A proposed standard set of problems to test finite element accuracy,” *Finite Elements in Analysis and Design*, 1, 3–20 (as named in the SAP2000 Example 5-001 text).
11. Computers and Structures, Inc., *SAP2000 Verification: Solids*, Examples 5-001, 5-003, 5-004 and 5-010.

**Note on provenance.** References 1, 2, 4 and 5 were read in full text. References 3, 6–10 are cited as listed in those papers or in the SAP2000 text and were not re-read. Equation numbers quoted above refer to refs. 1 and 2.

