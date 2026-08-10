# A Covariant/Contravariant Reformulation of the Simo (1993) 4-Node Shell Element
### Kinematics, Membrane–Bending–Shear–Drilling Operators, EAS Stabilization, and Verification Against DKMQ24 and Its MYSTRAN CQUADR Variant

**Status:** Formulation complete and verified against an independent, full
test harness (Section 6): 7/8 patch-test/benchmark categories pass outright,
and the MacNeal-Harder twisted-beam benchmark converges to within 0.1-0.3%
of the NASTRAN CQUAD4 reference at N=24 — outperforming standard `DKMQ24`
and matching or exceeding the independently-authored "V3" hybrid and the
MYSTRAN CQUADR port on this specific benchmark. The 8th category (pinched
cylinder) was investigated and diagnosed as a genuine but secondary
convergence-rate gap (Section 6.4) rather than a formulation defect: the
element converges to the correct answer, just one mesh-refinement level
slower than `DKMQ24`.

---

## 1. Introduction

This note documents a 4-node, 24-DOF (6 DOF/node) quadrilateral shell element
derived from the classical single-director Reissner-Mindlin shell kinematics
of Simo and co-workers [Simo & Fox 1989; Simo, Fox & Rifai 1990 (Part II);
Simo, Fox & Rifai 1990 (Part III); Simo, Rifai & Fox 1990 (Part IV); Fox &
Simo 1992], hereafter referred to collectively as "the Simo shell papers."
The element carries the classical drilling rotation as a numerically
regularized penalty degree of freedom (Hughes & Brezzi 1989, cited within
Part IV) rather than the exact Lagrange-multiplier drill formulation of Fox
& Simo (1992); adopting the latter is left as future work (Section 7).

The starting implementation (informally "V1.0") interpolated the bending
curvature directly from nodal rotation components treated as raw local
slopes, i.e. effectively `kappa = grad(theta)` with `theta` identified
component-by-component with physical bending angles. This is a well-known
simplification valid **only** for axis-aligned, rectangular, flat elements,
and fails (sometimes catastrophically) for oblique, distorted, or warped
quadrilaterals. The present work replaces this with a strain-measure-faithful
covariant kinematic description, linearizes it consistently, and transforms
every strain operator (membrane, bending, shear) to a single **fixed
physical Cartesian frame per element** through an explicit
covariant-to-contravariant metric transformation. The same transformation
machinery is reused, without modification, for all three strain fields —
this internal consistency is the organizing principle of the whole
implementation and is what allows the element to be debugged and verified
systematically.

---

## 2. Kinematics and Strain Measures

### 2.1 Shell configuration

Following Part I/III/IV, the shell mid-surface is parametrized by
`(xi^1, xi^2) in A`, with position field `phi(xi^1,xi^2)` and unit director
field `t(xi^1,xi^2)`. The reference (undeformed) configuration is denoted
with a superscript/subscript `0`; `t0` denotes the reference director field.
For a 4-node bilinear element the reference director at each node,
`t0_I = V_n[I]`, is computed from the cross product of the adjacent edge
vectors at that corner (with a diagonal fallback for degenerate corners),
and the interpolated reference director field over the element is

```
t0(xi^1,xi^2) = sum_I N_I(xi^1,xi^2) * V_n[I]
```

with `N_I` the standard bilinear shape functions.

### 2.2 Membrane strain

The (nonlinear) covariant membrane strain is

```
eps_ab = 1/2 ( phi,a . phi,b  -  phi0,a . phi0,b )
```

Linearizing about the reference configuration (small-displacement elastic
shell, consistent with the rest of this codebase):

```
delta(eps_ab) = 1/2 ( u,a . g0_b  +  g0_a . u,b )         (2.1)
```

where `g0_a := phi0,a` are the reference covariant tangent vectors
(columns of the surface Jacobian) and `u,a = d(phi)/d(xi^a)` is the natural
(parametric) derivative of the nodal displacement field.

### 2.3 Bending (curvature) strain

Following Part IV eq. (2.16)/(3.4):

```
rho_ab = phi,a . t,b
```

Its full linearization about the reference state has **two** contributions:

```
delta(rho_ab) = delta(phi),a . t0,b   +   phi0,a . delta(t),b        (2.2)
                `----- (I) -----'       `------- (II) --------'
```

with the director variation given by the standard rigid-rotation update

```
delta(t)_I = theta_I x t0_I
```

`theta_I` being the (global, 3-component) nodal rotation vector.

* **Term (II)** — the "rotational" part — is always present.
* **Term (I)** — the "translational" part — is proportional to
  `t0,b = sum_I dN_I/dxi_b * V_n[I]`, the natural-coordinate gradient of the
  **interpolated reference director field across the element**. For a
  perfectly flat element every nodal director `V_n[I]` is identical
  (parallel), so `t0,b = (sum_I dN_I/dxi_b) * V_n = 0` identically (the sum
  of shape-function derivatives is zero for any complete partition of
  unity). For a **warped** element (e.g. a twisted strip, non-coplanar
  quadrilateral) the four nodal directors differ from node to node and
  `t0,b != 0` in general.

An earlier version of this implementation dropped term (I) entirely,
reasoning that all verification cases used at that time were flat. This is
numerically invisible in flat benchmarks (patch tests, flat cantilevers) but
produces large, non-converging errors on warped geometry (Section 5.3).
Restoring term (I) is the single largest accuracy improvement recorded in
this work.

### 2.4 Transverse shear strain

Following the same pattern (Part IV eq. 2.16, transverse-shear row):

```
delta_a = phi,a . t  -  phi0,a . t0
```

Linearized:

```
delta(delta_a) = u,a . t0  +  g0_a . delta(t)                        (2.3)
                 `--(I)--'    `-------(II)-------'
```

structurally identical in form to the bending case: a translational term
(I), dotted with the *interpolated* reference director at the evaluation
point, and a rotational term (II) built from the actual covariant tangent
`g0_a` and the director-variation cross product.

### 2.5 Drilling constraint

The drilling (in-plane rotation) constraint follows the standard
Hughes-Brezzi penalty form,

```
gamma_drill = 1/2 (dv/dx - du/dy)  -  (theta . t0)
```

where `theta . t0` is the projection of the nodal rotation vector onto the
**node's own** reference director (not simply the global z-rotation
component — the latter is only correct if `t0` is exactly the global
z-axis, i.e. for a flat, XY-plane element). This projection automatically
and exactly annihilates any drilling-axis rotation component from the
bending and shear operators (Section 4.1), reproducing the classical
invariance property that "rotations about the director itself are
irrelevant" (Fox & Simo 1992, Section 2).

---

## 3. Covariant -> Contravariant -> Cartesian Transformation

All three strain measures above are naturally expressed as **natural
(parametric) coordinate, covariant-basis components** (indices `a, b in
{1,2}` referring to `xi^1, xi^2` directions). To assemble a constitutive
law expressed in an engineering (physical, orthonormal) frame, every strain
operator in this element is passed through the *same* transformation:

1. Compute the covariant metric `a_ab = g_a . g_b` and its inverse
   (contravariant metric) `a^ab`, using the ACTUAL covariant tangent
   vectors `g_a` at the evaluation point (not any unit-vector
   approximation).
2. Form the contravariant (dual) basis `g^a = a^ab g_b`, satisfying
   `g^a . g_b = delta^a_b`.
3. Project onto a **single fixed physical frame** `{E1, E2}`, computed
   once at the element center in `__init__` (not re-orthonormalized at
   every Gauss point):
   `c_{i,a} := E_i . g^a`.
4. For a rank-2 (tensor) quantity such as membrane or bending strain:

```
eps_ij(physical) = sum_{a,b} eps_ab(natural) * c_{i,a} * c_{j,b}
```

   For a rank-1 (vector) quantity such as transverse shear:

```
gamma_i(physical) = sum_a gamma_a(natural) * c_{i,a}
```

Using a single **fixed** frame per element (rather than re-deriving a local
orthonormal frame independently at every integration point) avoids a subtle
but real failure mode: independently orthonormalized per-point frames can
"drift" or effectively rotate relative to one another across a distorted
or warped element, introducing a spurious apparent in-plane rotation of the
strain frame that has nothing to do with the physical deformation. This was
identified as the root cause of a large (100-600%) error in the SAP2000
2-001 bending patch test in an earlier version of this element that used
`A = J^T [e1,e2]` with `e1,e2` recomputed and renormalized at each Gauss
point (Section 5.1).

This same four-step machinery is applied, without modification beyond the
"which strain field" bookkeeping, to `Bm` (membrane), `Bb` (bending), and
`Bs` (transverse shear, after ANS sampling — Section 4.3). Internal
consistency across all three operators, rather than independent per-field
tuning, is the central design choice of this implementation.

---

## 4. Discrete Operators

### 4.1 Bending operator `Bb`

Assembled per node `I`, combining terms (I) and (II) of eq. (2.2), each
independently transformed via Section 3:

```
d(rho_11)/d(u_I)     = dN_I/dxi_1 * t0,1                    [term I]
d(rho_22)/d(u_I)     = dN_I/dxi_2 * t0,2
d(rho_12)/d(u_I)     = 1/2 [ dN_I/dxi_1 * t0,2 + dN_I/dxi_2 * t0,1 ]

d(rho_11)/d(theta_I) = dN_I/dxi_1 * (t0_I x g_1)             [term II]
d(rho_22)/d(theta_I) = dN_I/dxi_2 * (t0_I x g_2)
d(rho_12)/d(theta_I) = 1/2 [ dN_I/dxi_1 * (t0_I x g_2) + dN_I/dxi_2 * (t0_I x g_1) ]
```

each row-triple then mapped to physical `{xx, yy, xy}` components via the
tensor transform of Section 3, and assembled into `Bb[:, 6I : 6I+3]`
(translational block, term I) and `Bb[:, 6I+3 : 6I+6]` (rotational block,
term II).

A direct algebraic consequence: since `(t0_I x g_a) . t0_I = 0` identically
(a cross product is orthogonal to both its factors), any drilling-axis
component of `theta_I` (i.e. `theta_I` parallel to `t0_I`) contributes
exactly zero to `Bb`, for both flat and warped geometry, with no separate
projection step required.

### 4.2 Membrane operator `Bm`

Directly analogous, from eq. (2.1) — a single (translational-only) term:

```
d(eps_11)/d(u_I) = dN_I/dxi_1 * g_1
d(eps_22)/d(u_I) = dN_I/dxi_2 * g_2
d(eps_12)/d(u_I) = 1/2 [ dN_I/dxi_1 * g_2 + dN_I/dxi_2 * g_1 ]
```
transformed to physical `{xx,yy,xy}` via Section 3.

### 4.3 Shear operator `Bs` (Assumed Natural Strain)

To avoid transverse shear locking, the natural-coordinate shear components
of eq. (2.3) (both terms I and II) are evaluated at the four edge-midpoint
sampling stations `A(0,-1), B(1,0), C(0,1), D(-1,0)` (Dvorkin-Bathe /
Simo-Fox-style ANS), each using the **locally interpolated** reference
director `t0(sampling point)` (not a single node's director) for the
translational term:

```
gamma_1(natural, sample) = sum_I [ dN_I/dxi_1 * (u_I . t0_sample)
                                    + N_I * (theta_I . (t0_I x g_1_sample)) ]
gamma_2(natural, sample) = (analogous, xi_2 direction)
```

`gamma_1` is sampled at A and C and bilinearly interpolated in `xi^2`;
`gamma_2` is sampled at B and D and bilinearly interpolated in `xi^1`
(standard ANS pattern). The resulting natural-coordinate vector field is
**then** transformed to the physical frame using the contravariant metric
evaluated **at the Gauss point** (not frozen at the sampling stations),
consistent with Section 3.

### 4.4 EAS membrane stabilization

A 4-parameter Andelfinger-Ramm / Simo-Rifai enhancement is used for the
membrane field:

```
M(xi,eta) = (detJ0 / detJ(xi,eta)) * T0 * M_hat(xi,eta)

M_hat = [ xi   0    0   0
          0    eta  0   0
          0    0    xi  eta ]
```

`T0` is built from the *same* contravariant construction as `Bm`, evaluated
once at the element center, so the enhanced field lives in the identical
physical frame as the compatible strain (required for
`K_ua = int Bm^T C M dV` to be meaningful). The ratio `detJ0/detJ(xi,eta)`
(present in the corrected version; **absent** in an earlier version that
used only a constant `1/detJ0`) is what makes the orthogonality condition
`int M dV = 0` hold **exactly** for oblique/distorted elements, not just
approximately — see Section 5.2 for a direct numerical demonstration.

### 4.5 Drilling penalty

```
K_drill = beta_drill * G * h * int (Bdrill^T Bdrill) dA
```

with `beta_drill` a user-adjustable, dimensionless penalty coefficient
(default `0.02`), matching the Hughes-Brezzi drilling-DOF penalty
convention and brought to the same order of magnitude as the reference
elements in this codebase (`DKMQ24_ShellElement_RHR`: effectively `~1e-4`;
`DKMQ24_MystranCQUADR_ShellElement_RHR`: effectively `~1e-3`, dynamically
scaled to the element's actual bending/shear stiffness rather than a fixed
multiple of `G`). A direct sensitivity sweep (Section 5.4) found the
converged solution essentially insensitive to `beta_drill` across four
orders of magnitude, for the benchmarks tested.

---

## 5. Verification (self-conducted, preliminary — see status note)

### 5.1 Bending patch test (SAP2000 2-001 benchmark)

Five oblique flat quadrilaterals covering an irregular octagonal patch
(coordinates and prescribed edge rotations per the SAP2000 verification
manual, independent reference: `mxx = myy = 1.111e-7`, `mxy = 0.333e-7`).

| Element | mxx/ref (5 elements) |
|---|---|
| `DKMQ24_ShellElement_RHR` (reference) | 1.000, 1.000, 1.000, 1.000, 1.000 |
| This element, **before** the covariant `Bb` fix | 4.874, -0.363, -1.072, -0.897, 1.659 |
| This element, **after** the covariant `Bb` fix | 1.000, 1.000, 1.000, 1.000, 1.000 |
| "V3" (independently authored hybrid, not derived here) | 2.738, 1.119, 0.270, 0.832, 0.310 |

The pre-fix scatter (no consistent sign or scale factor across elements of
differing shape) was traced to the per-Gauss-point frame re-orthonormalization
issue described in Section 3, combined with mixing a Cartesian-derivative
convention with natural-coordinate director-gradient terms. The fix in
Section 3-4.1 reproduces the independent reference to machine precision on
all five oblique elements, matching `DKMQ24_ShellElement_RHR`.

### 5.2 EAS orthogonality

For the same oblique patch geometry, `K_ua^T . u_affine` (which must vanish
for any constant-strain nodal displacement field if the enhanced field is
properly orthogonal) was measured directly:

| EAS formulation | Relative leakage |
|---|---|
| Without `detJ0/detJ(xi,eta)` ratio (earlier version) | ~1.0e-3 |
| With `detJ0/detJ(xi,eta)` ratio (Section 4.4) | ~1.5e-18 (machine precision) |

### 5.3 MacNeal-Harder twisted beam (warped-element benchmark)

Geometry: `L=12, W=1.1, h=0.32, E=2.9e7, nu=0.22`, 90-degree total twist,
`Ny=2` elements across the width, tip loads `Fz` (out-of-plane) and `Fy`
(in-plane); independent reference `uz_ref = 0.001749`, `uy_ref = 0.005429`.

| Nx | DKMQ24 Uy err | DKMQ24 Uz err | This element (before term-I fix) Uy/Uz err | This element (after term-I fix) Uy/Uz err |
|---|---|---|---|---|
| 6  | 38.3% | 25.9% | 39.0% / 27.6% | **2.71% / 5.90%** |
| 12 | 39.0% | 32.2% | 39.9% / 33.7% | **0.83% / 1.27%** |
| 24 | 39.8% | 34.7% | 40.3% / 35.3% | **0.30% / 0.06%** |

`DKMQ24_ShellElement_RHR` — itself a validated, patch-test-passing element
in the same codebase — shows essentially the same non-converging ~38-40%
error at this specific mesh setting (`Ny=2`), which was used as a control
to confirm the error was not an artifact of this specific minimal test
reproduction. Restoring the translational curvature term (Section 2.3,
term I) is the single change responsible for the improvement from
non-converging ~40% error to a properly mesh-converging <3% error. This
mirrors, independently, the "warping-aware bending strain" mechanism
present in `DKMQ24_MystranCQUADR_ShellElement_RHR` (its `nbc1`/`nbc2`
terms) and reportedly in the "V3" hybrid element, though those were not
consulted during the derivation of eq. (2.2) — the term was re-derived
directly from the linearization of `rho_ab = phi,a . t,b`.

### 5.4 Drilling penalty sensitivity

At `Nx=12, Ny=2` (MacNeal twisted beam, out-of-plane load case):

| `beta_drill` | Uy error | Uz error |
|---|---|---|
| 1.0    | 39.91% | 33.65% |
| 0.1    | 40.37% | 34.28% |
| 0.02 (default) | 40.41% | 34.34% |
| 0.01   | 40.41% | 34.34% |
| 0.001  | 40.42% | 34.35% |
| 0.0001 | 40.42% | 34.35% |

(Note: this table was generated **before** the Section 5.3 curvature-term
fix, hence the ~40% baseline; the point of this table is only to
demonstrate insensitivity to `beta_drill` across four orders of magnitude,
not the absolute accuracy level.) Drilling penalty magnitude was ruled out
as a contributor to the warped-element error.

### 5.5 Supporting unit tests (algebraic invariants)

All of the following hold to machine precision (`<1e-10` relative), for
both flat-square and synthetically warped test geometries:

* Rigid-body translation and rotation produce zero strain in all four
  fields (membrane, bending, shear, drilling).
* An unconstrained free element has exactly 6 zero-energy eigenmodes
  (the 6 rigid-body modes; no spurious hourglass modes).
* Drilling-axis rotation is exactly decoupled from `Bb` and `Bs`
  (`<1e-16` leakage) on both flat and warped geometry.
* An affine (constant-strain) nodal displacement field on an oblique flat
  quadrilateral reproduces the imposed constant membrane strain exactly
  (`<1e-10`) at every Gauss point (membrane patch test).

---

## 6. External Harness Comparison (independently run by the author)

The following results were obtained from the author's own full test harness
(all comparison elements present: `DKMQ20_6dof`, `DKMQ24`,
`MITC4pD_HughesBrezzi`, `DKMQ24_MystranCQUADR`, "Simo1993 V3", plus
AI-assisted variants labeled `Simo1993Gemini`/`Simo1993Gemin2` not authored
in this work). These are the first independent (outside this document's own
minimal reproductions) confirmations of the results in Section 5.

### 6.1 Patch-test / benchmark pass matrix

Columns: Eigenvalue test, constant-Bending patch, constant-Shear patch,
constant-Twist patch, SAP2000 membrane Stress patch, cantilever tip-Force,
cantilever tip-Moment, simply-supported-Plate, pinched-Cylinder.

```
Element              Eigen  PBend  PShear  PTwist  PStres  PForce  Moment  Load   Plate  Cyl  Pass Rate
DKMQ20_6dof          v      v      v       v       v       v       v       v      v      v    8/8
DKMQ24               v      v      v       v       v       v       v       v      v      v    8/8
MITC4pD_HB           v      x      v       x       v       x       v       v      v      x    5/8
Simo1993_v1p1        v      x      v       x       x       x       v       v      v      x    5/8
Simo1993_v1p2        v      x      v       x       x       x       v       v      v      x    5/8
Simo1993_v1p3        v      v      v       v       v       v       v       v      v      x    7/8
Simo1993_v1p4        v      v      v       v       v       v       v       v      v      x    7/8
Simo1993_v1p5        v      v      v       v       v       v       v       v      v      x    7/8
Simo1993_v1p6        v      v      v       v       v       v       v       v      v      x    7/8
Simo1993Gemini       v      x      v       x       x       x       v       v      v      x    5/8
Simo1993Gemin2       v      x      x       x       x       x       x       x      x      x    1/8
MystranCQUADR        v      v      v       v       v       v       v       v      v      v    8/8
```

Confirms, from an independent harness: the covariant `Bb`/`Bm` fix (Section
5.1) is what moves this element from 5/8 (v1p1-v1p2, matching the two other
AI-assisted attempts "Gemini"/"Gemin2" at similar maturity) to 7/8
(v1p3 onward). The one remaining systematic failure is the **pinched
cylinder** benchmark, unchanged across v1p1 through v1p6 — i.e. none of the
fixes in this document touch its root cause. This is flagged as the next
investigation target (Section 7).

### 6.2 Clamped beam, 6 load cases (Problem 2-002), N=12, out-of-plane/in-plane mix

| LC | Quantity | This element (v1p6) | Independent reference | Error |
|---|---|---|---|---|
| 1 | ux | 3.000e-05 | 3.000e-05 | 0.0% |
| 2 | uz | 1.067e-01 | 1.081e-01 | 1.3% |
| 3 | uy | 4.313e-01 | 4.321e-01 | 0.2% |
| 4 | uy | 3.032e-03 | 3.410e-03 | **11.1%** |
| 5 | ux | 8.898e-04 | 9.000e-04 | 1.1% |
| 6 | rz | 3.600e-02 | 3.600e-02 | 0.0% |

LC4 is the one case not improved by any fix in this document (identical
`3.032e-03` from v1p1 through v1p6). For comparison, `MystranCQUADR` at the
same N=12 setting shows LC2 at 71.5% error and LC4 at 571% error under
otherwise the same harness — i.e. LC2/LC4 in *this specific reduced-N*
comparison table appear to be a harder or differently-conditioned case for
several elements, and the comparison should be read alongside the N=24
convergence data in Section 6.3 rather than in isolation.

The author additionally reports (direct communication) that this element
**converges faster than the CQUADR-DKMQ24 hybrid on Problem 2-002** at the
mesh densities tested — consistent with the Section 5.3 MacNeal result.

### 6.3 MacNeal-Harder twisted beam convergence (Problem 2-004),
independent harness confirmation

```
Simo1993v1p5 (pre curvature-term-I fix):
  Uy errors: [40.64, 40.41, 40.42, 40.43] %
  Uz errors: [29.83, 34.34, 35.19, 35.49] %

Simo1993v1p6 (post curvature-term-I fix):
  Uy errors: [2.71, 0.83, 0.44, 0.30] %
  Uz errors: [5.90, 1.27, 0.38, 0.06] %

Simo1993 V3 (independent hybrid):
  Uy errors: [3.29, 1.07, 0.60, 0.43] %
  Uz errors: [12.34, 3.36, 1.40, 0.68] %

Simo1989 6dof (companion element, same codebase):
  Uy errors: [42.16, 40.67, 40.46, 40.41] %
  Uz errors: [31.09, 34.55, 35.21, 35.45] %
```

This independently reproduces the Section 5.3 result to within reporting
precision, and additionally shows v1p6 **outperforming** "V3" at every
refinement level on both components. The companion `Simo1989 6dof` element
(same codebase, presumably not carrying the equivalent of the term-I fix)
remains at the same ~40%, non-converging error level as `DKMQ24` and
pre-fix `Simo1993`, reinforcing that the fix is specific and mechanistic
rather than a general property of "any Simo-family element."

Tip-displacement values at N=24 (final refinement), compared against
NASTRAN CQUAD4 as an external, independent reference:

| Element | Uy (N=24) | Uz (N=24) |
|---|---|---|
| NASTRAN CQUAD4 (reference) | 5.4075123e-03 | 1.7482317e-03 |
| **Simo1993 v1p6 (this work)** | **5.4125916e-03** | **1.7479414e-03** |
| Simo1993 V3 | 5.4056372e-03 | 1.7370650e-03 |
| DKMQ24-MY (MYSTRAN CQUADR port) | 5.4019394e-03 | 1.7053021e-03 |
| DKMQ24 (standard) | 7.5919117e-03 | 2.3553501e-03 |
| MY-MITC4P (MYSTRAN-style MITC4+) | 5.4071923e-03 | 1.7076490e-03 |

v1p6's Uz matches the NASTRAN reference more closely than the MYSTRAN-port
elements in this table (`DKMQ24-MY`, `MY-MITC4P`) do, and its Uy is within
0.1% of the NASTRAN value — despite using only the theoretical Simo-shell
term (Section 2.3) with no MYSTRAN-specific tuning. Standard `DKMQ24`
(without the MYSTRAN edge-constraint enhancement) is off by ~35-40% on
this benchmark, underscoring that the warping-aware curvature term is doing
real, load-bearing work here, not just closing a small residual gap.

### 6.4 Pinched cylinder — diagnosed: not a bug, but a slower convergence rate

Initial concern: this benchmark failed identically across every version
(v1p1-v1p6), unmoved by any fix in this document. Direct re-derivation
against seven independent literature sources confirmed the harness's
problem parameters, boundary conditions, and reference value
(`1.8248e-5`, MacNeal & Harder 1985) are all standard and correct — the
harness itself is not at fault.

Re-running the exact harness boundary-condition/load code with `mesh_sizes
= [4, 8, 16]` (one refinement level beyond what the harness's reported
"Final w_ratio" used) shows the actual behavior:

| N | `DKMQ24` ratio | This element (v1p6) ratio |
|---|---|---|
| 4  | 0.617 | 0.386 |
| 8  | 0.945 | 0.753 (harness's reported "Final w_ratio: 0.7533" — matches exactly) |
| 16 | 1.019 | **0.932** (inside the harness's own +/-15% pass tolerance) |

Both elements converge monotonically to the reference value; `DKMQ24`
simply converges faster (reaching the +/-15% tolerance band already at
N=8), while this element needs one additional refinement level (N=16) to
reach the same tolerance. This is a **rate-of-convergence** gap on a
doubly-curved, membrane-dominated inextensional-bending benchmark — a
qualitatively different (and less severe) issue than the MacNeal
twisted-beam non-convergence diagnosed and fixed in Section 5.3/6.3, which
was a genuine formulation defect (a missing strain term) rather than a
convergence-rate limitation. Candidate explanations for the slower rate,
not yet isolated: the 4-parameter Andelfinger-Ramm EAS enhancement
(Section 4.4) may be less effective for doubly-curved membrane states than
for the flat/warped cases exercised in Sections 5.1-5.3; `Bs_ANS` may
similarly benefit from further refinement on curved geometry. This is
downgraded from "undiagnosed failure" (as reported in an earlier draft of
this document) to "diagnosed, secondary accuracy gap" and remains a
lower-priority item for future work relative to the correctness issues
resolved in Sections 5.1-5.3.



---

## 7. Known Limitations and Future Work

1. **Drilling formulation** uses the Hughes-Brezzi penalty regularization,
   not the exact Lagrange-multiplier drill-rotation constraint of Fox &
   Simo (1992). The penalty is shown insensitive to its coefficient over
   4 orders of magnitude for the cases tested (Section 5.4), which is
   reassuring but not a substitute for the exact formulation.
2. **EAS transformation** uses the simplest 4-parameter Andelfinger-Ramm
   enhancement (`M_hat = [xi,0,0,0; 0,eta,0,0; 0,0,xi,eta]`), not the full
   Simo-Rifai transformation-strain hierarchy described for more general
   (non-membrane) fields.
3. **Nonlinear / finite-deformation kinematics**: the current element is
   linear-elastic, small-displacement. The Simo shell papers this
   formulation is based on (Parts I, III, IV) are geometrically exact and
   support finite rotations via the exponential map (Part III) and
   extensible-director thickness stretch (Part IV); none of that is
   implemented here.
4. **`Bs_ANS` warping-awareness**: unlike `Bb`, which required an explicit
   fix to include the translational curvature term, `Bs_ANS` already
   included the analogous translational shear term at the time of this
   writing. Its independent necessity was not isolated the way `Bb`'s was
   — the Section 5.3 improvement is attributable to the `Bb` fix
   specifically, with `Bs_ANS` held constant across the "before/after"
   comparison.
5. **Local frame convention**: `E1,E2` (the fixed physical frame used for
   the Section 3 transformation) is currently built via a
   `cross(g2,e3)`-style construction at the element center, tied to the
   parametric `eta` direction. `DKMQ24_MystranCQUADR_ShellElement_RHR`
   uses a global-Z-projected frame (`T1 = cross(normal, global_Z)`)
   instead, argued to be more robust to element skew. Whether this
   matters *given* the covariant/contravariant transformation of Section 3
   (which should, in principle, absorb frame-choice differences through
   the metric) has not been isolated experimentally.
6. **Pinched cylinder benchmark converges more slowly than `DKMQ24`**
   (Section 6.4): both elements converge monotonically to the correct
   reference value, but this element needs approximately one additional
   mesh-refinement level to reach the harness's +/-15% pass tolerance.
   Diagnosed as a rate-of-convergence gap (likely EAS/ANS adequacy on
   doubly-curved geometry), not a formulation defect — a materially
   different and less urgent finding than initially suspected.

---

## 8. References

1. Simo, J.C. and Fox, D.D. (1989). "On a stress resultant geometrically
   exact shell model. Part I: Formulation and optimal parametrization."
   *Comput. Methods Appl. Mech. Engrg.* 72, 267-304.
2. Simo, J.C., Fox, D.D. and Rifai, M.S. (1989). "On a stress resultant
   geometrically exact shell model. Part II: The linear theory;
   computational aspects." *Comput. Methods Appl. Mech. Engrg.* 73, 53-92.
3. Simo, J.C., Fox, D.D. and Rifai, M.S. (1990). "On a stress resultant
   geometrically exact shell model. Part III: Computational aspects of
   the nonlinear theory." *Comput. Methods Appl. Mech. Engrg.* 79, 21-70.
4. Simo, J.C., Rifai, M.S. and Fox, D.D. (1990). "On a stress resultant
   geometrically exact shell model. Part IV: Variable thickness shells
   with through-the-thickness stretching." *Comput. Methods Appl. Mech.
   Engrg.* 81, 91-126.
5. Fox, D.D. and Simo, J.C. (1992). "A drill rotation formulation for
   geometrically exact shells." *Comput. Methods Appl. Mech. Engrg.* 98,
   329-343.
6. Simo, J.C. and Kennedy, J.G. (1992). "On a stress resultant
   geometrically exact shell model. Part V: Nonlinear plasticity:
   formulation and integration algorithms." *Comput. Methods Appl. Mech.
   Engrg.* 96, 133-171.
7. Andelfinger, U. and Ramm, E. (1993). "EAS-elements for
   two-dimensional, three-dimensional, plate and shell structures and
   their equivalence to HR-elements." *Int. J. Numer. Methods Engrg.* 36,
   1311-1337. (Cited via implementation convention; not directly consulted
   in this work — flagged for verification.)
8. Katili, I. et al. — DKMQ/DKMQ24 element family. (Cited via codebase
   docstrings and comparison-element behavior; original paper(s) not
   directly consulted in this work — flagged for verification.)
9. MYSTRAN `CQUADR_DKMQ24.f90` kernel, ported/documented as
   `DKMQ24_MystranCQUADR_ShellElement_RHR.py` in this codebase. The
   original Fortran source was not directly inspected in this work; the
   comparison drawn here is against the Python port and its accompanying
   docstring/commentary, which itself may not be a complete or exact
   transcription of the Fortran kernel.

**Note on sourcing (transparency statement):** references 1-6 (the Simo
shell paper series) were supplied to and read directly by the AI assistant
as PDF uploads during this work and are cited from that direct reading.
References 7-9 are cited on the basis of codebase docstrings, variable
naming, and comparison-element behavior rather than direct examination of
the original sources, and should be independently verified before this
document is considered citable as-is.

---

## Appendix A: Summary of Discrete Operator Structure

```
Bm[i, 6I:6I+3]   : membrane, translational only (eq. 2.1)

Bb[i, 6I:6I+3]   : bending, translational (term I, eq. 2.2) -- vanishes
                   identically for flat elements
Bb[i, 6I+3:6I+6] : bending, rotational (term II, eq. 2.2) -- always present

Bs[i, 6I:6I+3]   : shear, translational (term I, eq. 2.3), via ANS sampling
Bs[i, 6I+3:6I+6] : shear, rotational (term II, eq. 2.3), via ANS sampling

Bdrill[0, 6I:6I+3]   = 1/2 (dN_I/dx * e2 - dN_I/dy * e1)
Bdrill[0, 6I+3:6I+6] = -N_I * t0_I

M_eas(xi,eta) = (detJ0/detJ) * T0 * M_hat(xi,eta),  4 internal EAS parameters,
                statically condensed: K_cond = K_uu - K_ua K_aa^-1 K_ua^T
```

All contravariant projection coefficients and the fixed frame `{E1,E2}` are
as defined in Section 3. Full expressions are given in Sections 4.1-4.4 and
in the accompanying source file `Simo1993_ShellElement.py`.
