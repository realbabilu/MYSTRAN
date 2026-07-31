# CQUADR SIMO1993 formulation and validation notes

Date: July 31, 2026

This note documents the current MYSTRAN `CQUADR` Simo1993 branch:

- selector: `PARAM,QUADRTYP,SIMO`
- dispatcher: `Source/EMG/EMG1/EMG.f90`
- kernel: `Source/EMG/EMG4/CQUADR_SIMO1993.f90`
- Python reference: `D:\18a\bending_only\Shell\gemini2\\Simo1993_ShellElement_v1p6.py`

The implementation is intentionally treated as an isolated `CQUADR` kernel. It
does not call the `CQUAD4` MITC4/MITC4+ path and does not reuse the DKMQ24
kernel.

## Summary

`CQUADR_SIMO1993` is a 4-node, 6-DOF-per-node shell element based on the
Simo/Fox/Simo-Armero-Taylor style director kinematics. The current MYSTRAN port
tracks the Python `v1p6` behavior:

- membrane strain uses covariant-to-contravariant tensor transformation;
- bending curvature includes both rotational director curvature and the
  translational curvature term from nodal normal variation;
- transverse shear uses Simo/Dvorkin-Bathe style ANS sampling in natural
  coordinates, then transforms to the physical frame;
- membrane EAS uses a four-parameter Q1E4 enhanced field and static
  condensation;
- drilling stiffness is a small Hughes-Brezzi-style penalty with default
  scale `beta_drill = 0.02`.

The element currently passes the thick shell patch deck `problem_2_001` in
MYSTRAN F06 engineering-force output.

## Kinematics and DOF layout

Each node has the MYSTRAN shell slot:

```text
u_I = [ux, uy, uz, rx, ry, rz]
```

The reference director at node `I` is the nodal normal `t0_I`. In the Python
reference this is `self.V_n[I]`; in Fortran it is computed by
`CALC_NODAL_NORMALS`.

For a point `(xi, eta)`:

```text
g1 = x,xi
g2 = x,eta
a_ab = g_a . g_b
g^a = a^ab g_b
```

The physical frame is normally the element-center frame from `SIMO_FIXED_FRAME`.
For flat XY patch-test geometry the output/recovery frame is forced to global
`X-Y` so the F06 engineering forces are directly comparable with SAP/MSC/DKMQ
global reference values. For general 3D/warped geometry it remains the v1p6
center-fixed element frame.

## Membrane operator `Bm`

The membrane strain is formed as a tensor, not by simply differentiating in a
re-orthonormalized local frame at every Gauss point:

```text
eps_11 = u,1 . g1
eps_22 = u,2 . g2
eps_12 = 1/2 (u,1 . g2 + u,2 . g1)
```

Then the tensor is transformed to the fixed physical frame:

```text
eps_ij = eps_ab (E_i . g^a) (E_j . g^b)
```

For engineering shear output the third row is `2 eps_12`.

Fortran helper:

```text
T_STRAIN_FIXED_AT(...)
```

## Bending operator `Bb`

The current bending operator follows the `v1p6` director-gradient form:

```text
rho_ab = x,a . t,b
delta rho_ab = delta x,a . t0,b + x0,a . delta t,b
delta t_I = theta_I x t0_I
```

Two contributions are included:

1. rotational director curvature

```text
delta rho_11_rot = N_I,1 (t0_I x g1)
delta rho_22_rot = N_I,2 (t0_I x g2)
delta rho_12_rot = 1/2 [N_I,1 (t0_I x g2) + N_I,2 (t0_I x g1)]
```

2. translational curvature from nodal-normal variation

```text
t0,1 = sum_J N_J,1 t0_J
t0,2 = sum_J N_J,2 t0_J

delta rho_11_tr = N_I,1 t0,1
delta rho_22_tr = N_I,2 t0,2
delta rho_12_tr = 1/2 [N_I,1 t0,2 + N_I,2 t0,1]
```

Like membrane strain, bending curvature is transformed with the fixed-frame
contravariant tensor map.

This translational curvature term is zero for flat elements with constant nodal
normal, but it matters for warped/twisted geometry such as MacNeal twisted beam.

## Transverse shear `Bs_ANS`

The current shear is the important v1p6 correction. Shear is first sampled in
natural components at the ANS tying points:

```text
A = ( 0, -1)
B = ( 1,  0)
C = ( 0,  1)
D = (-1,  0)
```

At a sample point:

```text
gamma_a = u,a . t0 + g_a . (theta x t0)
```

Implemented nodally:

```text
gamma_1 <- N_I,1 t0_pt + N_I (t0_I x g1)
gamma_2 <- N_I,2 t0_pt + N_I (t0_I x g2)
```

Then ANS interpolation is done in natural components:

```text
gamma_1(xi,eta) = 1/2(1-eta) gamma_1(A) + 1/2(1+eta) gamma_1(C)
gamma_2(xi,eta) = 1/2(1-xi ) gamma_2(D) + 1/2(1+xi ) gamma_2(B)
```

Only after interpolation is shear transformed to physical components:

```text
gamma_i = gamma_a (E_i . g^a)
```

Fortran helper:

```text
T_SHEAR_FIXED_AT(...)
```

This avoids the older mistake of mixing Cartesian derivatives at the tying
points with unit local vectors. That older mix can look fine on square flat
elements but is not robust in distorted or twisted meshes.

## Membrane EAS

The membrane part uses a four-parameter Q1E4 enhanced strain mode:

```text
M_hat =
[ xi   0    0    0  ]
[ 0    eta  0    0  ]
[ 0    0    xi   eta]
```

The physical enhanced field is:

```text
M(xi,eta) = detJ0/detJ(xi,eta) * T0 * M_hat
```

where `T0` is the same fixed-frame contravariant tensor transform evaluated at
the element center. The stiffness is condensed:

```text
K_eff = K_uu - K_ua K_aa^-1 K_ua^T
```

For stress/force recovery, the effective membrane matrix must include the
condensed EAS contribution:

```text
B_eff = Bm - M K_aa^-1 K_ua^T
```

Without this recovery correction the solver displacement field can pass the
patch test while F06 engineering forces appear non-uniform.

## Drilling penalty

The current drilling strain is:

```text
B_drill,u = 1/2 (dN/dx1 E2 - dN/dx2 E1)
B_drill,theta = -N_I t0_I
```

The stiffness scale follows the Python default:

```text
k_drill = 0.02 * G * h
```

In Fortran this is implemented as:

```text
CDRILL = 2.0D-2 * SHELL_T(1,1) / (5.0D0/6.0D0)
```

because `SHELL_T(1,1)` contains the `5/6 * G * h` isotropic shear scale.

## Validation status

### Problem 2-001 thick patch test

Deck used:

```text
D:\18a\bending_only\Shell\gemini2\\working_mystran\prob_2_001_patch\prob_2_001_thick_my_cquadr_simo93.dat
```

Source deck:

```text
D:\18a\MYSTRAN_Validation-main\reference_msc\shell\prob_2_001_thick.dat
```

Variant:

```text
PARAM,QUADRTYP,SIMO
CQUADR cards
```

Last parsed MYSTRAN F06 result:

| Quantity | Result |
|---|---:|
| element count in subcase 2 | 5 |
| `Nxx,Nyy,Nxy` span | exact `[1.333333, 1.333333, 0.4]` |
| `Mxx,Myy,Mxy` span | exact `[1.111111e-07, 1.111111e-07, 3.333333e-08]` |
| max moment error vs exact | `0.000e+00` |
| max shear resultant | `3.670e-20` |

Interpretation:

```text
CQUADR/SIMO93 now passes the engineering-force patch output for problem 2-001.
```

### Problems 2-002 to 2-004

Existing comparison plots:

- `D:\18a\bending_only\Shell\gemini2\\prob_2_002_convergence.png`
- `D:\18a\bending_only\Shell\gemini2\\prob_2_003_convergence.png`
- `D:\18a\bending_only\Shell\gemini2\\prob_2_004_convergence.png`

The current documentation target is to compare:

- `MY-SIMO93` = `CQUADR + PARAM,QUADRTYP,SIMO`
- `MY-DKMQ24` = `CQUADR + PARAM,QUADRTYP,DKMQ24`
- `MY-MITCP+DHB` = `CQUADR + PARAM,QUADRTYP,MITC4PD`
- `MY-SIMO89` = `CQUAD4 + PARAM,QUAD4TYP,SIMO`
- `MY-MITC4P` = `CQUAD4 + PARAM,QUAD4TYP,MITC4+`
- `NASTRAN CQUAD4` = MSC Nastran deck with `CQUAD4`
- `NASTRAN CQUADR` = MSC Nastran deck mechanically converted from the same
  quad deck by replacing `CQUAD4` cards with `CQUADR`

Current regenerated plot files:

- `D:\18a\bending_only\Shell\gemini2\\value_convergence_2_002_thick.png`
- `D:\18a\bending_only\Shell\gemini2\\prob_2_003_convergence.png`
- `D:\18a\bending_only\Shell\gemini2\\macneal_twisted_beam_convergence.png`

MacNeal 2-004, N=24, after adding MSC Nastran `CQUADR`:

| Label | Uy | Uz |
|---|---:|---:|
| `MY-SIMO93` | `5.4056373e-03` | `1.7370647e-03` |
| `MY-DKMQ24` | `5.4019393e-03` | `1.7053023e-03` |
| `MY-MITC4P` | `5.4071923e-03` | `1.7076490e-03` |
| `NASTRAN CQUAD4` | `5.4075123e-03` | `1.7482317e-03` |
| `NASTRAN CQUADR` | `5.4028127e-03` | `1.7493023e-03` |

Problem 2-003 note: the current MSC Nastran F06 files for the mechanically
converted `CQUADR` decks are parsed and plotted, but the existing deck/output
combination gives `Uz = 0` for the out-of-plane subcase at the current target
node. Treat that curve as a deck/loading audit item before using it as a
physical benchmark.

## Open validation tasks

- Audit `problem_2_003` MSC Nastran out-of-plane loading/constraints before
  treating the Nastran `Uz=0` curve as a benchmark.
- Add final plot images directly into this note once the report layout is
  frozen.
