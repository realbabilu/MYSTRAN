# CQUAD4 DKMQ20 AU/K6ROT formulation and validation notes

Date: July 31, 2026

This note documents the current MYSTRAN `CQUAD4` DKMQ20 branch:

- selector: `PARAM,QUAD4TYP,DKMQ20`
- dispatcher: `Source/EMG/EMG1/EMG.f90`
- kernel: `Source/EMG/EMG4/CQUAD4_DKMQ20_RHR.f90`
- interface: `Source/Interfaces/CQUAD4_DKMQ20_RHR_Interface.f90`
- closest Python reference: `DKMQ20_ShellElement_RHR_6dof_AU_K6ROT.py`

The branch should not be described as the old 5-DOF/non-AU DKMQ20 prototype.
It is a MYSTRAN 24-DOF shell-slot implementation with an AU projection and
MYSTRAN K6ROT drilling stabilization.

## Summary

`CQUAD4_DKMQ20_RHR` is a 4-node quadrilateral shell branch for the `CQUAD4`
family.  Its physical formulation follows the DKMQ/DKMT-style assumed
bending/shear projection, but it is wired into MYSTRAN as a 6-DOF-per-node shell
element:

```text
u_I = [ux, uy, uz, rx, ry, rz]
```

The sixth DOF is not a physical DKMQ20 drilling strain mode.  It is stabilized
through MYSTRAN's K6ROT-compatible shell infrastructure so the element can live
inside MYSTRAN's standard shell assembly and constraint machinery.

The main implementation ingredients are:

- AU/ADELTA edge projection;
- bending operator with translational curvature from nodal-normal variation;
- shear operator through the `NGAM * AG * APHI * AINV_AU` path;
- 24x24 stiffness slot;
- translational lumped mass/pressure support consistent with the current
  MYSTRAN shell branch conventions;
- K6ROT-style handling for the drilling DOF.

## Difference from the old DKMQ20 5-DOF prototype

The old/simple Python `DKMQ20_ShellElement_RHR.py` is a useful DKMQ20 baseline,
but it is not the current MYSTRAN target.  Current MYSTRAN DKMQ20 is closer to:

```text
DKMQ20_ShellElement_RHR_6dof_AU_K6ROT.py
```

Important differences:

| Topic | Old DKMQ20 5-DOF prototype | Current MYSTRAN CQUAD4 DKMQ20 |
|---|---|---|
| Solver slot | physical 5 DOF/node | MYSTRAN 6 DOF/node shell slot |
| Drilling DOF | absent/fictitious outside the element | K6ROT-stabilized sixth DOF |
| Bending correction | simpler DKMQ-style bending path | AU/ADELTA projected bending correction |
| Nodal normals | simpler/local | nodal-normal based shell geometry |
| Curvature | no full normal-variation contribution | includes translational curvature from nodal-normal variation |
| Shear | older 5-DOF shear expression | `NGAM * AG * APHI * AINV_AU` projection |

Therefore, when comparing MYSTRAN results, the correct Python-side comparison
is the AU/K6ROT version, not the old non-AU file.

## AU/ADELTA projection

The current branch constructs edge auxiliary variables through:

```text
AU       = BUILD_AU(XYZ, NORMALS)
ADELTA   = BUILD_ADELTA(XYZ, thickness)
AINV_AU  = inv(ADELTA) * AU
```

For the current implementation `ADELTA` is diagonal, so the inverse is applied
row-wise:

```text
AINV_AU(i,:) = AU(i,:) / ADELTA(i,i)
```

`AU` maps nodal translations and rotations to edge variables.  The rotational
terms are built from the nodal normal/director system; the translational terms
come from edge-normal variation.

This AU layer is the main reason MYSTRAN `CQUAD4 DKMQ20` matches the
`DKMQ20-AU` Python experiment very closely on MacNeal 2-004.

## Bending operator

The bending matrix is assembled in the same structural pattern used by the
current DKMQ-family MYSTRAN ports:

```text
Bb = Bbeta + Bdel * AINV_AU
```

where:

- `Bbeta` carries the basic curvature terms;
- `Bdel * AINV_AU` carries the projected DKMQ correction;
- nodal-normal variation introduces translational curvature terms.

The translational curvature contribution is zero for flat elements with
constant normals, but it matters on warped/twisted geometries such as MacNeal
twisted beam.

## Shear operator

The transverse shear part is not an MITC tying-point shear.  It is the
DKMQ-family assumed shear projection:

```text
BSG = CO^T * NGAM * AG
Bs  = BSG * APHI * AINV_AU
```

where:

- `NGAM` is the assumed shear interpolation matrix;
- `AG` carries edge/sign scaling;
- `APHI` carries the thickness/shear-flexibility correction;
- `AINV_AU` maps the auxiliary edge variables back into nodal shell DOFs.

This path is powerful for some bending-dominated tests but can be slow to
converge in the current branch for particular shell verification problems.

## Drilling / K6ROT

The physical DKMQ20 formulation does not supply an independent robust drilling
strain.  MYSTRAN still needs a sixth nodal DOF for shell compatibility, so this
branch relies on K6ROT-style stabilization.

Practical interpretation:

```text
CQUAD4 DKMQ20 = DKMQ20 AU shell kernel + MYSTRAN 6-DOF shell slot + K6ROT drilling stabilization
```

It should not be documented as a Hughes-Brezzi drilling shell.  The
Hughes-Brezzi drilling experiments are separate branches.

## Validation status

The local comparison scripts used for the current validation are:

```text
D:\18a\bending_only\Shell\gemini2\shit\problem_2_001_mystran.py
D:\18a\bending_only\Shell\gemini2\shit\problem_2_002_mystran.py
D:\18a\bending_only\Shell\gemini2\shit\problem_2_003_mystran.py
D:\18a\bending_only\Shell\gemini2\shit\problem_2_004_mystran.py
```

SAP2000 thick reference values were taken from:

```text
C:\Program Files\Computers and Structures\SAP2000 20\Manuals\Verification\Analysis\Shells
```

### Problem 2-001 thick patch test

Current parsed F06 result:

| Quantity | MY-CQUAD4 DKMQ20 |
|---|---:|
| `Nxx` | `1.333333e+00` |
| `Nyy` | `1.333333e+00` |
| `Nxy` | `4.000000e-01` |
| `Mxx` | `1.111111e-07` |
| `Myy` | `1.111111e-07` |
| `Mxy` | `3.333333e-08` |
| `max |Q|` | `0.000e+00` |

Interpretation:

```text
CQUAD4 DKMQ20 passes the thick patch test in current MYSTRAN engineering-force output.
```

### Problem 2-002

SAP2000 thick reference:

| Load case | Reference |
|---:|---:|
| LC1 UX | `3.000e-05` |
| LC2 UZ | `1.072e-01` |
| LC3 UY | `4.321e-01` |
| LC4 abs(UY) | `2.240e-03` |
| LC5 abs(UX) | `8.990e-04` |
| LC6 RZ | `3.600e-02` |

Current last mesh (`N=24`) result:

| Load case | MY-CQUAD4 DKMQ20 | Relative error |
|---:|---:|---:|
| LC1 UX | `3.000000e-05` | `0.0%` |
| LC2 UZ | `6.358497e-02` | `40.7%` |
| LC3 UY | `4.320232e-01` | `0.018%` |
| LC4 abs(UY) | `7.723834e-03` | `244.8%` |
| LC5 abs(UX) | `5.294545e-04` | `41.1%` |
| LC6 RZ | `3.600000e-02` | `0.0%` |

Interpretation:

- LC1, LC3, and LC6 are excellent.
- LC2 and LC5 remain too stiff/low.
- LC4 converges very slowly and remains far from the SAP2000 reference at N=24.

In the current F06 set, `MY-CQUAD4 DKMQ20` and `MY-CQUADR DKMQ24` are
numerically identical for all Problem 2-002 load cases.  That does not prove two
independent formulations; it means both current branches share the same
effective AU/DKMQ-family behavior in this benchmark.

### Problem 2-003

SAP2000 thick reference:

```text
Uy = 0.0773
Uz = 0.4298
```

Current last mesh (`N=24`) result:

| Quantity | MY-CQUAD4 DKMQ20 | Relative error |
|---|---:|---:|
| `Uy` | `4.833892e-02` | `37.5%` |
| `Uz` | `7.520937e-01` | `75.0%` |

The `Uz` curve is the clearest slow-convergence symptom.  Coarse meshes are far
too flexible, and the value decreases only gradually as the mesh is refined:

| N | Uz |
|---:|---:|
| 2 | `6.034431e+00` |
| 3 | `7.024714e+00` |
| 4 | `6.462945e+00` |
| 6 | `4.298257e+00` |
| 8 | `2.838409e+00` |
| 12 | `1.569472e+00` |
| 16 | `1.094490e+00` |
| 24 | `7.520937e-01` |

In the same benchmark, CQUADR SIMO1993 and the DKMQ24N research branch stay
much closer to the SAP2000 `Uz` reference, while DKMQ20/DKMQ24 remain late.

### Problem 2-004 MacNeal twisted beam

SAP2000 thick reference:

```text
Uy = 0.005402
Uz = 0.001760
```

Current last mesh (`N=24`) result:

| Quantity | MY-CQUAD4 DKMQ20 | Relative error |
|---|---:|---:|
| `Uy` | `5.401939e-03` | `0.001%` |
| `Uz` | `1.705302e-03` | `3.11%` |

Problem 2-004 is a strong result for the AU/K6ROT DKMQ20 branch.  It also
matches the Python AU/K6ROT experiment to F06/parser rounding:

```text
Python DKMQ20-AU N=24 : Uy=5.4019394e-03, Uz=1.7053021e-03
MYSTRAN DKMQ20 N=24  : Uy=5.4019393e-03, Uz=1.7053023e-03
```

## Current interpretation

`CQUAD4 DKMQ20` is a useful and validated branch for:

- patch-test behavior;
- MacNeal twisted beam thick response;
- comparisons against the Python DKMQ20 AU/K6ROT implementation.

However, it is not uniformly strong across the SAP2000 Problem 2 set:

- it is late/poor in Problem 2-003 `Uz`;
- it is weak in Problem 2-002 LC2/LC4/LC5;
- it currently overlays `CQUADR DKMQ24` in Problems 2-002 and 2-003.

Recommended label in reports:

```text
MY-DKMQ20 = CQUAD4 + PARAM,QUAD4TYP,DKMQ20 = DKMQ20 AU/K6ROT branch
```

