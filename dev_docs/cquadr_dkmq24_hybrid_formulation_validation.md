# CQUADR DKMQ24 hybrid formulation and validation notes

Date: July 31, 2026

This note documents the current MYSTRAN `CQUADR` DKMQ24-family implementation.
It is intentionally written as an implementation note, not as a claim that the
code is a literal line-by-line transcription of Katili 2015.

Main paths:

- default/standard `CQUADR`: `Source/EMG/EMG4/CQUADR_DKMQ24.f90`
- normal/no-SNORM variant: `Source/EMG/EMG4/CQUADR_DKMQ24N.f90`
- selector: `PARAM,QUADRTYP,DKMQ24` or `PARAM,QUADRTYP,DKMQ24N`
- dispatcher: `Source/EMG/EMG1/EMG.f90`

## Summary

The current MYSTRAN DKMQ24 path is best described as a DKMQ24/DKMQ-family
hybrid shell implementation with MYSTRAN-specific 6-DOF shell infrastructure.

It contains recognizable DKMQ ingredients:

- 4-node quadrilateral DKMQ-style bending/shear interpolation;
- edge-based auxiliary variables;
- shear interpolation through `NGAM`, `AG`, and `APHI`;
- bending correction through `BDEL * AINV_AU`.

But it is not the same as the original standalone Python
`DKMQ24_ShellElement_RHR.py` path. In the current branch it also uses:

- nodal-normal handling, optionally SNORM-driven;
- AU/ADELTA projection machinery;
- translational curvature terms from nodal normal variation;
- MYSTRAN 24-DOF shell assembly slots;
- MYSTRAN drilling stabilization;
- MYSTRAN element transform and F06 recovery conventions.

This is why `MY-DKMQ24` can match the MYSTRAN-like AU family while differing
from the older Python DKMQ24 prototype.

## Nodal normals and SNORM

`CQUADR_DKMQ24` computes nodal normals from adjacent element edges. If GRID
SNORM data are present, the standard DKMQ24 path can replace the geometric
normal with the supplied SNORM after a compatibility dot-product check.

In the explicit `DKMQ24N` variant, SNORM substitution is disabled:

```text
IF ((QUADRTYP /= 'DKMQ24N ') .AND. ALLOCATED(GRID_SNORM)) THEN ...
```

Practical interpretation:

- `DKMQ24` = standard current CQUADR path, SNORM-aware;
- `DKMQ24N` = normal/geometric path, intended as a cleaner comparison against
  Python normal-based formulations.

## AU/ADELTA projection

The element constructs:

```text
AU       = BUILD_AU(XYZ, NORMALS)
ADELTA   = BUILD_ADELTA(XYZ, thickness)
AINV_AU  = inv(ADELTA) * AU
```

For the standard CQUADR DKMQ24 path, `ADELTA` is diagonal and the inverse is
applied row-wise:

```text
AINV_AU(i,:) = AU(i,:) / ADELTA(i,i)
```

This AU projection is the main reason the current MYSTRAN DKMQ24 behavior is
closer to the AU/K6ROT DKMQ-family Python experiments than to the older
standalone DKMQ24 Python file.

## Bending operator

The bending matrix is assembled as:

```text
Bb = Bbeta + Bdel * AINV_AU
```

where:

- `Bbeta` carries the basic curvature contribution;
- `Bdel * AINV_AU` carries the projected/edge-based DKMQ correction;
- `NBC1/NBC2` add the curvature contribution from variation of the nodal
  normal field.

The code path uses current geometry data:

```text
GEOMETRY_AT(...)
BB_AT(XYZ, NORMALS, XI, ETA, TV1, TV2, CO, BCMAT, AINV_AU)
```

Important implementation clue:

```text
NBC1 = DN(1,:) * BCM(1,1) + DN(2,:) * BCM(2,1)
NBC2 = DN(1,:) * BCM(1,2) + DN(2,:) * BCM(2,2)
```

Then translational curvature-like bending terms enter through the local
tangent directions:

```text
BBETA(1,trans) = T1 * NBC1
BBETA(2,trans) = T2 * NBC2
BBETA(3,trans) = T1 * NBC2 + T2 * NBC1
```

This is one of the practical differences versus a simpler flat-plate DKMQ
prototype.

## Shear operator

The shear operator is assembled through:

```text
NGAM
AG
APHI
AINV_AU
```

In simplified notation:

```text
BSG  = CO^T * NGAM * AG
Bs   = BSG * APHI * AINV_AU
```

where:

- `NGAM` is the assumed shear interpolation matrix;
- `AG` carries signed edge-length scaling;
- `APHI` carries the thickness/shear-flexibility correction;
- `AINV_AU` maps auxiliary edge variables back to nodal DOFs.

This is a DKMQ-family shear projection, not the same as MITC4+ tying-point
shear.

## Drilling / sixth DOF

The element runs in MYSTRAN's 24-DOF shell slot. The sixth DOF is stabilized by
the `DRILL_STIFFNESS` path in the CQUADR DKMQ implementation. This is part of
the production MYSTRAN behavior and should not be equated with a pure 5-DOF
DKMQ plate element.

Practical comparison rule:

```text
MY-DKMQ24 should be compared to a MYSTRAN-like DKMQ24/AU Python implementation,
not the older 5-DOF or non-AU DKMQ24 prototype.
```

## Relationship to Katili 2015

The current implementation is DKMQ/DKMT-inspired and shares several recognizable
ingredients with the literature family. However, based on the actual code, it
should be documented as:

```text
DKMQ24-family hybrid implementation in MYSTRAN shell infrastructure
```

not:

```text
pure Katili 2015 DKMQ24 transcription
```

Reasons:

- the actual code includes MYSTRAN-specific nodal-normal/SNORM behavior;
- bending includes nodal-normal variation terms;
- shear uses the AU/ADELTA projection path currently shared conceptually with
  the DKMQ20-AU experiments;
- drilling stabilization and 24-DOF shell transform are MYSTRAN production
  choices;
- MacNeal results show the current path aligns with `DKMQ20-AU`/MYSTRAN-like
  behavior, not the older Python `DKMQ24_ShellElement_RHR.py`.

## Validation status

### Problem 2-001 thick patch test

Current parsed MYSTRAN F06 results for standard DKMQ24 and DKMQ24N:

| Variant | N span | M span | max M error | max Q |
|---|---|---|---:|---:|
| `MY-DKMQ24` | exact `[1.333333, 1.333333, 0.4]` | exact `[1.111111e-07, 1.111111e-07, 3.333333e-08]` | `0.000e+00` | `5.331e-23` |
| `MY-DKMQ24N` | exact `[1.333333, 1.333333, 0.4]` | exact `[1.111111e-07, 1.111111e-07, 3.333333e-08]` | `0.000e+00` | `2.757e-23` |

Interpretation:

```text
Both DKMQ24 and DKMQ24N pass problem 2-001 thick patch output.
```

### MacNeal 2-004 thick reference point

From `dev_docs/shell_quad_formulation_notes.md`, last available N=24 values:

| Label | Uy | Uz |
|---|---:|---:|
| `MY-DKMQ24` | `5.4019393e-03` | `1.7053023e-03` |
| `MY-DKMQ20` | `5.4019393e-03` | `1.7053023e-03` |
| Python `DKMQ20-AU` | `5.4019394e-03` | `1.7053021e-03` |
| Python original `DKMQ24` | `7.5919117e-03` | `2.3553501e-03` |
| `MY-MITC4P` | `5.4071923e-03` | `1.7076490e-03` |
| `NASTRAN CQUAD4` | `5.4075123e-03` | `1.7482317e-03` |
| `NASTRAN CQUADR` | `5.4028127e-03` | `1.7493023e-03` |

Interpretation:

```text
MY-DKMQ24 follows the MYSTRAN/AU family in MacNeal 2-004, not the older
standalone Python DKMQ24 response.
```

## Planned validation plots

Final documentation should include refreshed plots for:

- `problem_2_001`: patch force/resultant table;
- `problem_2_002`: convergence plot;
- `problem_2_003`: convergence plot;
- `problem_2_004`: MacNeal convergence plot.

Primary head-to-head labels:

- `MY-DKMQ24` = `CQUADR + PARAM,QUADRTYP,DKMQ24`
- `MY-SIMO93` = `CQUADR + PARAM,QUADRTYP,SIMO`
- `MY-MITC4P` = `CQUAD4 + PARAM,QUAD4TYP,MITC4+`

Existing plot files in the Python workspace:

- `D:\18a\bending_only\Shell\gemini2\shit\value_convergence_2_002_thick.png`
- `D:\18a\bending_only\Shell\gemini2\shit\prob_2_003_convergence.png`
- `D:\18a\bending_only\Shell\gemini2\shit\macneal_twisted_beam_convergence.png`

These were regenerated after adding `NASTRAN CQUADR` to the solver-series
loader and plot whitelists. Problem 2-003 still needs a deck/loading audit for
the out-of-plane Nastran subcase because the current F06 target displacement is
parsed as zero for both CQUAD4 and CQUADR.
