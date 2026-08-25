# CQUADR DKMQ24N beta formulation and validation notes

Date: July 31, 2026

This note documents the current MYSTRAN `CQUADR` DKMQ24N branch:

- selector: `PARAM,QUADRTYP,DKMQ24N`
- dispatcher: `Source/EMG/EMG1/EMG.f90`
- kernel: `Source/EMG/EMG4/CQUADR_DKMQ24N.f90`
- interface: `Source/Interfaces/CQUADR_DKMQ24N_Interface.f90`

`DKMQ24N` is intentionally documented as a beta/research comparison branch.  It
is useful because it isolates the geometric-normal path from the standard
SNORM-aware `DKMQ24` branch, but the current validation results show that it is
not yet a production replacement for standard `CQUADR DKMQ24`.

## Purpose

The standard `CQUADR DKMQ24` branch is SNORM-aware.  If GRID SNORM data are
available, the element can substitute supplied nodal normals after a
compatibility check.

`DKMQ24N` was added to provide a cleaner comparison path:

```text
CQUADR DKMQ24N = DKMQ24-family AU shell branch using geometric normals only
```

The intended use is formulation research and debugging:

- compare SNORM-aware and geometric-normal behavior;
- isolate whether a result is caused by supplied nodal normals;
- compare MYSTRAN CQUADR behavior against Python normal-based DKMQ-family
  experiments;
- test whether geometric normals improve warped/cylindrical shell convergence.

## Relationship to DKMQ24 and DKMQ20

`DKMQ24N` shares the same broad DKMQ-family AU formulation style as the current
`DKMQ24` branch:

```text
AU       = BUILD_AU(XYZ, NORMALS)
ADELTA   = BUILD_ADELTA(XYZ, thickness)
AINV_AU  = inv(ADELTA) * AU
```

The important difference is how `NORMALS` are obtained and used.  In the
standard branch:

```text
DKMQ24 = SNORM-aware when GRID SNORM data are present
```

In the beta branch:

```text
DKMQ24N = geometric-normal/no-SNORM path
```

The `N` suffix should be read as "normal/geometric normal" rather than as a new
published element name.

Do not mix DKMQ24N F06 files into `MY-CQUADR DKMQ24` plots.  Because filenames
contain the substring `cquadr_dkmq24`, simple substring matching can accidentally
label `cquadr_dkmq24n` as `DKMQ24`.  This caused earlier misleading vertical or
zig-zag plot segments.  Validation scripts should check `cquadr_dkmq24n` before
`cquadr_dkmq24`, or filter DKMQ24N explicitly when the plot is intended to show
standard DKMQ24 only.

## Formulation skeleton

The branch follows the DKMQ-family implementation pattern already documented
for `CQUADR_DKMQ24`:

```text
Bb = Bbeta + Bdel * AINV_AU
```

and:

```text
BSG = CO^T * NGAM * AG
Bs  = BSG * APHI * AINV_AU
```

where:

- `AU` maps nodal shell DOFs to edge auxiliary variables;
- `ADELTA` scales those edge variables;
- `AINV_AU` maps projected edge quantities back to nodal DOFs;
- `Bbeta` carries basic bending curvature;
- `Bdel * AINV_AU` carries the DKMQ correction;
- `NGAM`, `AG`, and `APHI` build the assumed transverse shear projection.

As with the standard DKMQ24 branch, this is a MYSTRAN 24-DOF shell-slot element,
not a pure standalone 5-DOF plate prototype.

## Validation status

The focused comparison script used for the latest DKMQ24N check is:

```text
D:\18a\bending_only\Shell\gemini2\shit\problem_2_002_003_cquadr_focus.py
```

It compares:

```text
MY-CQUADR DKMQ24  = CQUADR + PARAM,QUADRTYP,DKMQ24
MY-CQUADR DKMQ24N = CQUADR + PARAM,QUADRTYP,DKMQ24N
MY-CQUADR SIMO1993 = CQUADR + PARAM,QUADRTYP,SIMO
```

against SAP2000 thick reference lines.

### Problem 2-001 thick patch test

Current parsed F06 status:

| Variant | Patch status |
|---|---|
| `MY-CQUADR DKMQ24` | pass |
| `MY-CQUADR DKMQ24N` | pass |

For DKMQ24N the current patch force/resultant output is exact within the
printed precision:

```text
Nxx = 1.333333
Nyy = 1.333333
Nxy = 0.4
Mxx = 1.111111e-07
Myy = 1.111111e-07
Mxy = 3.333333e-08
```

Patch pass is necessary, but not sufficient for production use.  The
convergence tests below expose beta-level issues.

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

Current last mesh (`N=24`) DKMQ24N result:

| Load case | MY-CQUADR DKMQ24N | Relative error |
|---:|---:|---:|
| LC1 UX | `3.000000e-05` | `0.0%` |
| LC2 UZ | `6.358497e-02` | `40.7%` |
| LC3 UY | `-2.722751e+09` | unstable/singular |
| LC4 abs(UY) | `7.494211e-04` | `66.5%` |
| LC5 abs(UX) | `5.294545e-04` | `41.1%` |
| LC6 RZ | `0.000000e+00` | `100.0%` |

Interpretation:

```text
DKMQ24N is not acceptable on Problem 2-002 in the current implementation.
```

The LC3 value is a clear instability/singularity symptom, not a normal
convergence miss.  LC6 also fails completely in the current run.  Therefore
DKMQ24N should stay beta/research despite its patch-test success.

### Problem 2-003

SAP2000 thick reference:

```text
Uy = 0.0773
Uz = 0.4298
```

Current last mesh (`N=24`) DKMQ24N result:

| Quantity | MY-CQUADR DKMQ24N | Relative error |
|---|---:|---:|
| `Uy` | `4.833892e-02` | `37.5%` |
| `Uz` | `4.697210e-01` | `9.3%` |

This is the main reason DKMQ24N remains interesting.  On Problem 2-003
out-of-plane response it is much better than standard DKMQ24:

| Variant | Uz at N=24 | Relative error vs SAP `0.4298` |
|---|---:|---:|
| `MY-CQUADR DKMQ24` | `7.520937e-01` | `75.0%` |
| `MY-CQUADR DKMQ24N` | `4.697210e-01` | `9.3%` |
| `MY-CQUADR SIMO1993` | `4.893872e-01` | `13.9%` |

However, the `Uy` response remains low and essentially follows the same poor
in-plane trend as DKMQ24/DKMQ20 in this problem.  DKMQ24N improves the
out-of-plane behavior, not the whole problem.

### Problem 2-004 MacNeal twisted beam

Older/partial DKMQ24N F06 files existed for only a subset of meshes before the
latest focused run.  They should not be used in standard DKMQ24 plots.  When
DKMQ24N is shown, it must be explicitly labeled as beta.

The current no-MIN production-style comparison should normally show:

```text
MY-CQUAD4 DKMQ20
MY-CQUAD4 SIMO1989
MY-CQUADR DKMQ24
MY-CQUADR SIMO1993
MY-CQUADR MITC4+DHB
NASTRAN CQUAD4/CQUADR
```

and exclude DKMQ24N unless the purpose is a beta/research comparison.

## Current recommendation

Recommended report label:

```text
MY-DKMQ24N = CQUADR + PARAM,QUADRTYP,DKMQ24N = geometric-normal DKMQ24 beta branch
```

Recommended status:

```text
beta / research only
```

Reasons:

- passes Problem 2-001 patch test;
- improves Problem 2-003 `Uz` significantly compared with standard DKMQ24;
- but fails Problem 2-002 LC3/LC6 badly in the current run;
- should not be mixed into standard `MY-DKMQ24` plots.

Practical use:

- keep it for diagnosing normal/SNORM sensitivity;
- keep it for studying why standard DKMQ24 is late in Problem 2-003 `Uz`;
- do not promote it as the default CQUADR formulation until Problem 2-002 is
  fixed.

