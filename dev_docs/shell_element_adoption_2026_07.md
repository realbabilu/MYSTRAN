# Shell element adoption notes, July 2026

This note records the current MYSTRAN shell element additions that should be
kept together as the next element-adoption commit.  It intentionally separates
production candidates from research/beta branches so later validation plots do
not hide formulation changes under similar names.

## Selector map

| MYSTRAN card family | PARAM selector | Kernel | Status |
|---|---|---|---|
| `CQUAD4` | `PARAM,QUAD4TYP,MIN4` | existing MIN4 path | baseline |
| `CQUAD4` | `PARAM,QUAD4TYP,MIN4T` | existing MIN4T path | baseline |
| `CQUAD4` | `PARAM,QUAD4TYP,MITC4` | existing MITC4 path | baseline |
| `CQUAD4` | `PARAM,QUAD4TYP,MITC4+` | existing MITC4+ path | baseline / strong MacNeal result |
| `CQUAD4` | `PARAM,QUAD4TYP,DKMQ20` | `CQUAD4_DKMQ20_RHR.f90` | new DKMQ20 AU/K6ROT branch |
| `CQUAD4` | `PARAM,QUAD4TYP,SIMO` | `CQUAD4_SIMO1989.f90` | new Simo1989-style 5-DOF shell with MYSTRAN K6ROT |
| `CQUADR` | `PARAM,QUADRTYP,DKMQ24` | `CQUADR_DKMQ24.f90` | standard DKMQ24-family branch |
| `CQUADR` | `PARAM,QUADRTYP,DKMQ24N` | `CQUADR_DKMQ24N.f90` | research beta, geometric-normal/no-SNORM comparison path |
| `CQUADR` | `PARAM,QUADRTYP,SIMO` | `CQUADR_SIMO1993.f90` | new Simo1993 branch, based on Python `Simo1993_ShellElement_v1p6.py` |
| `CQUADR` | `PARAM,QUADRTYP,MITC4PD` | `CQUADR_MITC4PHB.f90`, `CQUADR_MITC4PHB_B.f90` | beta; still not a patch-test-clean production recommendation |
| `CTRIA3` | `PARAM,TRIA3TYP,MITC3+` | `TPLT_MITC3P.f90` and MITC helpers | MITC3+ triangle branch using MYSTRAN shell 6-DOF infrastructure/K6ROT |
| `CTRIA3` | `PARAM,TRIA3TYP,MIN3` | existing MIN3 path | baseline triangle branch |
| `CTRIAR` | default DKMT18 path | existing DKMT18 path | baseline/recommended triangular `CTRIAR` branch |

## Formulation summary

### CQUAD4 DKMQ20 AU/K6ROT

The new `CQUAD4_DKMQ20_RHR` branch is not the old 5-DOF/non-AU Python DKMQ20.
It is the AU-compatible 24-DOF MYSTRAN shell-slot implementation:

```text
AU       = BUILD_AU(XYZ, NORMALS)
ADELTA   = BUILD_ADELTA(XYZ, thickness)
AINV_AU  = inv(ADELTA) * AU
```

The bending operator includes the translational curvature contribution from
nodal-normal variation.  Transverse shear uses the DKMQ-family
`NGAM * AG * APHI * AINV_AU` projection.  The sixth/drilling DOF is stabilized
through MYSTRAN K6ROT-style shell infrastructure rather than a Hughes-Brezzi
drilling field inside the DKMQ20 kernel.

### CQUAD4 SIMO1989

`CQUAD4_SIMO1989` is a quadrilateral Simo-style branch for the `CQUAD4` family.
It remains a 5-DOF physical shell formulation at the element level, with the
sixth DOF supplied by MYSTRAN K6ROT stabilization.  It should be compared
against the Simo1989/Simo-style Python shell family, not against the CQUADR
Simo1993 kernel.

### CQUADR DKMQ24 and DKMQ24N

`CQUADR_DKMQ24` is the standard DKMQ24-family MYSTRAN branch.  In the current
implementation it is best described as a DKMQ/DKMQ24 hybrid inside MYSTRAN's
6-DOF shell infrastructure, not as a literal pure Katili-paper transcription.
It is SNORM-aware when GRID SNORM data are supplied.

`CQUADR_DKMQ24N` is retained as a research/beta comparison path using geometric
normals only.  It should not be mixed into the `DKMQ24` validation label.  Some
old plot scripts accidentally matched `cquadr_dkmq24n` as `cquadr_dkmq24`; that
has to be filtered in validation tooling.

### CQUADR SIMO1993

`CQUADR_SIMO1993` is isolated from both DKMQ24 and CQUAD4 MITC4+.  It follows
the corrected Python `Simo1993_ShellElement_v1p6.py` structure:

- fixed-frame membrane tensor transformation;
- bending curvature with rotational director curvature plus translational
  curvature from nodal normal variation;
- Simo/Dvorkin-Bathe style ANS shear in natural components followed by
  physical-frame transformation;
- Q1E4 membrane EAS with static condensation;
- small drilling stabilization in the CQUADR Simo branch.

This branch passes the thick shell patch deck and performs well on the MacNeal
twisted-beam benchmark.

### CQUADR MITC4+D Hughes-Brezzi

`CQUADR_MITC4PHB` and its isolated B-matrix helper
`CQUADR_MITC4PHB_B` are retained as beta/research code.  The branch is useful
for comparing MITC4+ plus drilling/Hughes-Brezzi ideas in the CQUADR slot, but
it is not currently patch-test clean and should not be described as production
validated.

### Triangles: CTRIA3 MITC3+/MIN3 and CTRIAR DKMT18

`CTRIA3` has MIN3 and MITC3+ alternatives.  MITC3+ uses the MYSTRAN 6-DOF shell
slot and K6ROT-style stabilization; it should not be described as a
Hughes-Brezzi drilling formulation unless that branch is explicitly selected in
the Python-only experiments.

`CTRIAR` remains the DKMT18-family branch.  Python experiments for
`DKMT18`, `DKMT18_HughesBrezzi`, and `DKMT18MaknunSNORM` all pass the
triangular patch-test script, but the MYSTRAN production selector remains
`CTRIAR`/DKMT18.

## Quad validation summary: SAP2000 Problem 2-001 through 2-004, thick

Reference values were taken from the SAP2000 v20 verification PDFs under:

```text
[SAP2000 20\Manuals\Verification\Analysis\Shells]
https://docs.csiamerica.com/manuals/sap2000/Verification/Analysis/Shells/Problem%202-001.pdf
to
https://docs.csiamerica.com/manuals/sap2000/Verification/Analysis/Shells/Problem%202-017.pdf
```

The local comparison scripts are:

```text
problem_2_001_mystran.py
problem_2_002_mystran.py
problem_2_003_mystran.py
problem_2_004_mystran.py
```

### Problem 2-001 patch test

Patch target:

```text
Nxx = 1.333333333
Nyy = 1.333333333
Nxy = 0.4
Mxx = 1.111e-7
Myy = 1.111e-7
Mxy = 3.333e-8
```

Observed MYSTRAN status from the current generated F06 files:

| Element branch | Patch status |
|---|---|
| `MY-CQUAD4 DKMQ20` | pass |
| `MY-CQUAD4 SIMO1989` | pass |
| `MY-CQUADR DKMQ24` | pass |
| `MY-CQUADR SIMO1993` | pass |
| `MY-CQUAD4 MITC4+` / `MITC4` | not patch-perfect in current force recovery |
| `MY-CQUADR MITC4+DHB` | beta; not patch-test-clean |
| `MIN4`, `MIN4T` | baseline but not patch-perfect here |

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

Last mesh (`N=24`) from current MYSTRAN F06 parsing:

| Branch | LC1 | LC2 | LC3 | LC4 | LC5 | LC6 |
|---|---:|---:|---:|---:|---:|---:|
| `MY-CQUAD4 DKMQ20` | `3.000e-05` | `6.358e-02` | `4.320e-01` | `7.724e-03` | `5.295e-04` | `3.600e-02` |
| `MY-CQUADR DKMQ24` | `3.000e-05` | `6.358e-02` | `4.320e-01` | `7.724e-03` | `5.295e-04` | `3.600e-02` |
| `MY-CQUAD4 SIMO1989` | `3.000e-05` | `1.080e-01` | `4.319e-01` | `3.040e-03` | `9.000e-04` | `3.600e-02` |
| `MY-CQUADR SIMO1993` | `3.000e-05` | `1.077e-01` | `4.319e-01` | `3.040e-03` | `8.974e-04` | `3.600e-02` |
| `MY-CQUADR MITC4+DHB` | `3.000e-05` | `6.358e-02` | `4.319e-01` | `3.040e-03` | `5.295e-04` | `3.600e-02` |

Important note: in Problem 2-002 the current `MY-CQUADR DKMQ24` result is
numerically identical to `MY-CQUAD4 DKMQ20` for all six load cases.  Treat this
as a formulation equivalence/issue to revisit, not as independent confirmation
of two different kernels.

### Problem 2-003

SAP2000 thick reference:

```text
Uy = 0.0773
Uz = 0.4298
```

Last mesh (`N=24`) from current MYSTRAN F06 parsing:

| Branch | Uy | Uz |
|---|---:|---:|
| `MY-CQUAD4 MITC4+` | `4.834e-02` | `4.894e-01` |
| `MY-CQUAD4 DKMQ20` | `4.834e-02` | `7.521e-01` |
| `MY-CQUAD4 SIMO1989` | `8.841e-02` | `4.894e-01` |
| `MY-CQUADR DKMQ24` | `4.834e-02` | `7.521e-01` |
| `MY-CQUADR SIMO1993` | `8.808e-02` | `4.894e-01` |
| `MY-CQUADR MITC4+DHB` | `4.834e-02` | `4.894e-01` |

Problem 2-003 exposes the same issue as 2-002: `MY-CQUADR DKMQ24` and
`MY-CQUAD4 DKMQ20` are identical.  The `DKMQ24N` beta files must not be mixed
into the `DKMQ24` plot label; doing so creates misleading zig-zag or vertical
segments.

### Problem 2-004 MacNeal twisted beam

SAP2000 thick reference:

```text
Uy = 0.005402
Uz = 0.001760
```

Last mesh (`N=24`) from current MYSTRAN F06 parsing:

| Branch | Uy | Uz | Notes |
|---|---:|---:|---|
| `MY-CQUAD4 MITC4+` | `5.407e-03` | `1.708e-03` | close to SAP/NASTRAN |
| `MY-CQUAD4 MITC4` | `5.406e-03` | `1.706e-03` | close |
| `MY-CQUAD4 DKMQ20` | `5.402e-03` | `1.705e-03` | close; AU/K6ROT branch |
| `MY-CQUAD4 SIMO1989` | `5.413e-03` | `1.748e-03` | close |
| `MY-CQUADR DKMQ24` | `5.402e-03` | `1.705e-03` | close at N=24; only clean DKMQ24 points should be plotted |
| `MY-CQUADR SIMO1993` | `5.413e-03` | `1.748e-03` | close |
| `MY-CQUADR MITC4+DHB` | `3.401e-03` | `1.241e-03` | too stiff; beta |

The clean no-MIN plot should exclude `MIN4`/`MIN4T` because their MacNeal
values are extreme relative to the modern shell branches.  Also exclude
`DKMQ24N` unless explicitly making a beta/research plot.

## Triangle validation snapshot

Python triangular scripts run successfully for Problem 2-001 through 2-004:

```text
D:\18a\bending_only\Shell\gemini2\\prob_2_001_6dof_triangle.py
D:\18a\bending_only\Shell\gemini2\\prob_2_002_6dof_Triangle.py
D:\18a\bending_only\Shell\gemini2\\prob_2_003_6dof_Triangle.py
D:\18a\bending_only\Shell\gemini2\\prob_2_004_6dof_Triangle.py
```

Problem 2-001 triangular patch run status:

| Python triangle branch | Patch status |
|---|---|
| `DKMT18HB` | pass |
| `DKMT18` | pass |
| `DKMT18MaknunSNORM` | pass |
| `MITC3+_6D
MY-DKMT18
MY-MIN3
NASTRAN CTRIA3
```

Representative Problem 2-004 last-mesh MYSTRAN triangular values:

| Branch | Uy at N=24 | Uz at N=24 |
|---|---:|---:|
| `MY-MITC3+` | `6.562e-03` | `1.871e-03` |
| `MY-MIN3` | `7.121e-03` | `2.028e-03` |
| `MY-DKMT18` | `3.111e-03` | `1.313e-03` |

The triangle plots are useful for exploratory comparison, but their current
solver overlay parser still needs the same duplicate-file hygiene that was
added to the quad MYSTRAN comparison scripts.  Do not over-interpret duplicate
points in older triangle plots.

## Commit hygiene

Recommended staged set for this adoption commit:

- source kernels and interfaces for DKMQ20, Simo1989, DKMQ24N, Simo1993,
  MITC4PHB, and MITC4PHB_B;
- dispatcher/parameter/interface updates for `QUAD4TYP`, `QUADRTYP`, and
  triangle selectors;
- documentation in `dev_docs/` and `update/`;
- exclude transient `.bak`, `_smoke`, generated F06, and experimental plot
  artifacts unless a release note explicitly needs them.` | pass |
| `MITC3+D_HB` | pass |

The convergence scripts include MYSTRAN/NASTRAN overlays where matching F06
files exist.  Current MYSTRAN triangular labels of interest are:

```text
MY-MITC3+

