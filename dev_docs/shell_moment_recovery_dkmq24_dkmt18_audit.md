# DKMQ24 and DKMT18 Shell Moment-Recovery Audit

## Scope

This note isolates shell engineering-force recovery from shell stiffness.
It covers:

- MYSTRAN `CQUADR + DKMQ24`
- MYSTRAN `CQUAD4 + DKMQ20`
- MYSTRAN `CTRIAR + DKMT18`
- the corresponding Python DKMQ24 and DKMT18 implementations in
  `D:/18a/bending_only/Shell/gemini2`

The Simo and legacy MIN/MITC families are not redirected through these new
recovery paths. They retain their existing behavior.

## Main finding

DKMQ24 and DKMT18 do not have the same output-recovery contract.

| Family | MYSTRAN stress points | Kernel recovery | LINK9 behavior |
|---|---:|---|---|
| CQUADR/DKMQ24 | 5 | center plus four direct natural-corner operators for `Mxy`; other components remain Gauss samples | polynomial fit for normal components, direct overwrite for `Mxy` |
| CQUAD4/DKMQ20 | 5 | same isolated direct `Mxy` strategy | same isolated override |
| CTRIAR/DKMT18 | 1 | one operator at triangle centroid `(1/3,1/3)` | no corner polynomial extrapolation |

Therefore the DKMQ24 fix must not be copied mechanically to DKMT18.
`MODEL_STUF.NUM_SEi` is `5` for `QUADR` but only `1` for `TRIA3`, and the
current DKMT18 kernel fills only `BE1/BE2/BE3(:,:,1)`.

## August 2026 Python/Fortran parity audit

The Python port reproduces the CQUADR/DKMQ24 displacement essentially
exactly on problem 2-006:

```text
Python CQUADR/DKMQ24 Uz(49)  = -0.3233629
Fortran CQUADR/DKMQ24 Uz(49) = -0.3233629
```

The remaining discrepancy is therefore output recovery, not assembled shell
stiffness.

Two defects were found in the current Fortran recovery path:

1. `POLYNOM_FIT_STRE_STRN` only enters its quadrilateral branch for
   `TYPE(1:5) == 'QUAD4'` or `TYPE(1:5) == 'QUAD8'`.  `TYPE == 'QUADR   '`
   returns without fitting.  Consequently, the present CQUADR `Mxx/Myy`
   values are still raw Gauss samples placed in corner output slots, despite
   the intended policy in `OFP3_ELFE_2D`.
2. `CQUADR_DKMQ24` fills output slots 2:5 in natural-coordinate order
   `(+,+), (+,-), (-,+), (-,-)`.  `AGRID(1:4)` remains in element connectivity
   order, whose standard natural coordinates are `(-,-), (+,-), (+,+),
   (-,+)`.  Thus direct `Mxy` recovery is not attached to the corresponding
   GRID consistently.

The Python diagnostic path now implements the intended isolated policy:

- `Mxx/Myy`: default `Q4SURFIT=3`, i.e. least-squares `a + b*x + c*y` in
  actual element coordinates;
- `Mxy`: direct evaluation at the requested center/corner;
- scalar components are averaged as scalar MYSTRAN output values, rather than
  silently tensor-averaged in a different coordinate system.

Using a four-term natural-coordinate interpolation (`1,r,s,r*s`) was tested
and rejected because it does not reproduce MYSTRAN's default `Q4SURFIT=3`.

### Required isolated Fortran correction

Do not route all `QUADR` families through one generic change.  For
`QUADR/DKMQ24` only:

1. allow the DKMQ24 normal-moment rows to reach a `QUADR`-capable surface fit;
2. map direct corner `Mxy` to standard node order;
3. leave MIN, MITC, DKMT and Simo recovery unchanged;
4. rerun problem 2-001 (constant curvature) and problem 2-006 (curved shell)
   before accepting the patch.

### Simo decision

No Simo recovery change is justified by this audit.  On problem 2-006 the
existing Python and Fortran Simo trends are generally closer to the reference
than the current DKMQ24 recovery.  Simo must keep a separate writer/recovery
route until a matrix-level comparison demonstrates a Simo-specific defect.

## DKMQ24 result that motivated the isolated recovery

For the external-support twisting-moment samples in problem 2-006, the old
CQUADR/DKMQ24 output was approximately:

```text
133, 263, 510, 721, 860, 743, 581
```

The direct `Mxy` recovery path produced approximately:

```text
222, 372, 661, 907, 1095, 1112, 1113
```

The comparison reference was approximately:

```text
0, 370, 700, 990, 1210, 1310, 1280
```

This improves the trend substantially but is not an exact match. `Mxx` and
`Myy` remain on the original Gauss-sample/polynomial-fit path because the
direct nodal experiment did not improve them consistently.

The Python DKMQ24 implementation shows a similar qualitative recovery issue,
but not identical values. That does not yet prove a stiffness error. Different
local frames, tensor conventions, sampling positions, shared-grid averaging,
and output transformations can all change reported moments without changing
the assembled stiffness matrix.

## DKMT18 audit result

No Fortran DKMT18 recovery patch was applied in this audit. The existing path
already evaluates the bending operator directly at the centroid:

```fortran
XI0  = ONE/THREE
ETA0 = ONE/THREE
CALL BUILD_STAGE_BB_MAKNUN(..., BB)
BE2(1:3,1:18,1) = BB
```

The Python patch test was run with:

```powershell
$env:PYTHONIOENCODING='utf-8'
python prob_2_001_6dof_triangle.py
```

For all ten triangular elements, `DKMT18 Thick` recovered:

```text
Membrane: Sxx = 1333.333, Syy = 1333.333, Sxy = 400.000
Bending : Mxx = 1.1111e-7, Myy = 1.1111e-7, Mxy = 3.3333e-8
```

These match the independent patch-test target. DKMT18, MITC3+ 6-DOF, and
MITC3+D Hughes-Brezzi all passed this test.

This test validates constant membrane and constant-curvature recovery. It does
not validate curved-shell recovery, nonuniform moment gradients, shared-node
averaging, or local-to-global output transformation on distorted geometry.

## Python versus Fortran comparison requirements

Before changing DKMT18, compare exactly the same quantity at exactly the same
point:

1. Use centroid `(xi,eta)=(1/3,1/3)` in both implementations.
2. Compare the local `BB` matrix before multiplying by displacement.
3. Compare the local curvature vector before constitutive multiplication.
4. Compare local `Mxx, Myy, Mxy` before any coordinate rotation.
5. Apply the same tensor transformation exactly once.
6. Only then compare F06/OP2 engineering-force values.

Do not compare a Python centroid value with a MYSTRAN corner-extrapolated or
grid-averaged value.

## Questions for independent review

1. Does row 3 of every bending operator represent tensor curvature `kxy`,
   engineering twist `2*kxy`, or another sign convention?
2. Is the constitutive `D66` term consistent with that row-3 convention?
3. Is `Mxy` converted once, and only once, between tensor and engineering
   conventions?
4. Are triangle natural coordinates and node ordering identical between
   Python and Fortran?
5. Is the element-local basis right-handed and identically oriented?
6. Is the reported result element-centroid, element-corner, or averaged grid
   output?
7. If grid averaging is used, is it arithmetic, area weighted, angle filtered,
   property filtered, or surface filtered?
8. Is the local-to-basic/surface transformation being applied in both the
   kernel and output writer?

## Safe next experiment

Use a curved or twisted triangular benchmark where DKMT18 and the reference
have a clear nonuniform moment field. Instrument one element and print:

```text
BB_local
u_local
kappa_local
moment_local
moment_basic
```

Capture the same five objects in Python. If `BB_local` already differs, audit
the formulation. If only `moment_basic` differs, audit coordinate transforms.
If element values agree but grid values differ, audit averaging/output only.

Until that experiment identifies a mismatch, retain the current DKMT18
centroid recovery and keep the DKMQ24 direct-corner workaround isolated.
