# DSQK vs DKMQ20 Audit for Problem 2-001

Date: 2026-08-13

## Scope

This note audits why `CQUAD4_DSQK_RHR` fails patch test `2-001` in MYSTRAN while
`CQUAD4_DKMQ20_RHR` passes, and whether the discrepancy could be caused by
`GPSTRESS`/postprocessing instead of the element formulation.

## Bottom line

The earlier "DSQK kernel is wrong" diagnosis is no longer supported.

Current audit status:

- the Fortran `DSQK` local membrane recovery matches the upgraded Python
  reference when compared in the same element-local basis
- the Fortran DSQK element basis (`T56`, derived from `T1/T2/T3`) also matches
  the Python reference for the audited patch elements
- the remaining mismatch is downstream, in how `QUAD4/QUADR` shell stresses are
  interpreted or re-expressed in output space, especially in the `GPSTRESS`
  surface path

So the active problem is now best described as an output-space or
local-vs-surface/global handling issue, not a primary DSQK formulation-port
failure.

## Implemented output-path fix

The `GPSTRESS` surface path has now been patched in MYSTRAN so that shell
`QUAD4/QUADR` rows carry their shell basis (`TE`) through LINK9 and are rotated
from shell-local to surface-system before the OGS1 surface table is assembled.

Result on `prob_2_001_thick_cquad4_dsqk.dat`:

- the surface `GPSTRESS` table is now uniform at
  `[1333.333, 1333.333, 400.000]`
- this confirms that the earlier apparent DSQK patch failure in surface output
  was caused by postprocessing basis handling, not by the DSQK kernel

## Current remaining nuance

The standard shell stress block for `QUAD4` still writes direct element recovery
rows, then labels the table as basic-coordinate output.

After the latest patch:

- those rows are now rotated row-by-row toward basic output before formatting
- but the table still reflects the direct recovery points and not the
  surface-reconstructed/fitted field used by `GPSTRESS`

That means:

- `GPSTRESS` now shows the correct patch-uniform global/surface field
- the traditional shell table still preserves element-point recovery character
  and therefore does not yet collapse to the same uniform patch view

So the remaining work is no longer a basis bug. It is a product decision about
whether the standard shell table should continue to show direct recovery values
or should be rebuilt from the same postprocessed global/surface field used for
patch-style reporting.

## Updated evidence

### Local-to-local parity is good

For the MYSTRAN patch deck connectivity, direct Fortran-vs-Python comparison now
shows that DSQK center local membrane stresses agree element by element.

Representative local stress values from the Python reference for the MYSTRAN
deck ordering are:

- EID 1: `[1175.7656, 1490.9011, 367.6580]`
- EID 2: `[1359.5066, 1307.1600, 399.1428]`
- EID 3: `[1333.3333, 1333.3333, 400.0000]`
- EID 4: `[1333.3333, 1333.3333, 400.0000]`
- EID 5: `[1411.7798, 1254.8869, 392.2323]`

These match the corresponding Fortran `DEBUG(239)` recovery dump from
`ELEM_STRE_STRN_ARRAYS`.

### Element basis parity is good

`DEBUG(233)` output from `CQUAD4_DSQK_RHR.f90` shows that the element basis used
inside the Fortran DSQK kernel matches the Python reference.

Examples:

- EID 1 `T56` starts with `t1 = [ 0.9795777, -0.2010659, 0 ]`
- EID 2 `T56` starts with `t1 = [ 0.9994641,  0.0327342, 0 ]`
- EID 5 `T56` starts with `t1 = [ 0.9951333,  0.0985376, 0 ]`

These are the same directions previously extracted from the Python DSQK
reference for the same element ordering.

### Updated GPSTRESS diagnosis

The remaining issue now points to the shell output path:

- `SHELL_STRESS_OUTPUTS.f90` writes `QUAD4/QUADR` shell stress rows to `OGEL`
  directly from the current `STRESS` vector
- `WRITE_ELEM_STRESSES.f90` surface `GPSTRESS` reconstruction then consumes only
  `OGEL(:,2:4)` as if they were already surface-system stress components
- unlike the `QUAD8` path, the `QUAD4/QUADR` path does not apply an explicit
  `TRANSFORM_SHELL_STR` step before those surface rows are reused

This is consistent with the observed behavior:

- local DSQK recovery agrees with Python
- element-local stress can still look non-uniform across the patch
- the apparent patch failure shows up when those local stress components are
  treated as surface/global components in later output

## Evidence from F06

### DKMQ20 control case is self-consistent

From `prob_2_001_thick_cquad4_dkmq20.F06`:

- `GPSTRESS` surface stresses report:
  - `NORMAL-X = 1.333333E+03`
  - `NORMAL-Y = 1.333333E+03`
  - `SHEAR-XY = 4.000000E+02`
- engineering forces report:
  - `Nxx = 1.333333E+00`
  - `Nyy = 1.333333E+00`
  - `Nxy = 4.000000E-01`
  - `Mxx = 1.111111E-07`
  - `Myy = 1.111111E-07`
  - `Mxy = 3.333333E-08`

The automation parser returns the same values.

### DSQK failure is already in engineering forces

From `prob_2_001_thick_cquad4_dsqk.F06`:

- subcase-1 stress block already drifts:
  - `MAX* sxx = 1.41178E+03`
  - `MAX* syy = 1.49090E+03`
  - target should be `1.33333E+03`
- subcase-2 engineering forces also drift:
  - `Nxx = 1.411780E+00`
  - `Nyy = 1.490901E+00`
  - `Nxy = 4.000000E-01`
  - `Mxx = 7.447141E-07`
  - `Myy = 3.178054E-07`
  - `Mxy = 8.411761E-07`

This rules out a pure postprocessing bug.

## Implementation comparison

## 1. Membrane operator construction

`CQUAD4_DKMQ20_RHR` builds its membrane operator directly in 24-DOF shell form:

- `BMB = BM_AT(...)`
- `KMEM += BMB^T * A * BMB`

Reference:

- [CQUAD4_DKMQ20_RHR.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/CQUAD4_DKMQ20_RHR.f90#L144)
- [CQUAD4_DKMQ20_RHR.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/CQUAD4_DKMQ20_RHR.f90#L159)

`CQUAD4_DSQK_RHR` instead builds a 20-DOF reduced operator first, then lifts it
to 24 DOF with `T56`:

- `BM5 = BM_DSQK(...)`
- `KM5 += BM5^T * A * BM5`
- `KBASIC = T56^T * KTOT5 * T56`

Reference:

- [CQUAD4_DSQK_RHR.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/CQUAD4_DSQK_RHR.f90#L143)
- [CQUAD4_DSQK_RHR.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/CQUAD4_DSQK_RHR.f90#L156)
- [CQUAD4_DSQK_RHR.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/CQUAD4_DSQK_RHR.f90#L170)

Audit takeaway:

- `DKMQ20` never depends on a drilling-lifting map to define membrane strain.
- `DSQK` currently does.
- Since the patch-test membrane result is already wrong in DSQK, the first
  suspect is the `20 -> 24` lifting path, not `GPSTRESS`.

Direct compare:

- `DKMQ20` membrane strain is built from physical shell directions `T1/T2` and
  written straight into translational DOF columns of the 24-DOF element matrix.
- `DSQK` membrane strain is first built in a reduced 5-DOF/node system
  (`u1,u2,w,rot1,rot2`) and only later projected into MYSTRAN's 24-DOF shell
  through `T56`.
- That means any inconsistency in `BUILD_DSQK_BASIS` or `BUILD_T56` pollutes the
  membrane field before bending recovery is even reached.

## 2. Bending operator composition

`DKMQ20` constructs bending directly as:

- `BBB = BB_AT(...)`
- `KBEND += BBB^T * D * BBB`

Reference:

- [CQUAD4_DKMQ20_RHR.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/CQUAD4_DKMQ20_RHR.f90#L146)
- [CQUAD4_DKMQ20_RHR.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/CQUAD4_DKMQ20_RHR.f90#L160)

`DSQK` splits bending into:

- basic Kirchhoff part `BBU5`
- edge correction `BBD`
- assumed parameter map `AN4`
- recombined `BB5 = BBU5 + BBD * AN4`

Reference:

- [CQUAD4_DSQK_RHR.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/CQUAD4_DSQK_RHR.f90#L144)
- [CQUAD4_DSQK_RHR.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/CQUAD4_DSQK_RHR.f90#L146)
- [CQUAD4_DSQK_RHR.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/CQUAD4_DSQK_RHR.f90#L168)

Audit takeaway:

- This decomposition matches the intended DSQK family structure.
- However, because membrane is already off, bending cannot be trusted yet.
- The current DSQK failure is not a bending-only issue in MYSTRAN.

Direct compare:

- `DKMQ20` bending uses one consistent AU-family path:
  `BB_AT(..., AINV_AU)`.
- `DSQK` bending uses a different split path:
  `BBU_DSQK + BBDELTA_DSQK * AN_DSQK`.
- That difference is acceptable in principle, but it means `DSQK` is not a
  small variation of `DKMQ20`; it is a separate formulation and needs its own
  exact 6-DOF port, not an approximate family transplant.

## 3. Drilling stabilization policy differs sharply

`DKMQ20` currently uses no extra drilling penalty in this branch:

- `KDRILL = ZERO`

Reference:

- [CQUAD4_DKMQ20_RHR.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/CQUAD4_DKMQ20_RHR.f90#L171)

`DSQK` currently adds a nodal block penalty aligned only with the shell normal:

- `DRK = DRILL_SCALE * GVAL * THK * AREA`
- per-node addition on rotational block `DRK * e3 * e3^T`

Reference:

- [CQUAD4_DSQK_RHR.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/CQUAD4_DSQK_RHR.f90#L698)
- [CQUAD4_DSQK_RHR.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/CQUAD4_DSQK_RHR.f90#L705)
- [CQUAD4_DSQK_RHR.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/CQUAD4_DSQK_RHR.f90#L710)

Audit takeaway:

- This DSQK drill term is much more intrusive than DKMQ20.
- It may contaminate prescribed-displacement patch tests, especially when the
  `T56` lifting already couples reduced rotations into the full 24-DOF shell.
- Even so, because DSQK membrane stress is already wrong before any `Qabs`
  issue, drilling is likely a secondary amplifier, not the only root cause.

Direct compare:

- `DKMQ20` passes the patch test without any extra drill penalty in this branch.
- `DSQK` adds a normal-axis rotational penalty at every node.
- Since DSQK already shows spurious `T3/R1/R2` reactions while DKMQ20 stays near
  machine zero, this drill term is a high-value debug toggle candidate.

## 4. Recovery path is not the primary culprit

Both elements populate `BE1/BE2/BE3` from their Gauss-point `B` matrices and
then overwrite `BE2(3,:,gp)` at center/corners for twisting-moment recovery.

References:

- [CQUAD4_DKMQ20_RHR.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/CQUAD4_DKMQ20_RHR.f90#L178)
- [CQUAD4_DSQK_RHR.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/CQUAD4_DSQK_RHR.f90#L179)

Audit takeaway:

- DSQK engineering forces are already wrong.
- Therefore a bug only in `BE1/BE2` recovery cannot explain the full patch-test
  failure.

## 5. Additional symptom: DSQK creates larger out-of-plane SPC reactions

In the same patch deck:

- `DKMQ20` out-of-plane SPC reaction terms are around `1e-8` or smaller
- `DSQK` shows spurious `T3/R1/R2` reactions around `1e-5` to `1e-6`

Audit takeaway:

- DSQK currently introduces parasitic bending/drilling coupling under a test that
  should stay very clean.
- This supports the formulation-side diagnosis.

## Single-element parity result

A focused one-element compare was carried out against the Python reference for
element 1 at the first Gauss point (`xi = eta = -1/sqrt(3)`), using:

- Fortran debug dump from `DEBUG(233)` in
  [CQUAD4_DSQK_RHR.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/CQUAD4_DSQK_RHR.f90)
- Python reference script:
  [tmp_dsqk_element1_compare.py](D:/18a/tmp_dsqk_element1_compare.py)

The following objects match to normal numerical roundoff:

- local basis vectors `t1`, `t2`, `t3`
- assumed-parameter map `AN4`
- lift matrix `T56`
- membrane operator `BM5`
- basic bending operator `BBU5`
- edge correction operator `BBD`
- recombined bending operator `BB5`
- transverse shear operator `BS5`
- lifted operators `BM24` and `BB24`
- final reduced-to-6DOF stiffness `KBASIC`

Audit consequence:

- the DSQK formulation port itself is not the current root cause
- the `20 -> 24` lift strategy is not the current root cause
- the old suspicion that `T56` or `BM_DSQK` was the main failure source is no
  longer supported by the data

The only small observed difference is the separate drill penalty magnitude
(`0.072` in the current Fortran dump versus `0.0672` in the Python reference
for that one element), but the no-drill experiment already showed that this is
not what drives the patch-test membrane error.

## Working diagnosis

Primary suspect order after the direct Python-vs-Fortran parity check:

1. recovery/output path around `BE1/BE2/BE3`, `UEL`, and stress-resultant usage
2. element-axis transform usage through `TE` / `T24`
3. any downstream shell-force assembly that assumes a different local/global
   convention than the DSQK branch now uses
4. only after that, secondary cleanup such as drill scaling parity

## Important Python reference note

The current Python file `DSQK_ShellElement_RHR_6DOF.py` is not a native
`24 x 24` DSQK formulation in the same sense as `DKMQ20_ShellElement_RHR_6dof.py`.

It explicitly does:

- `K5 = super().k_local()`
- `T = _build_T_5to6()`
- `K24 = T.T @ K5 @ T`
- then adds a nodal drilling penalty about `t3`

So the present Python DSQK 6DOF reference is architecturally aligned with the
current MYSTRAN DSQK port:

- reduced 5-DOF kernel first
- then `20 -> 24` lifting
- then separate drilling stabilization

Audit consequence:

- Porting `DSQK_ShellElement_RHR_6DOF.py` "as is" will not remove the `T56`
  dependency, because that file also depends on the same idea.
- The key task is not "replace Fortran with native 24-DOF DSQK", because the
  present Python reference is not yet that object.
- The real task is to match the Python DSQK 6DOF wrapper exactly, especially:
  - basis construction
  - `T` map row meanings
  - drill penalty scaling and placement

## Drill-off experiment

A dedicated patch-test deck was run with:

- `PARAM,QUAD4TYP,DSQK`
- `DEBUG   238     1`

At that moment `DEBUG(238)` had been temporarily wired to disable
`BUILD_DSQK_DRILL`, but the run also revealed that `DEBUG(238)` is already used
by LK9 recovery/polyfit debugging (`POLYFIT238` / `RECOV238` output). The DSQK
no-drill switch was then moved to `DEBUG(250)` to avoid that collision in future
experiments.

Result:

- normal DSQK:
  - `Nxx = 1.411780`
  - `Nyy = 1.490901`
  - `Nxy = 0.400000`
  - `Mxx = 7.447141E-07`
  - `Myy = 3.178054E-07`
  - `Mxy = 8.411761E-07`
- DSQK with drill disabled:
  - `Nxx = 1.411780`
  - `Nyy = 1.490901`
  - `Nxy = 0.400000`
  - `Mxx = 7.447141E-07`
  - `Myy = 3.178054E-07`
  - `Mxy = 8.411761E-07`

Audit consequence:

- the current DSQK patch-test error is not caused by the added drilling term
- the dominant bug is upstream of drilling, inside the reduced DSQK kernel,
  its force recovery, or the way the reduced quantities are assembled back into
  MYSTRAN shell outputs

Additional recovery conclusion from the same run:

- `RECOV238` showed the wrong DSQK stresses already present at the element
  recovery points
- `POLYFIT238` then reproduced those same values essentially exactly

So the LK9 polynomial fit is not creating the DSQK error; it is only
propagating stresses that are already wrong upstream.

## Membrane trace refinement on August 13, 2026

After adding a dedicated `DEBUG(239)` trace in
[ELEM_STRE_STRN_ARRAYS.f90](D:/18a/MYSTRAN/Source/LK9/L92/ELEM_STRE_STRN_ARRAYS.f90),
the patch-test diagnosis tightened again.

For `prob_2_001_thick_cquad4_dsqk_dbg239.dat`:

- the membrane subcase global translations printed in the F06 match the DKMQ20
  control solution exactly
- element 1, recovery point 1 shows `UEB` equal to the exact patch field:
  - node 1: `(0, 0)`
  - node 2: `(5.0E-05, 4.0E-05)`
  - node 3: `(1.2E-04, 1.2E-04)`
  - node 4: `(6.0E-05, 1.2E-04)`
- but the same trace still gives
  - `STRAIN1 = [8.0304035E-04, 1.1969596E-03, 9.1914503E-04]`
  - `STRESS1 = [1.1757656E+03, 1.4909011E+03, 3.6765801E+02]`

This is a stronger result than the earlier LK9 audit:

- the solution field in basic coordinates is already correct
- the error appears before GPSTRESS and before LK9 polynomial fitting
- the bad membrane state is already encoded in the DSQK recovery operator that
  acts on `UEB`

Because `BE1` is stored from `BM24 * T24^T`, the recovery path should collapse
to `BM24 * UEB`. That means the current failure is no longer best described as a
generic `TE/T24` or LK9 output bug. The sharper hypothesis is:

- the Fortran DSQK membrane operator still reflects the older RHR-era kernel,
  while the upgraded Python DSQK path used as reference has changed in a way
  that restores the constant-strain patch behavior.

## Recommended next actions

1. Compare the patch-test membrane operator in the upgraded Python DSQK source
   against the Fortran `BM24` path in
   [CQUAD4_DSQK_RHR.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/CQUAD4_DSQK_RHR.f90),
   using the exact `UEB` field captured above.
2. Re-run the one-element Python/Fortran matrix audit for the exact distorted
   patch element and include `BM24 * UEB` directly, not only `K24`.
3. Keep the drill-scaling mismatch on the list, but treat it as a later parity
   cleanup item, not the current blocker for problem `2-001`.
