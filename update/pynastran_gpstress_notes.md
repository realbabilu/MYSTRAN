# pyNastran GPSTRESS / OGS1 Notes

Date: 2026-07-19

## Scope

This note records the current understanding of how `pyNastran` reads
`OGS1` surface grid-point stress tables (`GPSTRESS`) and what that means for
MYSTRAN validation against MSC/Nastran reference OP2 files.

Main local case used during this audit:

- `D:\18a\MYSTRAN_Validation-main\ctriar_cquadr\2_001_4_msc_gpstress.OP2`
- `D:\18a\MYSTRAN_Validation-main\ctriar_cquadr\2_001_4_mystran_gpstress.OP2`

## What pyNastran advertises

`pyNastran` exposes `OGS1` surface stress results through:

- `op2.grid_point_surface_stresses`

Its class currently reports the generic headers:

- `['nx', 'ny', 'txy', 'angle', 'majorP', 'minorP', 'tmax', 'ovm']`

Relevant local source paths:

- `C:\Users\arypr\AppData\Local\Programs\Python\Python312\Lib\site-packages\pyNastran\op2\tables\ogs_grid_point_stresses\ogs_surface_stresses.py`
- `C:\Users\arypr\AppData\Local\Programs\Python\Python312\Lib\site-packages\pyNastran\op2\tables\ogs_grid_point_stresses\ogs.py`

## What the MSC reference OP2 actually looks like

For the MSC reference deck above, `pyNastran` reads:

- same key structure as MYSTRAN
- shape `(1, 40, 8)` for each subcase
- repeated `(grid, element, fiber)` rows in the expected MSC order

However, the raw float payload does not behave like the generic header labels
after the first three fields.

Example from the MSC reference file, first `Z1` row:

- raw row:
  - `[1490.901, 1175.765625, 367.65802, 1733.333374, 933.333313, 1502.590332, 284.62851, 0.0]`

These values behave like:

1. `NX`
2. `NY`
3. `TXY`
4. principal/invariant value, not a literal angle
5. principal/invariant value
6. principal/invariant value
7. principal/invariant value
8. unused/zero in this case

So, while `pyNastran` labels the columns as if slot 4 were `angle`, the MSC
reference payload for this case clearly is not storing a raw angle there.

## What MYSTRAN was doing

The current MYSTRAN writer computes a classic F06-style surface row:

- `XX`
- `YY`
- `XY`
- `ANGLE`
- `MAJOR`
- `MINOR`
- `MAXSHEAR`
- `VONMISES`

That is internally reasonable for F06, but it does not match the observed MSC
reference OP2 payload convention for `OGS1`.

This is why an OP2 file can be:

- structurally valid
- readable by `pyNastran`
- correct in row count and record framing

yet still disagree numerically if a validator blindly trusts the generic
header names for the later `OGS1` slots.

## Current safe rule for validation

For `GRIDPOINTSURFACESTRESSES` in validation:

1. trust only the direct plane components from OP2:
   - `NX`
   - `NY`
   - `TXY`
2. derive these quantities from them:
   - `PRINCIPALANGLE`
   - `MAJOR`
   - `MINOR`
   - `MAXSHEAR`
   - `VONMISES`
3. do not assume the remaining raw `OGS1` slots have portable,
   solver-independent meanings unless verified against a known-good solver and
   case family

## Local parser hardening applied

`D:\18a\MYSTRAN_Validation-main\op2_query.py` was updated so that for
`GRIDPOINTSURFACESTRESSES` it now:

- reads `NX`, `NY`, `TXY` from the OP2 table
- recomputes:
  - `PRINCIPALANGLE`
  - `MAJOR`
  - `MINOR`
  - `MAXSHEAR`
  - `VONMISES`

This makes GPSTRESS validation stable across:

- MSC reference OP2
- MYSTRAN OP2

even when the raw invariant slots do not match the generic `pyNastran`
surface-table labels.

## Current MYSTRAN status

As of 2026-07-19, the local MYSTRAN `OGS1` path for the investigated shell
patch deck has these properties:

- OP2 table framing is valid
- `pyNastran` reads the file without crashing
- row count matches MSC for the audited deck
- row ordering matches MSC for the audited deck

So the remaining semantic caution is no longer about binary validity. It is
about interpretation of the later `OGS1` float fields.

## Recommended practice going forward

For GPSTRESS comparison workflows:

- use OP2 row structure for presence/order validation
- use recomputed invariants from `NX/NY/TXY` for value validation
- treat F06 and OP2 as related but not identical payload conventions

If exact raw-slot parity with MSC `OGS1` is required in the future, then
`WRITE_OGS1_SURFACE_STRESS` should be audited against more than one known-good
MSC/NX case to determine the true binary convention for the later invariant
slots, instead of assuming the generic `pyNastran` labels are exact.
