# OGS1 Writer Audit - 2026-07-19

Date: 2026-07-19

## Goal

Audit `WRITE_OGS1_SURFACE_STRESS` against available known-good MSC/NX-style
reference files and decide whether MYSTRAN should:

1. mimic the raw MSC `OGS1` payload more closely, or
2. keep the current structurally valid writer and rely on a hardened
   `op2_query.py` interpretation for validation.

## Writer under audit

- `D:\18a\MYSTRAN\Source\LK9\L91\WRITE_ELEM_STRESSES.f90`
- subroutine:
  - `WRITE_OGS1_SURFACE_STRESS`

## Reference files available in workspace

### Confirmed readable `OGS1` reference

- `D:\18a\MYSTRAN_Validation-main\ctriar_cquadr\2_001_4_msc_gpstress.OP2`

This file is the current known-good baseline in the local workspace for
surface-grid-point shell stress output.

### Additional nearby candidate files checked

- `D:\18a\MYSTRAN_Validation-main\ctriar_cquadr\2_001_4_msc_gpstress.op2.2`
- `D:\18a\MYSTRAN_Validation-main\ctriar_cquadr\2_001_4_msc_gpstress.op2.3`

Both of these read as ordinary shell-stress OP2 content in `pyNastran` and
did **not** expose a usable `grid_point_surface_stresses` container. In other
words:

- they are not additional `OGS1` references
- they cannot be used to confirm raw `OGS1` payload semantics

### NX references

No second or third confirmed NX `OGS1` shell-surface OP2 reference was found
in the current local workspace during this audit pass.

So, for raw `OGS1` semantics, the practical basis remains:

- one good MSC `OGS1` reference file
- one local MYSTRAN `OGS1` file produced from the same deck

## What is now confirmed good in MYSTRAN

For:

- `D:\18a\MYSTRAN_Validation-main\ctriar_cquadr\2_001_4_mystran_gpstress.OP2`

the following are now confirmed:

- `pyNastran` reads the file without crashing
- the old `struct.error: unpack requires a buffer of 4 bytes` issue is gone
- `grid_point_surface_stresses` exists
- row structure matches the MSC reference:
  - same subcase keys
  - same `(1, 40, 8)` shape
  - same repeated `(grid, element, fiber)` ordering

This means the **binary record framing and row-count problem is fixed**.

## What is not portable at the raw payload level

`pyNastran` advertises the generic `OGS1` surface headers as:

- `nx`
- `ny`
- `txy`
- `angle`
- `majorP`
- `minorP`
- `tmax`
- `ovm`

However, the MSC reference raw row does not behave like those labels after the
first three direct components.

### MSC reference example

First `Z1` row from:

- `2_001_4_msc_gpstress.OP2`

raw payload as read by `pyNastran`:

- `[1490.901, 1175.765625, 367.65802, 1733.333374, 933.333313, 1502.590332, 284.62851, 0.0]`

The first three values clearly match:

- `NX`
- `NY`
- `TXY`

But slot 4 is **not** a literal angle, despite the generic `pyNastran` label.

### MYSTRAN current writer example

First `Z1` row from:

- `2_001_4_mystran_gpstress.OP2`

raw payload:

- `[1333.333374, 1333.333374, 400.0, 45.0, 1733.333374, 933.333313, 400.0, 1502.590332]`

This row is internally consistent with MYSTRAN F06-style semantics:

- `NX`
- `NY`
- `TXY`
- `ANGLE`
- `MAJOR`
- `MINOR`
- `MAXSHEAR`
- `VONMISES`

So the current MYSTRAN writer is semantically self-consistent, but it does not
match the observed raw-slot behavior of the MSC reference OP2.

## Key finding

There is currently a split between:

1. **F06-style semantic ordering**
2. **observed raw MSC `OGS1` payload ordering**
3. **generic `pyNastran` surface-table labels**

The first three direct components are safe and portable:

- `NX`
- `NY`
- `TXY`

The later invariant/principal slots are **not yet proven portable** from the
available workspace evidence.

## Decision

### Recommended current direction

Do **not** change MYSTRAN writer semantics just to mimic the single available
MSC raw payload example.

Instead:

- keep the current MYSTRAN `OGS1` writer structurally valid
- keep the improved row-count / row-order behavior
- keep `op2_query.py` hardened so that for `GRIDPOINTSURFACESTRESSES` it:
  - trusts only `NX`, `NY`, `TXY`
  - recomputes:
    - `PRINCIPALANGLE`
    - `MAJOR`
    - `MINOR`
    - `MAXSHEAR`
    - `VONMISES`

### Why this is the safer choice

Because with the current evidence:

- there is only one confirmed MSC `OGS1` reference payload
- there is no second confirmed NX `OGS1` shell-surface reference in the
  workspace
- the generic `pyNastran` labels are not reliable proof of raw solver-slot
  semantics

Changing MYSTRAN writer now to imitate the raw MSC payload would risk:

- breaking the current internal consistency of MYSTRAN F06 vs OP2
- making the file less intuitive for future MYSTRAN-side maintenance
- overfitting to one observed solver payload without enough cross-checks

## When writer-side mimicry would become justified

Revisit writer-side payload mimicry only if at least one of these becomes true:

1. a second MSC `OGS1` shell-surface OP2 confirms the same raw-slot behavior
2. a known-good NX `OGS1` shell-surface OP2 confirms the same behavior
3. official binary-format documentation or DMAP-level evidence identifies the
   exact intended `OGS1` payload convention for those later slots

Until then, the best engineering choice is:

- stable writer
- robust validation parser

## Final conclusion

As of 2026-07-19:

- `WRITE_OGS1_SURFACE_STRESS` is now good enough at the binary-compatibility
  level
- the unresolved issue is interpretation of later invariant slots, not record
  validity
- MYSTRAN does **not** need to mimic the raw MSC payload yet
- `op2_query.py` should remain the compatibility layer for GPSTRESS value
  validation
