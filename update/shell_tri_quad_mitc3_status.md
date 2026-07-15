# Shell Triangle/Quadrilateral Element Status

## Scope

This note records the current `mystran18a` shell-element output status after the `CTRIAR`, `CQUADR`, `MITC3+`, global stress-coordinate, and baseline `OGS1` work.

The purpose is to avoid mixing three different topics:

- Bulk Data card compatibility.
- Element formulation used for stiffness/recovery.
- F06/OP2/NEU output support.

## Current element mapping

| Bulk Data card | Internal `TYPE` | Current formulation path | Current output path |
| --- | --- | --- | --- |
| `CTRIA3` | `TRIA3` | MIN3/MITC3+ path, depending branch settings | TRIA3 stress/strain/force writers |
| `CTRIAR` | `TRIA3` with DKMT sentinel | `CTRIAR_DKMT18` for DKMT-style shell stiffness | Reuses TRIA3 stress/strain/force writers |
| `CQUAD4` | `QUAD4` | MIN4/MITC4+ path, depending branch settings | QUAD4 stress/strain/force writers |
| `CQUADR` | `QUADR` | `CQUADR_DKMQ24` | QUAD4/CQUADR stress/strain/force writers |
| `CQUAD8` | `QUAD8` | MITC8 path | QUAD8 shell stress/strain/force writers |
| `CTRIA6` | not present | not implemented | not implemented |

## `CTRIAR`

`CTRIAR` is parsed by the existing `BD_CTRIA` path and stored internally as `TYPE='TRIA3'`.

The distinguishing marker is stored in `EDAT` as a DKMT18 sentinel:

```text
EDAT(... thickness key / CTRIAR sentinel) = -18
```

This lets the element-generation side dispatch to `CTRIAR_DKMT18` while still reusing the existing triangular output infrastructure.

Practical consequence:

- `CTRIAR` benefits from the same F06/OP2/NEU triangular stress/strain/force output path as `CTRIA3`.
- The new baseline `OGS1` writer includes `CTRIAR` because it checks internal `TYPE(1:5) == 'TRIA3'`.
- Output is currently corner-derived, not full MSC grid-point stress averaging.

## `CQUADR`

`CQUADR` is parsed as a separate internal shell type:

```text
TYPE = 'QUADR   '
```

Its stiffness path dispatches to:

```text
CQUADR_DKMQ24
```

The LINK9 output paths have been expanded so `QUADR` is treated with the QUAD4 family where appropriate:

- Element engineering forces.
- Element stresses.
- Element strains.
- FEMAP neutral element result vectors.
- Baseline `OGS1` grid-point surface stress output.

Practical consequence:

- `CQUADR` is not hidden as `QUAD4`; it has its own formulation path.
- Output compatibility intentionally reuses the QUAD4-family output layout where the result topology is four-corner shell data.
- The baseline `OGS1` writer includes `CQUADR` because it checks `TYPE == 'QUADR   '`.

## MITC3+

`MITC3+` is an alternate triangular shell formulation path for `CTRIA3`-family elements.

Current status:

- The active output type remains `TRIA3`.
- F06/OP2/NEU result writers do not need a separate `MITC3+` table name.
- Stress/strain/force values depend on the formulation/recovery arrays filled upstream, then reuse the standard TRIA3 output path.

Important distinction:

- `MITC3+` is a formulation/recovery choice.
- `CTRIA3`/`CTRIAR` are Bulk Data element names.
- OP2/F06 still see these as triangular shell result tables, not a new Nastran element table family.

## Baseline `OGS1` support

The current baseline `OGS1` OP2 writer supports:

- `CTRIA3`
- `CTRIAR`
- `CQUAD4`
- `CQUADR`
- `CQUAD8`

Current pyNastran readback on the local `GPSTRESS` patch-test deck confirms:

```text
['OUGV1', 'OES1X1', 'OGS1', 'OUGV1', 'OES1X1', 'OGS1']
grid_point_surface_stresses = 2
```

Limit:

This is corner-derived `OGS1` output. It is structurally OP2-compatible, but it is not yet the full MSC/Nastran `GPSTRESS` surface/volume averaging and grid-point stress recovery algorithm.

## `CTRIA6` status

`CTRIA6` is not currently implemented in this branch.

No active source support was found for:

- `CTRIA6`
- `TRIA6`
- `NCTRIA6`
- `MEDAT_CTRIA6`
- `BD_CTRIA6`

Therefore, adding only an OP2/F06 writer would be incorrect. There is no element stiffness, recovery, or stress array source for `CTRIA6` yet.

Minimum safe implementation path for real `CTRIA6` support:

1. Add parser/counting support:

```text
BD_CTRIA6 or extended BD_CTRIA
NCTRIA6
MEDAT_CTRIA6
ELMTYP entry 'TRIA6   '
NELGP = 6
```

2. Add EMG formulation:

```text
TRIA6 stiffness/mass/recovery path
material angle/offets/thickness handling
pressure/thermal/differential stiffness where applicable
```

3. Add LINK9 recovery sizing:

```text
NUM_SEi for TRIA6
MAXREQ_OGEL handling
ELEM_STRE_STRN_ARRAYS handling
CALC_ELEM_STRESSES / CALC_ELEM_STRAINS handling
```

4. Add output support:

```text
F06 stress/strain/force rows
OP2 OES/OEF/OSTR table mapping
baseline OGS1 rows
NEU element result vectors
```

5. Validate against MSC/NX:

```text
displacement
element force
center stress/strain
corner/grid stress/strain
OP2 readback through pyNastran
```

Until those steps exist, `CTRIA6` should stay documented as unsupported rather than partially accepted.

## Practical validation note

For current shell patch tests:

- Use `CTRIA3`/`CTRIAR` for 3-node triangle comparisons.
- Use `CQUAD4`/`CQUADR` for 4-node quadrilateral comparisons.
- Use `CQUAD8` only where the branch limitations are acceptable.
- Do not use `CTRIA6` as a validation target yet.

