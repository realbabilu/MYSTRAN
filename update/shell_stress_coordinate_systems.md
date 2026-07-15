# Shell Stress Coordinate Output

## Purpose

MYSTRAN now supports using `PARAM,STR_CID,0` to request shell stress and strain output in the basic/global coordinate system.

This is useful for verification decks whose reference values are expressed as global engineering components, such as SAP2000 patch-test tables. It should not be confused with Nastran-style element stress output, which is normally reported in element or material coordinates.

## Default behavior

The default remains unchanged:

```nastran
PARAM,STR_CID,-2
```

For shell elements this keeps output in the element/material-oriented system used by the existing MYSTRAN/Nastran-style stress paths.

This default is intentional because Nastran F06 shell stress tables are generally element-coordinate tables. MSC/Nastran prints shell tables with wording such as:

```text
STRESSES IN ELEMENT COORD SYSTEM
```

Therefore, default MYSTRAN validation against Nastran-style local component output should not use `PARAM,STR_CID,0`.

## Basic/global shell stress output

Use:

```nastran
PARAM,STR_CID,0
```

Effect:

- Shell membrane stress components are transformed to basic/global coordinates before F06/OP2/NEU output.
- Shell bending stress/curvature-style components are transformed consistently.
- Shell strain engineering shear terms are handled as tensor strain during transformation and converted back to engineering strain afterward.
- Principal stress, major/minor stress, and von Mises values are recomputed from the transformed components.

Supported shell element families in this pass:

- `CTRIA3`
- `CQUAD4`
- `CQUADR`
- `CQUAD8`

## Patch-test interpretation

For SAP2000 shell verification problem `2-001`, the reference membrane table is in global component form:

```text
Sxx = 1333
Syy = 1333
Sxy =  400
```

With local element output, only elements whose local x-y axes happen to align with the global x-y axes will show these same component values. Other elements can have different `Normal-X`, `Normal-Y`, and `Shear-XY` while preserving the same principal stresses.

With:

```nastran
PARAM,STR_CID,0
```

the same deck should report the global component values consistently for this patch test.

## Relation to GPSTRESS / GSTRESS

Altair OptiStruct documents `GPSTRESS` / `GSTRESS` as an I/O Options or Subcase Information request for grid point stress output. It is a nodal/grid-point stress result request, not a switch that changes element stress coordinates.

Reference:

```text
https://help.altair.com/hwsolvers/os/topics/solvers/os/gpstress_gstress_io_r.htm
```

Current MYSTRAN status:

- `GPSTRESS` is accepted as a compatibility alias.
- `GSTRESS` is accepted as a compatibility alias.
- `STRFIELD` is accepted as a compatibility Case Control request.
- `OUTPUT(POST)` is accepted without ending Case Control parsing.
- Post-processing `SURFACE` and `VOLUME` entries are accepted and counted.
- MSC-style post-set syntax such as `SET 1 ALL` is accepted after `OUTPUT(POST)`.
- `PARAM,POSTEXT,YES` is accepted as a compatibility no-op.
- The alias maps to `STRESS(CORNER)` and defaults `PARAM,STR_CID` to `0` when the user has not explicitly selected another stress coordinate system.
- A baseline MSC-readable `OGS1` OP2 table is written when `GPSTRESS/GSTRESS` is active and OP2 output is open.

This is intentional for the first Nastran-compatibility pass. `GPSTRESS/GSTRESS` should currently be read as "request accepted, global basic-system corner stress emitted, plus a baseline `OGS1` table" rather than full MSC grid-point stress recovery.

MSC/Nastran behavior note:

- MSC documentation says `GPSTRESS` calculates grid point stresses from adjoining plate/solid elements in a user-defined coordinate system.
- For meaningful grid-point stress recovery, MSC uses surfaces or volumes defined in the `OUTPUT(POST)` Case Control section.
- MSC also requires `STRESS` or `ELSTRESS` requests for the elements in the surfaces/volumes of interest.
- The `SURFACE`/`VOLUME` definitions are Case Control `OUTPUT(POST)` definitions, not Bulk Data entries.

The local test deck `2_001_4_msc_gpstress.dat` has been corrected to use the MSC-compatible `OUTPUT(POST)` surface path:

```nastran
GPSTRESS=ALL
STRFIELD=ALL
SUBCASE 1
  DISPLACEMENT = ALL
  STRESS(CENTER) = ALL
SUBCASE 2
  DISPLACEMENT = ALL
  STRESS(CENTER) = ALL
OUTPUT(POST)
SET 1 ALL
SURFACE 100 SET 1 NORMAL Z
BEGIN BULK
```

For this deck, `OUTPUT(POST)`, `SET`, and `SURFACE` were accepted by MSC when placed after the normal subcase definitions and before `BEGIN BULK`. Placing `OUTPUT(POST)` before the subcases caused MSC 2025.2 to reject the later `SURFACE` line as an illegal Case Control keyword.

The successful MSC run writes:

- `OGS1` for grid-point surface stresses.
- `OES1X1` for ordinary element stresses.
- `OUG1` for displacements.

`pyNastran` reads the resulting MSC OP2 with `grid_point_surface_stresses=2`, matching the two subcases in the deck.

Current MYSTRAN baseline readback for the same deck:

```text
['OUGV1', 'OES1X1', 'OGS1', 'OUGV1', 'OES1X1', 'OGS1']
grid_point_surface_stresses = 2
```

The MYSTRAN `OGS1` writer currently emits shell corner-derived Z1/Z2 values at grid IDs. It is structurally OP2-compatible with `pyNastran`, but it is not yet the full MSC grid-point stress averaging/recovery algorithm.

Example:

```nastran
GPSTRESS(PLOT) = ALL
```

is treated internally as:

```nastran
STRESS(CORNER,PLOT) = ALL
PARAM,STR_CID,0
```

unless the deck later overrides `PARAM,STR_CID`.

Remaining future implementation:

- Require or auto-diagnose the matching `STRESS/ELSTRESS` element stress request.
- Implement the MSC-compatible surface/volume grid-point stress recovery path separately from the current corner-stress alias.
- Replace the baseline corner-derived `OGS1` payload with the recovered MSC-style grid-point stress result.
- Keep `STRESS` as element-coordinate output by default; keep `GPSTRESS/GSTRESS` as the future MSC-compatible grid-point stress output.

## Practical validation rule

Use local/default output when comparing against Nastran element stress component tables.

Use `PARAM,STR_CID,0` when comparing against references that publish global stress components, such as many SAP2000 verification tables.

Use principal stresses as a coordinate-invariant comparison when validating solver behavior across different local coordinate conventions.
