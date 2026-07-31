# GPSTRESS Patch-Test Example

## Purpose

This note records the separate GPSTRESS/GSTRESS validation path.

`GPSTRESS` is not the same output as ordinary element stress. It is a grid-point stress request and is represented in OP2 as `OGS1`, while ordinary shell element stress is represented as `OES1X1`.

## Minimal MSC-style case-control shape

A GPSTRESS-style deck should keep ordinary element stress requested and add an `OUTPUT(POST)` surface definition:

```nastran
SET 10 = 100
GPSTRESS = 10
STRFIELD = 10

SUBCASE 1
  SPC = 1
  LOAD = 1
  DISPLACEMENT = ALL
  STRESS(CENTER,CORNER) = ALL

SUBCASE 2
  SPC = 2
  LOAD = 2
  DISPLACEMENT = ALL
  STRESS(CENTER,CORNER) = ALL

OUTPUT(POST)
SET 1 ALL
SURFACE 100 SET 1 NORMAL Z

BEGIN BULK
$PARAM,POST,-1
$PARAM,POSTEXT,YES
```

Notes:

- in the MSC/Nastran-style deck used for the shell patch test, `GPSTRESS` and `STRFIELD`
  are tied to the surface ID through `SET 10 = 100`
- if `OUTPUT(POST)` is used, `PARAM,POST,-1` and `PARAM,POSTEXT,YES` should stay commented out
  or omitted; enabling both paths together can cause a fatal conflict
- ordinary element stress tables and `GPSTRESS` surface tables are separate print paths, so
  repeated/cut-looking headers in F06 should not be interpreted as the same output family

For current MYSTRAN compatibility, `GPSTRESS/GSTRESS` also defaults the stress coordinate request to basic/global when the deck has not explicitly set `PARAM,STR_CID`.

## Current MYSTRAN behavior

The current implementation is a baseline compatibility writer:

- `GPSTRESS` and `GSTRESS` are accepted.
- `STRFIELD` is accepted.
- `OUTPUT(POST)`, `SET`, `SURFACE`, and `VOLUME` are accepted enough to keep MSC-style decks parseable.
- OP2 writes a readable `OGS1` table when OP2 output is enabled.
- The current `OGS1` payload is derived from shell corner stress values at grid IDs.

This is not yet full MSC grid-point stress recovery.

## OP2 mapping

Expected OP2 tables for a simple GPSTRESS shell patch deck:

```text
OUGV1   displacement
OES1X1  ordinary shell element stress
OGS1    grid-point surface stress
```

`pyNastran` exposes `OGS1` as:

```python
op2.grid_point_surface_stresses
```

The validator should not read `OGS1` through ordinary element stress paths such as:

```text
SC/#/SHELLSTRESSES/EID/#/...
```

If GPSTRESS is made part of automated validation, add an explicit path family such as:

```text
SC/#/GRIDPOINTSTRESSES/GID/#/SURFACE/#/Z1/XX
SC/#/GRIDPOINTSTRESSES/GID/#/SURFACE/#/Z2/XX
```

or another clearly named OGS/grid-point path. Do not overload `SHELLSTRESSES/EID`.

## Local test intent

Use GPSTRESS examples to verify:

- MYSTRAN accepts MSC-style `OUTPUT(POST)` decks.
- OP2 contains readable `OGS1`.
- `pyNastran` can read `grid_point_surface_stresses`.
- Ordinary `OES1X1` element stress remains available in the same OP2.

Use `PARAM,STR_CID,0` ordinary element stress tests for SAP2000 `2-001` global element-stress comparison.

Use `GPSTRESS` tests for future MSC-style grid-point stress recovery work.
