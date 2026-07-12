# NEU Writer Map for `v18.00.a`

This note maps the current FEMAP-neutral (`.neu`) output path in the `v18.00.a` branch.

It focuses on:

- file/block structure
- where nodal vectors are written
- where element vectors are written
- how `MODES` and `BUCKLING` are mapped
- what is and is not implemented today

## Scope

This branch currently implements a **NEU writer**, not a general FEMAP-neutral reader/importer.

There is no general-purpose parser that reads arbitrary `.neu` files back into MYSTRAN.

The neutral path is output-only:

- MYSTRAN computes results
- `LINK9` orchestrates set/vector emission
- helper routines assemble text lines
- the `.neu` file is written in FEMAP-style neutral blocks

## Core Files

Main orchestration:

- [LINK9.f90](D:/18a/MYSTRAN/Source/LK9/LINK9/LINK9.f90)

Neutral record helper layer:

- [FEMAP_NEU_WRITE_HELPERS.f90](D:/18a/MYSTRAN/Source/Modules/FEMAP_NEU_WRITE_HELPERS.f90)

Nodal vector writer:

- [WRITE_FEMAP_GRID_VECS.f90](D:/18a/MYSTRAN/Source/LK9/L91/WRITE_FEMAP_GRID_VECS.f90)

Element engineering force writer:

- [WRITE_FEMAP_ELFO_VECS.f90](D:/18a/MYSTRAN/Source/LK9/L91/WRITE_FEMAP_ELFO_VECS.f90)

Element stress writer:

- [WRITE_FEMAP_STRE_VECS.f90](D:/18a/MYSTRAN/Source/LK9/L91/WRITE_FEMAP_STRE_VECS.f90)

Element strain writer:

- [WRITE_FEMAP_STRN_VECS.f90](D:/18a/MYSTRAN/Source/LK9/L91/WRITE_FEMAP_STRN_VECS.f90)

Element-family preprocessors that populate FEMAP arrays:

- [OFP3_ELFE_1D.f90](D:/18a/MYSTRAN/Source/LK9/L92/OFP3_ELFE_1D.f90)
- [OFP3_ELFE_2D.f90](D:/18a/MYSTRAN/Source/LK9/L92/OFP3_ELFE_2D.f90)
- [OFP3_STRE_NO_PCOMP.f90](D:/18a/MYSTRAN/Source/LK9/L92/OFP3_STRE_NO_PCOMP.f90)
- [OFP3_STRE_PCOMP.f90](D:/18a/MYSTRAN/Source/LK9/L92/OFP3_STRE_PCOMP.f90)
- [OFP3_STRN_NO_PCOMP.f90](D:/18a/MYSTRAN/Source/LK9/L92/OFP3_STRN_NO_PCOMP.f90)
- [OFP3_STRN_PCOMP.f90](D:/18a/MYSTRAN/Source/LK9/L92/OFP3_STRN_PCOMP.f90)

Shared FEMAP staging arrays:

- `FEMAP_EL_NUMS`
- `FEMAP_EL_VECS`

These are allocated/deallocated through:

- [ALLOCATE_FEMAP_DATA.f90](D:/18a/MYSTRAN/Source/LK9/LINK9/ALLOCATE_FEMAP_DATA.f90)
- [DEALLOCATE_FEMAP_DATA.f90](D:/18a/MYSTRAN/Source/LK9/LINK9/DEALLOCATE_FEMAP_DATA.f90)

## File Header and Version

The FEMAP neutral version number written into block `100` comes from:

- [SCONTR.f90](D:/18a/MYSTRAN/Source/Modules/SCONTR.f90)

Current setting:

- `FEMAP_VERSION = 9.0`

That value is written in `LINK9` here:

- `WRITE(NEU_LINE,9013) FEMAP_VERSION`

in [LINK9.f90](D:/18a/MYSTRAN/Source/LK9/LINK9/LINK9.f90)

This was aligned to the neutral reference samples under:

- [reference_msc](D:/18a/MYSTRAN_Validation-main/reference_msc)

whose `.neu` block `100` headers use version `9.,`

## Block Structure

The helper layer emits a small set of recurring records:

- `NEU_WRITE_BLOCK_END`
  - writes `-1`
- `NEU_WRITE_SET_ID`
  - writes one set id line
- `NEU_WRITE_SET_VEC_HEADER`
  - writes `set_id, vec_id, 1`
- `NEU_WRITE_ANALYSIS_IDS`
  - writes `FROM_PROG, ANAL_TYPE`
- `NEU_WRITE_TRIPLE_REAL`
  - writes `min, max, abs`
- `NEU_WRITE_TEN_IDS`
  - writes a 10-entry vector-id/metadata block
- `NEU_WRITE_GRID_RANGE`
  - writes `grid_min, grid_max, ...`
- `NEU_WRITE_GRID_VALUE`
  - writes `id, value`
- `NEU_WRITE_ELEM_RANGE`
  - writes element range plus output/entity type
- `NEU_WRITE_ELEM_FLAGS`
  - writes calculation/centroid/component flags
- `NEU_WRITE_VECTOR_END`
  - writes vector terminator

At the file/set level the structure is:

1. block `100` file header
2. optional geometry snapshot:
   - block `450` set header
   - block `451` vectors
3. for each result set:
   - block `450` set header
   - block `451` vectors

## Analysis-Type Mapping

Analysis-type mapping is set in:

- `GET_FEMAP_ANAL_TYPE`

inside [LINK9.f90](D:/18a/MYSTRAN/Source/LK9/LINK9/LINK9.f90)

Current mapping:

- `STATICS -> 1`
- `NLSTATIC -> 10`
- `MODES -> 2`
- `MFREQ -> 5`
- `BUCKLING -> 7`
- `DFREQ -> 4`
- fallback/unknown -> `0`

`FEMAP_FROM_PROG` currently stays `0`.

## Set-ID Mapping

The result-set loop is the `JVEC` loop in `LINK9`.

Current `FEMAP_SET_ID` policy:

- `STATICS` / `NLSTATIC`
  - `FEMAP_SET_ID = SCNUM(JVEC)`
- `MODES`
  - `FEMAP_SET_ID = JVEC`
- `BUCKLING`, `LOAD_ISTEP = 1`
  - `FEMAP_SET_ID = SCNUM(JVEC)`
- `BUCKLING`, `LOAD_ISTEP = 2`
  - `FEMAP_SET_ID = JVEC`
- `GEN CB MODEL`
  - `FEMAP_SET_ID = JVEC`

For `MODES` and buckling eigen output, `INT_SC_NUM` is still traced back to the owning subcase when possible via:

- `MODE_SUBCASE`
- `NUM_EIGENS_SUB`
- `IS_BUCKLING_SUBCASE`

## Geometry Snapshot

Geometry snapshot is handled by:

- `WRITE_FEMAP_GEOM_GRID_SNAPSHOT`

in [LINK9.f90](D:/18a/MYSTRAN/Source/LK9/LINK9/LINK9.f90)

This is **phase-1 geometry only**, not a full FEMAP database export.

Current geometry vectors:

- `X basic coordinate`
- `Y basic coordinate`
- `Z basic coordinate`
- `Input coordinate system (CP)`
- `Output coordinate system (CD)`
- `Permanent SPC code`

Implementation helpers:

- `WRITE_ONE_GEOM_GRID_VEC`
- `WRITE_ONE_GEOM_GRID_INT_VEC`

Snapshot set id:

- `90000001`

## Nodal Result Output

Nodal vectors are written by:

- [WRITE_FEMAP_GRID_VECS.f90](D:/18a/MYSTRAN/Source/LK9/L91/WRITE_FEMAP_GRID_VECS.f90)

Families currently mapped:

- `DISP`
- `OLOA`
- `SPCF`
- `MPCF`

Dispatch comes from `LINK9` through:

- `OFP1` for `DISP` / `OLOAD`
- `OFP2` for `SPCF` / `MPCF`

For each family, the writer emits eight vectors:

- total translation
- `T1`
- `T2`
- `T3`
- total rotation
- `R1`
- `R2`
- `R3`

Special behavior for `DISP`:

- translations/rotations are transformed to basic coordinates if the grid has a nonzero output coordinate system

Vector-id offsets:

- `DISP` -> base `0`
- `OLOA` -> base `20000`
- `SPCF` -> base `30000`
- `MPCF` -> base `40000`

## Element Output Data Flow

Element result flow is:

1. `LINK9` decides the active set and requests output
2. `OFP3*` computes element quantities
3. `OFP3*` populates:
   - `FEMAP_EL_NUMS`
   - `FEMAP_EL_VECS`
4. family-specific writer emits one NEU vector per result component

The actual vector emission for element results is centralized in:

- `NEU_WRITE_ELEM_VECTOR`

in [FEMAP_NEU_WRITE_HELPERS.f90](D:/18a/MYSTRAN/Source/Modules/FEMAP_NEU_WRITE_HELPERS.f90)

That helper writes:

- set/vector header
- titles
- min/max/abs
- 20 ids
- element range
- element flags
- one `element_id, value` row per result row
- vector terminator

## Element Engineering Forces

Writer:

- [WRITE_FEMAP_ELFO_VECS.f90](D:/18a/MYSTRAN/Source/LK9/L91/WRITE_FEMAP_ELFO_VECS.f90)

Precompute paths:

- [OFP3_ELFE_1D.f90](D:/18a/MYSTRAN/Source/LK9/L92/OFP3_ELFE_1D.f90)
- [OFP3_ELFE_2D.f90](D:/18a/MYSTRAN/Source/LK9/L92/OFP3_ELFE_2D.f90)

Currently mapped families:

- `ROD`
- `BAR`
- `BEAM`
- `TRIA3K`
- `TRIA3`
- `QUAD4K`
- `QUAD4`
- `SHEAR`
- `ELAS1`
- `ELAS2`
- `ELAS3`
- `ELAS4`
- `BUSH`

Current `ELFO` vector-id offsets:

- `ROD` -> `50100`
- `BAR` -> `50200`
- `TRIA3K` -> `50300`
- `TRIA3` -> `50400`
- `QUAD4K` -> `50500`
- `QUAD4` -> `50600`
- `SHEAR` -> `50700`
- `ELAS1` -> `50800`
- `ELAS2` -> `50900`
- `ELAS3` -> `51000`
- `ELAS4` -> `51100`
- `BUSH` -> `51200`
- `BEAM` -> `51300`

Notes:

- `BEAM` is written as 12 component vectors for end/station-style force result channels
- `BAR`/`ROD` duplicate some component structure to match the neutral writer’s expected paired channels

## Element Stress

Writer:

- [WRITE_FEMAP_STRE_VECS.f90](D:/18a/MYSTRAN/Source/LK9/L91/WRITE_FEMAP_STRE_VECS.f90)

Precompute paths:

- [OFP3_STRE_NO_PCOMP.f90](D:/18a/MYSTRAN/Source/LK9/L92/OFP3_STRE_NO_PCOMP.f90)
- [OFP3_STRE_PCOMP.f90](D:/18a/MYSTRAN/Source/LK9/L92/OFP3_STRE_PCOMP.f90)

Currently mapped families:

- `ELAS`
- `ROD`
- `BAR`
- `BEAM`
- `TRIA3K`
- `TRIA3`
- `QUAD4K`
- `QUAD4`
- `HEXA8`
- `HEXA20`
- `PENTA6`
- `PENTA15`
- `TETRA4`
- `TETRA10`
- `SHEAR`

Current `STRE` vector-id offsets:

- `ELAS*` -> `60100`
- `ROD` -> `60200`
- `BAR` -> `60300`
- `TRIA3K` -> `60400`
- `TRIA3` -> `60500`
- `QUAD4K` -> `60600`
- `QUAD4` -> `60700`
- `HEXA8` -> `60800`
- `HEXA20` -> `60900`
- `PENTA6` -> `61000`
- `PENTA15` -> `61100`
- `TETRA4` -> `61200`
- `TETRA10` -> `61300`
- `SHEAR` -> `61400`
- `BEAM` -> `61500`

Notes:

- shells/non-PCOMP are expanded into top/bottom component vectors
- solids are written as 12 vectors
- if `STRE_OPT == 'VONMISES'`, the final scalar channel differs from octahedral mode

## Element Strain

Writer:

- [WRITE_FEMAP_STRN_VECS.f90](D:/18a/MYSTRAN/Source/LK9/L91/WRITE_FEMAP_STRN_VECS.f90)

Precompute paths:

- [OFP3_STRN_NO_PCOMP.f90](D:/18a/MYSTRAN/Source/LK9/L92/OFP3_STRN_NO_PCOMP.f90)
- [OFP3_STRN_PCOMP.f90](D:/18a/MYSTRAN/Source/LK9/L92/OFP3_STRN_PCOMP.f90)

Currently mapped families:

- `TRIA3K`
- `TRIA3`
- `QUAD4K`
- `QUAD4`
- `HEXA8`
- `HEXA20`
- `PENTA6`
- `PENTA15`
- `TETRA4`
- `TETRA10`
- `SHEAR`
- `BEAM`

Current `STRN` vector-id offsets:

- `TRIA3K` -> `70400`
- `TRIA3` -> `70500`
- `QUAD4K` -> `70600`
- `QUAD4` -> `70700`
- `HEXA8` -> `70800`
- `HEXA20` -> `70900`
- `PENTA6` -> `71000`
- `PENTA15` -> `71100`
- `TETRA4` -> `71200`
- `TETRA10` -> `71300`
- `SHEAR` -> `71400`
- `BEAM` -> `71500`

Notes:

- shell strain follows the same general component expansion pattern as shell stress
- `BEAM` strain gets 12 component vectors
- solids again use the 12-channel family

## Modes and Eigenvectors

`MODES` uses:

- `FEMAP_ANAL_TYPE = 2`
- `FEMAP_SET_ID = JVEC`

The displacement/eigenvector data is read from `UG_COL` in the `JVEC` loop and then pushed through the same nodal writer path as displacements.

So in NEU terms:

- there is no separate “eigen summary neutral block”
- the mode shape is represented as a result set with grid vectors

## Buckling

`BUCKLING` uses:

- `FEMAP_ANAL_TYPE = 7`

Two phases are handled differently:

1. `LOAD_ISTEP = 1`
   - preload/static part
   - `FEMAP_SET_ID = SCNUM(JVEC)`
2. `LOAD_ISTEP = 2`
   - eigen buckling part
   - `FEMAP_SET_ID = JVEC`

During the eigen buckling phase, `LINK9` backtracks the owning subcase using:

- `NUM_EIGENS_SUB`
- `IS_BUCKLING_SUBCASE`
- fallback `MODE_SUBCASE`

This lets result titles and output-family decisions still refer to the correct subcase context while writing per-mode NEU sets.

## Family Gating

`LINK9` currently derives these booleans:

- `WRITE_NEU_GEOM`
- `WRITE_NEU_DISP`
- `WRITE_NEU_OLOA`
- `WRITE_NEU_SPCF`
- `WRITE_NEU_MPCF`
- `WRITE_NEU_ELFO`
- `WRITE_NEU_STRE`
- `WRITE_NEU_STRN`

These are set from:

- `PRTNEU`
- whether any requests exist for the corresponding family

This is the current family-scoped routing base for NEU output in `OUTMODE=SMART`.

## Current Boundaries and Gaps

What is implemented:

- FEMAP-neutral file header and per-set/per-vector writing
- geometry snapshot phase-1
- nodal vectors for displacement/load/constraint-force families
- major element force/stress/strain families
- mode and buckling results as result sets

What is not claimed here:

- full FEMAP-neutral model import/export parity
- general-purpose NEU parsing into MYSTRAN
- full import-grade geometry/database blocks beyond the current snapshot approach

Known caveats:

- some `PCOMP`/composite paths still have partial or warning-only behavior
- several element writers still reflect historical MYSTRAN channel naming/ordering choices, not yet a full audited match to every reference neutral file
- `FEMAP_FROM_PROG` still writes as `0`

## Next Useful Audit Targets

If the goal is to move closer to the reference neutral files under `reference_msc`, the highest-value next checks are:

1. verify block `450/451` field semantics against the new FEMAP v9 references
2. audit `ID(20)` metadata per vector family
3. audit `OUT_TYPE`, `ENT_TYPE`, `CALC_WARN`, `COMP_DIR`, and `CENT_TOTAL`
4. compare shell stress/strain component ordering against the v9 reference files
5. compare buckling and modes set headers to confirm titles and set numbering are what FEMAP expects
