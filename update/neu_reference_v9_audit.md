## FEMAP Neutral v9 Reference Audit

This note maps the current `mystran18a` NEU writer against the text neutral files stored in:

- `D:\18a\MYSTRAN_Validation-main\reference_msc`
- `D:\18a\MYSTRAN_Validation-main\reference_msc\Benchmark`
- `D:\18a\MYSTRAN_Validation-main\reference_msc\vic\...`

The reference files appear to be FEMAP text neutral exports produced after importing BDF+OP2, not raw Nastran output files.

### Current branch status

- NEU header version has been updated to `9.0` in [SCONTR.f90](D:\18a\MYSTRAN\Source\Modules\SCONTR.f90:579).
- The writer remains an output-only path. There is no general NEU reader/parser in this branch.
- Core orchestration is still in [LINK9.f90](D:\18a\MYSTRAN\Source\LK9\LINK9\LINK9.f90).
- Low-level line assembly remains in [FEMAP_NEU_WRITE_HELPERS.f90](D:\18a\MYSTRAN\Source\Modules\FEMAP_NEU_WRITE_HELPERS.f90).

### Reference file structure observed

Typical reference files use:

1. Block `100`
   - Contains `<NULL>` and version line `9.,`

2. Block `450`
   - Defines output sets
   - Each set has:
     - set id
     - title text
     - analysis metadata pair such as `14, 2` or `14, 7`
     - one real value line
     - one integer line
     - `From:` source OP2 path
     - `Date : ...`
     - `<NULL>`
     - deck title / subtitle text

3. Block `451`
   - Defines vectors under the current set
   - Used for:
     - nodal displacement/load/constraint vectors
     - modal participation / effective-mass vectors
     - element force/stress/strain vectors

### Observed examples

#### Modes

From `103_one_subcase.neu`:

- block `450` contains repeated sets titled like:
  - `Mode 1, 1246.54 Hz`
  - `Mode 2, 1246.54 Hz`
- metadata line is:
  - `14,        2,`
- one real line contains frequency
- source/date/title trailer is present for every set

#### Buckling

From `105_one_statsub.neu`:

- preload/static set:
  - title `NASTRAN Case 1`
  - metadata `14,        1,`
- eigenvalue result sets:
  - titles like `Eigenvalue 1 1491.577`
  - metadata `14,        7,`

This matches the idea that buckling files contain a static preload-type set plus eigenvalue-style result sets.

### What the current MYSTRAN writer already does correctly

- Writes block `100`
- Writes block `450` and `451`
- Uses per-result set ids
- Writes per-vector min/max/abs statistics
- Writes sorted grid-based vectors
- Writes element vectors using family-specific vector ids
- Distinguishes analysis family in `GET_FEMAP_ANAL_TYPE`
- Supports statics, normal modes, MFREQ, DFREQ, buckling, and nonlinear statics at the set-metadata level

Recent local improvement on top of that:

- `MODES`
  - block `450` title now uses the pattern `Mode n, <freq> Hz`
  - block `450` scalar line now carries the modal frequency computed from the eigenvalue

- `BUCKLING` eigen step (`LOAD_ISTEP = 2`)
  - block `450` title now uses the pattern `Eigenvalue n <value>`
  - block `450` scalar line now carries the buckling eigenvalue

This removes the previous all-zero placeholder behavior for these solution families while leaving the rest of the writer design intact.

### What is still intentionally custom in MYSTRAN

#### Geometry snapshot

The current branch writes a phase-1 custom geometry snapshot:

- set id `90000001`
- vectors:
  - `91001` X basic
  - `91002` Y basic
  - `91003` Z basic
  - `91004` CP
  - `91005` CD
  - `91006` permanent SPC

This is emitted by `WRITE_FEMAP_GEOM_GRID_SNAPSHOT` in [LINK9.f90](D:\18a\MYSTRAN\Source\LK9\LINK9\LINK9.f90:1586).

This is useful and parseable, but it is not trying to replicate a full FEMAP model-neutral geometry export.

#### Set metadata trailer

Reference files include:

- `From: <path to op2>`
- `Date : ...`
- `<NULL>`
- one or more deck title lines

Current MYSTRAN set headers are much lighter. They do not yet attempt to mirror the reference export trailer.

Remaining gap in block `450`:

- no `From: <op2 path>`
- no `Date : ...`
- no trailing `<NULL>` text line per set
- no repeated deck title/subtitle trail after the primary set title

#### Program id

Current writer uses `FEMAP_FROM_PROG = 0` in `LINK9`.

If we want the file to identify itself more like Nastran-oriented neutral output, this field should be reviewed later. For now it is safer to keep it conservative than to claim a program identity that the file semantics do not fully match yet.

### Gap list for future work

1. Add optional richer block-450 metadata
   - source OP2 path
   - export date/time
   - title/subtitle trail

2. Audit whether current `FEMAP_ANAL_TYPE` values align with the reference export conventions for:
   - statics
   - modes
   - buckling preload sets
   - buckling eigenvalue sets

3. Compare block-451 vector headers line-by-line against reference files
   - nodal vectors
   - element force vectors
   - element stress vectors
   - element strain vectors
   - modal/effective-mass vectors

4. Decide whether MYSTRAN should keep the custom geometry snapshot permanently
   - recommended for now: keep it
   - it is useful, cheap, and does not block result import

### Practical conclusion

For the current branch, the important compatibility move was:

- switch the neutral header to version `9.0`
- stop writing a generic zero-only set title/scalar for `MODES` and buckling-eigen sets

That brings the file family in line with the reference neutral examples. The writer is still structurally simpler than the reference FEMAP exports, especially in block `450`, but it is now easier to continue toward closer text-level similarity without discarding the existing phase-1 geometry snapshot design.
