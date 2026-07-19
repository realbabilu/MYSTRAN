# GPSTRESS Recovery Design

## Scope

This note defines a staged implementation path to move local `GPSTRESS/GSTRESS`
from the current baseline compatibility writer toward MSC-style grid-point
surface stress recovery.

Current branch status:

- `GPSTRESS/GSTRESS` is accepted in Case Control.
- `OUTPUT(POST)`, `SET`, `SURFACE`, and `VOLUME` are accepted enough to keep
  MSC/NX-style decks parseable.
- `OGS1` is written in OP2.
- The present `OGS1` payload is derived from shell corner stress rows already
  stored in `OGEL`.
- `OUTPUT(POST)` phase-1 parsing now accepts MSC/NX-style post blocks such as
  `SET 1 ALL` and `SURFACE 100 SET 1 NORMAL Z` through a local GPSTRESS
  registry path, without routing those cards into the legacy general Case
  Control `SET` parser.

This is not yet true grid-point recovery/averaging.

## Current blocker

As of July 18, 2026, the `OUTPUT(POST)` parser path is no longer the limiting
issue for the shell GPSTRESS patch-test deck:

- `D:\18a\MYSTRAN_Validation-main\ctriar_cquadr\2_001_4_mystran_gpstress.dat`
  now runs to normal termination in MYSTRAN.
- The old `*ERROR 1405: SET ID 1 NOT FOUND` issue is resolved by the local
  post-set registry path.
- The remaining blocker is the OP2 `OGS1` table layout itself:
  `pyNastran` still aborts while entering `OGS1` with:
  `struct.error: unpack requires a buffer of 4 bytes`

This means the next work item is specifically:

- audit `WRITE_OGS1_SURFACE_STRESS` table-3/table-4 record structure against a
  known-good MSC/NX `OGS1`
- fix the binary record framing before spending more time on recovery weights
  or coordinate-system refinement

In short:

- parser path: working
- recovery patch collection: working enough to execute
- OP2 `OGS1` framing: still broken

## 2026-07-18 local follow-up: CQUAD4 + MITC4+ GPSTRESS

Local deck:

- `D:\18a\MYSTRAN_Validation-main\ctriar_cquadr\2_001_4_mitc4p_gpstress.dat`

Configuration:

- `CQUAD4`
- `PARAM,QUAD4TYP,MITC4+`
- `GPSTRESS=ALL`
- `OUTPUT(POST)` with `SURFACE 100 SET 1 NORMAL Z`
- `PARAM,STR_CID,` to recover in the basic/global system

Observed result:

- deck runs to normal termination
- F06 prints the expected:
  `S T R E S S E S   A T   G R I D   P O I N T S   - -     S U R F A C E     100`
- the patch-test membrane case is uniform at the grid-point surface output
- the plate-bending case is also recovered cleanly, with consistent `Z1/Z2/MID`
  rows and no obvious element-local scrambling

Practical conclusion:

- `CQUAD4 + MITC4+` now appears to support the current F06 `GPSTRESS` recovery
  path correctly for the local `2_001` patch-test deck
- this puts `CQUAD4(MITC4+)` in the same "F06 path looks healthy" bucket as the
  earlier `CQUADR` and `CTRIAR` checks
- this does **not** mean OP2 `OGS1` is fixed; the binary `OGS1` framing issue
  remains a separate blocker

## 2026-07-18 follow-up: `OGS1` table-3/table-4 framing recheck

Re-audit target:

- compare `WRITE_OGS1_SURFACE_STRESS` framing against a known-good MSC/NX
  `OGS1`
- then re-test `pyNastran` on:
  `D:\18a\MYSTRAN_Validation-main\ctriar_cquadr\2_001_4_mystran_gpstress.OP2`

Observed status on the current branch:

- the `OGS1` table is present in both MSC and MYSTRAN files
- `pyNastran` now opens both:
  - `2_001_4_msc_gpstress.OP2`
  - `2_001_4_mystran_gpstress.OP2`
- `pyNastran` reports `table_names` containing:
  - `OUGV1`
  - `OES1X1`
  - `OGS1`
  for both subcases
- `pyNastran` also populates `grid_point_surface_stresses` for the MYSTRAN
  file with two subcase entries

Low-level framing check around the `OGS1` open sequence:

- the marker/subtable pattern around `OGS1` matches the MSC structure at the
  level of:
  - table open marker
  - `WRITE_ITABLE(-3/-4/...)` style control records
  - table-3 header record start
- the visible binary difference in the early header chunk is only the run-date
  field, not a framing mismatch

Practical conclusion:

- for the currently generated local `2_001_4_mystran_gpstress.OP2`, the older
  `struct.error: unpack requires a buffer of 4 bytes` blocker is no longer
  reproduced
- the immediate table-3/table-4 framing path in
  `WRITE_OGS1_SURFACE_STRESS` is now good enough for `pyNastran` to step
  through `OGS1`
- remaining differences versus MSC are now at the result-content level, not at
  the gross binary-framing level

Important content difference still visible:

- MSC `OGS1` for this deck yields `(1, 40, 8)` and `(1, 40, 8)` data blocks
  for the two subcases
- MYSTRAN `OGS1` yields `(1, 16, 8)` and `(1, 16, 8)`

Interpretation:

- the current MYSTRAN `OGS1` path writes a smaller surface-recovery result set
  than MSC for this patch-test deck
- that is now a recovery/output-content parity issue, not a broken OP2 record
  framing issue

### Root cause of the current row-count mismatch

For `2_001_4_*_gpstress`:

- MSC `OGS1` writes `40` rows per subcase
- MYSTRAN writes `16` rows per subcase

The pattern shows why:

- MSC stores one `Z1` and one `Z2` row for every contributing
  `grid-element` corner pair on the surface
- MYSTRAN currently stores one `Z1` and one `Z2` row per unique grid on the
  surface recovery patch

That collapse happens explicitly in the current recovery kernel:

- `WRITE_ELEM_STRESSES.f90`
- `BUILD_SURFACE_PATCH_OUTPUT`
- current line of interest:
  `IF (FIND_INT(GID, OUT_GRIDS, NOUT) > 0) CYCLE`

So the present behavior is:

- repeated corner contributions for the same surface grid are merged into one
  averaged grid result

While the current MSC reference behavior for this deck is:

- keep the separate element-corner contributions in `OGS1`

This is why the binary result is now readable but still not parity-matched.

### Practical next decision

There are now two viable directions:

1. Keep the current MYSTRAN design:
   - `GPSTRESS` in `OGS1` means one recovered/averaged result per surface grid
   - cleaner as a recovery product
   - not byte/shape-compatible with MSC for this deck

2. Move toward MSC parity:
   - write one `OGS1` row per contributing `grid-element-fiber` corner pair
   - preserve repeated grids across adjacent elements
   - optionally keep a later averaged/F06-only view as a separate concept

If MSC compatibility is the priority, the next code change should be:

- remove the unique-grid collapse in `BUILD_SURFACE_PATCH_OUTPUT`
- expand the output arrays to carry one entry per contributing surface
  corner/sample rather than one entry per recovered grid

## Key assumption

For the first real recovery pass, assume the `SURFACE` request is explicitly
given by the deck and use that as the recovery patch definition.

That means:

- do not build free-form global patches from the whole mesh
- do not average all adjacent shell elements at a grid by default
- only average contributions from shell elements that belong to the requested
  `SURFACE`

This is closer to MSC intent and avoids accidental mixing across folds, sharp
corners, or unrelated panel groups.

## Behavioral target

For shell-only phase 1:

- `GPSTRESS/GSTRESS` requests grid-point surface stress recovery
- the recovery domain is the requested `SURFACE`
- each grid-point result is recovered from contributing shell elements on that
  surface
- all contributions are transformed into one target stress coordinate system
  before averaging
- `Z1` and `Z2` stay separate
- output is written to `OGS1`
- ordinary element stress output remains in `OES1X1`

## Phase 1 recovery rule

For each requested surface and each grid on that surface:

1. collect the shell elements on the surface that contain the grid
2. for each contributing element, recover the shell stress contribution at the
   matching corner
3. transform the contribution to the target stress system
4. average the contributions using nodal tributary area weighting
5. write one recovered `Z1` row and one recovered `Z2` row for that
   grid/surface pair

This is intentionally more sophisticated than a plain arithmetic mean, but much
 simpler than full multi-family MSC recovery.

## Why not use a plain mean

Plain averaging is too weak because it:

- gives the same influence to large and small adjacent elements
- overreacts to distorted mesh topology
- mixes good and poor patches equally
- is more likely to drift from MSC behavior

Phase 1 should therefore use:

- shell nodal tributary area weighting

and not:

- simple unweighted mean

## Weighting rule

For shell elements in a requested surface:

- `QUAD4/CQUADR`: each corner gets `element_area / 4`
- `TRIA3/CTRIAR`: each corner gets `element_area / 3`
- `QUAD8/TRIA6`: phase 1 may either:
  - use the same parent corner partition rule if only corner recovery is
    retained, or
  - be deferred until the surface membership and corner mapping path is stable

The weighting is applied independently for `Z1` and `Z2`.

Recovered component formula:

`sigma_gp = sum(w_i * sigma_i) / sum(w_i)`

for each requested component.

## Coordinate-system rule

All element contributions must be transformed before averaging.

Phase 1 target:

- honor explicit `PARAM,STR_CID`
- if the local GPSTRESS compatibility alias has defaulted `STR_CID` to `0`,
  average in the basic/global system

Do not average values that are still in per-element local coordinates when
element axes differ.

## Surface-driven grouping

The requested `SURFACE` should be treated as the primary recovery container.

Internal grouping for phase 1:

- group by `surface_id`
- group by `subcase`
- group by `grid_id`
- group by `fiber_side` (`Z1`, `Z2`)

Do not merge across different surfaces, even when they share the same grid.

## Later filters for phase 2

Once phase 1 is stable, add optional filters to reduce over-averaging:

- normal-angle break
- property break
- material break
- laminate/ply break

Suggested defaults for a later pass:

- only average shell contributions whose midsurface normals differ by less than
  a threshold such as `15 deg`

These filters are not required for the first working implementation.

## Data needed from Case Control

The current branch counts `SURFACE` and `VOLUME`, but does not yet retain
enough structured information to drive recovery.

Phase 1 should store:

- `SURFACE_ID`
- referenced `SET_ID`
- direction/orientation token from the `SURFACE` line
- resolved list of element IDs in the surface

At minimum, the parser must retain enough information to build:

- `surface -> element_ids`
- `surface -> grid_ids`

## Proposed source additions

### 1. Case Control storage

Primary files:

- `Source/Modules/CC_OUTPUT_DESCRIBERS.f90`
- `Source/LK1/L1A/LOADC.f90`

Add structured storage for parsed `OUTPUT(POST)` recovery definitions instead
of only counting surfaces/volumes.

Suggested new arrays/modules:

- `GP_SURFACE_IDS(:)`
- `GP_SURFACE_SETIDS(:)`
- `GP_SURFACE_NORMAL_MODE(:)`

If dynamic allocation is inconvenient in phase 1, use fixed-size arrays sized
from a conservative branch parameter and raise a clear fatal if exceeded.

### 2. Surface resolution helper

Add a new helper to resolve the requested `SET` into actual shell elements and
their corner grids.

Suggested new utility:

- `Source/LK9/L91/BUILD_GPSTRESS_SURFACE_PATCH.f90`

Responsibilities:

- read a surface definition
- resolve the referenced Case Control set
- retain only supported shell element families
- build the grid-to-contributing-element adjacency for that surface

### 3. Recovery kernel

Add a shell-specific recovery routine.

Suggested new utility:

- `Source/LK9/L91/RECOVER_GPSTRESS_SHELL_SURFACE.f90`

Responsibilities:

- input:
  - requested surface definition
  - shell family
  - `OGEL`, `EID_OUT_ARRAY`, `GID_OUT_ARRAY`
  - geometry needed for area weighting
- output:
  - recovered grid-point rows ready for F06/OP2 `OGS1`

Suggested recovered row fields:

- `grid_id`
- `element_or_surface_ref`
- `fiber_side`
- `oxx`
- `oyy`
- `txy`
- `angle`
- `major`
- `minor`
- `ovm_or_maxshear`

### 4. OGS1 writer hook

Primary file:

- `Source/LK9/L91/WRITE_ELEM_STRESSES.f90`

Change `WRITE_OGS1_SURFACE_STRESS` so it no longer writes direct element-corner
rows from `OGEL` for GPSTRESS.

Instead:

- call the recovery helper
- write recovered rows to `OGS1`

Important:

- keep `OES1X1` unchanged
- keep the current baseline path available behind a temporary fallback switch if
  needed during transition

## Reuse opportunities already present in branch

The branch already contains useful patterns:

- `CTETRA4S_SMOOTH_ASSEMBLY.f90`
  - shows a patch-based smoothing architecture
  - demonstrates weighted aggregation using topology and geometry
- `SURFACE_FIT.f90`
  - may be useful later if polynomial smoothing is needed
  - not required for phase 1 averaging

The shell GPSTRESS path should reuse the patch-assembly mindset from
`CTETRA4S_SMOOTH_ASSEMBLY`, but with shell-specific area weighting and
surface-defined patches.

## Phase 1 restrictions

To keep the implementation safe:

- support shells first
- support `TRIA3`, `CTRIAR`, `QUAD4`, `CQUADR` first
- keep solids out of phase 1
- keep composite averaging out of phase 1
- keep free-form global adjacency averaging out of phase 1

If unsupported element families are present in the requested surface:

- either skip them with warning, or
- issue a controlled warning that recovery falls back to the existing baseline

## Validation plan

### Stage A: no-change control

Use the current GPSTRESS patch deck and confirm:

- `OES1X1` remains unchanged
- `OGS1` still exists and is readable by `pyNastran`

### Stage B: flat shell patch

Use a flat multi-element shell patch where all contributions should match.

Expected:

- recovered GPSTRESS equals the common patch stress
- no drift from weighting scheme

### Stage C: mixed mesh density on one surface

Use adjacent large/small shell elements on one flat surface.

Expected:

- area-weighted average is stable
- better behavior than plain arithmetic mean

### Stage D: folded surface / two panels sharing a grid

Use two different surfaces that share a grid but have different normals.

Expected:

- no cross-surface averaging
- each surface writes its own recovered grid-point stress

## Future MSC-style extensions

After phase 1:

- normal-angle break inside a single surface
- property/material filters
- laminate-aware recovery
- solid `GPSTRESS` / `VOLUME` support
- richer F06 presentation mirroring MSC wording more closely

## Recommended implementation order

1. retain parsed `SURFACE` definitions, not just counts
2. build `surface -> element/grid` adjacency
3. implement shell area-weighted recovery
4. switch `OGS1` writer from direct-corner copy to recovered rows
5. validate against simple shell patch decks
6. add optional break filters later

## Bottom line

If the deck already provides the `SURFACE`, MYSTRAN should treat that surface
as the recovery patch and compute grid-point stress from surface-local shell
contributors using a weighted average. That is the cleanest first step toward
MSC-style GPSTRESS and avoids the main failure mode of naive global averaging.

## Current implementation status

As of July 18, 2026:

- `SURFACE id SET setid [NORMAL dir]` is now parsed and retained structurally.
- shell surface membership can now be resolved from the referenced Case Control
  `SET`.
- `WRITE_OGS1_SURFACE_STRESS` now contains a first compiled recovery path for
  shell GPSTRESS output.

Current phase-1 behavior in code:

- active only when `GPSTRESS/GSTRESS` is requested and at least one `SURFACE`
  definition exists
- limited to `TRIA3/CTRIAR` and `QUAD4/CQUADR`
- uses area-weighted accumulation of `XX/YY/XY`
- recomputes principal angle, major/minor, max shear, and von Mises from the
  averaged in-plane stress state
- preserves the old direct-corner `OGS1` path as fallback when no recovered
  rows are formed

Still not done:

- numerical audit against `reference_msc`
- advanced normal-based filtering inside one surface
- multi-surface duplicate-grid handling beyond the current conservative skip
- solid/volume GPSTRESS recovery
## 2026-07-18 status update

- `WRITE_ELEM_STRESSES.f90` now emits an F06 `GPSTRESS` surface block styled after MSC/Nastran:
  - `STRESSES AT GRID POINTS -- SURFACE ...`
  - `SURFACE X-AXIS X  NORMAL(Z-AXIS) Z ... CID 0`
  - per-grid `Z1`, `Z2`, and `MID` rows
- The OP2 `OGS1` path remains the primary binary output path and was left intact.
- Current numerical behavior is still a MYSTRAN-specific baseline recovery:
  - shell corner results are transformed to basic/global first
  - per-surface grid values are then area-weighted averaged from adjoining shell contributions
- This matches the intended coordinate-system direction for MSC-style `GPSTRESS`, but it is not yet the same recovery algorithm as MSC/Nastran.
- Patch-test comparison on `prob_2_001_*` shows:
  - F06 table shape/header is now much closer to MSC
  - values still differ from MSC because MSC uses a more refined grid-point surface recovery than the current simple weighted averaging
- Recommended next step:
  - keep the new F06 layout
  - refine the recovery kernel before adding `GPAVFORCE`, so the future force recovery follows the same surface/global philosophy with a better averaging scheme

## 2026-07-18 patch-test note

Further audit on the `prob_2_001_*` shell patch test showed a second issue in
the shell output chain:

- `QUAD4/CQUADR` stress/strain corner recovery was not following the same
  `SHELL_STR_ANGLE` post-fit transform path already used by `QUAD8`.
- This was patched in:
  - `Source/LK9/L92/OFP3_STRE_NO_PCOMP.f90`
  - `Source/LK9/L92/OFP3_STRN_NO_PCOMP.f90`

Effect:

- the recovered `GPSTRESS` field moves closer to a consistent surface/global
  field
- but the patch-test result is still not fully uniform like MSC/Nastran

Interpretation:

- one missing transform step did exist and is now fixed
- however, the remaining difference indicates that MSC-style `GPSTRESS` still
  needs a stronger surface-coordinate recovery stage, not only nodal
  accumulation

So the current implementation status is:

1. shell post-fit angle handling is less inconsistent than before
2. F06 `GPSTRESS` table formatting is now present
3. exact patch-test-uniform MSC-style grid-point stress recovery is still an
   open item

## 2026-07-18 deeper diagnosis on shell patch tests

Additional audit of the shell quad path shows that the main philosophical gap
with MSC/Nastran is not just the final nodal averaging.

### What MYSTRAN currently does

For `QUAD4/CQUADR`, the current LK9 shell output chain is effectively:

1. evaluate shell stress quantities at element sampling points
2. pass the rows into `POLYNOM_FIT_STRE_STRN`
3. extrapolate each scalar component independently to corner locations
4. build `OGEL` rows from those extrapolated scalar components
5. later average or fit those already-separated component rows for GPSTRESS

Important detail:

- `POLYNOM_FIT_STRE_STRN` fits `XX`, `YY`, `XY`, etc. as separate scalar
  fields over the element surface
- it does not recover one coupled in-plane stress tensor field as a tensor
  object over the patch

### Why this matters for patch tests

On distorted shell patches such as `prob_2_001`, a patch-test-consistent
surface recovery should aim to preserve one uniform physical stress tensor in a
common surface/global basis.

If `XX`, `YY`, and `XY` are extrapolated independently element by element
before patch recovery:

- the principal invariants can still stay nearly constant
- but the component split can drift from element to element
- later nodal averaging cannot fully repair that drift

That matches the observed MYSTRAN pattern:

- principal magnitudes are close
- `XX/YY/XY` vary around the surface
- MSC GPSTRESS remains much more uniform for the same patch test

### Practical interpretation

For MSC-style shell `GPSTRESS`, the desired pipeline is closer to:

1. transform each contributing shell sample to one common surface/global basis
2. recover a patch tensor field in that common basis
3. evaluate the recovered tensor at grid points
4. only then derive `ANGLE`, `MAJOR`, `MINOR`, `MAX SHEAR`, `VON MISES`

and not:

1. extrapolate scalar `XX`, scalar `YY`, scalar `XY` separately per element
2. average those component outputs afterward

### Recommended next technical direction

The next stronger implementation should therefore avoid treating shell
`GPSTRESS` as a post-average of already-finalized `OGEL(XX,YY,XY)` rows.

Instead, it should operate on a patch-level recovered in-plane tensor field.

For phase 2, candidate paths are:

- recover directly from per-sample shell stress tensors before
  `SHELL_STRESS_OUTPUTS`
- or reconstruct a patch-consistent strain/stress field from surface-grid
  displacements when the surface is flat and the request is explicitly
  `GPSTRESS`

For the current branch, this means:

- keep the new F06 GPSTRESS table layout
- keep the baseline OGS1 writer path
- recognize that exact MSC-like patch-test uniformity will likely require a
  patch-level tensor recovery stage earlier than the current `OGEL`-based
  averaging layer

## 2026-07-18 follow-up: `STR_CID` ordering diagnosis

Using a debug clone of `ctriar_cquadr/2_001_4_mystran_gpstress.dat`, the
behavior splits cleanly by `STR_CID`:

- with `PARAM,STR_CID,-1`, both the regular shell stress table and the
  GPSTRESS surface table become uniform and match the expected patch-test
  tensor for `prob_2_001`
- with `PARAM,STR_CID,0`, both tables become nonuniform

This is a strong indication that the remaining `prob_2_001` mismatch is not
coming from the `CQUADR` recovery itself.  The element-level probe in
`ELEM_STRE_STRN_ARRAYS` shows that the raw membrane/bending quantities are
already uniform before the later shell output transformation logic.

Working interpretation:

1. the shell recovery path for `CQUADR` is healthy enough for this patch test
2. the order/meaning of later shell stress transformations under active
   `STR_CID` is still inconsistent for this case
3. GPSTRESS inherits that transform problem because it currently consumes the
   already-finalized shell output rows

This narrows the next real fix to the shell output transform ordering, not to
further averaging tweaks.

## 2026-07-18 fix applied: suppress premature shell `STR_CID` rotation

The shell patch-test failure for `prob_2_001` under `PARAM,STR_CID,0` was
traced to a premature stress/strain coordinate transform in
`ELEM_STRE_STRN_ARRAYS`.

For shell families (`TRIA3`, `QUAD4`, `QUADR`, `QUAD8`):

- the raw recovered shell tensors are still consumed later by shell-specific
  extrapolation/output paths
- those later paths expect the unreoriented shell recovery basis
- rotating the raw shell tensors early with `STR_TENSOR_TRANSFORM` corrupts the
  data seen by corner recovery and GPSTRESS surface output

Applied fix:

- shell 2D families no longer apply `STR_CID` rotation inside
  `ELEM_STRE_STRN_ARRAYS`
- solid-family `STR_CID` handling is unchanged

Observed result on `tmp_2_001_4_mystran_gpstress_debug.dat`:

- subcase 1 shell stresses became uniform again
- `GPSTRESS` surface output became uniform again
- subcase 2 bending patch output also became uniform and symmetric in the
  expected way (`Z1`/`Z2` split, `MID` average)

This fix restores patch-test-consistent shell/basic output for the investigated
`CQUADR` case and confirms that the earlier mismatch was a transform-ordering
problem, not a `CQUADR` stress-recovery problem.

### Local follow-up in `ctriar_cquadr`: `CTRIAR` check

A local non-reference deck was added for the corresponding triangular patch:

- `D:\18a\MYSTRAN_Validation-main\ctriar_cquadr\2_001_3_mystran_gpstress.dat`

Using the same `GPSTRESS` + `SURFACE 100` + `PARAM,STR_CID,0` path, the
generated F06 showed:

- uniform subcase 1 grid-point stresses at all surface grids
- uniform subcase 2 `MID` values and the expected symmetric `Z1/Z2` split

So the transform-ordering fix is not limited to `CQUADR`; it also restores the
same patch-test behavior for the `CTRIAR` path in this local validation set.

## 2026-07-19 follow-up: OGS1 row-count/framing fixed, payload semantics split from F06

The current status for the local `GPSTRESS` OP2 path is now:

- `D:\18a\MYSTRAN_Validation-main\ctriar_cquadr\2_001_4_mystran_gpstress.dat`
  runs to normal termination
- `pyNastran` can read `2_001_4_mystran_gpstress.OP2` without the earlier
  `struct.error: unpack requires a buffer of 4 bytes`
- `OGS1` now matches the MSC reference row structure for this deck:
  - same subcase keys
  - same `(1, 40, 8)` data shape per subcase
  - same repeated `(grid, element, fiber)` ordering

So the old blocker was resolved at the binary table framing / row-count level.

### Important semantic note

While comparing the recovered `OGS1` payload against
`2_001_4_msc_gpstress.OP2`, a second issue became clear:

- MSC `OGS1` raw payload does **not** line up with the generic
  `pyNastran` header labels beyond the first three direct components
  (`NX`, `NY`, `TXY`)
- in practice, the MSC reference OP2 carries values that behave like:
  - direct components first
  - then principal/invariant-style values
  - without a reliable raw `ANGLE` slot in the same position suggested by the
    generic `GridPointSurfaceStressesArray.get_headers()` labels

For that reason, validation-side reading was hardened:

- `D:\18a\MYSTRAN_Validation-main\op2_query.py`
  now reconstructs `PRINCIPALANGLE`, `MAJOR`, `MINOR`, `MAXSHEAR`,
  and `VONMISES` for `GRIDPOINTSURFACESTRESSES` directly from the raw
  `NX/NY/TXY` values
- this makes MSC and MYSTRAN GPSTRESS interpretation consistent even when the
  raw `OGS1` invariant slots differ from the advertised generic labels

### Practical implication

For GPSTRESS validation work:

1. trust `NX`, `NY`, `TXY` from `OGS1`
2. derive principal quantities from those direct components
3. do not rely on the remaining `OGS1` float slots having one portable,
   solver-independent meaning unless the exact binary convention has been
   confirmed case-by-case
