# MYSTRAN 18a Source Backport Summary

Branch:

- `v18.00.a`

Main source commits prepared in this workspace:

- `ad5c292` `Backport v18 solver, output, and validation fixes`
- `14dfeb1` `Add v18 backport update notes`
- `7460d0c` `Backport validation, output, and solver updates to v18.00.a`

This note is intended to mirror the 2026 release progression for this branch, but only for source and build changes that belong in this repository.

It is intentionally limited to:

- `CMakeLists.txt`
- `BUILD.bat`
- Fortran source changes under `Source/`

It intentionally does not summarize:

- external validation workspace notes
- temporary compare artifacts
- KMSM-specific discussion not yet carried as a dedicated patch in this branch

## April 2026

This stage corresponds to the earlier `mystran18-cbeam-chase-feast` direction.

Main source themes:

- initial CBEAM-oriented source backports
- initial FEAST-related source hooks
- initial `CHASE`-related eigen extraction support
- early output and element-side source changes needed to support the expanded beam and modal paths

Relevant source areas:

- `Source/LK1/L1A-BD/BD_PARAM.F90`
- `Source/LK1/L1A-BD/BD_EIGRL.f90`
- `Source/LK4`
- `Source/LK9`

Current note for this branch:

- `CHASE` is not part of the retained June/July direction
- if it appears in parser or parameter handling text, it should not be treated as a featured July capability

## June 2026

This stage corresponds to the `mystran-18.00a-cbeam-mumps-feast` direction.

Main source themes:

- CBEAM-related source and print/output stability fixes
- MUMPS sparse-solver integration in the source and build system
- FEAST integration in the source and build system
- removal of `CHASE` from the active release direction

Main build-system changes:

- `CMakeLists.txt` updated so external dependency paths can be passed from the command line
- `USE_MUMPS` and `EXTERNAL_MUMPS` logic added or expanded
- `USE_FEAST` and `EXTERNAL_FEAST` logic added or expanded
- external/internal SuperLU selection improved
- BLAS/OpenBLAS linkage kept aligned across MYSTRAN and third-party solver paths
- `BUILD.bat` added for the intended MinGW rebuild flow

Main source areas touched:

- `Source/LK1/L1A-BD/BD_PARAM.F90`
- `Source/LK1/L1A-BD/BD_EIGR.f90`
- `Source/LK1/L1A-BD/BD_EIGRL.f90`
- `Source/LK1/L1D/RFORCE_PROC.f90`
- `Source/LK4/EIG_LANCZOS_ARPACK.f90`
- `Source/LK4/EIG_SUMMARY.f90`
- `Source/LK4/LINK4.f90`
- `Source/Modules/PARAMS.f90`
- `Source/Modules/SCONTR.f90`
- `Source/UTIL/REPORT_SOLVER_DISPATCH_POLICY.f90`

Source intent at this stage:

- allow `MUMPS` to be selected explicitly as sparse flavor
- allow `FEAST` to be selected explicitly for eigen extraction
- keep source behavior consistent with the build-time availability of the external libraries

## July 2026

This stage corresponds to the current `mystran18a` branch state prepared in this workspace.

Main source themes:

- continued solver-dispatch cleanup
- subcase and statsub source backports
- RFORCE and load-processing fixes
- K6ROT-related source additions and stabilization fixes
- `MEFFMASS/MPFACTOR` Case Control compatibility bridge plus backend corrections
- `GPSTRESS/GSTRESS` acceptance plus initial `OGS1` writer compatibility work
- thin-shell `PSHELL/CQUAD4` compatibility fix when `MID3` is blank
- `PBEAMZ` parser and beam-property expansion path for modifier/taper-style beam input
- shell stress, shell strain, and element force writer backports
- principal stress and principal strain helper backports
- `NEU` writer architecture cleanup for more consistent and faster text output
- initial fast text-writer integration for selected `F06` and `NEU` hot paths

Main new source files added:

- `Source/EMG/EMG4/CALC_K6ROT.f90`
- `Source/Interfaces/BUILD_KGGD_FROM_UG_Interface.f90`
- `Source/Interfaces/CALC_K6ROT_Interface.f90`
- `Source/Interfaces/CC_STATSUB_Interface.f90`
- `Source/Interfaces/PRINCIPAL_STRAIN_2D_Interface.f90`
- `Source/Interfaces/PRINCIPAL_STRESS_2D_Interface.f90`
- `Source/Interfaces/READ_L5A_UG_FOR_SUBCASE_Interface.f90`
- `Source/Interfaces/REBUILD_KLLD_FROM_KGGD_Interface.f90`
- `Source/LK1/L1A-CC/CC_STATSUB.f90`
- `Source/LK1/LINK1/BUILD_KGGD_FROM_UG.f90`
- `Source/LK2/REBUILD_KLLD_FROM_KGGD.f90`
- `Source/LK9/L91/PRINCIPAL_STRAIN_2D.f90`
- `Source/LK9/L91/PRINCIPAL_STRESS_2D.f90`
- `Source/Modules/FAST_OUTPUT_FORMATTERS.f90`
- `Source/Modules/FEMAP_NEU_WRITE_HELPERS.f90`
- `Source/UTIL/f06_fast.c`
- `Source/USE_IFs/BUILD_KGGD_FROM_UG_USE_IFs.f90`
- `Source/USE_IFs/CC_STATSUB_USE_IFs.f90`
- `Source/USE_IFs/READ_L5A_UG_FOR_SUBCASE_USE_IFs.f90`
- `Source/USE_IFs/REBUILD_KLLD_FROM_KGGD_USE_IFs.f90`

Fast-writer direction in the current workspace:

- `f06_fast.c` provides fast fixed-width numeric and ID formatting helpers used by Fortran through `ISO_C_BINDING`
- `FAST_OUTPUT_FORMATTERS.f90` is the thin wrapper layer used by Fortran writers
- `FEMAP_NEU_WRITE_HELPERS.f90` centralizes repeated neutral-file text assembly so the `NEU` writers no longer rely on many small internal formatted writes
- `WRT_REAL_TO_CHAR_VAR.f90` now routes its `1ES14.6` conversion through the fast formatter
- `WRITE_GRD_PRT_OUTPUTS.f90` uses the fast formatter path for grid-result rows and summary extrema text
- `WRITE_ELEM_ENGR_FORCE.f90` uses prebuilt line buffers for selected `F06` element-force rows
- `WRITE_ELEM_STRESSES.f90` and `WRITE_ELEM_STRAINS.f90` now use the same fast `I8 + 6 x E14.6` path for the simple `BUSH` and `USERIN` rows
- shell stress/strain `F06` rows now also have dedicated C-backed line builders for the mixed-format layouts:
  `1403/1404/1405/1406/1407` and `1703/1704/1706`
- `WRITE_PLY_STRESSES.f90` and `WRITE_PLY_STRAINS.f90` now route the per-ply data rows through buffer builders backed by
  the fast formatter layer instead of repeated `WRITE(F06,1405:1408)` row formatting
- the fast formatter layer now also covers the mixed-width numeric fields needed by layered shell output:
  `ES14.5`, `ES10.2`, and `F9.3`
- the shell helpers keep the legacy field widths (`ES11.3`, `ES13.5`, `F8.2`, `E9.1`) rather than coercing them to `E14.6`

Intent and scope:

- speed up the hottest plain-text output patterns without changing table semantics
- keep `F06` human-readable and layout-compatible
- keep `NEU` text semantics stable while reducing repeated formatted-write overhead
- defer the more complex mixed-format shell layouts such as `1ES13.5` and `F8.2` until a later isolated pass
- `Source/UTIL/READ_L5A_UG_FOR_SUBCASE.f90`

`MEFFMASS/MPFACTOR` compatibility and backend direction in the current workspace:

- added `Source/LK1/L1A-CC/CC_MPF_MEFM.f90` and `Source/Interfaces/CC_MPF_MEFM_Interface.f90`
- legacy forms such as `MEFFMASS = ALL` still work
- Nastran-style forms such as `MEFFMASS(ALL)=YES` and `MPFACTOR(ALL)=YES` are now accepted
- the current bridge explicitly maps these parenthesized descriptors:
  - `ALL`
  - `SUMMARY`
  - `PARTFAC`
  - `MEFFM`
  - `MEFFW`
  - `FRACSUM`
  - `GRID=`
- unsupported extended descriptors still fall back to current MYSTRAN behavior with transition warnings
- `Source/LK9/L92/OFP2.f90` now calculates `SOL MODES` participation factors and effective modal mass from
  `MGG * rigid-body displacement` at the active reference point instead of relying on the older SPC-force-style path
- `Source/LK9/L91/WRITE_MEFFMASS.f90` and `Source/LK9/L91/WRITE_MPFACTOR.f90` now honor per-subcase
  `MEFMLOC_SUB` / `MEFMGRID_SUB` requests during output grouping for normal modal runs
- this fixes the earlier case where free-free normal modes could print zero or reference-point-insensitive
  participation outputs

Current boundary of the `MEFFMASS/MPFACTOR` work:

- normal modal (`SOL MODES`) output now respects subcase reference-point changes in the backend and writer path
- Craig-Bampton (`GEN CB MODEL`) still uses the older global `TR6_MEFM` transform path in
  `Source/LK6/CALC_CB_MEFM_MPF.f90`
- because of that, the existing `LOADC` warning about conflicting per-subcase `MEFFMASS/MPFACTOR` reference settings
  remains valid for the Craig-Bampton branch and should not be removed yet

Main source areas heavily updated:

- `CMakeLists.txt`
- `Source/EMG`
- `Source/Interfaces`
- `Source/LK1`
- `Source/LK2`
- `Source/LK4`
- `Source/LK5`
- `Source/LK9`
- `Source/Modules`
- `Source/USE_IFs`
- `Source/UTIL`

Main retained July source changes:

- solver dispatch logic for sparse and eigen backends
- explicit support for `MUMPS` and `FEAST` source-side behavior
- subcase/state reconstruction support for multiple subcase workflows
- RFORCE corrections in the source path
- K6ROT helper introduction and wiring
- `GPSTRESS/GSTRESS` parser and OP2/F06 writer-side compatibility updates
- restored `PSHELL` thin-shell fallback so blank `MID3` inherits `MID2`
- `PBEAMZ` deck parsing, metadata storage, and section/station expansion
- writer-side shell stress/strain/force updates
- principal 2D stress/strain output updates

Main writer/output files updated:

- `Source/LK9/L91/WRITE_ELEM_STRESSES.f90`
- `Source/LK9/L91/WRITE_ELEM_STRAINS.f90`
- `Source/LK9/L91/WRITE_ELEM_ENGR_FORCE.f90`
- `Source/LK9/L91/PRINCIPAL_2D.f90`
- `Source/LK9/L91/PRINCIPAL_STRAIN_2D.f90`
- `Source/LK9/L91/PRINCIPAL_STRESS_2D.f90`
- `Source/LK9/L92/OFP3_ELFE_1D.f90`
- `Source/LK9/L92/OFP3_STRE_NO_PCOMP.f90`
- `Source/LK9/L92/SHELL_STRAIN_OUTPUTS.f90`
- `Source/LK9/L92/SHELL_STRESS_OUTPUTS.f90`

Main `NEU` architecture changes:

- added `Source/Modules/FEMAP_NEU_WRITE_HELPERS.f90` as the central text-output helper for FEMAP neutral records
- replaced scattered `WRITE(NEU,format)` calls with assembled text lines followed by `WRITE(NEU,'(A)')` inside the helper layer
- moved the phase-1 geometry snapshot path in `LINK9` to the same helper-based write pattern
- moved the main FEMAP vector writers to the same path:
  - `Source/LK9/L91/WRITE_FEMAP_GRID_VECS.f90`
  - `Source/LK9/L91/WRITE_FEMAP_ELFO_VECS.f90`
  - `Source/LK9/L91/WRITE_FEMAP_STRE_VECS.f90`
  - `Source/LK9/L91/WRITE_FEMAP_STRN_VECS.f90`
- removed the old pattern where each writer repeated its own block of `NEU` record formatting logic
- kept `F06` behavior separate; this cleanup is specifically for FEMAP-neutral text output
- `Source/LK9/LINK9/LINK9.f90` now gives `MODES` and buckling-eigen FEMAP set headers a more informative block-450
  title/scalar payload instead of the previous zero-only placeholder

Intent of the `NEU` writer cleanup:

- keep file semantics compatible with the existing FEMAP-neutral direction already present in the branch
- reduce repeated formatted file writes in the hot writer loops
- make future `NEU` result-family routing changes easier because record assembly is now centralized
- avoid introducing a custom neutral dialect; this is a writer-architecture cleanup, not a format redesign

Current boundary of this work:

- this is a phase-1 neutral cleanup only
- it covers geometry snapshot plus the main grid, element force, stress, and strain vector writers
- it does not yet claim full import-grade FEMAP geometry export beyond the current snapshot approach
- block `450` is now closer to FEMAP v9 reference exports for modal and buckling-eigen sets, but it still does not
  reproduce the full `From:/Date:/<NULL>/title-trailer` metadata used by FEMAP text exports

## GPSTRESS / GSTRESS update

Current July status:

- `GPSTRESS` and `GSTRESS` are accepted in the Case Control path.
- `STRFIELD` is accepted for the same workflow.
- MSC-style `OUTPUT(POST)` decks with `SET` and `SURFACE` are now parseable in the local MYSTRAN path.
- the OP2 writer now emits a readable `OGS1` table for the investigated shell patch-test workflow
  and the row/framing issues that previously caused `pyNastran` to abort were corrected in the current branch state
- local validation tooling was hardened so `GRIDPOINTSURFACESTRESSES` invariants are recomputed from
  `NX/NY/TXY` instead of trusting the later raw `OGS1` payload slots

Important scope note:

- this is an initial compatibility implementation, not full MSC-style weighted grid-point stress recovery yet
- ordinary shell element stress remains the main path for element-level validation
- `GPSTRESS` is the grid/surface-style path and should be treated as its own result family

Typical usage pattern:

```nastran
GPSTRESS = ALL
STRFIELD = ALL

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
PARAM,POST,-1
PARAM,POSTEXT,YES
```

Notes for current use:

- if the deck has not explicitly requested `PARAM,STR_CID`, the GPSTRESS/GSTRESS path currently defaults the
  stress-coordinate request toward the basic/global style needed by the patch-test comparisons
- when validating with OP2, `OGS1` should be treated separately from ordinary shell `OES1X1`
- useful local references:
  - `update/gpstress_patch_test_example.md`
  - `update/gpstress_recovery_design.md`
  - `update/pynastran_gpstress_notes.md`
  - `update/ogs1_writer_audit_2026-07-19.md`

## Thin-shell `PSHELL` / `CQUAD4` blank `MID3` compatibility fix

Current July status:

- `Source/LK1/L1A-BD/BD_PSHEL.f90` now restores the old compatibility fallback:
  when `MID3` is blank, it is reset to `MID2`
- this keeps thin-shell-style `PSHELL` input compatible with decks that specify membrane and bending material
  but leave the transverse-shear material field blank
- this fix matters especially for `CQUAD4` shell paths where older decks intended “thin-shell-ish” behavior
  but were losing the expected material carry-over when `MID3` was left blank

Typical usage pattern:

```nastran
PSHELL,1,1,0.001,1
```

Meaning in the current branch:

- `MID1 = 1` supplies membrane material
- `MID2 = 1` supplies bending material
- blank `MID3` now inherits `MID2` again for compatibility
- field 6 (`12I/TM^3`) still defaults to `1.0` when blank, so ordinary bending inertia remains active

Practical outcome:

- thin shell patch-test and compatibility decks no longer require a forced explicit `MID3` entry just to
  preserve the historic MYSTRAN/Nastran-style intent
- this is a compatibility fix, not a new shell formulation

## `PBEAMZ` beam-property deck support

Current July status:

- `PBEAMZ` was added as a beam-property input path built around `PBEAML`-style section definitions
  plus a small wizard-style metadata tail
- the current branch stores extra `PBEAMZ` metadata in the beam property arrays and expands the property
  into a beam-station representation for downstream beam use
- `PBEAMZ` is intended to support section-based beam input while separating:
  - physical section geometry
  - stiffness-only modifiers
  - optional taper metadata
  - optional rigid-offset metadata

Current recognized metadata tail:

- `STIFFMOD`
- `NSM`
- `TAPER`
- `STATIONS`
- `RIOFFSET` / `ROFSET`

Minimal non-tapered example:

```nastran
PBEAMZ,1,1,,BAR
+,DIM0A,12.0,12.0,
+,END
```

Modifier example:

```nastran
PBEAMZ,1,1,,BAR
+,DIM0A,12.0,12.0,
+,STIFFMOD,1000.,1.0,1.0,0.0,0.0
+,END
```

Station-control example:

```nastran
PBEAMZ,1,1,,BAR
+,DIM0A,12.0,12.0,
+,STATIONS,10,0.33,0.66
+,END
```

Current semantics:

- `DIM0A` defines the start/end section for a prismatic member
- a tapered member adds `TAPER` and an end-section block such as `DIM1A`
- default stationing currently expands the beam into 11 stations (`0.0, 0.1, ..., 1.0`) unless overridden
- `STIFFMOD` is stiffness-only; it is meant to avoid the older fake-geometry workaround where area/inertia were
  exaggerated just to emulate a modifier
- `NSM` applies to mass/selfweight intent, not stiffness

Local validation direction used in this workspace:

- `prob_001_inclined_frame_pbeaml.dat` is the main baseline deck
- `prob_001_inclined_frame_pbeamz_nomod.dat` is the no-modifier `PBEAMZ` comparison deck
- `prob_001_inclined_frame_pbeamz_mod.dat` is the modifier comparison deck
- broader notes live in:
  - `update/pbeamz_phase1.md`
  - `update/pbeamz_torsion_audit.md`

## Scope Note

This backport is not one isolated fix. It is a rolling 2026 maintenance line for `v18.00.a` covering:

- build configuration
- external solver integration
- modal and eigen extraction behavior
- subcase and statsub behavior
- RFORCE and load handling
- K6ROT and stabilization-related source changes
- shell and element output writer updates
