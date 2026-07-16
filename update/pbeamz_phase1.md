# PBEAMZ Phase 1

## Scope

- `PBEAMZ` is treated as a `PBEAML`-based beam property with extra metadata and beam-station expansion.
- `LOADB` and `LOADB0` now count `PBEAMZ` as a beam property entry.

## Stored metadata

`BD_PBEAML` recognizes the optional `PBEAMZ` tokens:

- `STIFFMOD`
- `RIOFFSET` / legacy `ROFSET`
- `TAPER`
- `AREAMOD`
- `I1MOD`
- `I2MOD`
- `K1MOD`
- `K2MOD`
- `JMOD`

The extra metadata is stored in unused `RPBEAM` columns `46:54`.

## Stationing behavior

- Input intent:
  - uniform `PBEAMZ` only needs the start section
  - tapered `PBEAMZ` only needs the start and end sections
- Default `PBEAMZ` sectioning now expands to 11 stations at `x/L = 0.0, 0.1, ..., 1.0`.
- That gives the beam a 10-segment base representation even when the input deck only supplies one or two end stations.
- If the deck contains interior concentrated `PLOAD1` points, those locations are inserted as extra stations when the beam runtime builds its active station list.
- This keeps the base 10-segment representation while preserving point-load breakpoints when they exist.

## Taper behavior

- `TAPER` accepts `LINEAR`, `PARABOLIC`, `CUBIC`, or numeric codes `1/2/3`.
- `NONE`, `NONTAPER`, `NON-TAPER`, `CONSTANT`, `PRISMATIC`, and `OFF` are accepted as explicit non-taper aliases.
- `TAPER_MODE <= 0` keeps the end-A section constant and avoids the taper interpolation path.
- `TAPER_MODE > 0` uses the taper path with linear area interpolation, linear torsion interpolation, and exponent-based inertia interpolation.
- The current interpolation rule is:
  - `I(x) = [ (I1^(1/n)) * (1 - x/L) + (I2^(1/n)) * (x/L) ]^n`
  - `n = 1` for linear, `n = 2` for parabolic, `n = 3` for cubic

## Notes

- The current pass is parser/storage plus beam station expansion.
- Solver-side use of the metadata is intentionally phased in through the existing beam path.
- Build completed successfully after the change.
- Probe runs on `reference_msc\cbeam\prob_010a_end_offsets` and `prob_023_space_frame_modal` still terminate with `*ERROR 1611` (`KGG` has zero terms). The same fatal appears on the baseline PBEAML deck, so this is not yet a clean `PBEAMZ`-specific validation case.
# PBEAMZ Phase 1 Notes

## Current state

- `PBEAML` baseline for `rob_001_inclined_frame_pbeaml.dat` now runs through `END OF JOB`.
- `PBEAMZ` nominal baseline `rob_001_inclined_frame_pbeamz_nomod.dat` also runs normally.
- `PBEAMZ` metadata is stored in `RPBEAM(46:54)` and is visible in `ELMDAT1` / `BEAM`.
- `PBEAMZ` default station expansion currently uses 11 stations when only end stations are supplied.

## Decks in reference_msc/cbeam

- `rob_001_inclined_frame_pbeaml.dat`
- `rob_001_inclined_frame_pbeamz_nomod.dat`
- `rob_001_inclined_frame_pbeamz_mod.dat`

## Follow-up

- Compare `PBEAMZ` nominal vs `PBEAML` baseline first.
- Then compare `PBEAMZ` modifier deck with `AREAMOD=1000.0` and `K1MOD/K2MOD=0.0`.
- If those are stable, add the axial-gravity / concentrated-load cases next.
