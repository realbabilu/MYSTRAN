# PBEAMZ Phase 1

## Scope

- `PBEAMZ` is treated as a `PBEAML`-based beam property with extra metadata and beam-station expansion.
- `LOADB` and `LOADB0` now count `PBEAMZ` as a beam property entry.
- The canonical compatibility reference remains:
  - `prob_001_inclined_frame_pbeaml.dat`
  - `prob_001_inclined_frame_pbeaml_msc.dat`
- `PBEAML` syntax must stay valid and unchanged; `PBEAMZ` should look like it first, then add wizard keywords after the section data.

## Stored metadata

`BD_PBEAML` recognizes the optional `PBEAMZ` tokens:

- `NSM`
- `STIFFMOD`
- `RIOFFSET` / legacy `ROFSET`
- `TAPER`
- `STATIONS`
- `AREAMOD`
- `I1MOD`
- `I2MOD`
- `K1MOD`
- `K2MOD`
- `JMOD`

The extra metadata is stored in unused `RPBEAM` columns `46:54`.

## Syntax shape

`PBEAMZ` follows the same section-body style as `PBEAML`, then adds a small wizard tail:

- the non-tapered form uses only the first section block
- the section dimensions are still written as raw values on the continuation line(s), just like `PBEAML`
- continuation labels such as `DIM0A` and `DIM1A` are accepted as ignored markers on those lines
- the section values are not written as `x/L` markers
- if taper is active, `TAPER` appears before the end-B section data
- the end-B section block is only present when taper is active
- `NSM` is optional and applies to mass/selfweight
- `STIFFMOD` is optional and applies to stiffness only
- `RIOFFSET` / `ROFSET` is optional and applies rigid offsets
- `STATIONS` is optional and uses the form `STATIONS,<base_segments>,<extra1>,<extra2>,<extra3>`
- the first value after `STATIONS` is the base segment count
- the optional extra values are extra breakpoints in `x/L` and are inserted into the beam station list
- `END` is mandatory for `PBEAMZ`

Suggested compact form:

```text
PBEAMZ, PID, MID, GROUP, TYPE/NAME
+ , DIM0A, <section A raw dimensions ...>
+ , TAPER, Tapertype
+ , DIM1A, <section B raw dimensions ...>    ! only when tapered
+ , NSM, realNSM1, realNSM2
+ , STIFFMOD, AREAmod, Imajmod, Iminmod, Ashear1mod, Ashear2mod, Torsionmod
+ , STATIONS, 10, 0.33, 0.66
+ , RIOFFSET, offset_i, offset_j
+ , END
```

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
- `PBEAMZ` smoke deck with explicit station control `prob_001_inclined_frame_pbeamz_stations.dat` also runs normally.
- `PBEAMZ` metadata is stored in `RPBEAM(46:54)` and is visible in `ELMDAT1` / `BEAM`.
- `PBEAMZ` default station expansion currently uses 11 stations when only end stations are supplied.
- `PBEAMZ` input syntax keeps the `DIM0A` continuation label in the deck and allows `STATIONS,10,0.33,0.66` as an optional station-control tail.

## Decks in reference_msc/cbeam

- `rob_001_inclined_frame_pbeaml.dat`
- `rob_001_inclined_frame_pbeamz_nomod.dat`
- `rob_001_inclined_frame_pbeamz_mod.dat`

## Follow-up

- Compare `PBEAMZ` nominal vs `PBEAML` baseline first.
- Then compare `PBEAMZ` modifier deck with `AREAMOD=1000.0` and `K1MOD/K2MOD=0.0`.
- Recheck `prob_001a` against the same syntax rule once the parser/doc are aligned.
- If those are stable, add the axial-gravity / concentrated-load cases next.
- Add a new radius-aware I-shape family, tentatively `I2`, instead of changing the meaning of the current `I` shape.
- `I2` should accept the same basic I-section dimensions plus root-radius data and compute `A`, `I1`, `I2`, and `J` automatically from the fuller section definition.
- Keep the current `I` shape as the simpler no-radius formulation so existing `PBEAML` / `PBEAMZ` decks stay stable.

## `prob_001_inclined_frame` comparison note

- `prob_001_inclined_frame_pbeamz_mod.dat` now matches `prob_001_inclined_frame.dat` on the parts that should match structurally:
  - subcase 2 to 7 displacements
  - SPC forces
  - grid-point force balance
  - beam engineering forces
- The remaining mismatch against `prob_001_inclined_frame.F06` is in the beam stress table, and that mismatch is expected for the current semantics:
  - the legacy reference `PBEAM` deck uses a fictitious geometric area (`A=144000`) to emulate an axial stiffness multiplier
  - `PBEAMZ` keeps the physical section geometry (`12 x 12`) and applies the multiplier only to stiffness
  - therefore local stress recovery is not numerically identical to the legacy fake-geometry deck even when the global response is matched
- This means `PBEAMZ` should be validated against the legacy deck primarily on response quantities, not on local stress magnitudes from the fake-area workaround deck.
