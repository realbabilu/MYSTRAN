# v18 Port Build Blockers (Current Branch Snapshot)

## Goal

Rebuild MYSTRAN on v18 branch with current CBEAM-port changes and Bernoulli-default experiments.

## Current state

- Build command used:
  - `cmake --build E:\mystran17\build --config Release -j 8`
- Build does **not** complete yet.

## Resolved blockers (already patched)

1. Legacy logging symbols (`WRT_LOG`, `F04`) missing in v18 `IOUNT1`:
   - Added compatibility members in `Source/Modules/IOUNT1.f90`.
2. Counter state symbols missing in v18 `SCONTR`:
   - Added `COUNTER_UPDATED`, `COUNTER_LIMITER` in `Source/Modules/SCONTR.f90`.
3. `PRESSURE_DATA_PROC` old call signatures:
   - Updated calls to `READERR`, `FILE_CLOSE`, `FILE_OPEN` to match v18 interfaces.

## Remaining blockers

Main remaining compile failures are in LK9 output/recovery files ported from v17 where interfaces drifted in v18:

- `Source/LK9/L91/WRITE_ELEM_ENGR_FORCE.f90`
  - Uses `ANS` logical unit from `IOUNT1`, but `ANS` no longer exists in v18 `IOUNT1`.
- `Source/LK9/L91/WRITE_ELEM_STRAINS.f90`
  - Same `ANS` unit mismatch.
- `Source/LK9/L91/WRITE_ELEM_STRESSES.f90`
  - Same `ANS` unit mismatch.
  - Additional procedure-call signature/type mismatches vs v18 helper interfaces (argument ordering and counts differ).

## Suggested next fix order

1. Align LK9 write/recovery files to v18 IO model:
   - replace/remove direct `ANS` unit writes with v18-supported outputs.
2. Reconcile changed helper subroutine signatures in `WRITE_ELEM_STRESSES` and `WRITE_ELEM_STRAINS`:
   - compare with upstream v18 versions and port only CBEAM-specific logic blocks.
3. Rebuild and then rerun static+eigen validation pack.

## Why this matters

Without full rebuild, validation currently uses existing binary in `Binaries`, so Bernoulli-default-at-source changes are not yet guaranteed to be active in executable behavior.
