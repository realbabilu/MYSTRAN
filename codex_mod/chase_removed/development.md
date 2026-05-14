# Chase Removal Development Note

## Goal

Remove CHASE from the active MYSTRAN project path so the current release surface is FEAST-centered.

## Why

FEAST already covers the eigenvalue workflow we are keeping for current validation and future response-spectrum use. CHASE was an extra backend with no unique requirement in the current direction.

## Live source files updated

- `CMakeLists.txt`
- `build_last.bat`
- `Source/LK1/L1A-BD/BD_PARAM.F90`
- `Source/LK1/L1A-BD/BD_EIGRL.f90`
- `Source/LK4/LINK4.f90`
- `Source/LK4/EIGRL_EXTRACT_SOLVERS.F90`

## Behavior change

- CHASE is no longer selectable from the active build path
- FEAST remains enabled
- the user-facing parser text no longer advertises CHASE as an active option

## Validation

- rebuild completed successfully after the wiring removal
- remaining warnings are legacy warnings unrelated to CHASE removal

## Notes

This archive is intentionally smaller than the older `chase_feast_add` snapshot.
It exists to document the removal decision after `optimization_rcm_v2`.
