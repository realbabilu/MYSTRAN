# Eigen8 Validation Proof (ARPACK vs FEAST vs CHASE)

## Scope
Deck family:
- `midas_eigen08_cbeam_axis_debug*.dat`

Compared runs:
- Baseline: `midas_eigen08_cbeam_axis_debug.dat`
- FEAST: `midas_eigen08_cbeam_axis_debug_feast.dat`
- CHASE: `midas_eigen08_cbeam_axis_debug_chase.dat`

## Result Summary
- All 3 runs ended normally (`MYSTRAN END` found in `.ERR`).
- Extracted eigenvalue count is identical: `5`.
- Mode 1–5 eigenvalue/radians/cycles are identical across ARPACK, FEAST, and CHASE.

## Numerical Evidence (from `.F06`)

| Mode | Eigenvalue | Radians | Cycles |
|---|---:|---:|---:|
| 1 | 5.301628E+05 | 7.281228E+02 | 1.158843E+02 |
| 2 | 5.873958E+05 | 7.664175E+02 | 1.219791E+02 |
| 3 | 8.103792E+05 | 9.002106E+02 | 1.432730E+02 |
| 4 | 1.993021E+06 | 1.411744E+03 | 2.246860E+02 |
| 5 | 7.023132E+06 | 2.650119E+03 | 4.217796E+02 |

These rows match for:
- `midas_eigen08_cbeam_axis_debug.F06`
- `midas_eigen08_cbeam_axis_debug_feast.F06`
- `midas_eigen08_cbeam_axis_debug_chase.F06`

## Files
How-to:
- `how_to_use.md`

Decks:
- `examples/midas_eigen08_cbeam_axis_debug.dat`
- `examples/midas_eigen08_cbeam_axis_debug_feast.dat`
- `examples/midas_eigen08_cbeam_axis_debug_chase.dat`

Outputs:
- `examples/midas_eigen08_cbeam_axis_debug.F06`
- `examples/midas_eigen08_cbeam_axis_debug_feast.F06`
- `examples/midas_eigen08_cbeam_axis_debug_chase.F06`
- `examples/midas_eigen08_cbeam_axis_debug.ERR`
- `examples/midas_eigen08_cbeam_axis_debug_feast.ERR`
- `examples/midas_eigen08_cbeam_axis_debug_chase.ERR`
