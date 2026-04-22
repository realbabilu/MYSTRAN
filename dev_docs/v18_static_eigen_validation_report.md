# v18 Static+Eigen Axis-Debug Validation (CBEAM/CBAR Subset)

Date: 2026-04-22 (rerun refreshed at 13:26 UTC)  
Executable: `E:\mystran17\mystran\Binaries\mystran.exe`

## Scope

Baseline rerun for 15 axis-debug decks:
- Eigen: 02, 03, 04, 05_1, 08, 09 modal, 12 (CBAR + CBEAM where available)
- Static: `midas_eigen09_pyramid_static_axis_debug.dat`

## Current Baseline Count (Original Decks, No Edits)

- PASS: 11
- FAIL: 4

Failing baseline decks:
1. `midas_eigen02_cbeam_axis_debug.dat`
   - `*ERROR 989` (`DPBSTF` KLL factorization failure)
2. `midas_eigen09_pyramid_modal_axis_debug.dat`
   - `*ERROR 9776` (NEV too high for `NDOFL=162`, `NEV=9`)
3. `midas_eigen09_pyramid_modal_cbeam_axis_debug.dat`
   - `*ERROR 9776` (NEV too high for `NDOFL=162`, `NEV=9`)
4. `midas_eigen12_cbeam_axis_debug.dat`
   - `*ERROR 989` (`DPBSTF` KLL factorization failure)

Static companion:
- `midas_eigen09_pyramid_static_axis_debug.dat`: PASS

## Low-Risk Fix Trial Refresh

Eigen09 modal:
- `EIGRL NEV 9 -> 6`:
  - `midas_eigen09_pyramid_modal_axis_debug_nev6.dat`: PASS
  - `midas_eigen09_pyramid_modal_cbeam_axis_debug_nev6.dat`: PASS

Eigen12 CBEAM:
- `EIGR MGIV -> INV`:
  - `midas_eigen12_cbeam_axis_debug_inv.dat`: PASS
- `SPC4-only`:
  - `midas_eigen12_cbeam_axis_debug_spc4.dat`: FAIL
- `MGIV + ART_MASS`:
  - `midas_eigen12_cbeam_axis_debug_mgiv_artmass.dat`: FAIL

Eigen02 CBEAM:
- `SPC4-only`:
  - `midas_eigen02_cbeam_axis_debug_spc4.dat`: FAIL
- `INV-only`:
  - `midas_eigen02_cbeam_axis_debug_inv.dat`: FAIL (`*ERROR 981`, KMSM SuperLU factorization)
- `MGIV + ART_MASS`:
  - `midas_eigen02_cbeam_axis_debug_mgiv_artmass.dat`: FAIL
- `INV + ART_MASS`:
  - `midas_eigen02_cbeam_axis_debug_inv_artmass.dat`: PASS
  - `midas_eigen02_cbeam_axis_debug_inv_artmass_spc4.dat`: PASS

## Exactly What Changed This Refresh

1. Re-ran all 15 baseline subset decks and regenerated status from `.ERR`.
2. Re-ran low-risk fix trial decks for eigen02/eigen09/eigen12.
3. Updated:
   - `dev_docs/v18_static_eigen_validation_summary.csv`
   - `dev_docs/v18_static_eigen_validation_report.md`

## Decks Still Failing (Baseline Only)

- `midas_eigen02_cbeam_axis_debug.dat`
- `midas_eigen09_pyramid_modal_axis_debug.dat`
- `midas_eigen09_pyramid_modal_cbeam_axis_debug.dat`
- `midas_eigen12_cbeam_axis_debug.dat`
