# FEAST/CHASE Dispatch Proof (LANCZOS path)

## Why previous MGIV-only evidence was insufficient
If deck uses `EIGR,...,MGIV`, the run stays on MGIV path and `PARAM,LANCMETH` is not used.

So dispatch proof must use `EIGRL` (LANCZOS branch in LINK4).

## Test decks used
- `midas_eigen08_cbeam_axis_debug_lanc_arpack.dat`
- `midas_eigen08_cbeam_axis_debug_lanc_chase.dat`
- `midas_eigen08_cbeam_axis_debug_lanc_feast.dat`
- `midas_eigen08_cbeam_axis_debug_lanc_feast_interval.dat`

All were generated from Eigen8 CBEAM deck with `EIGRL` and `PARAM,LANCMETH,...`.

## Runtime evidence (`.ERR`) - revalidated 2026-04-24
- `midas_eigen08_cbeam_axis_debug_lanc_arpack.dat`
  - no 49xx warnings (direct ARPACK baseline path)
- `midas_eigen08_cbeam_axis_debug_lanc_chase.dat`
  - `*WARNING 4916: CHASE NATIVE REQUIRES EIGRL FREQUENCY INTERVAL (V1/V2). USING ARPACK LANCZOS FALLBACK.`
- `midas_eigen08_cbeam_axis_debug_lanc_chase_interval_native.dat`
  - `*WARNING 4918: CHASE EXTERNAL BACKEND LINKED. RUNNING NATIVE CHASE PATH FOR MODES.`
- `midas_eigen08_cbeam_axis_debug_lanc_feast.dat`
  - `*WARNING 4915: FEAST NATIVE REQUIRES EIGRL FREQUENCY INTERVAL (V1/V2). USING ARPACK LANCZOS FALLBACK.`
- `midas_eigen08_cbeam_axis_debug_lanc_feast_interval.dat`
  - `*WARNING 4918: FEAST EXTERNAL BACKEND LINKED. RUNNING NATIVE FEAST PATH FOR MODES.`
  - `*WARNING 4919: ... REGULARIZING ... FOR NATIVE FEAST ATTEMPT.`
  - `*WARNING 4916: FEAST NATIVE FAILED (INFO=3). USING ARPACK LANCZOS FALLBACK.`

These warnings are emitted only inside:
- `EIG_LANCZOS_CHASE.f90`
- `EIG_LANCZOS_FEAST.f90`

So this confirms LINK4 dispatch to CHASE/FEAST wrappers is active on LANCZOS path.

## F06 evidence
All LANC decks show LANCZOS summary line:
- `E I G E N V A L U E   A N A L Y S I S   S U M M A R Y   (LANCZOS ... )`

This is expected for EIGRL branch.

## Quick modal comparison (first extracted cycles)
- ARPACK baseline / CHASE-no-interval / FEAST-no-interval / FEAST-interval-fallback:
  - same first modes (e.g., ~115.8844, 121.9792, 143.2731, 224.6861 Hz)
- CHASE interval-native:
  - different modal set (e.g., ~320.4222, 432.7091, 568.0889, 759.1902 Hz)
  - indicates native CHASE is honoring interval targeting differently than fallback path.

## Conclusion
- Dispatch proof: **confirmed** for FEAST/CHASE entry points via LANCZOS (`EIGRL`).
- Native backend link: **confirmed active** for both CHASE and FEAST (`4918` appears).
- FEAST native robustness: **still not converged** on Eigen8 interval deck (falls back with `INFO=3`).
- FEAST fallback policy: now **ARPACK Lanczos** (no MGIV fallback in Lanczos route).
- CHASE native: runs, but interval-native modal set currently differs from ARPACK baseline/fallback; needs follow-up validation against target interval expectations.
