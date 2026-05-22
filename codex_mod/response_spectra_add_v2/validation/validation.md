! --- response_spectrum_mystran_add begin --- !
# Validation Summary (RS Writer)

## Scope
- Target case: Example 1-024
- Output package: geometry + modal + RS X/Y + linear sign/100-30 combos
- Delivery artifact: `artifacts/model-033-1-024-rs-writer.neu`

## What Was Validated
1. Geometry block imported in FEMAP without errors.
2. Modal sets recognized (Mode 1..4).
3. RS directional sets recognized: RS X Direction, RS Y Direction.
4. Combo sets recognized (sign variants + 100/30 variants).
5. Vector labels normalized to FEMAP style: Total Translation/T1/T2/T3/Total Rotation/R1/R2/R3.
6. P1 output profile package added and documented:
   - `FAST` (PLOT-only, SOL112-safe)
   - `CHECK` (PRINT,PLOT, SOL112-safe)
   - `FEMAP` (FAST + optional NEU post-step)
7. Wrapper execution smoke for `FAST` and `CHECK` on `Example-1-024` completed with normal solver termination.

## Acceptance Checks
- Constrained base nodes remain zero/near-zero in modal sets.
- RS combinations computed by per-node linear superposition (no repeated solver run needed).
- Package is external post-processing/writer; solver core unchanged.
- `run_femap_mode.ps1` now fails hard on solver non-zero exit and checks F06/ERR tail for abnormal termination markers.

## Runtime Stabilization Fixes (P1)
1. Disabled always-on high-frequency ERR debug spam for CQUAD4 output loops unless debug flag enabled:
   - debug print now behind `DEBUG(200) > 0`.
2. Added guard in `WRITE_GRD_OP2_OUTPUTS`:
   - fallback for invalid `ISUBCASE_INDEX`
   - fallback for undefined `ANALYSIS_CODE` in MFREQ/SOL112 path
3. Result: `Example-1-024_sol112_srss_femap_fast.bdf` and `_femap_check.bdf` now complete as `Finished solution`.

## SOL111 Compatibility Check (Example 1-024)
1. Converted `Example-1-024_sol112_srss.bdf` using `tools/sol112_to_sol111_compat.py`.
2. Ran generated deck in MYSTRAN (`SOL111`).
3. Status: `Finished solution` (normal termination), with non-fatal warnings only.
4. Optional `--clean-sol111` mode available to reduce non-essential log noise.
5. NEU numerical comparison (`tools/compare_neu_sets.py`) between SOL112 and SOL111 clean output:
   - `common_vectors = 8`
   - `max_abs_diff = 0.0`
   - `max_rel_diff = 0.0`
   - status: `same`

## Added Benchmark Packs

### Problem 1-020
- folder:
  - `validation/problem_1_020`
- purpose:
  - validated modal basis plus deterministic SRSS reconstruction reference
- key files:
  - `problem_1_020_validated.md`
  - `problem_1_020_modal.dat`
  - `benchmark_1_020_modal.py`
  - `benchmark_1_020_validate.py`
- use:
  - keep `1-020` as the lightweight modal/SRSS reference while `rsa_nastran` migrates input semantics

### Problem 1-024
- folder:
  - `validation/problem_1_024`
- purpose:
  - validated modal/static benchmark plus direct-reference RSA decks
- key files:
  - `problem_1_024_validated.md`
  - `problem_1_024_mystran_dense_modal_beam_swap.dat`
  - `problem_1_024_static_unit_load.dat`
  - `problem_1_024_sol112_srss_directref_beamz.bdf`
  - `problem_1_024_sol112_cqc_directref_beamz.bdf`
  - `problem_1_024_sol112_abs_directref_beamz.bdf`
  - `problem_1_024_sol112_10pct_directref_beamz.bdf`
- use:
  - `1-024` is the direct-reference RSA checkpoint for the corrected beam-axis and parser-safe benchmark set

### Problem 1-025
- primary working benchmark remains:
  - `docs/problem_1_025_stage1_note.md`
- use:
  - `1-025` is the active migration benchmark for `rsa_nastran` features:
    - `SDAMP`
    - `TABDMP1`
    - `DTI,SPECSEL`
    - `PARAM,OPTION`
    - active mass participation

! --- response_spectrum_mystran_add end --- !
