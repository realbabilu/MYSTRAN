# FEAST/CHASE Surrogate Validation (v18)

Date: 2026-04-22  
Executable: `E:\mystran17\mystran\Binaries\mystran.exe`

## Context

External FEAST/ChASE libraries are now linked in the v18 build
(`MYSTRAN_USE_EXTERNAL_FEAST=ON`, `MYSTRAN_USE_EXTERNAL_CHASE=ON`).
Current status:

- `LANCMETH=FEAST` now has a **native FEAST generalized sparse path** in `SOL 103` (modes),
  with controlled fallback to MGIV for unsupported/unsafe cases.
- `LANCMETH=CHASE` is still **MGIV surrogate**.
- `LANCMETH=ARPACK` remains native ARPACK path.

## Deck Used

- `midas_eigen09_pyramid_modal_axis_debug_nev6.dat` (ARPACK baseline)
- `midas_eigen09_pyramid_modal_axis_debug_nev6_feast.dat`
- `midas_eigen09_pyramid_modal_axis_debug_nev6_chase.dat`
- `feast_modal_chain_fullrank.dat` (native FEAST smoke test)

## Result

| Deck | Status | Runtime (s) | Notes |
|---|---|---:|---|
| `..._nev6.dat` | PASS | 0.0625 | native ARPACK |
| `..._nev6_feast.dat` | PASS | ~0.09 | FEAST native entered (`4918`), then fallback `4919` because `MLL` mass diagonal is not full-rank (12/162 positive) |
| `..._nev6_chase.dat` | PASS | 0.109 | `*WARNING 4914` CHASE external linked, generalized bridge still MGIV surrogate |
| `feast_modal_chain_fullrank.dat` | PASS | ~0.00 | FEAST native executed successfully (only `4918`, no `4916/4919`) |

Beam pyramid runs still extract eigenvalues via MGIV fallback by design (mass-singular generalized system).
Full-rank mini deck confirms FEAST native backend path is functional.

## Interpretation

- FEAST is now integrated as a real native option for generalized sparse modal solve when
  `MLL` is suitable (full-rank positive diagonal in current gate).
- For common beam modal decks with lumped translational mass only, FEAST currently auto-falls-back
  to MGIV using explicit diagnostics (`4919`) instead of ambiguous FEAST error (`INFO=2`).
- CHASE remains surrogate and should be promoted in a separate step.

## Next Step

1. Extend FEAST path to handle mass-singular generalized systems robustly (or add reduction strategy to finite-mass subspace).
2. Promote CHASE from surrogate to native path.
3. Keep warning codes split by state:
   - `4911`: FEAST backend missing
   - `4915`: FEAST interval missing
   - `4916`: FEAST failed and fallback used
   - `4918`: FEAST native entered
   - `4919`: FEAST skipped due to mass-rank gate
