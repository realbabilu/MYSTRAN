# Eigensolver Integration Plan (ARPACK + FEAST + CHASE)

Date: 2026-04-22

## Current v18 Status (Confirmed)

1. `LANCMETH` now accepts `ARPACK`, `FEAST`, and `CHASE`:
   - [BD_PARAM.f90](/E:/mystran17/mystran/Source/LK1/L1A-BD/BD_PARAM.f90:1138)
   - [PARAMS.f90](/E:/mystran17/mystran/Source/Modules/PARAMS.f90:207)
2. LINK4 now dispatches FEAST/CHASE to dedicated routines:
   - [LINK4.f90](/E:/mystran17/mystran/Source/LK4/LINK4.f90:273)
   - [EIG_LANCZOS_FEAST.f90](/E:/mystran17/mystran/Source/LK4/EIG_LANCZOS_FEAST.f90:1)
   - [EIG_LANCZOS_CHASE.f90](/E:/mystran17/mystran/Source/LK4/EIG_LANCZOS_CHASE.f90:1)
3. Current FEAST/CHASE backend behavior is controlled surrogate (non-ARPACK) with warning:
   - FEAST -> MGIV surrogate (`*WARNING 4911`)
   - CHASE -> MGIV surrogate (`*WARNING 4912`)
4. In sparse mode, ARPACK shift-invert already uses SuperLU factor+solve (`KMSM`):
   - [ARPACK_LANCZOS_EIG.f](/E:/mystran17/mystran/Source/Modules/ARPACK/ARPACK_LANCZOS_EIG.f:667)
   - [ARPACK_LANCZOS_EIG.f](/E:/mystran17/mystran/Source/Modules/ARPACK/ARPACK_LANCZOS_EIG.f:852)

So your observation is correct: v18 is not LAPACK-only for ARPACK when `SOLLIB=SPARSE`.

## Why FEAST/CHASE Is Valuable

- FEAST: strong for interval/subspace extraction, can be robust for clustered spectra.
- CHASE: good performance for large Hermitian problems with block filtering.
- For Clement-like and large sparse benchmarks, ARPACK+SuperLU can be very fast; FEAST/CHASE should be tested as additional options, not immediate replacement.

## Latest Local Benchmark Snapshot (eigenbenchmark folder)

From `E:\mystran17\eigenbenchmark\arpack_chase_feast_superlu.exe` on 2026-04-22:

- ChASE: ~2.25 s
- ARPACK DENSE: ~3.92 s
- FEAST DENSE: ~1.38 s
- ARPACK + SUPERLU: ~0.40 s

Lowest and first five eigenvalues matched numerically across methods (Clement matrix case), supporting solver consistency in that benchmark.

## Recommended Integration Strategy

### Phase 1 (Low risk): External benchmark harness first

Use your `E:\\mystran17\\eigenbenchmark` tools to create a reproducible matrix benchmark report:
- ARPACK (dense)
- ARPACK + SuperLU (shift-invert)
- FEAST
- CHASE

Keep this outside MYSTRAN first to lock API, memory, and failure behavior.

### Phase 2 (Done): New LANCMETH options + controlled surrogate

Implemented values:
- `PARAM,LANCMETH,ARPACK`
- `PARAM,LANCMETH,FEAST`
- `PARAM,LANCMETH,CHASE`

Non-ARPACK methods currently route through dedicated wrappers and execute MGIV surrogate path.

### Phase 3 (Higher risk): In-core solver adapter layer

Introduce one interface module in LINK4 level:
- `EIG_LANCZOS_ARPACK` (existing)
- `EIG_LANCZOS_FEAST` (new)
- `EIG_LANCZOS_CHASE` (new)

And centralize common responsibilities:
- matrix/operator build (`KMSM` / generalized form)
- factorization reuse policy
- convergence/error mapping to MYSTRAN F06/ERR style

### Phase 4: Robustness and fallback

If FEAST/CHASE fails:
- auto fallback to ARPACK for same run (optional param-controlled)
- print clear reason and solver stats

## Decision Guidance

- Keep ARPACK+SuperLU as default sparse eigensolver for now.
- Add FEAST/CHASE as experimental selectable methods.
- Do not replace ARPACK path until validation matrix shows consistent gains across:
  - beam axis-debug
  - shell/solid modal cases
  - ill-conditioned and mechanism-prone models.

## Build Hooks Added (External FEAST/ChASE)

`CMakeLists.txt` now has optional hooks for external libraries:

- `MYSTRAN_USE_EXTERNAL_FEAST` (default `OFF`)
- `MYSTRAN_USE_EXTERNAL_CHASE` (default `OFF`)
- `MYSTRAN_FEAST_INCLUDE_DIR`, `MYSTRAN_FEAST_LIBRARY`
- `MYSTRAN_CHASE_INCLUDE_DIR`, `MYSTRAN_CHASE_LIBRARY`

Current local sandbox does not expose `C:\gcc\feast32\libfeast.a` and `C:\gcc\chase32\nompi\libchase_f.a`, so native linking was not activated in this run.
