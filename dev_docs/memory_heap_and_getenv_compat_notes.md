# Memory Heap Safety + Intel getenv Compatibility (v18 port)

## Why this update was needed
- Large models can trigger stack/heap pressure when local automatic arrays scale with `NDOFG`.
- Intel toolchains (`ifx/ifort`) are more sensitive to non-standard `GETENV` usage across build environments.

## Changes applied
- Replaced `GETENV` with standard `GET_ENVIRONMENT_VARIABLE` in:
  - `Source/MAIN/GET_MYSTRAN_DIR.f90`

- Converted large local work arrays to `ALLOCATABLE` (heap-backed) with explicit deallocation in:
  - `Source/LK1/L1E/SPARSE_KGG.f90`
  - `Source/LK1/L1E/SPARSE_KGGD.f90`
  - `Source/LK2/LINK2.f90`
  - `Source/UTIL/PRT_MATS_ON_RESTART.f90`

## Expected impact
- Lower risk of runtime failure from large local arrays for big node/DOF models.
- Better portability/compatibility for Intel compiler builds while preserving behavior.

## Notes
- `WINAMEM` logic still exists and is still active in current code paths.
- This patch does not remove `WINAMEM`; it only reduces stack-risk hot spots and modernizes environment variable retrieval.
