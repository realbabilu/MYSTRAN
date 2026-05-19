# RBE2 `dy/dz` Swap Bug

## Status

Historical hypothesis only. This package is now superseded by the proven deck-generation bug in [bug_rbe2_card_continuation](C:/mystran3/wt_optimization_rcm_v2/codex_mod/bug_rbe2_card_continuation/bug_rbe2_card_continuation.md).

The old `dy/dz` suspicion was raised before the generated `RBE2` cards were checked against the actual MYSTRAN `BD_RBE2` field contract. Once the `RBE2` decks were rewritten with proper continuation cards, `RBE2 126` and explicit diaphragm `MPC` became numerically identical for `Problem 1-024`.

That means this package should be read as a record of an intermediate debugging hypothesis, not as the current root cause.

## Bug Marker

Snippet marker used for this issue:

```fortran
! --- codex_mod\\bug_rbe2_dy_dz ---
```

See [rbe2_delta_bug_snippet.f90](C:/mystran3/wt_optimization_rcm_v2/codex_mod/bug_rbe2_dy_dz/rbe2_delta_bug_snippet.f90).

## Why It Is Superseded

The stronger evidence now points elsewhere:

- [BD_RBE2.f90](C:/mystran3/wt_optimization_rcm_v2/Source/LK1/L1A-BD/BD_RBE2.f90) only reads dependent grids from fields `5-9` on the first card and `2-9` on continuation cards
- the original generated benchmark decks overflowed that format and silently truncated slave grids
- after fixing deck continuation, `RBE2` and explicit `MPC` matched exactly

Current numerical evidence after the deck fix:

- `RBE2 126`: `0.274964, 0.268458, 0.087143, 0.085224 s`
- explicit `MPC`: `0.274964, 0.268458, 0.087143, 0.085224 s`

Static comparison after the deck fix:

- joint `28`, `Ux = 1.234744000e-04`, `Rz = -1.176715000e-06` for both
- joint `29`, `Ux = 2.265431000e-04`, `Rz = -1.136647000e-06` for both

## Current Source Cross-Check

The active `RBE2_PROC` convention now matches the rigid-body reference builder in [RB_DISP_MATRIX_PROC.f90](C:/mystran3/wt_optimization_rcm_v2/Source/LK1/L1C/RB_DISP_MATRIX_PROC.f90):

- `ux` couples to `+dz * ry - dy * rz`
- `uy` couples to `-dz * rx + dx * rz`
- `uz` couples to `+dy * rx - dx * ry`

So the present source state is not the main proven issue for `1-024`.

## Historical Wrong vs Correct Snippet

Pre-patch:

```fortran
DELTA_0(1,2) =  (RGRID(GRID_ID_ROW_NUM_D,3) - RGRID(GRID_ID_ROW_NUM_I,3))
DELTA_0(1,3) = -(RGRID(GRID_ID_ROW_NUM_D,2) - RGRID(GRID_ID_ROW_NUM_I,2))
```

Post-patch:

```fortran
DELTA_0(1,2) =  (RGRID(GRID_ID_ROW_NUM_D,2) - RGRID(GRID_ID_ROW_NUM_I,2))
DELTA_0(1,3) = -(RGRID(GRID_ID_ROW_NUM_D,3) - RGRID(GRID_ID_ROW_NUM_I,3))
```

The rest of the block stays:

```fortran
DELTA_0(2,1) = -DELTA_0(1,2)
DELTA_0(2,3) =  (x_dep - x_ind)
DELTA_0(3,1) = -DELTA_0(1,3)
DELTA_0(3,2) = -DELTA_0(2,3)
```

## Included Files

- [rbe2_delta_bug_snippet.f90](C:/mystran3/wt_optimization_rcm_v2/codex_mod/bug_rbe2_dy_dz/rbe2_delta_bug_snippet.f90)
- [mystran/Source/LK1/L1D/RBE2_PROC.f90](C:/mystran3/wt_optimization_rcm_v2/codex_mod/bug_rbe2_dy_dz/mystran/Source/LK1/L1D/RBE2_PROC.f90)
