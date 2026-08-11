# MYSTRAN Shell Two-Step Overhaul

## Goal

Make shell normal handling explicit and two-step:

1. Preprocess shell geometry and build `GRID_SNORM`.
2. Let shell element routines consume `GRID_SNORM` when present.

## Current State

- `LINK0` already contains `UPDATE_GENERATED_SHELL_SNORM`.
- `GRID_SNORM` already exists in `MODEL_STUF`.
- `CTRIA3_T3FF` already consumes `GRID_SNORM`.
- `CQUAD8_SIMOEAS1` and `CTRIA6_SIMO1993` now also consume `GRID_SNORM`.

## What Changed

- `CTRIA6_SIMO1993` now falls back to geometric normals, but prefers `GRID_SNORM` when available.
- `CQUAD8_SIMOEAS1` now does the same.
- `LINK0` now generates `GRID_SNORM` for `TRIA6` as well, so quadratic triangles participate in the same preprocessing path.

## Intent

- Keep the preprocessing phase geometry-only.
- Keep the element phase focused on stiffness/mass/load assembly.
- Make crease-aware nodal normals available consistently across the supported shell formulations.

## Notes

- `TRIA6` is currently always included in generated nodal-normal preprocessing.
- Manual `SNORM` data still has priority over generated normals.
- If a node is flagged as a crease/junction by the preprocessing pass, the generated normal is suppressed and the element falls back to its local geometric normal.

