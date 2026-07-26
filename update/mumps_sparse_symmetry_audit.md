# MUMPS Sparse Symmetry Audit

Date:

- `2026-07-26`

Scope:

- `Source/LK1/L1A-BD/BD_PARAM.F90`
- `Source/Modules/DMUMPS_STUF.f90`
- `Source/LK3/LINK3.f90`
- `Source/LK2/REDUCE_KFF_TO_KAA.f90`
- `Source/LK4/EIG_INV_PWR.f90`

## Why this was added

The `18a` sparse path had two separate behaviors:

- `SUPERLU` factored sparse stiffness-style matrices as general sparse matrices
- `MUMPS` could still be forced down a symmetric-only path that kept only one triangle

This mismatch was exposed by `CQUADR` smoke decks where:

- `SUPERLU` solved normally
- `MUMPS` failed during factorization

The practical issue was not that `CQUADR` itself was unsupported. The problem was
that a matrix stored in `NONSYM` sparse form could still be handed to `MUMPS` as
if it were guaranteed symmetric.

## Parser update

`PARAM,SPARSE_FLAVOR,<value>` now works as a standalone 3-word parameter:

```nastran
PARAM,SOLLIB,SPARSE
PARAM,SPARSE_FLAVOR,SUPERLU
```

or

```nastran
PARAM,SOLLIB,SPARSE
PARAM,SPARSE_FLAVOR,MUMPS
```

This keeps `PARAM,SOLLIB,SPARSE` simple while still allowing the sparse backend
to be selected explicitly.

## New MUMPS-only audit policy

The new audit is intentionally narrow:

- it runs only on the `MUMPS` sparse path
- it runs once per matrix before factorization
- it checks both pair-pattern symmetry and numeric symmetry
- it does not change `SUPERLU`, banded, or dense paths

For stiffness-like sparse matrices handled here:

- `KLL` in `LINK3`
- `KOO` in `REDUCE_KFF_TO_KAA`
- `KMSM` in `EIG_INV_PWR`

the runtime behavior is now:

1. If `SPARSTOR='SYM'`, treat the matrix as symmetric for `MUMPS`.
2. If `SPARSTOR='NONSYM'`, audit the actual CRS matrix contents.
3. If the matrix is pairwise and numerically symmetric, use `MUMPS` symmetric mode.
4. Otherwise, use `MUMPS` general sparse mode.

This gives `MUMPS` a matrix-mode decision based on the matrix that MYSTRAN
actually assembled, not just on the storage policy alone.

## F06 reporting

When the audit runs, MYSTRAN writes an informational line to the `F06`
identifying which mode was chosen for `MUMPS`.

Meaning:

- mode `Y` = factor as symmetric
- mode `N` = factor as general sparse

This was added so deck triage can see whether a sparse failure is due to the
matrix itself or due to the symmetry assumption used for factorization.

## Investigated smoke deck

The investigated workspace deck was:

- `D:\18a\bending_only\Shell\gemini2\working_mystran\prob_2_002_mesh_02_quadr.dat`

Observed behavior after the patch set:

- `SUPERLU` path terminates normally
- `MUMPS` path also terminates normally after the audit-based dispatch

## Cost / performance note

This audit does add some overhead, but it is intentionally limited:

- only on `MUMPS`
- only once per matrix
- only for the targeted stiffness-like matrices above

So it is a safety/compatibility guard, not a global sparse slowdown.
