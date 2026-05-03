# Development Notes

## Scope

This archive covers the sparse `MUMPS` path added to MYSTRAN around COO/triplet-style transfer.

There are really two connected pieces:

1. `LINK3` direct sparse factor/solve with `MUMPS`
2. `LK4` / `ARPACK` sparse shift-invert using `MUMPS`

## Files included

- `CMakeLists.txt`
- `BD_PARAM.F90`
- `PARAMS.f90`
- `LINK3.f90`
- `DMUMPS_STUF.f90`
- `DSBAND_PREFAC.f`
- `EIG_INV_PWR.f90`
- `EIG_LANCZOS_ARPACK.f90`
- `EIG_LANCZOS_ARPACK_ADAPTIVE.f90`
- `ARPACK_LANCZOS_EIG.f`
- `REPORT_SOLVER_DISPATCH_POLICY.f90`

## Design summary

### LINK3

`LINK3` already had a sparse `SUPERLU` path. The new work added:

- `SPARSE_FLAVOR=MUMPS`
- direct factorization of `KLL`
- direct solve on the sparse factorization
- cleanup/free path after solve

The local MUMPS factorization in `LINK3` builds triplet arrays from MYSTRAN CRS data and feeds those directly to MUMPS.

### Shared helper module

`DMUMPS_STUF.f90` was added so the eigen side would not duplicate MUMPS runtime setup and solve code repeatedly.

It centralizes:

- compile-time guard with `DMUMPS_Solver`
- MUMPS runtime/init state
- CRS to MUMPS triplet conversion
- factor
- solve
- free/finalize

### ARPACK / LK4

Before this work, sparse ARPACK effectively meant:

- factor with `SUPERLU`
- solve with `SUPERLU`

The new work extends that to:

- factor with `MUMPS`
- solve with `MUMPS`

The touched areas were:

- inverse power helper
- ARPACK adaptive path
- ARPACK reverse-communication sparse solve calls
- backend-report diagnostics

## Important implementation note

The working local library family on this machine was:

- `C:\\gcc\\mumps32_nonmpi\\libdmumps.a`
- `C:\\gcc\\mumps32_nonmpi\\libmumps_common.a`
- `C:\\gcc\\mumps32_nonmpi\\libpord.a`
- `C:\\gcc\\mumps32_nonmpi\\libmpiseq.a`

That matched the earlier battle-solver environment better than the other guessed MUMPS tree.

## Commits

- `9ef5536` — `Add direct DMUMPS LINK3 sparse path`
- `ce7d7f8` — `Add ARPACK sparse MUMPS backend`
