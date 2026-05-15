# reduced_warning

Archive bundle for the warning-reduction pass in the fresh `mystran3` tree.

## Scope

- reduce actionable `maybe-uninitialized` warnings in active code paths
- reduce `realloc-lhs` warnings in the active eigen solver path
- improve human readability in dense matrix operations by splitting nested expressions into smaller steps
- keep ARPACK changes minimal while the broader overhaul is still in progress

## Included files

- `Source/EMG/EMG4/QPLT3.f90`
- `Source/LK4/EIGRL_EXTRACT_SOLVERS.F90`
- `Source/LK9/L92/OFP3_ELFE_1D.f90`
- `Source/Modules/ARPACK/ARPACK_LANCZOS_EIG.f`

## Notes

- `QPLT3.f90` initializes temporary triangle-side data and debug indices to silence false-positive uninitialized warnings.
- `OFP3_ELFE_1D.f90` initializes temporary engineering-force storage used by the LINK9 element-output path.
- `EIGRL_EXTRACT_SOLVERS.F90` now uses smaller, more readable matrix-operation steps and exports the buckling helper entry points so the compiler no longer flags them as unused private procedures.
- `ARPACK_LANCZOS_EIG.f` only received a minimal initialization tweak; deeper ARPACK cleanup is intentionally deferred.

## Snippet marker

`! --- reduced_warning begin --- !`

`! --- reduced_warning end --- !`
