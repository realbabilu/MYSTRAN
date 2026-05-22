# Wrapper Sources

This directory is source-only. It is meant to hold small bridge files used when linking MYSTRAN against external third-party libraries.

Included here:

- `superlu/c_fortran_dgssv.c` : Fortran/C bridge source used by the external SuperLU path.

Not present as a separate shim:

- FEAST currently does not use a separate `wrapper/feast/*` bridge file in this public branch.
- FEAST integration is implemented directly in `Source/LK4/EIGRL_EXTRACT_SOLVERS.F90`.

Not duplicated here:

- MUMPS integration logic lives in `Source/Modules/DMUMPS_STUF.f90`
- FEAST and RSA solver selection logic lives in the main MYSTRAN source tree

This public branch does not bundle third-party binary libraries in `wrapper/`.

