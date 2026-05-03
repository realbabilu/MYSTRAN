# MUMPS COO Add

This folder archives the direct `MUMPS` sparse work in two stages:

1. `LINK3` static solve support with direct CRS/COO-style handoff into MUMPS
2. `ARPACK` sparse shift-invert support using the same MUMPS backend

Copied source files here are snapshots from the live tree after the work was finished.
The copied files are annotated with:

- `! --- MUMPS_COO add begin --- !`
- `! --- MUMPS_COO add end --- !`

Those markers are for archive reading only, so the modifications can be scanned quickly.

Main outcomes:
- `PARAM,SOLLIB,SPARSE,MUMPS` works in `LINK3`
- `ARPACK` sparse modal solve can use `MUMPS` instead of `SUPERLU`
- diagnostics now report `ARPACK + MUMPS` honestly when that backend is active

See also:
- [development.md](D:/fortran/mystran2/codex_mod/mumps_coo_add/development.md)
- [validation.md](D:/fortran/mystran2/codex_mod/mumps_coo_add/validation.md)
