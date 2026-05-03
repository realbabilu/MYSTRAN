# Validation

## 1. LINK3 static solve

Earlier smoke validation used the bar deck from the direct MUMPS LINK3 trial.

Observed result:

- `SUPERLU` and `MUMPS` both terminated normally
- solved baseline displacement vector matched exactly
- reported comparison:
  - `max_abs_diff = 0.000000000000000E+00`

Interpretation:

- the direct `LINK3` sparse `MUMPS` path reproduced the existing sparse `SUPERLU` result on the checked static case

## 2. ARPACK sparse modal solve

Validation decks copied into this folder:

- [cbar_modal_103_superlu.dat](D:/fortran/mystran2/codex_mod/mumps_coo_add/cbar_modal_103_superlu.dat)
- [cbar_modal_103_mumps.dat](D:/fortran/mystran2/codex_mod/mumps_coo_add/cbar_modal_103_mumps.dat)

Configuration:

- `SOL 103`
- `ARPACK`
- `PARAM,SOLLIB,SPARSE,SUPERLU` or `PARAM,SOLLIB,SPARSE,MUMPS`
- same geometry/material otherwise

Final checked runs:

- both terminated normally
- both extracted `2` eigenvalues

Eigenvalue comparison:

- `SUPERLU`
  - `1.337580E+05`
  - `1.337580E+05`
- `MUMPS`
  - `1.337580E+05`
  - `1.337580E+05`

PowerShell comparison on the generated `.F06` files gave:

- `max_abs_eig_diff = 2.775558E-17`

## 3. Diagnostic banner check

The `MUMPS` run `.ERR` showed:

- `SOLLIB=SPARSE WILL USE SPARSE_FLAVOR=MUMPS`
- `ARPACK LINEAR BACKEND ... : ARPACK + MUMPS`

That confirms the backend reporting was updated along with the solver path itself.
