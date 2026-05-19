# shell_renovation Fortran Change Snapshot

This folder keeps a path-preserving snapshot of the Fortran files touched or
used as direct oracles for the shell renovation work.  The files under
`Source/...` are copied from the live `MYSTRANSolver-18.0.0/Source/...` tree and
were hash-checked against that tree after the copy.

Primary source markers are:

```fortran
! --- shell_renovation begin --- !
! --- shell_renovation end --- !
```

Older MITC3+ parameter/interface additions still use:

```fortran
! --- mitc3plus_add begin --- !
! --- mitc3plus_add end --- !
```

Main element kernels:

- `Source/EMG/EMG4/CQUADR_DKMQ24.f90`
- `Source/EMG/EMG4/CTRIAR_DKMT18.f90`
- `Source/EMG/EMG4/TPLT_MITC3P.f90`
- `Source/EMG/EMG4/TREL1.f90`
- `Source/EMG/EMG4/TMEM1.f90`
- `Source/EMG/EMG4/MITC_INITIALIZE.f90`

Dispatch, card parsing, parameter, and interface support:

- `Source/EMG/EMG1/EMG.f90`
- `Source/EMG/EMG1/SHELL_ABD_MATRICES.f90`
- `Source/LK1/L1A/LOADB.f90`
- `Source/LK1/L1A/LOADB0.f90`
- `Source/LK1/L1A-BD/BD_CQUAD.f90`
- `Source/LK1/L1A-BD/BD_CTRIA.f90`
- `Source/LK1/L1A-BD/BD_PARAM.F90`
- `Source/LK1/L1C/ELEM_PROP_MATL_IIDS.f90`
- `Source/Interfaces/*DKMQ24*`
- `Source/Interfaces/*DKMT18*`
- `Source/Interfaces/*MITC3P*`
- `Source/USE_IFs/EMG_USE_IFs.f90`
- `Source/USE_IFs/TREL1_USE_IFs.f90`
- `Source/Modules/PARAMS.f90`

Output/recovery support:

- `Source/LK9/L91/WRITE_ELEM_STRAINS.f90`
- `Source/LK9/L91/WRITE_ELEM_STRESSES.f90`
- `Source/LK9/L92/CALC_ELEM_STRAINS.f90`
- `Source/LK9/L92/CALC_ELEM_STRESSES.f90`
- `Source/LK9/L92/SHELL_STRAIN_OUTPUTS.f90`
- `Source/LK9/L92/SHELL_STRESS_OUTPUTS.f90`

Debug/oracle files copied because they were used during reconciliation:

- `Source/EMG/EMG4/MITC4.f90`
- `Source/EMG/EMG4/QMEM1.f90`

Related notes and run results are in:

- `codex_mod/shell_renovation/buckling06_notes.md`

Convenience copies of the main kernels are also refreshed at the
`codex_mod/shell_renovation` folder root.

## SKIP_K6ROT note

- The shell snapshot includes EMG.f90 logic that uses SKIP_K6ROT.
- This flag is **not** declared inside the archived ortran_changes bundle.
- When applying the shell changes into a live tree, add the following declaration to Source/Modules/MODEL_STUF.f90 near the other shell state flags:

`ortran
CHARACTER(1*BYTE) :: SKIP_K6ROT = 'N'
`

- Without that declaration, EMG.f90 will fail to compile with SKIP_K6ROT not found in module MODEL_STUF.
- Scope: this flag is only toggled inside EMG.f90 for selected shell paths and is intended to suppress K6ROT tweaks there.
