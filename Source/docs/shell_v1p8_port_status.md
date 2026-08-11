# MYSTRAN Shell v1p8 Port Status

## Scope

This note summarizes the current MYSTRAN-side Fortran port status for the upgraded Python shell work, with emphasis on:

- default shell selectors
- two-step `GRID_SNORM` preprocessing
- thermal `alpha` / `PTE` support
- pressure `PPE` support

## Selector Status

### CQUAD8 default

`CQUAD8` now defaults to the Simo Q8 branch instead of `MITC8`.

Relevant implementation points:

- `PARAM,QUAD8TYP` default is `SIMOEAS1`
- `SIMOQ8` is accepted as an alias and mapped to `SIMOEAS1`
- the dispatcher routes `SIMOEAS1` and `SIMOQ8` to `CQUAD8_SIMOEAS1`

Files:

- [PARAMS.f90](D:/18a/MYSTRAN/Source/Modules/PARAMS.f90)
- [BD_PARAM.F90](D:/18a/MYSTRAN/Source/LK1/L1A-BD/BD_PARAM.F90)
- [EMG.f90](D:/18a/MYSTRAN/Source/EMG/EMG1/EMG.f90)

## Two-Step SNORM Status

### Architecture

The intended shell-normal flow is:

1. preprocess shell geometry in `LINK0`
2. build or update `GRID_SNORM`
3. let shell elements consume `GRID_SNORM` when present

This keeps crease-aware nodal-normal generation outside the element stiffness routines.

### Current participants

The generated `GRID_SNORM` preprocessing path now covers:

- `CTRIA3_T3FF`
- `CTRIA6_SIMO1993`
- `CQUADR_Q4RS`
- `CQUAD8_SIMOEAS1`
- `CTRIAR_DKMT18`
- `CQUADR_DKMQ24`
- `CQUADR_DKMQ24N`

### Notes

- manual `SNORM` data still takes priority over generated normals
- crease/junction suppression remains a preprocessing concern
- `CTRIA6` has been added to the generated normal path via `UPDATE_GENERATED_SHELL_SNORM`

Files:

- [LINK0.f90](D:/18a/MYSTRAN/Source/LK1/LINK1/LINK0.f90)
- [shell_two_step_overhaul.md](D:/18a/MYSTRAN/Source/docs/shell_two_step_overhaul.md)

## Thermal Alpha Status

### Shell elements

Thermal support here means the element accepts shell thermal expansion data through `ALPVEC`, temperature through `DT`, and produces `PTE`.

Current status:

- `MITC4`: supported
- `MITC8`: now ported from thermal stub to active `OPT(2)` path
- `CQUAD8_SIMOEAS1`: now ported from thermal stub to active `OPT(2)` path
- `CTRIA6_SIMO1993`: now ported from thermal stub to active `OPT(2)` path

Implementation note:

- the new `CQUAD8_SIMOEAS1` and `CTRIA6_SIMO1993` thermal paths currently use shell resultants assembled from `SHELL_A`, `SHELL_D`, and `SHELL_T`
- this is intended to activate `alpha` support consistently with the current shell stiffness formulation

Files:

- [MITC8.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/MITC8.f90)
- [CQUAD8_SIMOEAS1.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/CQUAD8_SIMOEAS1.f90)
- [CTRIA6_SIMO1993.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/CTRIA6_SIMO1993.f90)
- [MITC4.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/MITC4.f90)

### New solid elements

The `NEWSOLID` family already has active thermal support through `ALPVEC`, `DT`, `TREF`, and `PTE` in the existing solid element routines.

Confirmed files:

- [HEXA.f90](D:/18a/MYSTRAN/Source/EMG/EMG5/HEXA.f90)
- [PENTA.f90](D:/18a/MYSTRAN/Source/EMG/EMG5/PENTA.f90)
- [PYRAM.f90](D:/18a/MYSTRAN/Source/EMG/EMG5/PYRAM.f90)
- [TETRA.f90](D:/18a/MYSTRAN/Source/EMG/EMG5/TETRA.f90)

## Pressure Status

Pressure support here means the element has an active `OPT(5)` path and writes `PPE`.

Confirmed shell coverage includes:

- `MITC4`
- `MITC8`
- `CTRIA3_T3FF`
- `CTRIA6_SIMO1993`
- `CQUADR_Q4RS`
- `CQUADR_SIMO1993`
- `CQUADR_Q4EASANS`
- `CQUADR_DKMQ24`
- `CQUADR_DKMQ24N`
- `CTRIAR_DKMT18`
- `CQUAD8_SIMOEAS1`

This means the audited linear and quadratic shell branches already have pressure assembly in place.

## Existing CQUAD4 and CTRIA3 Normal Handling

The older `CQUAD4` / `CTRIA3` paths are not identical to the new two-step crease-aware flow.

Practical distinction:

- existing paths may already consume `SNORM` or compatible nodal-normal data
- the new path is more explicit because `LINK0` first builds `GRID_SNORM`, including crease-aware suppression logic, and only then do the element routines consume the result

So the existing implementations are not necessarily wrong; they are simply less centralized than the newer two-step pipeline used by the upgraded shell formulations.

## Remaining Caution

This document records source-level port status only.

Recommended next step after code changes:

1. compile MYSTRAN
2. run shell thermal regression cases for `MITC8`, `CQUAD8_SIMOEAS1`, and `CTRIA6_SIMO1993`
3. compare pressure and thermal results against the Python reference models where available
