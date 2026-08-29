# MYSTRAN Shell v1p8 Port Status

## Scope

This note summarizes the current MYSTRAN-side Fortran port status for the upgraded Python shell work, with emphasis on:

- default shell selectors
- two-step `GRID_SNORM` preprocessing
- thermal `alpha` / `PTE` support
- pressure `PPE` support

## Selector Status

### CQUAD8 default

`CQUAD8` now defaults to the Python-aligned Simo Q8 branch instead of `MITC8`.

Relevant implementation points:

- `PARAM,QUAD8TYP` default is `SIMOQ8`
- accepted Simo aliases normalize onto `SIMOQ8`
- the dispatcher routes `SIMOQ8` to `CQUAD8_SIMOQ8`

Files:

- [PARAMS.f90](D:/18a/MYSTRAN/Source/Modules/PARAMS.f90)
- [BD_PARAM.F90](D:/18a/MYSTRAN/Source/LK1/L1A-BD/BD_PARAM.F90)
- [EMG.f90](D:/18a/MYSTRAN/Source/EMG/EMG1/EMG.f90)

### CTRIA6 default

`CTRIA6` now has an explicit formulation selector through `PARAM,TRIA6TYP`.

Current implementation points:

- `PARAM,TRIA6TYP` default is `SIMOT6`
- accepted Simo aliases are normalized onto `SIMOT6`
- the dispatcher routes `SIMOT6` to `CTRIA6_SIMO1993`
- the dispatcher routes `MITC6` to `CTRIA6_MITC6`
- the dispatcher routes `MH6T` to `CTRIA6_MH6T`
- the dispatcher routes `REZAIEE` to `CTRIA6_REZAIEE`

Current caution:

- `SIMOT6` remains the default branch
- `MITC6`, `MH6T`, and `REZAIEE` are now active selectable branches rather than fallback placeholders

Files:

- [PARAMS.f90](D:/18a/MYSTRAN/Source/Modules/PARAMS.f90)
- [BD_PARAM.F90](D:/18a/MYSTRAN/Source/LK1/L1A-BD/BD_PARAM.F90)
- [EMG.f90](D:/18a/MYSTRAN/Source/EMG/EMG1/EMG.f90)
- [CTRIA6_SIMO1993.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/CTRIA6_SIMO1993.f90)
- [CTRIA6_MITC6.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/CTRIA6_MITC6.f90)
- [CTRIA6_MH6T.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/CTRIA6_MH6T.f90)
- [CTRIA6_REZAIEE.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/CTRIA6_REZAIEE.f90)

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
- `CQUAD8_SIMOQ8`
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
- `CQUAD8_SIMOQ8`: now ported from thermal stub to active `OPT(2)` path
- `CTRIA6_SIMO1993`: now ported from thermal stub to active `OPT(2)` path
- `CTRIA6_MITC6`: ported with the same active `OPT(2)` thermal path while replacing only membrane/shear with MITC tying
- `CTRIA6_MH6T`: ported with the same active `OPT(2)` thermal path while replacing membrane/shear with the MacNeal assumed-strain construction
- `CTRIA6_REZAIEE`: ported with the same active `OPT(2)` thermal path while replacing membrane/shear with the Rezaiee-based formulation

Implementation note:

- the new `CQUAD8_SIMOQ8` and `CTRIA6_SIMO1993` thermal paths currently use shell resultants assembled from `SHELL_A`, `SHELL_D`, and `SHELL_T`
- this is intended to activate `alpha` support consistently with the current shell stiffness formulation

Files:

- [MITC8.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/MITC8.f90)
- [CQUAD8_SIMOQ8.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/CQUAD8_SIMOQ8.f90)
- [CTRIA6_SIMO1993.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/CTRIA6_SIMO1993.f90)
- [CTRIA6_MITC6.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/CTRIA6_MITC6.f90)
- [CTRIA6_MH6T.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/CTRIA6_MH6T.f90)
- [CTRIA6_REZAIEE.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/CTRIA6_REZAIEE.f90)
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
- `CTRIA6_MITC6`
- `CTRIA6_MH6T`
- `CTRIA6_REZAIEE`
- `CQUADR_Q4RS`
- `CQUADR_SIMO1993`
- `CQUADR_Q4EASANS`
- `CQUADR_DKMQ24`
- `CQUADR_DKMQ24N`
- `CTRIAR_DKMT18`
- `CQUAD8_SIMOQ8`

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
2. run shell thermal regression cases for `MITC8`, `CQUAD8_SIMOQ8`, and `CTRIA6_SIMO1993`
3. compare pressure and thermal results against the Python reference models where available

## CTRIA6 MITC6 Smoke Validation

Quick smoke validation on August 15, 2026 with the rebuilt binary:

- `D:\18a\MYSTRAN_Validation-main\working\prob_2_002_nx01_ctria6_simot6.dat` terminated normally
- `D:\18a\MYSTRAN_Validation-main\working\prob_2_002_nx01_ctria6_mitc6.dat` terminated normally
- `D:\18a\MYSTRAN_Validation-main\working\prob_2_002_nx04_ctria6_mitc6.dat` terminated normally
- `D:\18a\MYSTRAN_Validation-main\working\prob_2_003_nx02_ctria6_mitc6.dat` terminated normally

Initial comparison note:

- for the `prob_2_002_nx01` cantilever-style check, `SIMOT6` and `MITC6` produce very close global displacement fields while remaining distinct formulations
- this confirms the new `PARAM,TRIA6TYP,MITC6` path is active through assembly, solve, and output, not just compile-time reachable

## CTRIA6 MH6T / REZAIEE Smoke Validation

Quick smoke validation on August 15, 2026 with the rebuilt binary:

- `D:\18a\MYSTRAN_Validation-main\working\prob_2_002_nx01_ctria6_mh6t.dat` terminated normally
- `D:\18a\MYSTRAN_Validation-main\working\prob_2_003_nx02_ctria6_mh6t.dat` terminated normally
- `D:\18a\MYSTRAN_Validation-main\working\prob_2_002_nx01_ctria6_rezaiee.dat` terminated normally
- `D:\18a\MYSTRAN_Validation-main\working\prob_2_003_nx02_ctria6_rezaiee.dat` terminated normally

Initial comparison note:

- these runs confirm that `PARAM,TRIA6TYP,MH6T` and `PARAM,TRIA6TYP,REZAIEE` are both live through assembly, solve, and output
- `SIMOT6` remains the default selector, but all four T6 branches are now directly callable from bulk data

## Shell Output Audit Status

### Family coverage

The current shell output path audit now covers the main shell families:

- `CQUAD4`
- `CQUADR`
- `CTRIA3`
- `CTRIAR`
- `CTRIA6`
- `CQUAD8`

Implementation notes:

- `CTRIAR` follows the shared `TRIA3` output family path
- `CTRIA6` is now included in the same shell stress/strain dispatch and writer routing used by the triangular shell family
- `CQUAD8` remains on its own quadrilateral shell output family path together with `CQUAD4` and `CQUADR`

Relevant files:

- [CALC_ELEM_STRESSES.f90](D:/18a/MYSTRAN/Source/LK9/L92/CALC_ELEM_STRESSES.f90)
- [CALC_ELEM_STRAINS.f90](D:/18a/MYSTRAN/Source/LK9/L92/CALC_ELEM_STRAINS.f90)
- [SHELL_STRESS_OUTPUTS.f90](D:/18a/MYSTRAN/Source/LK9/L92/SHELL_STRESS_OUTPUTS.f90)
- [SHELL_STRAIN_OUTPUTS.f90](D:/18a/MYSTRAN/Source/LK9/L92/SHELL_STRAIN_OUTPUTS.f90)
- [OFP3_STRE_NO_PCOMP.f90](D:/18a/MYSTRAN/Source/LK9/L92/OFP3_STRE_NO_PCOMP.f90)
- [OFP3_STRN_NO_PCOMP.f90](D:/18a/MYSTRAN/Source/LK9/L92/OFP3_STRN_NO_PCOMP.f90)
- [WRITE_ELEM_STRESSES.f90](D:/18a/MYSTRAN/Source/LK9/L91/WRITE_ELEM_STRESSES.f90)
- [WRITE_ELEM_STRAINS.f90](D:/18a/MYSTRAN/Source/LK9/L91/WRITE_ELEM_STRAINS.f90)

### GPSTRESS / reported-basis cleanup

The shell stress output cleanup now does two things:

- surface `GPSTRESS` output rotates shell-local stress rows into the requested surface basis before patch reconstruction
- the main quadrilateral shell table labels now refer to the reported coordinate system instead of the old element-coordinate wording

This is an output-basis cleanup, not a change to the shell stiffness kernels.

### CTRIA6 generated-normal fix

During verification of the `CTRIA6` path, a `LINK0` generated-normal bug was found and fixed:

- the preprocessing pass still stored shell grid rows in arrays sized for only 4 nodes
- `CTRIA6` therefore overran the temporary `BGRD` storage during generated `GRID_SNORM` assembly
- the generated-normal storage now accepts 6-node shell connectivity, while the flat triangle normal still uses the corner-triangle geometry

Relevant file:

- [LINK0.f90](D:/18a/MYSTRAN/Source/LK1/LINK1/LINK0.f90)

### Verification notes

Quick verification on August 13, 2026:

- `prob_2_002_nx01_ctria6.dat` runs to normal termination with the patched binary
- the original Python-generated `prob_2_002_nx01_cquad8.dat` still needs a CQUAD8 continuation card for MYSTRAN bulk-data parsing
- a MYSTRAN-compatible corrected deck confirms that the `CQUAD8` path also runs to normal termination once the required continuation and `PSHELL` convention are respected
