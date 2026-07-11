# MYSTRAN 18a Backport Summary

Branch:

- `v18.00.a`

Latest local backport commit prepared in this workspace:

- `ad5c292` `Backport v18 solver, output, and validation fixes`

## Main Work Areas

### 1. Build and dependency selection

- Reworked `CMakeLists.txt` so external dependency paths can be passed in from the command line instead of hardcoding `C:/gcc/...`.
- Improved logic around:
  - `USE_MUMPS`
  - `USE_FEAST`
  - `EXTERNAL_MUMPS`
  - `EXTERNAL_FEAST`
  - external or internal SuperLU
- Kept BLAS/OpenBLAS usage aligned so MYSTRAN and dependent libraries point to the same BLAS import library.
- Added `BUILD.bat` to make the intended Windows MinGW rebuild flow easier to repeat.

### 2. Solver dispatch and sparse-library behavior

- Fixed or improved routing so decks requesting alternate solver libraries do not silently fall back to the wrong sparse path.
- Added reporting and support code for solver dispatch policy.
- Updated sparse and modal flow code touching:
  - MUMPS
  - FEAST
  - SuperLU
  - ARPACK/Lanczos handling

### 3. Subcase and statsub support

- Backported subcase-related fixes so multi-subcase behavior is closer to later MYSTRAN work and to the validation reference set.
- Added or updated support files for:
  - `CC_STATSUB`
  - reading prior subcase displacement data
  - rebuilding reduced stiffness data from prior state where needed

### 4. RFORCE and load-processing fixes

- Updated RFORCE parsing and processing paths.
- Adjusted load and case-control handling in the LK1 flow so requested solver or recovery behavior is honored more consistently.

### 5. K6ROT and stabilization work

- Added `CALC_K6ROT.f90` and related interfaces/use-files.
- Integrated K6ROT-related fixes into shell and element workflows where needed.
- Improved several shell/problem cases that were previously sensitive to stabilization details.

### 6. Stress, strain, force, and shell output backports

- Backported writer-side improvements so shell stress/strain/force output is closer to MYSTRAN 19 and closer to MSC/Nastran expectations.
- Added separate principal stress/strain helper routines:
  - `PRINCIPAL_STRAIN_2D.f90`
  - `PRINCIPAL_STRESS_2D.f90`
- Updated major output paths including:
  - `WRITE_ELEM_STRESSES.f90`
  - `WRITE_ELEM_STRAINS.f90`
  - `WRITE_ELEM_ENGR_FORCE.f90`
  - `OFP3_ELFE_1D.f90`
  - `OFP3_STRE_NO_PCOMP.f90`
  - shell output helpers in `LK9/L92`

### 7. Element and matrix/data-structure support

- Updated several EMG, LK1, LK2, LK4, LK5, and LK9 routines used by:
  - shell recovery
  - element offsets
  - eigen setup
  - reduced matrix recovery
  - output descriptor plumbing
- Expanded supporting interfaces and `USE_IFs` so the backported routines are wired correctly in 18a.

## Validation and reference work done alongside the code

The code backport was paired with a large amount of validation and reference cleanup in `D:\18a\MYSTRAN_Validation-main`.

Major supporting work included:

- comparing MYSTRAN 18a against:
  - MYSTRAN dev
  - MYSTRAN 19
  - MSC/Nastran
  - NX/Nastran in selected cases
- repairing reference decks that had:
  - missing OP2 output
  - duplicate/load-ID conflicts
  - bailout/autospc requirements
  - parser-path mismatches
- improving the OP2-vs-F06 validation workflow
- separating true parser gaps from reference-deck stabilization issues

Two useful generated validation maps are:

- `D:\18a\MYSTRAN_Validation-main\reference_msc_bailout_map.csv`
- `D:\18a\MYSTRAN_Validation-main\reference_msc_needs_bailout.csv`

These identify which MSC reference decks explicitly need `BAILOUT` or related stabilization handling and help keep that separate from actual MYSTRAN defects.

## Scope note

This backport is not a single isolated bugfix. It is a combined 18a maintenance pass covering:

- build configuration
- external solver integration
- modal/subcase behavior
- RFORCE and load handling
- shell output formatting and recovery
- stabilization and K6ROT-related behavior
- validation and parser support around MSC/NX reference outputs

