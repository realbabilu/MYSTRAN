# Shell Quad Formulation Matrix

Date: July 20, 2026

This note summarizes how the current `MYSTRAN 18.0.A` branch treats the main
4-node shell/plate quadrilateral paths:

- `QUAD4TYP=MIN4`
- `QUAD4TYP=MIN4T`
- `QUAD4TYP=MITC4`
- `QUAD4TYP=MITC4+`
- `CQUADR`

The goal is not theory completeness. The goal is to make clear:

- where membrane comes from
- where bending comes from
- where shear comes from
- whether the path can run through shell composite (`PCOMP`) logic
- whether the path tends to behave better as thin or thick

## Executive Summary

For `CQUAD4`, the formulation is selected by:

- `PARAM,QUAD4TYP,MIN4`
- `PARAM,QUAD4TYP,MIN4T`
- `PARAM,QUAD4TYP,MITC4`
- `PARAM,QUAD4TYP,MITC4+`

Source:

- [BD_PARAM.F90](D:/18a/MYSTRAN/Source/LK1/L1A-BD/BD_PARAM.F90)
- [PARAMS.f90](D:/18a/MYSTRAN/Source/Modules/PARAMS.f90)
- [EMG.f90](D:/18a/MYSTRAN/Source/EMG/EMG1/EMG.f90)

`CQUADR` is a separate element path in this branch and currently routes to the
DKMQ24 implementation:

- [CQUADR_DKMQ24.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/CQUADR_DKMQ24.f90)
- [EMG.f90](D:/18a/MYSTRAN/Source/EMG/EMG1/EMG.f90)

## Matrix

| Formulation | Membrane | Bending | Shear | Composite path | Thin/thick tendency | Main source path |
|---|---|---|---|---|---|---|
| `MIN4` | Q4 shell membrane path | Tessler MIN4 plate block | MIN4 constrained shear block | Yes, but shell coupling matrix `SHELL_B` path is limited | Better as general Mindlin, can be used thin-ish with care | `QDEL1 + QPLT2` |
| `MIN4T` | Q4 shell membrane path | 4 MIN3 sub-triangles, reduced to quad | MIN3/MIN4T shear from sub-trias | Yes, but more fragile, especially with orthotropy/composite coupling | Intended as thick Mindlin family, can be sensitive in distorted/thin cases | `QDEL1 + QPLT3 + TPLT2` |
| `MITC4` | MITC4 shell membrane part | MITC4 bending part | MITC covariant/direct-interpolated shear | Yes, with shell ply loop in `MITC4` | Best all-round shell path; usually preferred for thin to moderately thick | `MITC4 + MITC4_B` |
| `MITC4+` | Same overall shell path as `MITC4` | MITC4+ bending/shear form inside `MITC4_B` | MITC4+ covariant/direct-interpolated shear | Yes, with same shell ply loop path as `MITC4` | Intended improvement over MITC4, especially on distorted shells | `MITC4 + MITC4_B` |
| `CQUADR` | DKMQ24 membrane block | DKMQ24 bending block | DKMQ24 shear block | Currently best treated as non-composite/isotropic shell path unless separately audited | Very strong candidate for distorted quad plate/shell bending behavior | `CQUADR_DKMQ24` |

## Per-Formulation Notes

### 1. `MIN4`

`MIN4` is not a pure plate-only element at usage level. In MYSTRAN it is used
as the `QUAD4` shell formulation option.

Behavior split:

- membrane comes from the shell quad membrane path
- bending and transverse shear come from the MIN4 plate block

Main files:

- [QDEL1.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/QDEL1.f90)
- [QPLT2.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/QPLT2.f90)
- [BBMIN4.f90](D:/18a/MYSTRAN/Source/EMG/EMG6/BBMIN4.f90)
- [BSMIN4.f90](D:/18a/MYSTRAN/Source/EMG/EMG6/BSMIN4.f90)
- [MIN4SH.f90](D:/18a/MYSTRAN/Source/EMG/EMG7/MIN4SH.f90)

Practical reading:

- membrane is shell-like
- bending/shear are Mindlin MIN4-like
- this is a shell usage path with a plate-derived bending/shear core

### 2. `MIN4T`

`MIN4T` is the alternate `QUAD4` shell formulation assembled from four
sub-triangles.

Main files:

- [QDEL1.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/QDEL1.f90)
- [QPLT3.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/QPLT3.f90)
- [TPLT2.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/TPLT2.f90)

Reduction mode:

- controlled by `PARAM,MIN4TRED,B54` or `STC`

Source:

- [BD_PARAM.F90](D:/18a/MYSTRAN/Source/LK1/L1A-BD/BD_PARAM.F90)
- [PARAMS.f90](D:/18a/MYSTRAN/Source/Modules/PARAMS.f90)

Important caution:

- the code itself still contains notes about MIN4T orthotropic issues
- use extra care before treating `MIN4T` as the first-choice composite quad

Relevant source note:

- [MATERIAL_PROPS_2D.f90](D:/18a/MYSTRAN/Source/EMG/EMG8/MATERIAL_PROPS_2D.f90)

### 3. `MITC4`

`MITC4` runs through the dedicated shell routine:

- [MITC4.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/MITC4.f90)
- [MITC4_B.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/MITC4_B.f90)

This is the most direct “shell” path of the four `QUAD4TYP` options in this
branch.

Behavior:

- membrane: shell membrane part in `MITC4`
- bending: MITC4 bending block
- shear: MITC covariant/direct-interpolated shear

Composite handling:

- `MITC4.f90` explicitly has the shell composite/ply loop path
- this makes it the natural quad candidate when shell laminate/composite work
  is needed

### 4. `MITC4+`

`MITC4+` is not a separate top-level element routine. It uses the same top
level routine as `MITC4`, but switches internals in `MITC4_B`.

Main files:

- [MITC4.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/MITC4.f90)
- [MITC4_B.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/MITC4_B.f90)
- [MITC4_COVARIANT_STRAIN_DIRECT_INTERPOLATION.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/MITC4_COVARIANT_STRAIN_DIRECT_INTERPOLATION.f90)

Important code clue:

- `MITC4_B.f90` has explicit branching on `QUAD4TYP == 'MITC4+'`

Practical reading:

- same shell framework as `MITC4`
- different bending/shear interpolation details
- should be thought of as a shell formulation variant, not a separate deck
  element name

### 5. `CQUADR`

In this branch, `CQUADR` is routed to:

- [CQUADR_DKMQ24.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/CQUADR_DKMQ24.f90)

Dispatcher:

- [EMG.f90](D:/18a/MYSTRAN/Source/EMG/EMG1/EMG.f90)

Outputs:

- [WRITE_ELEM_STRESSES.f90](D:/18a/MYSTRAN/Source/LK9/L91/WRITE_ELEM_STRESSES.f90)
- [WRITE_ELEM_STRAINS.f90](D:/18a/MYSTRAN/Source/LK9/L91/WRITE_ELEM_STRAINS.f90)

Practical reading:

- this is the clean DKMQ-style path in the current branch
- if the user wants a non-`QUAD4TYP` quadrilateral to compare against
  `CQUAD4`, `CQUADR` is that path

## Composite Path Notes

Composite detection in shell paths is done through:

- [IS_ELEM_PCOMP_PROPS.f90](D:/18a/MYSTRAN/Source/EMG/EMG1/IS_ELEM_PCOMP_PROPS.f90)
- [SHELL_ABD_MATRICES.f90](D:/18a/MYSTRAN/Source/EMG/EMG1/SHELL_ABD_MATRICES.f90)

General rule in this branch:

- shell-type quads can enter `PCOMP` logic
- but not every stiffness coupling path is equally mature for every
  formulation

Important limitation:

- nonzero `SHELL_B` coupling for some shell paths still has explicit “code not
  written yet” guards in the older shell combination routines

See:

- [QDEL1.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/QDEL1.f90)
- [TREL1.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/TREL1.f90)

So “composite path = yes” here means:

- the formulation can enter shell laminate logic
- not that every possible unsymmetric laminate coupling case is equally mature

## Thin/Thick Reading

This is a practical engineering summary, not an absolute theorem.

- `MIN4`
  - general Mindlin family
  - usable for thick and moderate thin
  - may need care on locking/distortion comparisons
- `MIN4T`
  - also Mindlin family, but more special and more fragile
  - use when specifically studying that formulation
- `MITC4`
  - safest all-round shell default in many practical studies
- `MITC4+`
  - intended improvement over `MITC4`, especially for distortion-sensitive
    cases
- `CQUADR`
  - strong comparison path when plate/shell bending behavior on quads is the
    main focus

## Recommendation Snapshot

If the goal is:

- general shell laminate/composite work:
  - try `CQUAD4 + QUAD4TYP=MITC4` first
  - then compare `MITC4+`
- compare older Mindlin quad families:
  - use `MIN4` and `MIN4T`
- compare against a separate quadrilateral family:
  - use `CQUADR`
- distorted quad bending studies:
  - compare `MITC4+` and `CQUADR`

## Suggested Validation Table Extension

For future validation spreadsheets/CSV, use these columns:

- `Deck`
- `Element family`
- `Formulation`
- `Membrane path`
- `Bending path`
- `Shear path`
- `Composite requested`
- `Composite actually exercised`
- `Thin/Thick label`
- `Reference solver`
- `Pass/Fail`
- `Notes`

