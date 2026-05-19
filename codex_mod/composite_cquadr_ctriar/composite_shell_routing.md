# Composite Shell Routing

## Purpose

This note defines the intended MYSTRAN routing for laminated/composite shell
elements using the current shell renovation basis.

The Python benchmark reference is:

- [composite_buckling_v5_operator.py](D:/python/t3ff/composite_buckling_v5_operator.py)

That reference uses laminated shell elements directly:

- `CompDKMQ24`
- `CompDKMT18`

and does **not** reinterpret a composite request as a generic legacy shell.

The same design should be followed in MYSTRAN.

## Routing Rule

If the user requests laminated/composite shell behavior through `PCOMP`:

- `CQUADR` should route to `DKMQ24`
- `CTRIAR` should route to `DKMT18`

This applies to:

- static
- modal
- buckling

It should **not** fall back to a generic quad/tri shell branch that ignores
the intended laminated DKMQ24/DKMT18 basis.

## Why This Is The Right Mapping

The current MYSTRAN shell pipeline already has the needed laminate foundation:

- [SHELL_ABD_MATRICES.f90](D:/fortran/mystran3/MYSTRANSolver-18.0.0/Source/EMG/EMG1/SHELL_ABD_MATRICES.f90)

That routine computes shell laminate matrices and related scaling for `PCOMP`:

- `SHELL_A`
- `SHELL_D`
- `SHELL_T`
- `FCONV`

The renovated shell elements already consume these shell matrices:

- [CQUADR_DKMQ24.f90](D:/fortran/mystran3/MYSTRANSolver-18.0.0/Source/EMG/EMG4/CQUADR_DKMQ24.f90)
- [CTRIAR_DKMT18.f90](D:/fortran/mystran3/MYSTRANSolver-18.0.0/Source/EMG/EMG4/CTRIAR_DKMT18.f90)

So the clean architecture is:

- laminate definition and ABD assembly live in the common shell property path
- laminated shell kinematics live in DKMQ24/DKMT18

## Python Reference Behavior

The Python v5 operator uses:

- a `Laminate`
- laminate prebuckling resultants from prescribed membrane/curvature state
- shell element operators based on:
  - `CompDKMQ24`
  - `CompDKMT18`

Core relation:

```text
N = A eps0 + B kappa
M = B eps0 + D kappa
```

Then the benchmark buckling operator is formed from that laminate state.

This means the Python reference is fundamentally:

- laminated DKMQ24
- laminated DKMT18

not:

- generic shell + separate composite post-processing

## Recommended MYSTRAN Implementation Stages

### Stage 1: Laminated Static / Modal

Goal:

- make `PCOMP + CQUADR` use laminated `DKMQ24`
- make `PCOMP + CTRIAR` use laminated `DKMT18`

Expected behavior:

- use `SHELL_A`, `SHELL_D`, `SHELL_T`
- use `FCONV`
- preserve current shell-local SNORM-aware behavior where applicable

This stage should be implemented first because it is the lowest-risk and most
direct extension of the shell renovation work already completed.

### Stage 2: Laminated Buckling Benchmark Parity

Goal:

- match the benchmark intent of `composite_buckling_v5_operator.py`

Important note:

The Python v5 script is a benchmark-level operator path, not yet a full
production prebuckling static solve with element-by-element stress recovery.

So this stage should first target:

- prescribed `eps0`
- prescribed `kappa`
- laminate resultants `N, M`
- buckling operator parity

This is the cleanest way to validate the laminated buckling formulation before
connecting it to the full MYSTRAN prebuckling recovery workflow.

### Stage 3: Production Composite Buckling Path

Goal:

- use MYSTRAN prebuckling static recovery
- recover laminate resultants/stresses consistently
- assemble `KGGD` for laminated DKMQ24/DKMT18 through the full solver path

This should only be attempted after Stage 2 parity is stable.

## Files Most Likely To Be Touched

Primary shell routing and property path:

- [EMG.f90](D:/fortran/mystran3/MYSTRANSolver-18.0.0/Source/EMG/EMG1/EMG.f90)
- [SHELL_ABD_MATRICES.f90](D:/fortran/mystran3/MYSTRANSolver-18.0.0/Source/EMG/EMG1/SHELL_ABD_MATRICES.f90)
- [ELEM_PROP_MATL_IIDS.f90](D:/fortran/mystran3/MYSTRANSolver-18.0.0/Source/LK1/L1C/ELEM_PROP_MATL_IIDS.f90)
- [BD_CQUAD.f90](D:/fortran/mystran3/MYSTRANSolver-18.0.0/Source/LK1/L1A-BD/BD_CQUAD.f90)
- [BD_CTRIA.f90](D:/fortran/mystran3/MYSTRANSolver-18.0.0/Source/LK1/L1A-BD/BD_CTRIA.f90)

Element kernels:

- [CQUADR_DKMQ24.f90](D:/fortran/mystran3/MYSTRANSolver-18.0.0/Source/EMG/EMG4/CQUADR_DKMQ24.f90)
- [CTRIAR_DKMT18.f90](D:/fortran/mystran3/MYSTRANSolver-18.0.0/Source/EMG/EMG4/CTRIAR_DKMT18.f90)

Buckling-related follow-up:

- element `OPT(6)` / `KED` paths in the shell kernels
- stress/resultant recovery bridge feeding laminated buckling operators

## Practical Recommendation

For the next implementation pass:

1. lock `PCOMP -> DKMQ24/DKMT18` routing first
2. validate laminated static/modal first
3. add benchmark laminated buckling operator parity second
4. postpone full production composite buckling recovery until the benchmark
   operator path is trusted

That order keeps the work incremental and makes debugging much easier.
