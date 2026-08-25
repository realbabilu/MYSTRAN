# DSQK CQUAD4 Adoption Plan

## Decision

Adopt DSQK as a new `CQUAD4` family branch, not as a `CQUADR` branch.

Recommended selector:

- `PARAM,QUAD4TYP,DSQK`

Reference implementation to port first:

- [DSQK_ShellElement_RHR_6DOF.py](D:/18a/bending_only/Shell/gemini2/shit/claudeBuckling/got/DSQK_ShellElement_RHR_6DOF.py)

Underlying formulation note:

- the physical shell formulation is still the 5-DOF DSQK kernel from
  [DSQK_ShellElement_RHR.py](D:/18a/bending_only/Shell/gemini2/shit/claudeBuckling/got/DSQK_ShellElement_RHR.py)
- the `6DOF` file is a solver-interface wrapper that expands the kernel to
  standard shell `6 dof/node`

## Why CQUAD4

DSQK is closer to the current `CQUAD4` family than to `CQUADR`.

Main reasons:

- 4-node bilinear quadrilateral shell
- base kinematics are 5-DOF with added drilling stabilization
- transverse shear uses Bathe-Dvorkin / MITC4-style assumed shear
- bending is DSQK / discrete-Kirchhoff flavored, not Simo1993 EAS-membrane

Closest MYSTRAN relatives by formulation and architecture:

1. `MITC4`
2. `CQUAD4_SIMO1989`
3. `CQUADR_DKMQ24`
4. `CQUADR_SIMO1993`

## K6ROT Versus DSQK Drilling

These are not the same mechanism.

### DSQK Python drilling

In [DSQK_ShellElement_RHR_6DOF.py](D:/18a/bending_only/Shell/gemini2/shit/claudeBuckling/got/DSQK_ShellElement_RHR_6DOF.py):

- `K24 = T^T K5 T`
- then drilling stabilization is added as a nodal penalty about the actual
  shell normal `t3`
- scale:
  `drill_k = 0.03 * G * h * area`
- applied to rotational DOFs with `outer(t3,t3)` at each node

This is a local rotational penalty attached directly to the drilling direction.

### MYSTRAN CALC_K6ROT

In [CALC_K6ROT.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/CALC_K6ROT.f90):

- stiffness scale is approximately
  `1e-6 * K6ROT * SHELL_A(3,3) * area`
- it builds a drilling constraint vector `B`
- that vector couples the drilling DOF to adjacent-node translations and
  rotations
- it is not just a diagonal penalty on the drilling rotation

So even if both are "drilling stabilizers", they are structurally different:

- DSQK Python drilling: direct penalty on drilling rotation about shell normal
- MYSTRAN `CALC_K6ROT`: compatibility-style constraint spring assembled through
  a `B^T B` construction

### Practical consequence

Do not assume `CALC_K6ROT` is equivalent to the DSQK Python drilling term.

If exact Python parity is the short-term goal, the first DSQK MYSTRAN branch
should carry its own local drilling stabilization inside the element routine.

After that, a second experiment can test whether a "kernel without drilling +
external `CALC_K6ROT`" version is close enough for production use.

## Recommended MYSTRAN Strategy

### Phase 1

Create a new `CQUAD4` branch:

- `PARAM,QUAD4TYP,DSQK`
- new routine name suggestion:
  `CQUAD4_DSQK_RHR.f90`

Implementation policy:

- port from the Python `6DOF` wrapper
- keep DSQK drilling stabilization inside the element branch
- skip `CALC_K6ROT` for this branch

This matches how `QUADR` branches with built-in drilling already behave in
[EMG.f90](D:/18a/MYSTRAN/Source/EMG/EMG1/EMG.f90).

### Phase 2

Optional comparison branch:

- DSQK kernel without internal drilling
- rely on `CALC_K6ROT`
- compare results against Phase 1

This will tell whether common MYSTRAN drilling can replace the Python DSQK
drilling without materially changing benchmark behavior.

## Benchmark Priority

The requested problem set spans membrane, bending, warping, shell curvature,
eigen, and load-transfer behavior. The best order is not arbitrary.

### First-wave acceptance tests

Use these first because they are most sensitive to 4-node shell kinematics and
drilling compatibility:

1. `problem_2_001` patch test
2. `problem_2_003_quad.py`
3. `problem_2_002_cantilever.py`
4. `problem_2_004_macnealtwist.py`
5. `problem_2_005_clamped.py`
6. `problem_2_008_eigen.py`

Why:

- `problem_2_001` patch test is the first gate for membrane/bending consistency;
  if DSQK does not pass this cleanly, later shell benchmarks are not trustworthy
- `problem_2_003_quad.py` is a cheap early gate for basic quadrilateral shell
  behavior before moving to the more expensive curved or convergence studies
- cantilever checks basic bending/membrane coupling
- MacNeal twist is highly sensitive to warped/twisted shell behavior
- clamped plate is a strong quad-family discriminator already used in the repo
- eigen checks whether drilling and mass choices create spurious modes

### Second-wave shell benchmarks

Use next for curved-shell and load-transfer behavior:

7. `problem_2_006_scordelisloroof.py`
8. `problem_2_007_hemisphere.py`
9. `problem_2_010_pipe.py`

Why:

- Scordelis-Lo checks cylindrical shell bending/shear behavior
- hemisphere checks curved shell response under stronger geometric coupling
- pipe adds pressure/load-transfer sensitivity

### Specialized or later

10. `problem_2_005_triangle_simple.py`
11. `problem_2_009_elastic_foundation.py`

These are still useful, but less critical for deciding the initial DSQK
placement and drilling policy.

## Comparison Matrix

For each benchmark, compare:

1. Python DSQK 6DOF reference
2. MYSTRAN `CQUAD4 DSQK` with built-in DSQK drilling
3. MYSTRAN `CQUAD4 DSQK` with external `CALC_K6ROT` only
4. MYSTRAN `CQUAD4 MITC4`
5. MYSTRAN `CQUAD4 SIMO`
6. MYSTRAN `CQUADR DKMQ24`
7. MYSTRAN `CQUADR SIMO`
8. NASTRAN `CQUAD4`
9. NASTRAN `CQUADR`

Primary observables:

- static tip displacements
- shell stress/resultant trends
- eigenvalue convergence
- sensitivity to warped geometry
- pressure/load-transfer behavior where applicable

## Implementation Recommendation

Short version:

- classify DSQK as `CQUAD4`
- port the Python `6DOF` wrapper first
- keep built-in DSQK drilling for the first parity branch
- do not substitute `CALC_K6ROT` until benchmark evidence says the two are
  interchangeable enough

## Current Repo Status

The selector scaffolding is now wired into the MYSTRAN source:

- `PARAM,QUAD4TYP,DSQK` is reserved in the parser/dispatcher path
- `EMG` dispatch now calls `CQUAD4_DSQK_RHR`
- the `QUAD4/DSQK` branch now skips external `CALC_K6ROT`
- `CQUAD4_DSQK_RHR.f90` is no longer a fail-fast stub

Current DSQK phase-1 kernel scope in Fortran:

- standalone element-local helper set kept inside `CQUAD4_DSQK_RHR.f90`
- constant DSQK mean-plane basis using diagonal construction
- averaged nodal normal support with `GRID_SNORM` override
- membrane `Bm`
- bending `Bbu + BbDelta * An`
- MITC4-style shear tying
- internal drilling penalty following the Python `6DOF` wrapper strategy
- thermal membrane load vector using `ALPVEC` and nodal `DT`
- pressure load vector
- diagonal translational and rotational inertia
- shell recovery matrices `BE1/BE2/BE3`
- geometric stiffness scaffold using DSQK membrane stress-resultant recovery

Still to validate or refine against the Python reference and benchmarks:

- exact thickness/material scaling in the drilling and inertia paths
- patch-test parity and warped-shell parity
- buckling sensitivity versus Python DSQK and versus `CALC_K6ROT`
- whether any pressure or thermal sign/convention adjustments are needed after
  first benchmark runs

Relevant files:

- [PARAMS.f90](D:/18a/MYSTRAN/Source/Modules/PARAMS.f90)
- [BD_PARAM.F90](D:/18a/MYSTRAN/Source/LK1/L1A-BD/BD_PARAM.F90)
- [EMG.f90](D:/18a/MYSTRAN/Source/EMG/EMG1/EMG.f90)
- [EMG_USE_IFs.f90](D:/18a/MYSTRAN/Source/USE_IFs/EMG_USE_IFs.f90)
- [CQUAD4_DSQK_RHR.f90](D:/18a/MYSTRAN/Source/EMG/EMG4/CQUAD4_DSQK_RHR.f90)
- [CQUAD4_DSQK_RHR_Interface.f90](D:/18a/MYSTRAN/Source/Interfaces/CQUAD4_DSQK_RHR_Interface.f90)
