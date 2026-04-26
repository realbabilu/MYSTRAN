<!-- ! --- chase_feast_add --- begin ! -->
# EIGRL Extract Methods V1

## Summary

This branch adds a new `EIGRL`-driven extract-method framework for `SOL 103`.

Supported method values on the first `EIGRL` continuation are:

- `ARPACK`
- `CHASE`
- `FEAST`
- `SUBSP`
- `DENSE`

`ARPACK` remains the default when no extract-method continuation is present.

## Input Format

Primary `EIGRL` line stays unchanged:

```text
EIGRL,SID,V1,V2,ND,MSGLVL,NCVFACL,SIGMA,NORM
```

New positional continuation:

```text
,EXTRACT_METHOD,METHOD_MODE,OPT1,OPT2,OPT3,OPT4,OPT5,OPT6
```

Examples:

```text
EIGRL,1,,,6
,DENSE,,64
```

```text
EIGRL,1,,,6
,SUBSP,,24,1.0E-6,40
```

```text
EIGRL,1,,,6
,CHASE,,64,1.0E-10,80
```

```text
EIGRL,1,0.,700.,6
,FEAST,,48,8,60,8,1.10
```

## Method Mapping

### ARPACK

- `OPT1` = `NCVFACL`
- `OPT2` = `SIGMA`
- `OPT3` = `NEV_DELT`
- `OPT4` = `MODE`
- `OPT5` = `LAP_MAT_TYPE`

### CHASE

- `OPT1` = `NEX`
- `OPT2` = `TOL`
- `OPT3` = `MAX_ITER`
- `OPT4` = `DEG`

### FEAST

- `OPT1` = `M0`
- `OPT2` = `TOL_DIGITS`
- `OPT3` = `MAX_LOOP`
- `OPT4` = `N_CONTOUR`
- `OPT5` = `SEARCH_SCALE`

### SUBSP

- `OPT1` = `NSUB`
- `OPT2` = `TOL`
- `OPT3` = `MAX_ITER`

### DENSE

- `OPT1` = `NEX`

## Defaults

If the continuation method line is omitted:

- `EIGRL` defaults to `ARPACK`

If the method is present but method-specific options are blank:

- `ARPACK`
  - `NEV_DELT = 2`
  - `MODE = 3`
  - `LAP_MAT_TYPE = DGB`
- `CHASE`
  - `NEX = max(64, 4*ND)`
  - `TOL = 1.0D-10`
  - `MAX_ITER = 80`
- `FEAST`
  - `M0 = max(48, 2*(ND + 16))`
  - `TOL_DIGITS = 8`
  - `MAX_LOOP = 60`
  - `N_CONTOUR = 8`
  - `SEARCH_SCALE = 1.10`
- `SUBSP`
  - `NSUB = max(ND + 16, 24)`
  - `TOL = 1.0D-6`
  - `MAX_ITER = 40`
- `DENSE`
  - `NEX = max(64, 4*ND)`

## Alias / Back Compatibility

`PARAM,LANCMETH,...` is still accepted as a deprecated alias.

Precedence is:

1. `EIGRL` extract-method continuation
2. `PARAM,LANCMETH,...`
3. default `ARPACK`

## Implementation Notes

### Shared Condensed Core

`CHASE`, `FEAST`, `SUBSP`, and `DENSE` all build the same condensed modal problem:

1. detect active DOF from positive diagonal mass entries
2. partition active vs zero-mass DOF
3. condense zero-mass DOF through Schur complement
4. build:
   - `Kcond`
   - diagonal active mass `Mcond`
5. solve:
   - `DENSE`, `SUBSP`, `CHASE` on standardized `Astd = M^(-1/2) Kcond M^(-1/2)`
   - `FEAST` on generalized condensed problem `Kcond x = lambda Mcond x`

Eigenvectors are expanded back to the full `L`-set before the normal MYSTRAN modal post-processing continues.

### DENSE

`DENSE` is the validated reference path:

- condensed problem
- standardized operator
- `DSYEV`

This is intended as a gold-standard / validation backend, not the large-model production default.

### SUBSP

`SUBSP` uses inverse subspace iteration on the standardized condensed operator.

### CHASE

`CHASE` uses native ChASE on the standardized condensed operator when the external library is enabled at CMake configure time.

If native ChASE is unavailable or unsupported for the request, the implementation falls back to `DENSE`.

### FEAST

`FEAST` uses native dense generalized FEAST on the condensed generalized problem when the external library is enabled at CMake configure time.

If native FEAST is unavailable or unsupported for the request, the implementation falls back to `DENSE`.

## Current Scope

- `SOL 103` only
- consistent-mass modal workflow
- modes path only

Buckling still falls back to ARPACK/Lanczos.

## Validation

Validated locally with the shell benchmark model `eigen13`.

Observed first 6 frequencies:

- `DENSE`
  - `79.91010`
  - `134.1940`
  - `239.8078`
  - `327.4830`
  - `364.4808`
  - `520.4631`
- `SUBSP`
  - matched `DENSE`
- `CHASE`
  - matched `DENSE` within last printed digits
- `FEAST`
  - matched `DENSE`

The post-processing path also completed normally:

- generalized mass
- mass renormalization
- LINK5/LINK9 modal output flow

## Files Touched

- `Source/LK1/L1A-BD/BD_EIGRL.f90`
- `Source/LK1/L1A-BD/BD_PARAM.F90`
- `Source/Modules/MODEL_STUF.f90`
- `Source/Modules/PARAMS.f90`
- `Source/UTIL/WRITE_L1M.f90`
- `Source/UTIL/READ_L1M.f90`
- `Source/LK4/LINK4.f90`
- `Source/LK4/EIGRL_EXTRACT_SOLVERS.F90`
<!-- ! --- chase_feast_add --- end ! -->
- matching interface files
- `CMakeLists.txt`
