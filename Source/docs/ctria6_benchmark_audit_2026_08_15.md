# CTRIA6 Benchmark Audit

Date: August 15, 2026

## Scope

This note records two additional Python-side benchmark checks for the active
quadratic triangular shell family:

- `SimoT6`
- `MITC6`
- `MH6T`
- `REZAIEE`

The goal here is to decide which benchmark cases are currently trustworthy as
targets for MYSTRAN `CTRIA6` port audit, and which cases are still dominated by
benchmark-harness issues.

## Files exercised

- `D:\18a\bending_only\Shell\gemini2\shit\validation\q8\battle_2_004_twisted.py`
- `D:\18a\bending_only\Shell\gemini2\shit\validation\q8\battle_2_006_scordelis.py`

## Problem 2-004 Twisted Beam

### Status

This benchmark is currently usable as a formulation discriminator.

### Python results at `nx=24`

Fixed diagonal:

| Element | `Uz` error | `Uy` error |
| --- | ---: | ---: |
| `Q8` | `+1.742%` | `+7.437%` |
| `SimoT6` | `+0.546%` | `+2.974%` |
| `MITC6` | `+0.328%` | `+2.425%` |
| `MH6T` | `-25.001%` | `-45.978%` |

Alternating diagonal:

| Element | `Uz` error | `Uy` error |
| --- | ---: | ---: |
| `SimoT6` | `+1.187%` | `+6.021%` |
| `MITC6` | `+1.137%` | `+5.520%` |
| `MH6T` | `-24.703%` | `-44.067%` |

### Reading

- `SimoT6` is good on this benchmark.
- `MITC6` is also good on this benchmark.
- `MH6T` is distinctly bad here and is much too soft.
- the diagonal pattern changes the numbers a bit, but does not change the
  overall ranking
- this means `battle_2_004_twisted.py` is a useful case for separating
  `SIMOT6` and `MITC6` quality from `MH6T`

## Problem 2-006 Scordelis-Lo Roof

### Status

This benchmark is **not** currently usable as a trustworthy `CTRIA6` port
target.

### Harness repair

The original `battle_2_006_scordelis.py` imported symbols that were not
exported by `problem_2_006_scordelis_q8.py`. On August 15, 2026 the harness was
patched so it now:

- carries its own Q8 shape functions
- carries its own Q8 consistent gravity integration
- runs all four T6 branches:
  - `SimoT6`
  - `MITC6`
  - `MH6T`
  - `REZAIEE`

### Python results after repair

Reference:

- `Uz(node49) = -0.3086`

Measured:

| Element | `Uz(node49)` | error |
| --- | ---: | ---: |
| `Q8` | `-0.062859` | `-79.63%` |
| `SimoT6` | `-0.026842` | `-91.30%` |
| `MITC6` | `-0.027014` | `-91.25%` |
| `MH6T` | `-0.034232` | `-88.91%` |
| `REZAIEE` | `-0.027675` | `-91.03%` |

### Reading

- this is not a `SimoT6`-only problem
- all T6 branches are far too stiff in magnitude
- even the Q8 branch remains far from the published reference
- because both Q8 and T6 are still poor here, `battle_2_006_scordelis.py`
  should currently be treated as a benchmark-harness audit, not a MYSTRAN port
  equivalence benchmark

## Practical conclusion

For the immediate MYSTRAN `CTRIA6` audit:

1. keep using `prob_2_003` as the direct Python-vs-MYSTRAN mismatch trigger
2. add `battle_2_004_twisted.py` as a trusted Python-side quality discriminator
3. do **not** use `battle_2_006_scordelis.py` yet to judge whether the Fortran
   port is right or wrong
4. treat `MH6T` as lower-confidence for default or parity targeting until its
   behavior is better understood

## MYSTRAN Twisted-Beam Follow-Up

Two MYSTRAN `CTRIA6` decks were created from the existing thick SAP Example
`2-004` geometry by splitting each quadrilateral patch into two fixed-diagonal
`CTRIA6` elements:

- `D:\18a\MYSTRAN_Validation-main\working\prob_2_004_thick_ctria6_simot6.dat`
- `D:\18a\MYSTRAN_Validation-main\working\prob_2_004_thick_ctria6_mitc6.dat`

These decks are mesh-equivalent to the existing `2-004` validation geometry,
which corresponds to the Python twisted-beam mesh with `nx=12`.

### Python reference at `nx=12`

| Formulation | `Uz` | `Uy` |
| --- | ---: | ---: |
| `SIMOT6` | `1.760041E-03` | `5.608056E-03` |
| `MITC6` | `1.710267E-03` | `5.413650E-03` |

### MYSTRAN tip average over `GRID 13, 26, 39`

Subcase mapping:

- subcase 1: load case `OUT`, compare `T3 = Uz`
- subcase 2: load case `IN`, compare `T2 = Uy`

| Formulation | MYSTRAN `Uz` | MY/PY | MYSTRAN `Uy` | MY/PY |
| --- | ---: | ---: | ---: | ---: |
| `SIMOT6` | `1.023854E-04` | `0.058172` | `7.002635E-05` | `0.012487` |
| `MITC6` | `1.121423E-04` | `0.065570` | `7.693630E-05` | `0.014212` |

### Reading

- both `SIMOT6` and `MITC6` are still far too stiff in MYSTRAN on the twisted
  beam
- this is now visible on a benchmark where the Python `SimoT6` and `MITC6`
  references are healthy
- the twisted-beam mismatch therefore strengthens the earlier conclusion from
  `prob_2_003`: the remaining problem is inside the Fortran `CTRIA6`
  formulation path, not just in a questionable Python benchmark harness

## Recommended next step

Use `SimoT6` and `MITC6` as the primary parity targets for the next Fortran
audit pass, because:

- both behave well in `battle_2_004_twisted.py`
- both are already live in MYSTRAN through `PARAM,TRIA6TYP`
- both are more credible short-term references than `MH6T`
