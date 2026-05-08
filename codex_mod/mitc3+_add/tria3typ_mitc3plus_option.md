# TRIA3TYP=MITC3+ Documentation

## Purpose

`PARAM,TRIA3TYP,MITC3+` adds `MITC3+` as a new **CTRIA3-class shell triangle option** in MYSTRAN.

This is intentionally documented as:

- `CTRIA3` shell topology
- `6 dof/node` shell behavior
- not a refined drilling family under `CTRIAR`
- analogous to `MITC4` / `MITC4+` living under `CQUAD4`

## User-Facing Syntax

```text
PARAM,TRIA3TYP,MIN3
PARAM,TRIA3TYP,MITC3+
```

Current meanings:

- `MIN3`
  - legacy Tessler/Hughes `TPLT2` triangle shell
- `MITC3+`
  - Lee/Lee/Bathe `MITC3+` triangle shell

## Routing

When the card is:

- `CTRIA3`

and the parameter is:

- `PARAM,TRIA3TYP,MITC3+`

the triangle shell bending/shear branch routes to:

- [TPLT_MITC3P.f90](D:/mystran2/MYSTRANSolver-18.0.0/Source/EMG/EMG4/TPLT_MITC3P.f90:1)

Legacy fallback remains:

- [TPLT2.f90](D:/mystran2/MYSTRANSolver-18.0.0/Source/EMG/EMG4/TPLT2.f90:1)

## Files Changed

All source edits for this addition were marked with the snippet guard:

- `! --- mitc3plus_add begin --- !`
- `! --- mitc3plus_add end --- !`

Primary files:

- [PARAMS.f90](D:/mystran2/MYSTRANSolver-18.0.0/Source/Modules/PARAMS.f90:406)
- [BD_PARAM.F90](D:/mystran2/MYSTRANSolver-18.0.0/Source/LK1/L1A-BD/BD_PARAM.F90:2604)
- [TREL1_USE_IFs.f90](D:/mystran2/MYSTRANSolver-18.0.0/Source/USE_IFs/TREL1_USE_IFs.f90:59)
- [TREL1.f90](D:/mystran2/MYSTRANSolver-18.0.0/Source/EMG/EMG4/TREL1.f90:199)
- [TPLT_MITC3P.f90](D:/mystran2/MYSTRANSolver-18.0.0/Source/EMG/EMG4/TPLT_MITC3P.f90:1)
- [TPLT_MITC3P_Interface.f90](D:/mystran2/MYSTRANSolver-18.0.0/Source/Interfaces/TPLT_MITC3P_Interface.f90:1)

## Validation Summary

### One-Element Closure

`CTRIA3 + PARAM,TRIA3TYP,MITC3+` now matches the uploaded Python `MITC3+` one-element shell stiffness very closely:

- full relative Frobenius difference: `5.331447947E-07`
- membrane block: `8.813135389E-08`
- plate block: `2.027922243E-06`

The remaining visible difference is only the small `rz` drilling penalty definition.

### Shell Benchmark Closure

Using MYSTRAN-generated shell decks derived from the legacy `MITC4` benchmark decks:

| Load Case | Theory | MYSTRAN MITC3+ | Error % |
|---|---:|---:|---:|
| `Static-22` | `-3.09000000E-01` | `-2.45034300E-01` | `20.701` |
| `Static-23` | `4.51970000E-04` | `4.56968000E-04` | `1.106` |
| `Static-24` | `9.40000000E-02` | `8.66592200E-02` | `7.809` |
| `Static-33 In-plane` | `-5.42400000E-03` | `-5.32377200E-03` | `1.848` |
| `Static-33 Out-of-plane` | `-1.75400000E-03` | `-1.87920400E-03` | `7.138` |

## MYSTRAN MIN3 vs MITC3+

At the current benchmark state:

- `MITC3+` is stronger on `Static-23` and especially `Static-33`
- `MIN3` is stronger on `Static-22` and `Static-24`

So the new option should be documented as:

- a valid alternative `CTRIA3` shell kernel
- not a blanket replacement for `MIN3`

## Journal / Reference

Primary formulation reference:

- `D:\mystran2\mitc3+\The_MITC3+_shell_element_and_its_performance.pdf`

Supporting local references:

- `D:\mystran2\mitc3+\katili2018.pdf`
- `D:\mystran2\mitc3+\katili2019.pdf`
- `D:\mystran2\mitc3+\rahmawati2018.pdf`

Python implementation used for closure:

- [mitc3plus.py](D:/mystran2/mitc3+/mitc3plus.py:1)
