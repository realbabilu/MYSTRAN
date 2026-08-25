# DKMQ24EAS Problem 2-003 Audit

Date: 2026-08-25

Scope: `PARAM,QUADRTYP,DKM24EA` compared with `D:\18a\python\linear\DKMQ24_EAS4_ShellElement_RHR_standalone.py` using `battle_2_003_combined.py`.

## Finding

The large `UZ` drift in Problem 2-003 was caused by a formulation mismatch, not by the deck translation. The Python class is `DKMQ24 + EAS4`, while MYSTRAN's `DKM24EA` path was using the legacy AU DKMQ24 branch plus EAS4. Since EAS4 only modifies the membrane block, MYSTRAN's out-of-plane `UZ` response stayed identical to `DKM24AU`.

The `SNORM` generator was also tested. Removing `DKM24EA` from generated-SNORM requirements did not change this flat quarter-annulus result, so SNORM is not the cause of the 2-003 drift.

The remaining large `UY` error was in the flat EAS membrane frame. Python's `DKMQ24_EAS4` computes the flat-element membrane/EAS frame from the surface basis (`e1 = cross(g2,e3)`) and forms the EAS local Jacobian from that same frame. The MYSTRAN EAS4 path was still using `GEOMETRY_AT`, which chooses an arbitrary in-plane frame from a global reference vector. That is harmless for some patch cases, but not for the curved in-plane `UY` loading in Problem 2-003.

## Implemented

- `Source/EMG/EMG4/CQUADR_DKMQ24N.f90`
  - `DKM24EA` now uses the native DKMQ/DKMQ20-style base branch before EAS4 condensation, instead of the AU base branch.
  - The EAS4 membrane path has a fixed centroid-frame `Bm` helper, matching the Python implementation's flat-element policy.
  - The EAS4 membrane `Bm` and EAS `G` helpers now use a Python-style surface basis/local Jacobian for flat `DKM24EA`, without changing the older `GEOMETRY_AT` behavior used by other CQUADR variants.
- `Source/Modules/PARAMS.f90`
  - `DKM24EA` no longer requests generated SNORM by default. `DKM24AU` and `Q4RS` still do.

## Result After Change

Against Python `DKMQ24EAS`:

| nx | UY error | UZ error |
| ---: | ---: | ---: |
| 4 | -0.055% | 1.543% |
| 8 | -0.000% | 1.527% |
| 12 | 0.000% | 1.515% |
| 16 | 0.000% | 1.503% |

The AU-like `UZ` problem is fixed, and the `UY` EAS membrane behavior now matches Python. The remaining observed difference in Problem 2-003 is the roughly 1.5% `UZ` offset, which should be audited separately on the bending/shear side if tighter parity is required.
