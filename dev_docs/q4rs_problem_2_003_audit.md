# Q4RS Problem 2-003 Audit

Date: 2026-08-25

Scope: `PARAM,QUADRTYP,Q4RS` compared with `D:\18a\python\linear\Q4RS_ShellElement_Std_standalone.py` using `battle_2_003_combined.py`.

## Finding

The Q4RS `UY` response already matched Python closely, so the Problem 2-003 discrepancy was isolated to the out-of-plane `UZ` response and the bending/shear side.

The mismatch came from the Q4RS shear stabilization length. Python initializes the element diameter from zero and then uses the maximum distance from node 1 to nodes 2, 3, and 4. MYSTRAN initialized the same diameter with a lower bound of `1.0`. In the refined quarter-annulus meshes, the true element diameter can be less than `1.0`, so MYSTRAN used a smaller stabilization factor than Python and made Q4RS too flexible in `UZ`.

## Implemented

- `Source/EMG/EMG4/CQUADR_Q4RS.f90`
  - `QUAD_DIAMETER` now initializes `D` to zero, matching the Python reference.
  - A fallback `IF (D <= ZERO) D = ONE` preserves the degenerate-element guard without imposing a unit lower bound on valid small elements.

## Result After Change

Against Python `Q4RS`:

| nx | UY error | UZ error |
| ---: | ---: | ---: |
| 4 | -0.014% | -1.162% |
| 8 | -0.002% | 0.038% |
| 12 | -0.000% | 0.520% |
| 16 | -0.000% | 0.771% |

The old `UZ` plateau around `1.0` is removed. MYSTRAN now follows the Python Q4RS convergence curve for both `UY` and `UZ`.
