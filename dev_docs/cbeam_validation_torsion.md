# CBEAM Validation 03: Torsion

This note records the torsion check for the current `CBEAM` implementation.

Input deck:

- [cbeam_torsion_validation.dat](E:/mystran17/mystran/Binaries/cbeam_torsion_validation.dat)

Output file:

- [cbeam_torsion_validation.F06](E:/mystran17/mystran/Binaries/cbeam_torsion_validation.F06)

## Model data

- `T = 1`
- `L = 1`
- `G = 4.1666666667E+6`
- `J = 2.6041666667E-4`

## Closed-form target

For a prismatic beam in pure Saint-Venant torsion:

- `theta_x = T L / (G J)`
- `theta_x = 1 / (4.1666666667E+6 * 2.6041666667E-4)`
- `theta_x = 9.216000E-04`

Expected other displacement/rotation components:

- `T1 = 0`
- `T2 = 0`
- `T3 = 0`
- `R2 = 0`
- `R3 = 0`

Expected reaction at the fixed node:

- `Mx = -1`

## MYSTRAN result

From [cbeam_torsion_validation.F06](E:/mystran17/mystran/Binaries/cbeam_torsion_validation.F06):

- free node `1002`: `R1 = 9.216000E-04`
- all other free-node components: `0.0`
- fixed node `1001` reaction: `Mx = -1.000000E+00`
- element engineering force table: torque `= 1.000000E+00`

## Comparison

- target `R1 = 9.216000E-04`
- MYSTRAN `R1 = 9.216000E-04`
- relative error `= 0%` within printed precision

Reaction check:

- expected fixed-end reaction `Mx = -1`
- MYSTRAN fixed-end reaction `Mx = -1.000000E+00`
- result: exact match in this test

Force recovery check:

- expected element torque `= 1`
- MYSTRAN element torque `= 1.000000E+00`
- result: exact match in this test

## Conclusion

This check validates the torsion part of the 12-DOF `CBEAM` formulation independently from axial and bending behavior.

For this one-element static case, the current `CBEAM` gives:

- exact torsional rotation match to `TL/(GJ)`
- exact fixed-end torque reaction
- exact recovered element torque

So the torsion block is behaving correctly in this benchmark.
