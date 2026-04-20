# CBEAM Validation 02: Axial Extension

This note records the axial-extension check for the current `CBEAM` implementation.

Input deck:

- [cbeam_axial_validation.dat](E:/mystran17/mystran/Binaries/cbeam_axial_validation.dat)

Output file:

- [cbeam_axial_validation.F06](E:/mystran17/mystran/Binaries/cbeam_axial_validation.F06)

## Model data

- `P = 1`
- `L = 1`
- `E = 1.0E+7`
- `A = 0.025`

## Closed-form target

For a prismatic bar/beam in pure axial extension:

- `u = P L / (E A)`
- `u = 1 / (1.0E+7 * 0.025)`
- `u = 4.000000E-06`

Expected other displacement/rotation components:

- `T2 = 0`
- `T3 = 0`
- `R1 = 0`
- `R2 = 0`
- `R3 = 0`

Expected reaction at the fixed node:

- `Fx = -1`

## MYSTRAN result

From [cbeam_axial_validation.F06](E:/mystran17/mystran/Binaries/cbeam_axial_validation.F06):

- free node `1002`: `T1 = 4.000000E-06`
- all other free-node components: `0.0`
- fixed node `1001` reaction: `Fx = -1.000000E+00`
- element engineering force table: axial force `= 1.000000E+00`

## Comparison

- target `u = 4.000000E-06`
- MYSTRAN `T1 = 4.000000E-06`
- relative error `= 0%` within printed precision

Reaction check:

- expected fixed-end reaction `Fx = -1`
- MYSTRAN fixed-end reaction `Fx = -1.000000E+00`
- result: exact match in this test

Force recovery check:

- expected element axial force `= 1`
- MYSTRAN element axial force `= 1.000000E+00`
- result: exact match in this test

## Conclusion

This check validates the axial part of the 12-DOF `CBEAM` formulation independently from bending and torsion.

For this one-element static case, the current `CBEAM` gives:

- exact axial displacement match to `PL/(EA)`
- exact fixed-end reaction
- exact recovered axial element force

So the axial block is behaving correctly in this benchmark.
