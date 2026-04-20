# CBEAM Validation: PLOAD1 Axial And Torsion

This note validates the extended `PLOAD1` beam load path for the new full-span uniform local load types:

- `FXE` = local axial distributed force
- `MXE` = local torsional distributed moment

Validation input:

- [cbeam_pload1_axial_torsion_validation.dat](E:/mystran17/mystran/Binaries/cbeam_pload1_axial_torsion_validation.dat)

Geometry/material used:

- `L = 1.0`
- `A = 0.025`
- `J = 2.6041666667E-4`
- `E = 1.0E+7`
- `G = 4.1666666667E+6`

## Exact targets

### Subcase 1: uniform local axial load

For a cantilever beam with a uniform axial load `qx` over the full span:

- fixed-end reaction: `R = qx L`
- free-end axial displacement: `u(L) = qx L^2 / (2 E A)`

With `qx = 1.0`:

- reaction target = `1.000000E+00`
- displacement target = `2.000000E-06`

### Subcase 2: uniform local torsional moment

For a cantilever beam with a uniform distributed torsional moment `mx` over the full span:

- fixed-end reaction torque: `T = mx L`
- free-end twist: `theta(L) = mx L^2 / (2 G J)`

With `mx = 1.0`:

- reaction target = `1.000000E+00`
- rotation target = `4.608000E-04`

## MYSTRAN results

From:

- [cbeam_pload1_axial_torsion_validation.F06](E:/mystran17/mystran/Binaries/cbeam_pload1_axial_torsion_validation.F06)

### Subcase 1: `FXE`

- node 2 axial displacement `T1 = 2.000000E-06`
- fixed-end reaction `Fx = -1.000000E+00`
- recovered element axial force `= 1.000000E+00`
- recovered element torque `= 0.000000E+00`

This matches the exact axial target exactly within printed precision.

### Subcase 2: `MXE`

- node 2 torsional rotation `R1 = 4.608000E-04`
- fixed-end reaction `Mx = -1.000000E+00`
- recovered element torque `= 1.000000E+00`
- recovered element axial force `= 0.000000E+00`

This matches the exact torsion target exactly within printed precision.

## Conclusion

The current full-span uniform local `PLOAD1` path is now validated for:

- transverse local `y/z` distributed force
- axial local `x` distributed force
- torsional local `x` distributed moment

The validation run also exposed and fixed a pre-existing `PLOAD1` storage collision between subcases. `PDATA` is now sized and filled so different subcases can carry different beam/bar element loads without overwriting each other.
