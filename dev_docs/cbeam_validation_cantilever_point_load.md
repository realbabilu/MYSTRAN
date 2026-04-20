# CBEAM Validation 01: Cantilever With End Point Load

This note records the first simple validation of the current `CBEAM` implementation in this repository.

Case used:

- 2-node cantilever beam
- point load at free end
- one element
- load acts in the local bending plane associated with `T2/R3`

Input deck:

- [katili_beam_single.dat](E:/mystran17/mystran/Binaries/katili_beam_single.dat)

Output file:

- [katili_beam_single.F06](E:/mystran17/mystran/Binaries/katili_beam_single.F06)

## Model data

- `P = 1`
- `L = 1`
- `E = 1.0E+7`
- `G = 4.1666666667E+6`
- `A = 0.025`
- `I = 1.3020833333E-4`
- `k = 5/6 = 0.8333333333`

## Closed-form targets

For a cantilever beam with a point load at the free end:

Bernoulli-Euler:

- free-end rotation
  - `theta_B = P L^2 / (2 E I)`
  - `theta_B = 3.840000E-04`
- free-end deflection
  - `delta_B = P L^3 / (3 E I)`
  - `delta_B = 2.560000E-04`

Timoshenko:

- bending rotation stays the same for this comparison
  - `theta_T = P L^2 / (2 E I)`
  - `theta_T = 3.840000E-04`
- total free-end deflection
  - `delta_T = P L^3 / (3 E I) + P L / (k G A)`
  - shear term `delta_s = P L / (k G A) = 1.152000E-05`
  - `delta_T = 2.675200E-04`

## MYSTRAN result

From [katili_beam_single.F06](E:/mystran17/mystran/Binaries/katili_beam_single.F06), free node `1002`:

- `T2 = -2.675200E-04`
- `R3 = -3.840000E-04`

The negative sign is due to the downward load direction. For magnitude comparison:

- `|delta_MYSTRAN| = 2.675200E-04`
- `|theta_MYSTRAN| = 3.840000E-04`

## Comparison

Rotation:

- Bernoulli target: `3.840000E-04`
- Timoshenko target: `3.840000E-04`
- MYSTRAN: `3.840000E-04`
- Result: exact match in this test

Deflection:

- Bernoulli target: `2.560000E-04`
- Timoshenko target: `2.675200E-04`
- MYSTRAN: `2.675200E-04`
- Result: matches Timoshenko target, not pure Bernoulli target

Relative differences:

- vs Bernoulli deflection:
  - `(2.675200E-04 - 2.560000E-04) / 2.560000E-04 = 4.5%`
- vs Timoshenko deflection:
  - `0%` within printed precision

## Conclusion

For this one-element cantilever point-load case, the current `CBEAM` behaves as a shear-deformable beam:

- rotation matches the classical beam target
- deflection matches the Timoshenko target
- deflection is larger than the Bernoulli target by the expected shear contribution

This is a good first sign that the current formulation is acting like the intended DSB/Timoshenko-type beam rather than a pure Bernoulli beam.

## Next recommended validation

The next clean checks are:

1. axial extension: compare against `u = P L / (E A)`
2. torsion: compare against `theta_x = T L / (G J)`
3. cantilever thin beam: verify convergence toward Bernoulli as shear effect becomes small
4. fixed-simple or simple-fixed-roll benchmark using equivalent nodal loads for uniform loading
