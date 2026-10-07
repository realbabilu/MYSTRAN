# SIMOQ8 V13 alignment and MITC8 V4 upgrade

`PARAM,QUAD8TYP,SIMOQ8` retains its existing physical-gradient, three-translation,
pointwise-prestress 4x4 geometric stiffness. The V13 alignment documents that
behavior and removes duplicate gradient intermediates; static stiffness is
unchanged.

`PARAM,QUAD8TYP,MITC8` now follows the default supplied Python V4 configuration:

- Kikuchi field with tensor6/BDG4 shear blending:
  `w=max(0,1-15/s)`, `s=sqrt(element area)/thickness`.
- Automatic membrane tying on curved/warped geometry; direct membrane on flat
  straight-sided elements. ILS is disabled; drilling coefficient is 0.01.
- Consistent/HRZ mass with three rotary inertia components. Material density
  contributes rotary inertia; PSHELL NSM contributes midsurface translation mass.
- Recovery, thermal forces and prestress recovery use the upgraded operators.
  The physical TEMPP1 connectivity-gradient sign convention is retained.

Validation was completed in the accompanying workspace Python harnesses before
commit. The original source, executable, decks and solver outputs were saved
separately from the upgraded results; all seven before/after deck hashes match.

| Matched benchmark | Before | After |
|---|---:|---:|
| 2-005 16x16 CQUAD8, center Uz (inch) | -13.51560497 | -12.97469330 |
| 2-005 signed error vs -12.97 | +4.2067% | +0.0362% |
| 2-016 shear, 16x16, t=.01, stress error | +8.0733% | -0.1817% |
| 2-016 bending, 16x16, t=.01, stress error | +7.2814% | -0.3242% |
| 2-017 2x12, first buckling factor | 126.5028 | 126.5028 |
| 2-017 4x24, first buckling factor | 126.4854 | 126.4854 |

The two 2-016 t=1 first factors are unchanged at F06 precision. Thin-plate
critical stress is factor/t for unit reference line load; references are
0.042501034 (shear) and 0.116491057 (bending).

Checks completed:

- CMake build succeeded.
- SIMOQ8: 12 compiled production Kg comparisons against standalone V13;
  maximum relative matrix difference 1.465e-15. Static operator regression,
  five SOL105 cases and the unchanged 2-005 SOL101 response passed.
- MITC8: 270 B-operator points across 15 geometry/thickness configurations;
  maximum absolute difference against Python V4 3.109e-15. Consistent/HRZ mass
  matrices agree within 6.662e-16.
- MITC8 pressure/displacement and equivalent FORCE comparisons passed.
- Six native SOL105 comparisons against Python V4 passed; maximum relative
  first-three-factor difference 4.114e-7, internal residual below 2.004e-8,
  translation MAC above .999. Float32 OP2 modes can have much larger
  reconstructed residuals, so the gate uses internal residual and mode MAC.
- Six SOL103 density/NSM/HRZ checks passed; maximum relative five-frequency
  difference 3.413e-7.
- Five thermal cases passed: free/restrained expansion, gradient patches in
  both windings and the retained annulus benchmark.

The workspace evidence is under `test/mitc8_upgrade/v4_checkpoint` and
`test/q8_upgrade`; those directories are outside this repository checkout.
Executables and generated solver outputs are not tracked here.

Scope remains isotropic shells; geometric stiffness is the flat membrane
initial-stress tangent without director/rotation stress terms. The separate
Enhanced pressure report's 3D geometry-reproduction correction is not part of
the supplied V4 and was not ported. General curved modified-field pressure
retains its previously identified first-moment limitation.
