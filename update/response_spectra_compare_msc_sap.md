# Response Spectrum Comparison: MYSTRAN vs MSC vs SAP2000

Date: July 14, 2026

## Scope

This note compares the current MYSTRAN RSA path against:

- MSC/Nastran 2025.2 reference runs
- available SAP2000 benchmark references

Primary files:

- MYSTRAN:
  - `D:\18a\femap_RSA\problem_1_020_rsa_srss_new.F06`
  - `D:\18a\femap_RSA\problem_1_024_rsa_srss_new.F06`
- MSC:
  - `D:\18a\femap_RSA\problem_1_020_rsa_srss_msc.f06`
  - `D:\18a\femap_RSA\problem_1_024_rsa_srss_msc.f06`
- SAP references:
  - `D:\18a\femap_RSA\Problem 1-020.pdf`
  - `D:\18a\femap_RSA\Problem 1-024.pdf`
  - `D:\18a\femap_RSA\Example 1-024.s2k.txt`

## Problem 1-020

### Independent benchmark target

From the SAP2000 / Chopra comparison pages for Example 1-020:

- periods:
  - mode 1 period = `1.562 s`
  - mode 2 period = `0.5868 s`
- response-spectrum displacement:
  - `Ux(joint 2) = 7.566 in` in the independent column
  - `Ux(joint 3) = 18.81 in` in the independent column
- bending moments `M33` from the independent column:
  - elem 1 @ joint 1 = `12624`
  - elem 1 @ joint 2 = `6792`
  - elem 2 @ joint 2 = `6024`
  - elem 2 @ joint 3 = `5220`
  - elem 5 @ joint 2 = `9792`
  - elem 6 @ joint 5 = `9792`
  - elem 7 @ joint 3 = `5220`
  - elem 8 @ joint 6 = `5220`

### Modal check

MYSTRAN and MSC match in modes:

- mode 1 cycles:
  - MYSTRAN `6.401494E-01`
  - MSC `6.401511E-01`
- mode 2 cycles:
  - MYSTRAN `1.704106E+00`
  - MSC `1.704109E+00`

These also match the older SAP-validated note:

- `T1 = 1.562131 s`
- `T2 = 0.586817 s`

Current MYSTRAN RSA periods also match the independent target:

- mode 1:
  - MYSTRAN `1 / 6.401494E-01 = 1.562135 s`
  - independent `1.562 s`
- mode 2:
  - MYSTRAN `1 / 1.704106E+00 = 0.586818 s`
  - independent `0.5868 s`

### Response displacement

Key X translations:

- grid 2:
  - MYSTRAN `7.582045`
  - MSC `4.205205E-02`
  - old SAP-style validated value `7.576084`
  - independent target `7.566`
- grid 3:
  - MYSTRAN `1.885334E+01`
  - MSC `1.068353E-01`
  - old SAP-style validated value `1.883852E+01`
  - independent target `1.881000E+01`

Relative to the SAP-validated MYSTRAN baseline:

- grid 2 error = `+0.079%`
- grid 3 error = `+0.079%`

Relative to the independent benchmark:

- grid 2 error = `+0.212%`
- grid 3 error = `+0.230%`

Conclusion:

- the current MYSTRAN RSA result is now tightly aligned with the older direct SAP validation
- the remaining drift is small and is within the same range as the SAP vs independent comparison shown in the benchmark pages
- the MSC-normalized comparison deck is still not numerically equivalent in RSA scaling, even though the modal solution matches

### Internal forces

Representative bending results compared to the older SAP-validated MYSTRAN benchmark:

- element 1:
  - MYSTRAN `M2A = 1.264577E+04`
  - SAP-validated reference `1.263583E+04`
  - error = `+0.079%`
  - MYSTRAN `M2B = 6.798488E+03`
  - SAP-validated reference `6.793146E+03`
  - error = `+0.079%`
- element 2:
  - MYSTRAN `M2A = 6.027772E+03`
  - SAP-validated reference `6.023038E+03`
  - error = `+0.079%`
  - MYSTRAN `M2B = 5.226572E+03`
  - SAP-validated reference `5.222469E+03`
  - error = `+0.079%`
- element 5:
  - MYSTRAN `M2A = 9.818062E+03`
  - SAP-validated reference `9.810580E+03`
  - error = `+0.076%`

Relative to the independent benchmark table:

- elem 1 @ joint 1:
  - MYSTRAN `12645.77`
  - independent `12624`
  - error = `+0.172%`
- elem 1 @ joint 2:
  - MYSTRAN `6798.488`
  - independent `6792`
  - error = `+0.096%`
- elem 2 @ joint 2:
  - MYSTRAN `6027.772`
  - independent `6024`
  - error = `+0.063%`
- elem 2 @ joint 3:
  - MYSTRAN `5226.572`
  - independent `5220`
  - error = `+0.126%`
- elem 5 @ joint 2:
  - MYSTRAN `9818.062`
  - independent `9792`
  - error = `+0.266%`
- elem 6 @ joint 5:
  - MYSTRAN `9818.062`
  - independent `9792`
  - error = `+0.266%`
- elem 7 @ joint 3:
  - MYSTRAN `5226.572`
  - independent `5220`
  - error = `+0.126%`
- elem 8 @ joint 6:
  - MYSTRAN `5226.572`
  - independent `5220`
  - error = `+0.126%`

MSC force levels are much smaller:

- MSC element 1 `M2A = 6.861236E+01`
- MSC element 1 `M2B = 3.345469E+01`

Conclusion:

- `problem_1_020` now tracks the prior SAP-style validation path closely in both displacement and bending moments
- the integrated RSA path still gets periods right, and after the TABLED1 fix the displacement and force levels are no longer materially weaker than the older validated path
- the commercial MSC baseline is useful for parser/format comparison here, but not as an equivalent RSA scaling reference

### Root cause that was fixed

The main break in `problem_1_020` was not the SRSS combiner itself. The issue was in `BD_TABLED1`:

- the prior parser only consumed two `(x,y)` spectrum pairs per continuation line
- the benchmark deck uses up to four `(x,y)` pairs per continuation line
- this dropped the plateau points around the second-mode period
- the result was an under-read spectral displacement for mode 2 and therefore underpredicted SRSS displacement and moments

After extending `BD_TABLED1` to read all continuation-line pairs, the `SCRSPEC` summary in the F06 now shows:

- mode 1 `Sd = 1.375697E+01`
- mode 2 `Sd = 4.566921E+00`

That corrected the SRSS result back to the expected range.

## Problem 1-024

### Modal check

SAP model units in the `.s2k` file are:

- `CurrUnits="Kip, ft, F"`

SAP benchmark periods from the PDF:

- mode 1 `0.2271 s`
- mode 2 `0.2156 s`
- mode 3 `0.0733 s`
- mode 4 `0.0720 s`

MSC modes:

- `4.404087 Hz`
- `4.637502 Hz`
- `13.63416 Hz`
- `13.88791 Hz`

MYSTRAN modes:

- `4.404087 Hz`
- `4.637502 Hz`
- `13.63416 Hz`

Important difference:

- MSC extracts 4 modes
- current MYSTRAN run extracts only 3 modes

So the modal side is only partially matched:

- first 3 modes match exactly
- mode 4 is missing in the current MYSTRAN RSA run

### Response displacement

Target SAP quantity from the benchmark PDF:

- joint/grid 29 X displacement = `0.02012 ft`

Current results:

- MYSTRAN grid 29 `T1 = 2.414246E-01`
- MSC grid 29 `T1 = 4.957481E-04`

Ratios:

- MYSTRAN / SAP = `11.999`
- MSC / SAP = `0.02464`
- MYSTRAN / MSC = `486.99`

Conclusion:

- `problem_1_024` is not yet correct in the new RSA path
- this is not a small drift; it is a large scaling mismatch
- the first three modes are correct, so the problem is in RSA application/scaling, not in the basic modal solution

### Velocity and acceleration

Grid 29 combined values:

- MYSTRAN:
  - velocity `T1 = 6.688920E+00`
  - acceleration `T1 = 1.872786E+02`
- MSC:
  - velocity `T1 = 1.376000E-02`
  - acceleration `T1 = 3.916757E-01`

These show the same large scaling break seen in displacement.

### Internal forces

Representative element 1 comparison:

- MYSTRAN:
  - `M1A = 1.346909E+03`
  - `M2A = 6.143293E+01`
  - `FX  = 1.139921E+02`
- MSC:
  - `M1A = 2.781558E+00`
  - `M2A = 1.364042E-01`
  - `FX  = 2.343958E-01`

Ratios:

- `M1A ratio = 484.23`
- `M2A ratio = 450.37`
- `FX ratio  = 486.32`

Conclusion:

- internal-force scaling is broken in the same way as displacement/velocity/acceleration
- this is consistent with an RSA application factor problem, not a local force-recovery bug

## Practical status

### Reliable today

- `problem_1_020`:
  - modal solution is correct
  - RSA displacement and main bending-force levels remain close to the older SAP-validated path

### Not reliable yet

- `problem_1_024`:
  - first 3 modes are correct
  - mode 4 is missing in the current MYSTRAN RSA run
  - combined displacement, velocity, acceleration, and bar forces are over-scaled by about `450x` to `487x` vs MSC
  - grid 29 X displacement is `12x` above the SAP target

## Immediate next debug targets

- ensure `problem_1_024` extracts the fourth mode in the active RSA run
- audit the RSA application scale path after modal extraction:
  - spectrum ordinate scaling
  - gravity conversion
  - support-excitation participation scaling
  - modal combination accumulation
- compare `problem_1_024` against the older modal-only / postprocessed benchmark path before changing the modal solver
