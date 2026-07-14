# Response Spectrum Beam Validation

Date: July 14, 2026

## Scope

This note records the current beam/frame response-spectrum validation status for the
commercial-Nastran-compatible path based on:

- `TABLED1`
- `DTI,SPECSEL`
- `SUPORT`
- `PARAM,SCRSPEC`
- `PARAM,OPTION`

Primary working folder:

- `D:\18a\femap_RSA`

Reference file:

- `D:\18a\femap_RSA\beam_rs_msc.f06`

## Recommended validation decks

Use these as the current beam/frame RSA validation set:

- `beam_RS.dat`
  - primary SRSS validation deck
  - `SOL SEMODES`
  - requests combined displacement, velocity, acceleration, stress, SPC force, and ESE
- `beam_RS_plot.dat`
  - same as `beam_RS.dat` but with `DISPLACEMENT(PRINT,PLOT)=ALL`
  - use this to confirm the OP2/NEU plot-output route stays alive
- `beam_RS_cqc.dat`
  - local combination check for `PARAM,OPTION,CQC`
  - requires its own commercial reference if numerical comparison is needed
- `beam_RS_cqc_plot.dat`
  - same as `beam_RS_cqc.dat` with plot-route output enabled
- `beam_RS_abs.dat`
  - local combination check for `PARAM,OPTION,ABS`
  - requires its own commercial reference if numerical comparison is needed
- `beam_RS_easy.dat`
  - narrow smoke/check deck for the simplified compatibility path
  - not the right deck for direct comparison against `beam_rs_msc.f06`

## Current numerical baseline vs MSC

The useful direct comparison is:

- reference: `beam_rs_msc.f06`
- test: `beam_RS.F06`
- test: `beam_RS_plot.F06`

Current beam SRSS comparison summary:

- mode 1 eigenvalue:
  - MSC: `1.017487E+03`
  - MYSTRAN: `1.017683E+03`
  - relative gap: about `1.93e-4`
- mode 2 eigenvalue:
  - MSC: `6.011450E+04`
  - MYSTRAN: `6.082800E+04`
  - relative gap: about `1.19e-2`
- combined displacement at `GRID 2`:
  - `TY`: MSC `4.173022E-01`, MYSTRAN `4.172221E-01`
  - `RZ`: MSC `1.041883E-01`, MYSTRAN `1.041884E-01`
- combined velocity at `GRID 2`:
  - `TY`: MSC `1.331114E+01`, MYSTRAN `1.330986E+01`
  - `RZ`: MSC `3.323408E+00`, MYSTRAN `3.323730E+00`
- combined acceleration at `GRID 2`:
  - `TY`: MSC `4.246000E+02`, MYSTRAN `4.246000E+02`
  - `RZ`: MSC `1.060104E+02`, MYSTRAN `1.060308E+02`

The largest reported relative errors from the generic comparator are dominated by
near-zero terms where MSC prints tiny nonzero values or exact zeroes. For this beam
deck, the meaningful active response DOFs are the ones above.

`beam_RS_plot.F06` currently matches `beam_RS.F06` numerically for the beam RSA vector
outputs, which is the expected result. The difference is only that the plot deck also
routes the combined displacement payload to OP2/NEU.

## Validation-parser fix applied

To make the beam RSA validation usable again, `D:\18a\MYSTRAN_Validation-main\f06_query.py`
was updated to recognize MYSTRAN's combined response-spectrum vector blocks:

- `V E L O C I T I E S`
- `A C C E L E R A T I O N S`

It also now avoids trying to parse the summary line
`RESPONSE SPECTRUM ... VECTOR` as if it were the numeric table body, and waits for the
actual detailed MYSTRAN table header instead.

Without that parser fix, beam RSA velocity and acceleration looked "missing" in the
validation tooling even though MYSTRAN had written them to F06 correctly.

## Practical use

For frame/beam regression work, use this order:

1. `beam_RS.dat` for the main commercial-style SRSS check
2. `beam_RS_plot.dat` for the same check plus OP2/NEU routing
3. `beam_RS_cqc.dat` and `beam_RS_abs.dat` for combination-method sanity checks
4. `beam_RS_easy.dat` only as a narrow smoke deck, not as the primary commercial comparison

## July 14, 2026 frame RSA status

The commercial-style `SOL SEMODES + SCRSPEC` path now also emits a combined
`ELFORCE(ENGR)` F06 summary for 1D members:

- heading:
  - `RESPONSE SPECTRUM SRSS COMBINED ELEMENT ENGINEERING FORCE SUMMARY (1D)`
- current scope:
  - `BAR`
  - `BEAM`
  - `ROD`
  - `BUSH`
  - `ELAS`
- current combination support:
  - `SRSS`
  - `ABS`
  - `CQC` still reports unsupported for this 1D force summary

The active frame decks:

- `D:\18a\femap_RSA\problem_1_020_rsa_srss_new.dat`
- `D:\18a\femap_RSA\problem_1_024_rsa_srss_new.dat`
- `D:\18a\femap_RSA\problem_1_025_rsa_srss_new.dat`

were updated to request:

- `ELFORCE(PRINT) = ALL`

Observed status after the latest scaling fix:

- `problem_1_020_rsa_srss_new.F06`
  - acceleration-spectrum scaling now follows the SAP2000/S2K interpretation again
    for these inch-second-pound decks:
    - `PARAM,GRAV = 386.098`
    - `PARAM,RSTYPE = PERG`
    - `TABLED1` ordinates treated as spectral acceleration in `g`
  - combined displacement now lands on the expected SAP2000 scale
  - `GRID 2 / T1 = 7.478772E+00`
  - `GRID 3 / T1 = 1.882860E+01`
  - representative force rows now also scale correctly:
    - `BAR 5 M1b = 9.810184E+03`
    - `BAR 6 M2b = 9.810184E+03`
  - the BAR 1D summary labels were corrected to match the legacy writer column order:
    - col 2 = `M2a`
    - col 3 = `M1b`
  - `DEBUG,56,1` audit shows the residual force mismatch is not coming from the RSA
    combiner itself:
    - `UEL/PEL` are internally consistent for the modal BAR recovery path
    - the current RSA deck uses `CBAR/PBAR`
    - the older `1-020` package that was previously marked validated uses
      `CBEAM/PBEAM`
  - because of that, the remaining moment gap should be treated first as a
    formulation/deck-parity issue, not immediately as an RSA algorithm defect
  - residual gaps remain on some member-end moments, so this path is now usable for
    engineering comparison but not yet a final exact-match validation replacement for
    the older deterministic post-processing package

Additional deck prepared for apples-to-apples comparison with the older validation
package:

- `D:\18a\femap_RSA\problem_1_020_rsa_srss_cbeam_ref.dat`
  - same commercial-style RSA control path
  - switches the frame members to `CBEAM/PBEAM` so force comparison can be made
    against the older validated `1-020` package without mixing BAR and BEAM
    formulations
- `problem_1_024_rsa_srss_new.F06`
  - writes both combined 1D force summary and combined displacement block
- `problem_1_025_rsa_srss_new.F06`
  - writes both combined 1D force summary and combined displacement block

## July 14, 2026 `1-025` CBEAM parity trial

An isolated parity deck was added:

- `D:\18a\femap_RSA\problem_1_025_rsa_srss_cbeam_ref.dat`

Purpose:

- try the `1-025` RSA path with `CBEAM/PBEAM`
- avoid mixing the current `CROD` surrogate deck with a beam-style reference model

Observed result:

- the deck runs to completion
- modal frequencies are healthy and match the expected older surrogate values:
  - mode 1: `3.059159 Hz`
  - mode 2: `3.118760 Hz`
- but the deck still carries some cards from the older axial-surrogate package that are
  not yet part of the clean RSA compatibility path:
  - `PARAM,GRAV` is warned as ignored in this deck form
  - `DAREA` entries are warned as not processed

Interpretation:

- this confirms the `CBEAM/PBEAM` structure itself is numerically stable for the modal
  part
- it does **not** yet provide a clean RSA comparison baseline, because the excitation
  path is still contaminated by legacy surrogate cards
- therefore the next safe comparison basis for `1-025` should be:
  - keep the new RSA loading path
  - rebuild the element topology as `CBEAM/PBEAM`
  - remove the legacy `DAREA/RLOAD1` surrogate remnants that are not needed by the
    current MYSTRAN commercial-style RSA route

## July 14, 2026 `1-025` CBEAM deck cleanup

`D:\18a\femap_RSA\problem_1_025_rsa_srss_cbeam_ref.dat` was then cleaned to use the
same narrow commercial-style RSA input route as the active `SCRSPEC` beam decks:

- `SDAMP = 4`
- `DTI,SPECSEL -> TABLED1 1`
- `DLOAD,30,1.0,1.0,1`
- removed legacy `RLOAD1/DAREA` excitation remnants

Observed result after cleanup:

- modal basis remains the validated surrogate basis:
  - mode 1: `3.059159 Hz`
  - mode 2: `3.118760 Hz`
- `SCRSPEC` summary now reports:
  - mode 1 `Sd = 8.471552E-01`
  - mode 2 `Sd = 7.992520E-01`
- combined displacement at `GRID 29` becomes:
  - `T1 = 2.157005E-01`
  - `T2 = 2.573753E-01`

Interpretation:

- the huge `~8.33E+01` displacement level from the mixed legacy deck was an input-path
  contamination issue, not a modal extraction issue
- the cleaned `CBEAM` surrogate is now numerically sane enough to compare against SAP or
  MSC/NX references
- it is still a surrogate beam/release model, so any remaining gap against the original
  benchmark should be treated as model-interpretation drift first, not immediately as an
  RSA combiner bug

## July 14, 2026 `1-025` SRSS benchmark check against the SAP PDF

The local PDF reference:

- `D:\18a\femap_RSA\Problem 1-025.pdf`

states the following `SRSS` benchmark targets for the axial-only braced-frame
interpretation:

- modal frequencies:
  - mode 1: `3.0592 Hz`
  - mode 2: `3.1188 Hz`
- roof center-of-mass displacement at joint/grid `51`:
  - `X = 0.7372 in`
  - `Y = 0.7372 in`
  - `RZ = 0.000252 rad`
- representative frame axial forces:
  - element `1`: `200.55 kips`
  - element `4`: `139.57 kips`
  - element `6`: `86.48 kips`

Current MYSTRAN result from the cleaned
`D:\18a\femap_RSA\problem_1_025_rsa_srss_cbeam_ref.dat`:

- modal frequencies:
  - mode 1: `3.059159 Hz`
  - mode 2: `3.118760 Hz`
- joint/grid `51` combined displacement:
  - `T1 = 7.371697E-01 in`
  - `T2 = 7.371697E-01 in`
  - `R3 = 2.518462E-04 rad`
- combined 1D engineering-force summary:
  - `BEAM 1 FX = 2.005469E+02 kips`
  - `BEAM 4 FX = 1.395693E+02 kips`
  - `BEAM 6 FX = 8.648336E+01 kips`

Practical conclusion:

- the cleaned `1-025` `CBEAM/PBEAM` `SRSS` path now reproduces the SAP benchmark to
  engineering precision
- this validates the commercial-style `TABLED1 + DTI,SPECSEL + SUPORT + SCRSPEC`
  application path on this benchmark, provided the model uses the axial-only beam
  surrogate that matches the published benchmark assumptions

## July 14, 2026 `1-025` CQC and ABS check against the SAP PDF

Additional MYSTRAN runs were executed on:

- `D:\18a\femap_RSA\problem_1_025_rsa_cqc_cbeam_ref.dat`
- `D:\18a\femap_RSA\problem_1_025_rsa_abs_cbeam_ref.dat`

For `CQC`, the roof center-of-mass result at `GRID 51` matches the SAP benchmark:

- displacement:
  - `T1 = 1.032876E+00 in`
  - `T2 = 1.414383E-01 in`
  - `R3 = 2.518462E-04 rad`

For `ABS`, the roof center-of-mass result and representative axial forces also match:

- displacement at `GRID 51`:
  - `T1 = 1.042279E+00 in`
  - `T2 = 1.042279E+00 in`
  - `R3 = 2.518462E-04 rad`
- representative axial forces:
  - `BEAM 1 FX = 2.819903E+02 kips`
  - `BEAM 4 FX = 1.962493E+02 kips`
  - `BEAM 6 FX = 1.216048E+02 kips`

Practical conclusion:

- `SRSS`, `CQC`, and `ABS` now all reproduce the published `1-025` SAP benchmark on the
  cleaned MYSTRAN `CBEAM/PBEAM` surrogate path

## July 14, 2026 MSC comparison attempt for `1-025`

MSC-normalized deck copies were created:

- `D:\18a\femap_RSA\problem_1_025_rsa_srss_msc.dat`
- `D:\18a\femap_RSA\problem_1_025_rsa_cqc_msc.dat`
- `D:\18a\femap_RSA\problem_1_025_rsa_abs_msc.dat`

The following parser-level fixes were applied only to those copies:

- removed `MPFACTOR = YES`
- normalized `DTI,SPECSEL`
- normalized `TABDMP1` terminator syntax
- shifted `CBEAM` release codes so MSC does not read `45` as `OFFT`

Observed MSC behavior on `problem_1_025_rsa_srss_msc.dat`:

- case-control syntax can be normalized successfully
- bulk parsing can be normalized successfully
- MSC then aborts in element generation with repeated:
  - `USER FATAL MESSAGE 2040 (EMG) SINGULAR MATRIX FOR ELEMENT ...`

Interpretation:

- the current MYSTRAN axial-only `CBEAM/PBEAM` surrogate uses near-zero bending/torsion
  terms that MYSTRAN accepts for this benchmark path
- MSC does not accept that surrogate as a valid beam property set, so the run fails
  before any response-spectrum `DISPLACEMENT/VELOCITY/ACCELERATION` blocks are produced
- this is a surrogate-element compatibility limit, not evidence of an RSA combiner error

Next safe MSC comparison path:

- use the already validated MYSTRAN `1-025` results against the SAP PDF for
  `SRSS/CQC/ABS`
- if an MSC comparison is still required, build a separate MSC-safe surrogate
  (`CROD/PROD` or a non-singular frame property set) and compare global response
  quantities only
