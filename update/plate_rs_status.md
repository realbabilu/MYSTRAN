# Plate RS Status

Date: July 14, 2026

## Active files

- RSA deck:
  - `D:\18a\femap_RSA\plate_RS.dat`
- Matching modal-only deck derived from the same structure:
  - `D:\18a\femap_RSA\plate_RS_modal_from_rsa_dense.dat`

## Modal-basis check

The modal-only deck was created from the active RSA deck by removing:

- `PARAM,SCRSPEC`
- `PARAM,OPTION`
- `DLOAD`
- `DTI,SPECSEL`
- `TABLED1`
- `TABDMP1`
- `SUPORT`
- RSA-only velocity/acceleration output requests

Observed result:

- `plate_RS.dat` and `plate_RS_modal_from_rsa_dense.dat` produce the same extracted modal basis in the active `1..200 Hz` window.
- Extracted modes:
  - mode 1: `9.851300 Hz`
  - mode 2: `60.98133 Hz`
  - mode 3: `121.8857 Hz`
  - mode 4: `164.0251 Hz`

This means the current `MODES/SEMODES + SCRSPEC` isolation patch is behaving correctly on the plate model too.

## Current RSA observation

The current RSA summary in `plate_RS.F06` reports:

- `SUPORT grid=28155, component=2`
- fallback global-component weighting
- `max(|UG|)=1.556642E-12`

The retained-mode active participation values are also tiny:

- mode 1 gamma: `6.623102E-12`
- mode 2 gamma: `1.061706E-12`
- mode 3 gamma: `3.079759E-13`
- mode 4 gamma: `-5.913904E-13`

Interpretation:

- the model is no longer suffering from the earlier `SUPORT -> structural R-set` corruption
- the present near-zero RSA response is instead tied to the actual excitation-direction/participation content of this plate setup

## Immediate consequence

For `plate_RS`, the next debugging target is not eigen extraction.

The next target is:

- why the chosen `SUPORT` direction on this model produces near-zero modal participation in the current MYSTRAN RSA path
- and whether that is expected from the deck definition or indicates a remaining support-direction mapping issue
