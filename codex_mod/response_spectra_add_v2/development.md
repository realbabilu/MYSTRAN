# development.md

## Scope

This archive captures the current response-spectrum implementation state in `mystran3`.

The main active areas are:

- parser and storage support for response-spectrum bulk data
- `PARAM,RSTYPE,...` handling
- `MFREQ` solution-family routing
- modal participation / effective-mass table output in F06
- one-direction benchmark work for the SAP verification problems

## What is known to work

### Problem 1-020

- Modal basis was matched after bringing the SAP axial-modification effect into the MYSTRAN model.
- One-direction response spectrum in X was made to produce the target displacement results for the benchmark.
- This is the cleanest stage-1 confirmation that the `MFREQ` path can run and produce valid RS results.

### Problem 1-021

- This is eigenvalue-only, not response spectrum.
- The modal results were matched after correcting the restraint representation in the MYSTRAN deck.

### Problem 1-024

- This benchmark is one-direction response spectrum in global X only.
- The compared methods are modal-combination rules, not directional combinations:
  - `SRSS`
  - `CQC`
  - `ABS`
  - `NRC 10 Percent`
- `SRSS`, `CQC`, and `ABS` are active in the current implementation.
- `10 Percent` logic was also added, but the benchmark case itself may collapse to the same result as `SRSS` when no mode pair lies within the 10 percent frequency cluster rule.
- The model-equivalence work for periods and story-mass representation remains delicate.

## Problem 1-025 status

`Problem 1-025` is not yet closed.

Important points:

- It is still a one-direction response-spectrum benchmark, not a true two-direction combination benchmark.
- It adds torsional coupling through translational mass plus `MMI3` at the diaphragm master points.
- The current MYSTRAN experiments reproduced:
  - rigid diaphragm constraints,
  - master concentrated mass,
  - release-aware frame trials,
  but the model still behaved like a mechanism / near-zero-mode system.

### Why it is still open

The issue is no longer just one missing parser or one missing RS routine.

The remaining uncertainty is in the structural equivalence of the SAP frame model to the simplified MYSTRAN frame deck, especially the interaction of:

- rigid diaphragms,
- concentrated story mass / inertia,
- frame releases,
- and the chosen bar/beam representation.

### Release trial note

A release-aware deck was generated using SAP `M2/M3` release data mapped to MYSTRAN pin flags (`56` at the corresponding end). That trial ran, but it did not remove the near-zero modal behavior.

So the release mapping was worth testing, but it was not the final missing piece.

## Included design notes

- `mystran_response_spectrum_load2mass_spec.md`
- `sol111_vs_sol103_vs_sol112.md`
- `problem_1_025_stage1_note.md`

These three notes are the main conceptual breadcrumbs for the current implementation and stage-1 verification path.
