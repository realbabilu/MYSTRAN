# RESPONSE_SPECTRA_ADD_V2

This folder consolidates the current response-spectrum support work from:

- `response_spectra_add`
- `response_spectrum_add`
- `response_spectrum_mystran_add`

## Purpose

Use this bundle as the single codex archive for response-spectrum work in `mystran3`.

It combines:

- solver-side response-spectrum source snapshots
- FEMAP/NEU postprocessing workflow
- SOL112-to-SOL111 compatibility tooling
- validation notes and status files
- stage-1 mass/load2mass notes that affected RS verification

## Layout

- `source_snapshot/`
  Solver-side Fortran files tied to response-spectrum support.
- `scripts/`
  Postprocessing and wrapper scripts.
- `tools/`
  Helper tools and compatibility scripts.
- `validation/`
  Validation notes and machine-readable status files.
- `docs/`
  Merged documentation from prior packs.
- `artifacts/`
  Sample generated output.

## Notes

- `response_spectrum_mystran_add` was the richer NEU/FEMAP workflow pack and is the main source for `scripts/`, `tools/`, and `validation/`.
- `response_spectra_add` remains the source for the solver `source_snapshot/` and stage-1 notes.
- Duplicate docs were preserved with suffixed filenames where needed.
