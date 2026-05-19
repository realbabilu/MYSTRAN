# RESPONSE_SPECTRA_ADD

This folder bundles the current response-spectrum work from `mystran3`.

It is intended as an internal codex archive focused on:

- bulk-data intake for response-spectrum cards
- `RSTYPE` handling (`PERG`, `PERA`, `FRQG`, `FRQA`)
- `SOL103` / `SOL111` / `SOL112` modal-to-RS workflow
- `MFREQ` runtime plumbing in `LINK4`, `LINK5`, and `LINK9`
- F06 / OP2 participation and effective-mass output additions
- stage-1 mass-path work that affects RS verification models

## Main subfolders

- `source_snapshot/`
  Current Fortran files most directly tied to the response-spectrum implementation.
- `docs/`
  Current design and benchmark notes.

## Key source files

1. `RESPONSE_SPECTRA_STUF.f90`
2. `BD_DAREA.f90`
3. `BD_DLOAD.f90`
4. `BD_FREQ1.f90`
5. `BD_RLOAD1.f90`
6. `BD_TABLED1.f90`
7. `LOADE.f90`
8. `LINK4.f90`
9. `LINK5.f90`
10. `LINK9.f90`
11. `OFP2.f90`
12. `WRITE_GRD_OP2_OUTPUTS.f90`
13. `WRITE_MPFACTOR.f90`
14. `WRITE_MEFFMASS.f90`
15. `OUTPUT2_WRITE_TABLE.f90`
16. `MGGS_MASS_MATRIX.f90`

## Notes

- This archive intentionally focuses on response-spectrum and stage-1 verification plumbing.
- It does not try to capture the separate warning-reduction sweep.
- `Problem 1-025` is still open; the included note documents what has already been tried.
