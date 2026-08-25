# Shell Bending Patch Sign Audit - 2026-08-25

Scope: `python/linear/battle_2_001_patch.py` against the generated MYSTRAN F06 files in
`python/linear/working_mystran_2_001`, using MSC/Nastran as the sign convention reference.

## Decisions

- `MITC4PD` in the patch decks is the `CQUADR_MITC4PHB` implementation path. `PARAM,QUADRTYP,MITC4PD`
  dispatches to `CQUADR_MITC4PHB`.
- `Q4RS` is acceptable for the GPSTRESS comparison target. Some local force-table `MXY` rows are mixed by
  element coordinate orientation, but GPSTRESS is uniform and matches the Nastran-style patch result.
- `DSQK` is not treated as a simple sign flip. Its center/element-force behavior is in the same broad class
  as the Python DSQK reference, but its GPSTRESS bending shear was much worse when the distorted DSQK corner
  bending recovery was averaged as nodal stress.
- Original/older shell paths such as `MIN3`, `MIN4`, `MIN4T`, `MITC4`, and `MITC4+` are intentionally left
  alone for this pass.

## Fix Applied

- `battle_2_001_patch.py` now writes LC3 constant-shear patch decks with components `45` constrained at
  every distorted-patch node, matching the Python reference field. The boundary nodes still receive the
  component `3` SPC/SPCD transverse displacement values. This removes the artificial LC3 flexibility that
  made `CQUADR/SIMO` report `w=0.005780553` instead of the Python/reference `w=0.0055`.
- `TPLT_MITC3P` no longer flips the rotational recovery columns after forming `BB_REC` and `BS_REC`.
  The stiffness path still applies its physical rotation convention through `KPHYS`; only the stress/force
  recovery double flip was removed.
- `CQUAD4_DSQK_RHR` now feeds DSQK center bending recovery into the bending corner recovery rows for GPSTRESS
  use, avoiding the inflated corner `MXY` values. `OFP3_STRE_NO_PCOMP` also treats DSQK recovery rows like
  direct output rows instead of running them through the generic quad polynomial fit.

## Verification Snapshot

- `MITC3+` membrane GPSTRESS: `sxx=1333.333`, `syy=1333.333`, `sxy=400.0`, about `0.0167%` error.
- `MITC3+` LC2 GPSTRESS comparison: `mxx=1.111111166e-07`, `myy=1.111111166e-07`,
  `mxy=3.333333333e-08`, about `0.0400%` error in the Python absolute-valued comparison.
- Raw MYSTRAN `MITC3+` force output follows the Nastran sign convention at the matching reference rows:
  `MX=-1.111111E-07`, `MY=-1.111111E-07`, `MXY=+3.333333E-08`.
- `DSQK` membrane GPSTRESS remains `sxx=1333.333`, `syy=1333.333`, `sxy=400.0`, about `0.0167%` error.
- `DSQK` LC2 GPSTRESS comparison improved from the bad `273%` case to about `22.17%`:
  `mxx=1.046709479e-07`, `myy=1.030220062e-07`, `mxy=3.616032083e-08`.
- Regenerated `prob_2_001_patchshear_cquadr_simo.dat` with the LC3 all-node rotation constraints and reran
  MYSTRAN: `SIMO93` LC3 now gives `w_internal=0.0055`, `ratio=1.0`, `0.0%` error. LC4/LC5 were already exact.
