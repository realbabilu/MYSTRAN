# dev.cons - CHASE/FEAST Add-on

This note captures developer considerations specific to the CHASE/FEAST addition package.

## 1) Scope boundary

- This package is for Lanczos-family routing only (`EIGRL` path).
- `EIGR/MGIV` must not be treated as FEAST/CHASE execution.
- Keep warning guard in LINK4 so users know when `LANCMETH` is ignored.

## 2) Safety rule for FEAST

- FEAST is enabled as experimental-native path with fallback.
- On native failure, fallback target should remain ARPACK Lanczos for compatibility and predictability.
- Do not silently switch to MGIV in this package.

## 3) Parameter stability

- FEAST controls are sensitive; aggressive tuning can trigger backend error codes.
- Keep FEAST control values conservative unless validated against regression set.
- Any tuning change should be validated at least on:
  - Eigen8 quick-pass set
  - Eigen13 modal shell case

## 4) Evidence-first updates

- Every functional change should update proof files in this folder:
  - `validation_eigen8_proof.md`
  - `validation_eigen13_proof.md`
  - `lancmeth_compare_eigen13.csv`

## 5) Snapshot sync rule

- If source changes in core files, mirror updated files into:
  - `howbuild/chase_feast_add/mystran/...`
- Keep snippet log updated in:
  - `snippets_chase_feast_add.md`

This keeps package reproducible when shared outside the main tree.
