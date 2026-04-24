# Eigen13 Validation Proof (ARPACK vs CHASE vs FEAST)

## Scope

Compared decks:
- `examples/eigen13_arpack.dat`
- `examples/eigen13_chase.dat`
- `examples/eigen13_feast_true.dat`

Compared outputs:
- `examples/eigen13_arpack.F06/.ERR`
- `examples/eigen13_chase.F06/.ERR`
- `examples/eigen13_feast_true.F06/.ERR`

## Dispatch status

- ARPACK deck: baseline Lanczos.
- CHASE deck: fallback to ARPACK (ERR shows `WARNING 4916: CHASE NATIVE REQUIRES EIGRL FREQUENCY INTERVAL`).
- FEAST true deck: native FEAST active (ERR shows `WARNING 4918`, no `4916`).

## Frequency comparison

See full table in:
- `lancmeth_compare_eigen13.csv`

Summary:
- CHASE (current deck) matches ARPACK exactly due to fallback.
- FEAST native tuned is very close to ARPACK (max delta mode 1 = `0.015380%`).

## FEAST caution (important)

FEAST behavior is sensitive to:
- interval definition (`EIGRL V1/V2` must be valid),
- matrix storage flag passed to backend,
- mass regularization floor for massless rotational DOF.

Current stable setting in this package:
- FEAST uses full sparse storage path (`UPLO='F'` in source snapshot).
- conservative regularization floor (`MREG_FLOOR = max(EPS1, avg_mass_diag*1e-10)`).

If FEAST starts falling back again (`WARNING 4916`), first check:
1. deck interval is defined and realistic,
2. native FEAST library is actually linked,
3. no accidental changes on FEAST `FPM` controls.
