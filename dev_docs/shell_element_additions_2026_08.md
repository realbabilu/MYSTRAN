# Shell Element Additions - August 2026

## Commit

The main implementation was committed as:

```text
88c09188 Add shell element mass and pressure support
```

It was later pushed to branch `v18.00.a` through merge commit:

```text
03b67dd3 Merge origin/v18.00.a
```

## Scope

This package adds and wires several experimental shell element families, with
the immediate production focus on mass support for SOL 103 and pressure-load
support for SOL 101.

Stress/force recovery is only partially touched where needed by existing
output paths. It should not be treated as fully validated for all new element
families yet.

## Added or exposed shell formulations

Quadrilateral families:

- `CQUAD4_DKMT20`
- `CQUAD8_SIMOEAS1`
- `CQUADR_DKM24AU`
- `CQUADR_DKM24EA`
- `CQUADR_Q4EASANS`
- `CQUADR_Q4RS`

Triangular families:

- `CTRIA3_T3FF`
- `CTRIA6_SIMO1993`
- `CTRIAR_MITC3PHB`
- `CTRIAR_T3FFD`

Related interface modules and bulk-data readers were added for the new element
entry points, including `CTRIA6` sizing/read support.

## Mass support

The package extends shell mass handling so the new quadratic and selector-based
shell elements can participate in SOL 103. The intended control remains through
the existing mass mode path, including the `COUPMASS` parameter for consistent
versus lumped behavior where supported.

Validation run during development included:

- `CQUAD8` / `MITC8`, consistent mass with `PARAM,COUPMASS,1`
- `CQUAD8` / `MITC8`, lumped/default mass
- `CQUAD8` / `SIMOEAS1`, consistent mass with `PARAM,COUPMASS,1`
- `CTRIA6` / `SIMO1993`, consistent mass with `PARAM,COUPMASS,1`

The HW20 path was intentionally not included in the staged package.

## Pressure support

Pressure-load support was added or fixed for the quadratic shell paths:

- `CQUAD8` / `MITC8`
- `CQUAD8` / `SIMOEAS1`
- `CTRIA6` / `SIMO1993`

Important pressure-loader fixes:

- `PLOAD4` first pressure value now reads from field 4 directly.
- `PLOAD4` shell type classification recognizes `CTRIA*` and `CQUAD*`
  prefixes.
- `ELMDAT2` recognizes quadratic shell pressure data through both `TYPE` and
  `ETYPE`.
- `EMG` now calls `ELMDAT2` for `QUAD8`, so `PRESS` is populated before the
  element pressure vector is generated.

Post-merge pressure verification:

```text
pressure_my_cquad8_mitc8.dat     -> PN=8, PL=4, normal termination
pressure_my_cquad8_simoeas1.dat  -> PN=8, PL=4, normal termination
pressure_my_ctria6_simo.dat      -> PN=6, PL=3, normal termination
```

Reaction sanity checks before push:

```text
CQUAD8 unit square pressure: SPC total Z = -1.000000E+00
CTRIA6 half-square pressure: SPC total Z = -5.000000E-01
```

## Stress and force recovery status

Stress/force recovery is not globally complete for every newly added element.
Current notes:

- Existing validated paths such as `CQUAD4`/Simo and `CTRIA3`/T3FF remain the
  reference behavior for their families.
- `GPSTRESS` and element force recovery still require element-by-element audit
  for the new quadratic and selector formulations.
- DKMQ24/DKMT18 moment recovery has a separate audit note:
  `dev_docs/shell_moment_recovery_dkmq24_dkmt18_audit.md`.

## Recommended next validation

Before declaring the whole element package complete:

- rerun SOL 103 rigid-body/eigen checks for Q8/T6 in consistent and lumped mass
  modes;
- rerun SOL 101 pressure patch tests for CQUAD8 and CTRIA6;
- audit SOL 101 element force recovery and GPSTRESS per formulation;
- audit SOL 105 KED/geometric stiffness support separately, because it is not
  covered by the mass/pressure verification above.
