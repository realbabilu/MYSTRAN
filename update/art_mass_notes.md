# ART_MASS Notes

This note records the intended `mystran18a` behavior for `PARAM,ART_MASS`.

## Goal

Keep `ART_MASS` as an explicit compatibility trigger, but make the common case simpler and closer to typical Nastran usage.

## Current Intended Semantics

Primary usage:

- `PARAM,ART_MASS,Y`

Meaning:

- enable artificial mass insertion
- reset both fallback values to `1.0E-6`
- use:
  - translational artificial mass = `1.0E-6`
  - rotational artificial mass = `1.0E-6`

Advanced usage:

- `PARAM,ART_MASS,Y,<trans_mass>,<rot_mass>`

Meaning:

- enable artificial mass insertion
- override the default fallback values explicitly

Disable:

- `PARAM,ART_MASS,N`

Meaning:

- do not add artificial mass terms

## Why This Behavior

The branch already stores:

- `ART_TRAN_MASS = 1.0E-6`
- `ART_ROT_MASS  = 1.0E-6`

but these values only matter when `ART_MASS='Y'`.

The important cleanup is:

- `PARAM,ART_MASS,Y` should not accidentally inherit stale override values from an earlier `PARAM,ART_MASS,Y,x,y`
- plain `Y` should always mean a clean fallback to `1.0E-6`

## Scope of Use

Artificial mass is inserted in the concentrated/grid mass assembly path to prevent zero diagonal mass terms on grid DOF where needed.

This is a numerical compatibility/stabilization feature, not a silent always-on physics change.

## Compatibility Position

Recommended user-facing rule:

1. Keep `PARAM,ART_MASS,Y` as the explicit trigger.
2. Treat fields 4 and 5 as advanced override fields.
3. Do not automatically force artificial mass for every run without a trigger.

That keeps behavior predictable for:

- normal modes
- buckling
- modal participation / effective mass reporting
- decks that intentionally want to expose singular or massless behavior
