# PBEAM Pretension Status In MYSTRAN 17

This note clarifies whether `PBEAM` in this codebase has a direct pretension input like some commercial Nastran variants.

## Current status

There is no dedicated `PBEAM` pretension/preload field in the current MYSTRAN source.

Code basis:

- `PBEAM` parsing is in [BD_PBEAM.f90](E:/mystran17/mystran/Source/LK1/L1A-BD/BD_PBEAM.f90)
- `RPBEAM` fields are geometry/stress/shear/warping/offset-related (A, I1, I2, I12, J, NSM, C/D/E/F, K1, K2, CW, etc.)
- no field is parsed or stored as beam axial pretension
- no `PRETENS`-style bulk card path was found in `Source/LK1/L1A-BD`

## Practical meaning

For MYSTRAN 17, `PBEAM` defines section properties; it does not define an initial axial force state by itself.

If an axial force state is needed in `SOL 101`, use load/path options such as:

- nodal `FORCE` in beam axis direction
- element distributed axial load `PLOAD1` with `TYPE=FXE`
- thermal loading path (for thermal strain-driven axial force)

These produce axial force through the load vector, not through a dedicated `PBEAM` pretension property.
