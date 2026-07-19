# PSHELLZ Concept

## Purpose

`PSHELLZ` is a proposed shell-property extension for per-mode stiffness modifiers, similar in intent to MIDAS-style shell stiffness scaling.

This is a concept note only. It is not ready for implementation or production use yet.

## Proposed syntax

```text
PSHELLZ,PID,MID1,T,MID2,12I/T**3,MID3,TS/T,MID4
+,STIFFMOD,FXX,FYY,FXY,MXX,MYY,MXY,VXZ,VYZ
```

## Modifier meaning

- `FXX` = in-plane axial stiffness scale in shell local `x`
- `FYY` = in-plane axial stiffness scale in shell local `y`
- `FXY` = in-plane in-plane shear stiffness scale
- `MXX` = bending stiffness scale associated with local `x`
- `MYY` = bending stiffness scale associated with local `y`
- `MXY` = twisting / coupled bending stiffness scale
- `VXZ` = transverse shear stiffness scale in local `xz`
- `VYZ` = transverse shear stiffness scale in local `yz`

Default intent:

- every modifier defaults to `1.0`
- `PSHELL` remains unchanged
- `PSHELLZ` affects stiffness only, not mass, density, NSM, or gravity

## Design rules

- Scale constitutive stiffness terms, not geometry fields directly.
- Do not reinterpret thickness or material density as part of stiffness modification.
- Keep local-axis semantics explicit; all factors apply in shell local coordinates.
- Keep composite / laminate support out of the first implementation phase.

## Suggested rollout

1. Add parser and metadata storage only.
2. Implement stiffness scaling for:
   - `CQUAD4`
   - `CQUADR`
   - `CTRIA3`
   - `CTRIAR`
3. Validate on small scalar benchmarks before extending to other shell families.

## Minimum validation before coding

- membrane-only patch test
- bending-only plate test
- thick-shell transverse shear test
- local-axis rotation / orthotropic sanity test
- check that mass and self-weight are unchanged

## Status

Deferred. No validated benchmark set exists yet, so `PSHELLZ` should remain a design note until test coverage is ready.
