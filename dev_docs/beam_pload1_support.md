# PLOAD1 Beam/Bar Support

`PLOAD1` is now supported in a limited first implementation for `CBAR`, `BART`, and `CBEAM`.

This v1 path converts the element line load into an equivalent nodal load vector `PPE`, then sends it through the normal `EPTL -> SYS_LOAD -> PG` load assembly path.

## Supported syntax

Use the standard small-field `PLOAD1` layout:

```text
PLOAD1  SID  EID  TYPE  SCALE  X1  P1  X2  P2
```

Supported values in this v1 implementation:

- `TYPE = FYE`, `FY`, or `Y`
  - uniform distributed load in the beam local `+y` direction
- `TYPE = FZE`, `FZ`, or `Z`
  - uniform distributed load in the beam local `+z` direction
- `TYPE = FXE`
  - uniform distributed load in the beam local `+x` direction
- `TYPE = MXE`
  - uniform distributed torsional moment about the beam local `+x` axis
- `TYPE = MYE`
  - distributed bending moment about the beam local `+y` axis
- `TYPE = MZE`
  - distributed bending moment about the beam local `+z` axis
- `SCALE = FR`
  - fractional position along the element
- `0.0 <= X1 <= 1.0`
- `0.0 <= X2 <= 1.0`
- `X2 >= X1`
- `P1` and `P2` may be different (linearly varying load over the span segment)

This now supports:

- full-span uniform load (`X1=0, X2=1, P1=P2`)
- partial-span uniform load (`X1<X2`, `P1=P2`)
- partial-span linearly varying load (`X1<X2`, `P1/=P2`)
- concentrated-in-element load as the degenerate case `X1 = X2`

## Not supported yet

- `LE` scale
- global-direction beam load types `FX/FY/FZ/MX/MY/MZ`
- projected-length scales `LEPR` and `FRPR`
- multiple `PLOAD1` entries for the same `(subcase, element, component)` combination

## Equivalent nodal load used

For each active component in `PLOAD1`, the code builds equivalent nodal loads in local beam axes and then uses the normal element-load transform/assembly path.

- distributed case (`X2 > X1`): numerical integration of shape functions over `[X1, X2]`
- concentrated case (`X2 = X1`): shape functions evaluated at the point location

For transverse load components, cubic Hermite beam functions are used for the `(U, R, U, R)` bending DOF sets.  
For axial force and torsional moment components, linear axial/torsion interpolation functions are used.

For a full-span uniform transverse line load `q` over beam length `L`, this reduces to the standard consistent beam load vector:

- local `y` load on DOFs `(UY, RZ, UY, RZ)`

```text
q * [ L/2, L^2/12, L/2, -L^2/12 ]
```

- local `z` load on DOFs `(UZ, RY, UZ, RY)`

```text
q * [ L/2, -L^2/12, L/2, L^2/12 ]
```

These are inserted into the 12-DOF beam/bar element vector and then transformed to basic/global coordinates by the usual MYSTRAN element-load path.

For the added full-span uniform local load types:

- local `x` force on DOFs `(UX, UX)`

```text
qx * [ L/2, L/2 ]
```

- local `x` torsional moment on DOFs `(RX, RX)`

```text
mx * [ L/2, L/2 ]
```

- local `y` distributed bending moment on DOFs `(RY, RY)`

```text
my * [ L/2, L/2 ]
```

- local `z` distributed bending moment on DOFs `(RZ, RZ)`

```text
mz * [ L/2, L/2 ]
```

## Notes

- In the current implementation, one entry per component is allowed per element and per subcase. A duplicate for the same component triggers an input error.
- Local-axis interpretation follows the existing MYSTRAN right-handed beam element system.
- For compatibility with earlier local testing, `FY/FZ` are still accepted as aliases for local beam-axis loads in this first implementation, even though commercial Nastran reserves `FYE/FZE` for element-axis loading.
- The current implementation stores beam/bar `PLOAD1` data per internal subcase, so different subcases may now carry different beam element loads without overwriting each other.
